# =============================================================================
# 10_gap_fill_bd.R
# Gap-fill bulk density using Random Forest trained on measured values
# =============================================================================

library(terra)
library(dplyr)
library(randomForest)
library(caret)
library(ggplot2)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

INPUT_FILE <- "./Global_MRT_code/outputs/09_with_landcover.rds"
OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR <- "./Global_MRT_code/plots/10_bd_gapfill"

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# RF predictors (>90% coverage in target population)
RF_PREDICTORS <- c(
  # Terrain
  "terrain_ruggedness", "terrain_elev_mean", "terrain_slope_mean", 
  "terrain_northness", "terrain_eastness",
  # Climate
  "koppen_value", "koppen_main_num",
  # Land cover
  "lc_trees", "lc_grassland", "lc_shrubs", "lc_cropland", "lc_built",
  "lc_bare", "lc_snow", "lc_water", "lc_wetland", "lc_mangroves", "lc_moss",
  # Soil class
  "soilclass_wrb_value", "soilclass_func_code",
  # SoilGrids properties (not BD)
  "sg_clay", "sg_sand", "sg_cfvo", "sg_nitrogen", "sg_cec", "sg_phh2o",
  # Productivity
  "npp_mean", "npp_cv",
  # Disturbance
  "burn_count_total"
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# -----------------------------------------------------------------------------
# Cross-validation comparison: RF vs SoilGrids
# -----------------------------------------------------------------------------

run_cv_comparison <- function(sites, k = 10) {
  log_step(sprintf("Running %d-fold CV: RF vs SoilGrids", k))
  
  # Data with both measured BD and SoilGrids
  cv_data <- sites %>%
    filter(!is.na(bulk_density_oven_dry), !is.na(sg_bdod)) %>%
    select(bulk_density_oven_dry, sg_bdod, all_of(RF_PREDICTORS)) %>%
    na.omit()
  
  cat(sprintf("  CV samples: %s\n", format(nrow(cv_data), big.mark = ",")))
  
  set.seed(42)
  folds <- createFolds(cv_data$bulk_density_oven_dry, k = k)
  rf_preds <- rep(NA, nrow(cv_data))
  
  t0 <- Sys.time()
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    
    rf_model <- randomForest(
      bulk_density_oven_dry ~ . - sg_bdod,
      data = cv_data[-test_idx, ],
      ntree = 200
    )
    rf_preds[test_idx] <- predict(rf_model, cv_data[test_idx, ])
    
    cat(sprintf("\r  Fold %d/%d (%.1f min)", i, k, 
                difftime(Sys.time(), t0, units = "mins")))
  }
  cat("\n")
  
  # Calculate metrics
  obs <- cv_data$bulk_density_oven_dry
  sg <- cv_data$sg_bdod
  
  results <- data.frame(
    model = c("SoilGrids", "RF"),
    n = nrow(cv_data),
    R2 = c(cor(obs, sg)^2, cor(obs, rf_preds)^2),
    RMSE = c(sqrt(mean((obs - sg)^2)), sqrt(mean((obs - rf_preds)^2))),
    MAE = c(mean(abs(obs - sg)), mean(abs(obs - rf_preds))),
    bias = c(mean(sg - obs), mean(rf_preds - obs))
  )
  
  cat("\n  === CV Results ===\n")
  print(results)
  cat(sprintf("\n  RF improvement: ΔR² = +%.3f, ΔRMSE = %.3f g/cm³\n",
              results$R2[2] - results$R2[1],
              results$RMSE[2] - results$RMSE[1]))
  
  # Store predictions for plotting
  cv_predictions <- data.frame(
    observed = obs,
    soilgrids = sg,
    rf_cv = rf_preds
  )
  
  return(list(metrics = results, predictions = cv_predictions))
}

# -----------------------------------------------------------------------------
# Train final RF model
# -----------------------------------------------------------------------------

train_bd_model <- function(sites) {
  log_step("Training final RF model for bulk density")
  
  train_data <- sites %>%
    filter(!is.na(bulk_density_oven_dry)) %>%
    select(bulk_density_oven_dry, all_of(RF_PREDICTORS)) %>%
    na.omit()
  
  cat(sprintf("  Training samples: %s\n", format(nrow(train_data), big.mark = ",")))
  
  set.seed(42)
  model <- randomForest(
    bulk_density_oven_dry ~ .,
    data = train_data,
    ntree = 200,
    importance = TRUE,
    do.trace = 50
  )
  
  return(model)
}

# -----------------------------------------------------------------------------
# Apply gap-filling
# -----------------------------------------------------------------------------

apply_gapfill <- function(sites, model) {
  log_step("Applying gap-fill")
  
  # Initialize
  sites$bd_filled <- sites$bulk_density_oven_dry
  sites$bd_source <- ifelse(!is.na(sites$bulk_density_oven_dry), "measured", NA_character_)
  
  # RF predictions where possible
  needs_fill <- which(is.na(sites$bd_filled))
  can_predict <- complete.cases(sites[needs_fill, RF_PREDICTORS])
  fill_rf_idx <- needs_fill[can_predict]
  
  cat(sprintf("  Measured: %s\n", format(sum(!is.na(sites$bulk_density_oven_dry)), big.mark = ",")))
  cat(sprintf("  RF gap-fill: %s\n", format(length(fill_rf_idx), big.mark = ",")))
  
  if (length(fill_rf_idx) > 0) {
    sites$bd_filled[fill_rf_idx] <- predict(model, sites[fill_rf_idx, ])
    sites$bd_source[fill_rf_idx] <- "rf_model"
  }
  
  # SoilGrids fallback
  still_missing <- which(is.na(sites$bd_filled) & !is.na(sites$sg_bdod))
  cat(sprintf("  SoilGrids fallback: %s\n", format(length(still_missing), big.mark = ",")))
  
  if (length(still_missing) > 0) {
    sites$bd_filled[still_missing] <- sites$sg_bdod[still_missing]
    sites$bd_source[still_missing] <- "soilgrids"
  }
  
  return(sites)
}

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

create_bd_diagnostics <- function(sites, model, cv_results) {
  log_step("Creating diagnostics")
  
  # --- Save CV metrics (for paper) ---
  write.csv(cv_results$metrics, file.path(PLOT_DIR, "cv_comparison_metrics.csv"), row.names = FALSE)
  
  # --- Source distribution ---
  source_dist <- sites %>%
    count(bd_source) %>%
    mutate(pct = 100 * n / nrow(sites))
  
  write.csv(source_dist, file.path(PLOT_DIR, "bd_source_distribution.csv"), row.names = FALSE)
  cat("\n  BD source distribution:\n")
  print(source_dist)
  
  # --- Variable importance ---
  imp <- importance(model)
  imp_df <- data.frame(
    variable = rownames(imp),
    IncMSE = imp[, "%IncMSE"]
  ) %>% arrange(desc(IncMSE))
  
  write.csv(imp_df, file.path(PLOT_DIR, "bd_variable_importance.csv"), row.names = FALSE)
  
  cat("\n  Top 10 predictors:\n")
  print(head(imp_df, 10))
  
  # --- Importance plot ---
  p_imp <- ggplot(head(imp_df, 15), aes(x = reorder(variable, IncMSE), y = IncMSE)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
      title = "RF Variable Importance for Bulk Density",
      x = NULL, y = "% Increase in MSE"
    ) +
    theme_minimal()
  
  ggsave(file.path(PLOT_DIR, "bd_importance_plot.png"), p_imp, width = 8, height = 6, dpi = 200)
  
  # --- CV comparison scatterplots (for paper) ---
  cv_preds <- cv_results$predictions
  
  png(file.path(PLOT_DIR, "cv_comparison_scatterplots.png"), width = 2400, height = 1200, res = 200)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # SoilGrids
  sg_r2 <- cv_results$metrics$R2[1]
  sg_rmse <- cv_results$metrics$RMSE[1]
  plot(cv_preds$soilgrids, cv_preds$observed,
       pch = ".", col = rgb(0, 0, 0, 0.15),
       xlab = "SoilGrids BD (g/cm³)", ylab = "Measured BD (g/cm³)",
       main = sprintf("SoilGrids (R² = %.3f, RMSE = %.3f)", sg_r2, sg_rmse),
       xlim = c(0, 2.2), ylim = c(0, 2.2))
  abline(0, 1, col = "red", lwd = 2)
  
  # RF
  rf_r2 <- cv_results$metrics$R2[2]
  rf_rmse <- cv_results$metrics$RMSE[2]
  plot(cv_preds$rf_cv, cv_preds$observed,
       pch = ".", col = rgb(0, 0, 0, 0.15),
       xlab = "RF Predicted BD (g/cm³)", ylab = "Measured BD (g/cm³)",
       main = sprintf("Random Forest CV (R² = %.3f, RMSE = %.3f)", rf_r2, rf_rmse),
       xlim = c(0, 2.2), ylim = c(0, 2.2))
  abline(0, 1, col = "red", lwd = 2)
  
  dev.off()
  
  # --- Residual comparison ---
  cv_preds$resid_sg <- cv_preds$observed - cv_preds$soilgrids
  cv_preds$resid_rf <- cv_preds$observed - cv_preds$rf_cv
  
  png(file.path(PLOT_DIR, "cv_residual_histograms.png"), width = 2400, height = 1200, res = 200)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  hist(cv_preds$resid_sg, breaks = 100, col = "gray70",
       main = sprintf("SoilGrids Residuals (bias = %.3f)", cv_results$metrics$bias[1]),
       xlab = "Residual (Measured - Predicted)", xlim = c(-1.5, 1.5))
  abline(v = 0, col = "red", lwd = 2)
  
  hist(cv_preds$resid_rf, breaks = 100, col = "steelblue",
       main = sprintf("RF Residuals (bias = %.3f)", cv_results$metrics$bias[2]),
       xlab = "Residual (Measured - Predicted)", xlim = c(-1.5, 1.5))
  abline(v = 0, col = "red", lwd = 2)
  
  dev.off()
  
  cat(sprintf("\n  Diagnostics saved to: %s\n", PLOT_DIR))
}

# =============================================================================
# Main execution
# =============================================================================

log_step("Starting bulk density gap-filling (Step 10)")

# Load data
sites <- readRDS(INPUT_FILE)
cat(sprintf("  Loaded: %s sites, %d variables\n", 
            format(nrow(sites), big.mark = ","), ncol(sites)))

# CV comparison (RF vs SoilGrids)
cv_results <- run_cv_comparison(sites)

# Train final model
bd_model <- train_bd_model(sites)

# Apply gap-fill
sites <- apply_gapfill(sites, bd_model)

# Coverage summary
has_mrt_vars <- !is.na(sites$organic_carbon) & !is.na(sites$bd_filled) & !is.na(sites$npp_mean)
cat(sprintf("\nComplete cases for MRT: %s (%.1f%%)\n",
            format(sum(has_mrt_vars), big.mark = ","),
            100 * sum(has_mrt_vars) / nrow(sites)))

# Diagnostics
create_bd_diagnostics(sites, bd_model, cv_results)

# Save
log_step("Saving outputs")
saveRDS(bd_model, file.path(OUTPUT_DIR, "10_bd_rf_model.rds"))
saveRDS(sites, file.path(OUTPUT_DIR, "10_with_bd_filled.rds"))
saveRDS(cv_results, file.path(OUTPUT_DIR, "10_cv_results.rds"))
cat("  Saved: 10_bd_rf_model.rds\n")
cat("  Saved: 10_with_bd_filled.rds\n")
cat("  Saved: 10_cv_results.rds\n")

log_step("Bulk density gap-filling complete!")

