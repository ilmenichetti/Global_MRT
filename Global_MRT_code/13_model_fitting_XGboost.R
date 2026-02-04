################################################################################
# 13b_fit_MRT_models_xgboost.R
#
# ALTERNATIVE to 13_fit_MRT_models.R - Uses XGBoost instead of Random Forest
# Often outperforms RF by 5-15% due to better gradient handling
#
# Fit XGBoost models to predict Mean Residence Time (MRT)
# Using nested model structure to decompose variance by predictor class
#
# Model Structure:
#   M1: Climate only                          - Baseline
#   M2: Climate + Physical                    - Physical protection
#   M3: Climate + Chemical                    - Chemical stabilization
#   M4: Climate + Biological                  - Biological regulation
#   M5: Climate + Physical + Chemical         - Abiotic combined
#   M6: Full (Climate + Physical + Chemical + Biological)
#
# Cross-validation: 2° spatial blocks
# Variable importance: Gain-based
#
# References:
#   Chen T, Guestrin C (2016) XGBoost: A Scalable Tree Boosting System.
#     Proceedings of the 22nd ACM SIGKDD. DOI: 10.1145/2939672.2939785
#
# Input:  ./Global_MRT_code/outputs/12b_model_ready.rds
# Output: ./Global_MRT_code/outputs/13b_xgb_models.rds
#         ./Global_MRT_code/outputs/13b_model_comparison.csv
#         ./Global_MRT_code/outputs/13b_variable_importance.csv
#
# Author: Lorenzo
# Date: 2026-01-08
################################################################################

library(dplyr)
library(tidyr)
library(xgboost)
library(ggplot2)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots")

INPUT_FILE  <- file.path(OUTPUT_DIR, "12b_model_ready.rds")

# Create directories
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# XGBoost parameters (can be tuned)
XGB_PARAMS <- list(
  objective = "reg:squarederror",
  eta = 0.1,               # Learning rate (lower = more trees needed, but more precise)
  max_depth = 6,           # Tree depth (6-10 typical)
  min_child_weight = 10,   # Minimum samples per leaf
  subsample = 0.8,         # Row sampling
  colsample_bytree = 0.8,  # Column sampling per tree
  gamma = 0                # Minimum loss reduction for split
)

XGB_NROUNDS <- 500         # Maximum boosting rounds
XGB_EARLY_STOP <- 50       # Early stopping patience
XGB_NTHREAD <- parallel::detectCores() - 1

# Cross-validation parameters
CV_METHOD <- "spatial_grid"
CV_FOLDS <- 10
BLOCK_SIZE <- 2  # degrees

# Set seed for reproducibility
set.seed(42)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  XGBOOST MODEL FITTING (Step 13b - Alternative)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Using", XGB_NTHREAD, "threads\n\n")

# =============================================================================
# DEFINE PREDICTOR VARIABLE GROUPS
# =============================================================================

CLIMATE_VARS <- c(
  "temperature_seasonality",
  "soil_temperature_0_20cm",
  "soil_moisture_0_20cm",
  "snow_cover_mean",
  "potential_evaporation",
  "lai_high_veg",
  "lai_low_veg"
)

PHYSICAL_VARS <- c(
  "sg_clay", "sg_sand", "sg_bdod", "sg_cfvo",
  "terrain_elev_mean", "terrain_slope_mean",
  "terrain_northness", "terrain_eastness", "terrain_ruggedness",
  "lc_trees", "lc_grassland", "lc_shrubs", "lc_cropland", "lc_bare", "lc_wetland",
  "hansen_any_loss", "burn_count_before_obs"
)

CHEMICAL_VARS <- c(
  "sg_phh2o", "sg_cec", "sg_nitrogen", "soilclass_func_code"
)

BIOLOGICAL_VARS <- c(
  "fungal_proportion",
  "AM_roots_colonized", "EcM_roots_colonized", "EcM_AM_root_ratio",
  "AM_richness", "EcM_richness", "AM_endemism", "EcM_endemism",
  "EcM_AM_richness_ratio"
)

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

cat("Loading data...\n")
soil_data <- readRDS(INPUT_FILE)
cat("  Loaded", format(nrow(soil_data), big.mark = ","), "observations\n\n")

# Filter to valid MRT values
cat("Filtering to valid MRT observations...\n")
model_data <- soil_data %>%
  filter(MRT_QC == "valid") %>%
  filter(!is.na(MRT_years)) %>%
  filter(MRT_years > 0 & MRT_years < Inf)

cat("  Valid observations:", format(nrow(model_data), big.mark = ","), "\n\n")

# Remove outliers (extreme 1% tails)
cat("Removing MRT outliers (1st and 99th percentile)...\n")
q01 <- quantile(model_data$MRT_years, 0.01, na.rm = TRUE)
q99 <- quantile(model_data$MRT_years, 0.99, na.rm = TRUE)
n_before <- nrow(model_data)
model_data <- model_data %>% filter(MRT_years > q01 & MRT_years < q99)
n_removed <- n_before - nrow(model_data)
cat("  Removed", n_removed, "outliers (", round(100 * n_removed / n_before, 1), "%)\n")
cat("  Kept MRT range:", round(q01, 1), "-", round(q99, 1), "years\n")
cat("  Remaining observations:", format(nrow(model_data), big.mark = ","), "\n\n")

# Log-transform response
cat("Log-transforming MRT response variable...\n")
model_data$log_MRT <- log(model_data$MRT_years)
cat("  log_MRT range:", round(min(model_data$log_MRT), 2), "to", 
    round(max(model_data$log_MRT), 2), "\n\n")

# =============================================================================
# CHECK VARIABLE AVAILABILITY
# =============================================================================

check_vars <- function(var_list, data) {
  available <- var_list[var_list %in% names(data)]
  missing <- var_list[!var_list %in% names(data)]
  if (length(missing) > 0) {
    cat("  Missing:", paste(missing, collapse = ", "), "\n")
  }
  return(available)
}

cat("Checking variable availability...\n")
CLIMATE_VARS <- check_vars(CLIMATE_VARS, model_data)
cat("  Climate:", length(CLIMATE_VARS), "\n")
PHYSICAL_VARS <- check_vars(PHYSICAL_VARS, model_data)
cat("  Physical:", length(PHYSICAL_VARS), "\n")
CHEMICAL_VARS <- check_vars(CHEMICAL_VARS, model_data)
cat("  Chemical:", length(CHEMICAL_VARS), "\n")
BIOLOGICAL_VARS <- check_vars(BIOLOGICAL_VARS, model_data)
cat("  Biological:", length(BIOLOGICAL_VARS), "\n\n")

# =============================================================================
# DEFINE MODEL CONFIGURATIONS
# =============================================================================

model_configs <- list(
  M1_climate = list(
    name = "M1: Climate",
    vars = CLIMATE_VARS
  ),
  M2_climate_physical = list(
    name = "M2: Climate + Physical",
    vars = c(CLIMATE_VARS, PHYSICAL_VARS)
  ),
  M3_climate_chemical = list(
    name = "M3: Climate + Chemical",
    vars = c(CLIMATE_VARS, CHEMICAL_VARS)
  ),
  M4_climate_biological = list(
    name = "M4: Climate + Biological",
    vars = c(CLIMATE_VARS, BIOLOGICAL_VARS)
  ),
  M5_climate_phys_chem = list(
    name = "M5: Climate + Physical + Chemical",
    vars = c(CLIMATE_VARS, PHYSICAL_VARS, CHEMICAL_VARS)
  ),
  M6_full = list(
    name = "M6: Full Model",
    vars = c(CLIMATE_VARS, PHYSICAL_VARS, CHEMICAL_VARS, BIOLOGICAL_VARS)
  )
)

cat("Model configurations:\n")
for (m in names(model_configs)) {
  cat("  ", model_configs[[m]]$name, ":", length(model_configs[[m]]$vars), "predictors\n")
}
cat("\n")

# =============================================================================
# CREATE SPATIAL BLOCK CV FOLDS
# =============================================================================

cat("Creating cross-validation folds...\n")
cat("  Method:", CV_METHOD, "\n")
cat("  Block size:", BLOCK_SIZE, "°\n")

model_data <- model_data %>%
  mutate(
    grid_x = floor(longitude_decimal_degrees / BLOCK_SIZE),
    grid_y = floor(latitude_decimal_degrees / BLOCK_SIZE),
    grid_id = paste(grid_x, grid_y, sep = "_")
  )

unique_grids <- unique(model_data$grid_id)
n_grids <- length(unique_grids)
cat("  Unique grid cells:", n_grids, "\n")

set.seed(42)
grid_folds <- data.frame(
  grid_id = unique_grids,
  cv_fold = sample(1:CV_FOLDS, n_grids, replace = TRUE)
)

model_data <- model_data %>%
  left_join(grid_folds, by = "grid_id")

cat("  Assigned to", CV_FOLDS, "folds\n\n")

# =============================================================================
# HELPER FUNCTION: Fit XGBoost with CV
# =============================================================================

fit_xgb_spatial_cv <- function(data, predictors, response = "log_MRT", 
                               fold_col = "cv_fold", params = XGB_PARAMS) {
  
  # Remove rows with missing predictors
  complete_idx <- complete.cases(data[, c(response, predictors, fold_col)])
  data_complete <- data[complete_idx, ]
  
  n_complete <- nrow(data_complete)
  cat("    Complete cases:", n_complete, "/", nrow(data), "\n")
  
  if (n_complete < 100) {
    warning("Very few complete cases!")
    return(NULL)
  }
  
  # Prepare data
  X <- as.matrix(data_complete[, predictors])
  y <- data_complete[[response]]
  folds <- data_complete[[fold_col]]
  
  unique_folds <- sort(unique(folds))
  n_folds <- length(unique_folds)
  cat("    Folds:", n_folds, "\n")
  
  # Storage
  predictions <- rep(NA, n_complete)
  best_nrounds <- c()
  
  # CV loop
  cat("    CV folds: ")
  for (fold in unique_folds) {
    cat(fold, "")
    
    test_idx <- which(folds == fold)
    train_idx <- which(folds != fold)
    
    # Create DMatrix
    dtrain <- xgb.DMatrix(data = X[train_idx, , drop = FALSE], label = y[train_idx])
    dtest <- xgb.DMatrix(data = X[test_idx, , drop = FALSE], label = y[test_idx])
    
    # Train with early stopping
    watchlist <- list(train = dtrain, eval = dtest)
    
    model <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = XGB_NROUNDS,
      watchlist = watchlist,
      early_stopping_rounds = XGB_EARLY_STOP,
      verbose = 0,
      nthread = XGB_NTHREAD
    )
    
    best_nrounds <- c(best_nrounds, model$best_iteration)
    
    # Predict
    predictions[test_idx] <- predict(model, dtest)
  }
  cat("done\n")
  
  # Compute CV metrics
  observed <- y
  residuals <- observed - predictions
  
  SS_res <- sum(residuals^2, na.rm = TRUE)
  SS_tot <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  
  R2_cv <- 1 - SS_res / SS_tot
  RMSE_cv <- sqrt(mean(residuals^2, na.rm = TRUE))
  MAE_cv <- mean(abs(residuals), na.rm = TRUE)
  
  cat("    Avg best rounds:", round(mean(best_nrounds)), "\n")
  
  # Fit final model on all data
  dtrain_full <- xgb.DMatrix(data = X, label = y)
  
  final_model <- xgb.train(
    params = params,
    data = dtrain_full,
    nrounds = round(mean(best_nrounds) * 1.1),  # Slightly more rounds for full data
    verbose = 0,
    nthread = XGB_NTHREAD
  )
  
  # Get importance
  importance <- xgb.importance(feature_names = predictors, model = final_model)
  
  return(list(
    model = final_model,
    metrics = data.frame(
      n_obs = n_complete,
      n_predictors = length(predictors),
      R2_cv = R2_cv,
      RMSE_cv = RMSE_cv,
      MAE_cv = MAE_cv,
      avg_nrounds = mean(best_nrounds)
    ),
    predictions_cv = data.frame(
      observed = observed,
      predicted = predictions,
      fold = folds
    ),
    importance = importance
  ))
}

# =============================================================================
# FIT ALL MODELS
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  FITTING MODELS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

results <- list()
total_start <- Sys.time()

for (model_id in names(model_configs)) {
  
  config <- model_configs[[model_id]]
  cat("Fitting", config$name, "...\n")
  model_start <- Sys.time()
  
  predictors <- config$vars[config$vars %in% names(model_data)]
  cat("  Predictors:", length(predictors), "\n")
  
  if (length(predictors) < 1) {
    cat("  SKIPPED - no predictors available\n\n")
    next
  }
  
  result <- fit_xgb_spatial_cv(
    data = model_data,
    predictors = predictors,
    response = "log_MRT",
    fold_col = "cv_fold"
  )
  
  if (!is.null(result)) {
    result$config <- config
    result$predictors <- predictors
    results[[model_id]] <- result
    
    elapsed <- round(difftime(Sys.time(), model_start, units = "mins"), 1)
    cat("  R² (CV):", round(result$metrics$R2_cv, 3), "\n")
    cat("  RMSE (CV):", round(result$metrics$RMSE_cv, 3), "(log scale)\n")
    cat("  Time:", elapsed, "minutes\n\n")
  }
}

total_elapsed <- round(difftime(Sys.time(), total_start, units = "mins"), 1)
cat("Total fitting time:", total_elapsed, "minutes\n\n")

# =============================================================================
# COMPILE RESULTS
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  MODEL COMPARISON\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

comparison_df <- do.call(rbind, lapply(names(results), function(m) {
  data.frame(
    model_id = m,
    model_name = results[[m]]$config$name,
    n_predictors = results[[m]]$metrics$n_predictors,
    n_obs = results[[m]]$metrics$n_obs,
    R2_cv = results[[m]]$metrics$R2_cv,
    RMSE_cv = results[[m]]$metrics$RMSE_cv,
    MAE_cv = results[[m]]$metrics$MAE_cv,
    avg_nrounds = results[[m]]$metrics$avg_nrounds
  )
}))

# Calculate delta R² vs baseline
if ("M1_climate" %in% names(results)) {
  baseline_R2 <- results$M1_climate$metrics$R2_cv
  comparison_df$delta_R2 <- comparison_df$R2_cv - baseline_R2
} else {
  comparison_df$delta_R2 <- NA
}

print(comparison_df %>% 
        select(model_name, n_predictors, R2_cv, RMSE_cv, delta_R2) %>%
        mutate(across(where(is.numeric), ~round(., 3))),
      row.names = FALSE)

# =============================================================================
# VARIABLE IMPORTANCE
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  VARIABLE IMPORTANCE (Full Model - Gain)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

if ("M6_full" %in% names(results)) {
  
  importance_df <- results$M6_full$importance %>%
    mutate(
      var_group = case_when(
        Feature %in% CLIMATE_VARS ~ "Climate",
        Feature %in% PHYSICAL_VARS ~ "Physical",
        Feature %in% CHEMICAL_VARS ~ "Chemical",
        Feature %in% BIOLOGICAL_VARS ~ "Biological",
        TRUE ~ "Other"
      )
    ) %>%
    rename(variable = Feature, importance = Gain) %>%
    arrange(desc(importance))
  
  cat("Top 20 predictors (Gain-based importance):\n")
  print(head(importance_df %>% select(variable, importance, var_group), 20), row.names = FALSE)
  
  cat("\nImportance by variable group:\n")
  group_importance <- importance_df %>%
    group_by(var_group) %>%
    summarise(
      n_vars = n(),
      total_importance = sum(importance),
      mean_importance = mean(importance),
      .groups = "drop"
    ) %>%
    arrange(desc(total_importance))
  
  print(group_importance, row.names = FALSE)
}

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  SAVING OUTPUTS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Save models
models_output <- file.path(OUTPUT_DIR, "13b_xgb_models.rds")
saveRDS(results, models_output)
cat("✓ Models saved:", models_output, "\n")

# Save comparison
comparison_output <- file.path(OUTPUT_DIR, "13b_model_comparison.csv")
write.csv(comparison_df, comparison_output, row.names = FALSE)
cat("✓ Model comparison saved:", comparison_output, "\n")

# Save importance
if ("M6_full" %in% names(results)) {
  importance_output <- file.path(OUTPUT_DIR, "13b_variable_importance.csv")
  write.csv(importance_df, importance_output, row.names = FALSE)
  cat("✓ Variable importance saved:", importance_output, "\n")
}

# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================

cat("\nGenerating diagnostic plots...\n")

if ("M6_full" %in% names(results)) {
  
  cv_preds <- results$M6_full$predictions_cv
  cv_preds$observed_MRT <- exp(cv_preds$observed)
  cv_preds$predicted_MRT <- exp(cv_preds$predicted)
  
  p1 <- ggplot(cv_preds, aes(x = observed_MRT, y = predicted_MRT)) +
    geom_hex(bins = 50) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    scale_fill_viridis_c(trans = "log10") +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "XGBoost M6 Full Model: Observed vs Predicted MRT",
      subtitle = paste0("R² (log-scale) = ", round(results$M6_full$metrics$R2_cv, 3),
                        ", RMSE (log) = ", round(results$M6_full$metrics$RMSE_cv, 3)),
      x = "Observed MRT (years)",
      y = "Predicted MRT (years)"
    ) +
    theme_bw()
  
  ggsave(file.path(PLOT_DIR, "13b_xgb_obs_vs_pred_full.png"), p1, 
         width = 8, height = 7, dpi = 150)
  cat("✓ Saved:", file.path(PLOT_DIR, "13b_xgb_obs_vs_pred_full.png"), "\n")
  
  # Model comparison
  p2 <- ggplot(comparison_df, aes(x = reorder(model_name, R2_cv), y = R2_cv)) +
    geom_col(fill = "darkgreen") +
    geom_text(aes(label = round(R2_cv, 3)), hjust = -0.1, size = 3) +
    coord_flip() +
    ylim(0, max(comparison_df$R2_cv) * 1.15) +
    labs(
      title = "XGBoost Model Comparison: Cross-validated R² (log-scale)",
      x = "",
      y = "R² (spatial CV)"
    ) +
    theme_bw()
  
  ggsave(file.path(PLOT_DIR, "13b_xgb_model_comparison.png"), p2,
         width = 8, height = 5, dpi = 150)
  cat("✓ Saved:", file.path(PLOT_DIR, "13b_xgb_model_comparison.png"), "\n")
  
  # Variable importance
  p3 <- ggplot(head(importance_df, 20), 
               aes(x = reorder(variable, importance), y = importance, fill = var_group)) +
    geom_col() +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "XGBoost Variable Importance (Full Model)",
      subtitle = "Gain-based importance - top 20 predictors",
      x = "",
      y = "Importance (Gain)",
      fill = "Variable Group"
    ) +
    theme_bw()
  
  ggsave(file.path(PLOT_DIR, "13b_xgb_variable_importance.png"), p3,
         width = 10, height = 8, dpi = 150)
  cat("✓ Saved:", file.path(PLOT_DIR, "13b_xgb_variable_importance.png"), "\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  STEP 13b (XGBoost) COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("To use these models for prediction, run:\n")
cat("  source('14b_predict_global_MRT_xgboost.R')\n")
cat("\nOr compare with RF results in 13_model_comparison.csv\n")

# =============================================================================
# METHODS TEXT (for manuscript)
# =============================================================================
#
# Gradient boosted tree models (XGBoost; Chen & Guestrin 2016) were fitted to
# predict log-transformed MRT. XGBoost was chosen as an alternative to Random
# Forest due to its ability to capture complex interactions and its typically
# superior performance on regression tasks. Models were trained with the
# following hyperparameters: learning rate (eta) = 0.1, maximum tree depth = 6,
# subsample ratio = 0.8, column sample ratio = 0.8. Early stopping was used
# with a patience of 50 rounds to prevent overfitting.
#
# Six nested model configurations were compared to decompose variance explained
# by different predictor classes. Model performance was evaluated using spatial
# block cross-validation with 2° × 2° grid cells randomly assigned to k=10 folds.
#
# Variable importance was assessed using the gain metric, which measures the
# improvement in accuracy contributed by each feature across all splits where
# it is used.
# =============================================================================