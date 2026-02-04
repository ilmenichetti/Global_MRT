################################################################################
# 13_fit_MRT_models.R
#
# Fit Random Forest models to predict Mean Residence Time (MRT)
# Using nested model structure to decompose variance by predictor class
#
# Model Structure (4 mechanistic groups beyond climate baseline):
#   M1: Climate only                          - Baseline
#   M2: Climate + Edaphic                     - Intrinsic soil properties
#   M3: Climate + LandUse                     - Land cover & disturbance
#   M4: Climate + Biological                  - Soil ecology
#   M5: Climate + Edaphic + LandUse           - Abiotic combined
#   M6: Climate + Edaphic + Biological        - Soil-focused
#   M7: Full (Climate + Edaphic + LandUse + Biological)
#
# Variable Groups:
#   CLIMATE   - Temperature, moisture controls (baseline decomposition drivers)
#   EDAPHIC   - Intrinsic soil properties (texture, chemistry, terrain)
#   LANDUSE   - Land cover fractions & disturbance history
#   BIOLOGICAL - Soil ecology (mycorrhizal associations, fungal proportion)
#
# Cross-validation: Spatial block CV (2° grid cells)
# Variable importance: Permutation-based
#
# Input:  ./Global_MRT_code/outputs/12b_model_ready.rds
# Output: ./Global_MRT_code/outputs/13_rf_models.rds
#         ./Global_MRT_code/outputs/13_model_comparison.csv
#         ./Global_MRT_code/outputs/13_variable_importance.csv
#         ./Global_MRT_code/outputs/13_model_config.rds
#
# Author: Lorenzo
# Date: 2026-01-12
################################################################################

library(dplyr)
library(tidyr)
library(ranger)
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

# Random Forest parameters
RF_NTREES <- 500
RF_MTRY   <- NULL  # Will be set to sqrt(p) by default
RF_MIN_NODE_SIZE <- 5
RF_NUM_THREADS <- parallel::detectCores() - 1  # Leave one core free

# Cross-validation parameters
CV_FOLDS <- 10

# Set seed for reproducibility
set.seed(42)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  RANDOM FOREST MODEL FITTING (Step 13)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Using", RF_NUM_THREADS, "threads for parallel processing\n\n")

# =============================================================================
# DEFINE PREDICTOR VARIABLE GROUPS (4-GROUP STRUCTURE)
# =============================================================================
#
# Conceptual framework:
#
# 1. CLIMATE: Temperature and moisture controls on decomposition rates
#    - Direct effects on microbial activity and enzyme kinetics
#    - Baseline model - these are the primary drivers ESMs already capture
#    - NOTE: LAI excluded as it correlates with productivity used in MRT calc
#
# 2. EDAPHIC: Intrinsic soil properties (slow-changing)
#    - Physical: Texture, bulk density, coarse fragments
#    - Chemical: pH, CEC, nitrogen, soil classification
#    - Topographic: Elevation, slope, aspect (affect microclimate/drainage)
#    - These represent "what the soil IS" - mineral matrix properties
#
# 3. LANDUSE: Land cover and disturbance regime (extrinsic/anthropogenic)
#    - Current land cover composition
#    - Historical disturbance (fire, forest loss)
#    - These represent "what happens TO the soil" - perturbation regime
#
# 4. BIOLOGICAL: Soil ecological community
#    - Fungal:bacterial ratio
#    - Mycorrhizal associations (AM vs EcM colonization, richness)
#    - These represent the decomposer community composition
# =============================================================================

CLIMATE_VARS <- c(
  "temperature_seasonality",   # Temperature seasonality (SD) - ERA5
  "soil_temperature_0_20cm",   # Soil temperature 0-20cm (°C) - ERA5
  "soil_moisture_0_20cm",      # Soil moisture 0-20cm (m³/m³) - ERA5
  "snow_cover_mean",           # Mean snow cover fraction - ERA5
  "potential_evaporation",      # Potential evaporation (m) - ERA5
  "koppen_value"               # Köppen-Geiger climate zone 
  )

EDAPHIC_VARS <- c(
  # Texture - from SoilGrids
  "sg_clay",              # Clay content (%) - SoilGrids
  "sg_sand",              # Sand content (%) - SoilGrids
  
  # Bulk density & structure - from SoilGrids
  "sg_bdod",              # Bulk density (cg/cm³) - SoilGrids
  "sg_cfvo",              # Coarse fragments volume (%) - SoilGrids
  
  # Chemistry - from SoilGrids
  "sg_phh2o",             # Soil pH - SoilGrids
  "sg_cec",               # Cation exchange capacity (cmol/kg) - SoilGrids
  "sg_nitrogen",          # Total nitrogen - SoilGrids
  
  # Soil classification
  "soilclass_func_code",  # WRB functional group code (1-6)
  
  # Topography (affects microclimate, drainage, erosion)
  "terrain_elev_mean",    # Elevation (m)
  "terrain_slope_mean",   # Slope (degrees)
  "terrain_northness",    # Aspect northness (-1 to 1)
  "terrain_eastness",     # Aspect eastness (-1 to 1)
  "terrain_ruggedness"    # Terrain ruggedness index
)

LANDUSE_VARS <- c(
  # Landcover fractions - ESA WorldCover 2021
  "lc_trees",             # Tree cover fraction
  "lc_grassland",         # Grassland fraction
  "lc_shrubs",            # Shrubland fraction
  "lc_cropland",          # Cropland fraction
  "lc_bare",              # Bare/sparse fraction
  "lc_wetland",           # Wetland fraction
  
  # Disturbance history
  "hansen_any_loss",      # Binary: any forest loss before observation
  "burn_count_before_obs" # Number of fires before observation
)

BIOLOGICAL_VARS <- c(
  "fungal_proportion",      # F:B ratio (Yu et al. 2022)
  "AM_roots_colonized",     # AM root colonization (Barceló et al. 2023)
  "EcM_roots_colonized",    # EcM root colonization (Barceló et al. 2023)
  "EcM_AM_root_ratio",      # EcM:AM ratio (derived)
  "AM_richness",            # AM fungal richness (SPUN 2025)
  "EcM_richness",           # EcM fungal richness (SPUN 2025)
  "AM_endemism",            # AM rarity-weighted richness (SPUN 2025)
  "EcM_endemism",           # EcM rarity-weighted richness (SPUN 2025)
  "EcM_AM_richness_ratio"   # EcM:AM richness ratio (derived)
)

# Store for export (used in step 14 and 15)
VAR_GROUPS <- list(
  CLIMATE = CLIMATE_VARS,
  EDAPHIC = EDAPHIC_VARS,
  LANDUSE = LANDUSE_VARS,
  BIOLOGICAL = BIOLOGICAL_VARS
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

# =============================================================================
# LOG-TRANSFORM RESPONSE VARIABLE
# =============================================================================

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

cat("Log-transforming MRT response variable...\n")
model_data$log_MRT <- log(model_data$MRT_years)

cat("  MRT_years summary:\n")
print(summary(model_data$MRT_years))
cat("\n  log_MRT summary:\n")
print(summary(model_data$log_MRT))
cat("\n")

# Check which variables are actually available
check_vars <- function(var_list, data) {
  available <- var_list[var_list %in% names(data)]
  missing <- var_list[!var_list %in% names(data)]
  if (length(missing) > 0) {
    cat("  Missing variables:", paste(missing, collapse = ", "), "\n")
  }
  return(available)
}

cat("Checking variable availability...\n")
cat("Climate:\n")
CLIMATE_VARS <- check_vars(CLIMATE_VARS, model_data)
cat("  Available:", length(CLIMATE_VARS), "\n")

cat("Edaphic:\n")
EDAPHIC_VARS <- check_vars(EDAPHIC_VARS, model_data)
cat("  Available:", length(EDAPHIC_VARS), "\n")

cat("LandUse:\n")
LANDUSE_VARS <- check_vars(LANDUSE_VARS, model_data)
cat("  Available:", length(LANDUSE_VARS), "\n")

cat("Biological:\n")
BIOLOGICAL_VARS <- check_vars(BIOLOGICAL_VARS, model_data)
cat("  Available:", length(BIOLOGICAL_VARS), "\n\n")

# Update VAR_GROUPS with available variables
VAR_GROUPS <- list(
  CLIMATE = CLIMATE_VARS,
  EDAPHIC = EDAPHIC_VARS,
  LANDUSE = LANDUSE_VARS,
  BIOLOGICAL = BIOLOGICAL_VARS
)

# =============================================================================
# DEFINE MODEL CONFIGURATIONS
# =============================================================================

model_configs <- list(
  M1_climate = list(
    name = "M1: Climate",
    vars = CLIMATE_VARS,
    description = "Climate variables only - baseline model"
  ),
  M2_climate_edaphic = list(
    name = "M2: Climate + Edaphic",
    vars = c(CLIMATE_VARS, EDAPHIC_VARS),
    description = "Climate + intrinsic soil properties"
  ),
  M3_climate_landuse = list(
    name = "M3: Climate + LandUse",
    vars = c(CLIMATE_VARS, LANDUSE_VARS),
    description = "Climate + land cover & disturbance"
  ),
  M4_climate_biological = list(
    name = "M4: Climate + Biological",
    vars = c(CLIMATE_VARS, BIOLOGICAL_VARS),
    description = "Climate + soil ecology"
  ),
  M5_climate_edaphic_landuse = list(
    name = "M5: Climate + Edaphic + LandUse",
    vars = c(CLIMATE_VARS, EDAPHIC_VARS, LANDUSE_VARS),
    description = "All abiotic variables"
  ),
  M6_climate_edaphic_bio = list(
    name = "M6: Climate + Edaphic + Bio",
    vars = c(CLIMATE_VARS, EDAPHIC_VARS, BIOLOGICAL_VARS),
    description = "Soil-focused (no land use)"
  ),
  M7_full = list(
    name = "M7: Full Model",
    vars = c(CLIMATE_VARS, EDAPHIC_VARS, LANDUSE_VARS, BIOLOGICAL_VARS),
    description = "All predictor variables"
  )
)

cat("Model configurations:\n")
for (m in names(model_configs)) {
  cat("  ", model_configs[[m]]$name, ":", length(model_configs[[m]]$vars), "predictors\n")
}
cat("\n")

# =============================================================================
# CROSS-VALIDATION STRATEGY
# =============================================================================

CV_METHOD <- "spatial_grid"
BLOCK_SIZE <- 2  # Degrees (2° ≈ 200km at equator)

cat("Creating cross-validation folds...\n")
cat("  Method:", CV_METHOD, "\n")

if (CV_METHOD == "spatial_grid") {
  # Create grid cell IDs
  model_data <- model_data %>%
    mutate(
      grid_x = floor(longitude_decimal_degrees / BLOCK_SIZE),
      grid_y = floor(latitude_decimal_degrees / BLOCK_SIZE),
      grid_id = paste(grid_x, grid_y, sep = "_")
    )
  
  # Get unique grid cells
  unique_grids <- unique(model_data$grid_id)
  n_grids <- length(unique_grids)
  cat("  Block size:", BLOCK_SIZE, "°\n")
  cat("  Unique grid cells:", n_grids, "\n")
  
  # Randomly assign grid cells to folds
  set.seed(42)
  grid_folds <- data.frame(
    grid_id = unique_grids,
    cv_fold = sample(1:CV_FOLDS, n_grids, replace = TRUE)
  )
  
  # Merge back to data
  model_data <- model_data %>%
    left_join(grid_folds, by = "grid_id")
  
  cat("  Assigned to", CV_FOLDS, "folds\n")
}

# Report fold distribution
cat("\n  Fold sizes:\n")
fold_table <- table(model_data$cv_fold, useNA = "ifany")
print(fold_table)
cat("\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to fit RF and compute CV metrics
fit_rf_spatial_cv <- function(data, predictors, response = "log_MRT", 
                              fold_col = "cv_fold") {
  
  # Remove rows with missing predictors or fold assignment
  complete_idx <- complete.cases(data[, c(response, predictors, fold_col)])
  data_complete <- data[complete_idx, ]
  
  n_complete <- nrow(data_complete)
  cat("    Complete cases:", n_complete, "/", nrow(data), "\n")
  
  if (n_complete < 100) {
    warning("Very few complete cases!")
    return(NULL)
  }
  
  # Get actual folds present in complete data
  folds <- sort(unique(data_complete[[fold_col]]))
  n_folds <- length(folds)
  cat("    Folds:", n_folds, "\n")
  
  # Storage for predictions
  predictions <- rep(NA, n_complete)
  observed <- data_complete[[response]]
  
  # Spatial CV loop
  cat("    CV folds: ")
  for (fold in folds) {
    cat(fold, "")
    test_idx <- data_complete[[fold_col]] == fold
    train_idx <- !test_idx
    
    # Fit model on training data
    train_data <- data_complete[train_idx, ]
    test_data <- data_complete[test_idx, ]
    
    # Fit Random Forest
    rf_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    
    rf_model <- ranger(
      formula = rf_formula,
      data = train_data,
      num.trees = RF_NTREES,
      min.node.size = RF_MIN_NODE_SIZE,
      num.threads = RF_NUM_THREADS,
      importance = "none",
      seed = 42
    )
    
    # Predict on test fold
    predictions[test_idx] <- predict(rf_model, data = test_data)$predictions
  }
  cat("done\n")
  
  # Compute CV metrics
  residuals <- observed - predictions
  
  SS_res <- sum(residuals^2, na.rm = TRUE)
  SS_tot <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  
  R2_cv <- 1 - SS_res / SS_tot
  RMSE_cv <- sqrt(mean(residuals^2, na.rm = TRUE))
  MAE_cv <- mean(abs(residuals), na.rm = TRUE)
  
  # Fit final model on all data (for variable importance and prediction)
  rf_final <- ranger(
    formula = rf_formula,
    data = data_complete,
    num.trees = RF_NTREES,
    min.node.size = RF_MIN_NODE_SIZE,
    num.threads = RF_NUM_THREADS,
    importance = "permutation",
    seed = 42
  )
  
  return(list(
    model = rf_final,
    metrics = data.frame(
      n_obs = n_complete,
      n_predictors = length(predictors),
      R2_cv = R2_cv,
      RMSE_cv = RMSE_cv,
      MAE_cv = MAE_cv,
      R2_oob = rf_final$r.squared,
      RMSE_oob = sqrt(rf_final$prediction.error)
    ),
    predictions_cv = data.frame(
      observed = observed,
      predicted = predictions,
      fold = data_complete[[fold_col]]
    ),
    importance = rf_final$variable.importance
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
  
  # Get available predictors for this model
  predictors <- config$vars[config$vars %in% names(model_data)]
  cat("  Predictors:", length(predictors), "\n")
  
  if (length(predictors) < 1) {
    cat("  SKIPPED - no predictors available\n\n")
    next
  }
  
  # Fit model with spatial CV
  result <- fit_rf_spatial_cv(
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
    cat("  RMSE (CV):", round(result$metrics$RMSE_cv, 3), "(log-scale)\n")
    cat("  R² (OOB):", round(result$metrics$R2_oob, 3), "\n")
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

# Create comparison table
comparison_df <- do.call(rbind, lapply(names(results), function(m) {
  data.frame(
    model_id = m,
    model_name = results[[m]]$config$name,
    n_predictors = results[[m]]$metrics$n_predictors,
    n_obs = results[[m]]$metrics$n_obs,
    R2_cv = results[[m]]$metrics$R2_cv,
    RMSE_cv = results[[m]]$metrics$RMSE_cv,
    MAE_cv = results[[m]]$metrics$MAE_cv,
    R2_oob = results[[m]]$metrics$R2_oob
  )
}))

# Calculate incremental R² compared to climate-only baseline
if ("M1_climate" %in% names(results)) {
  baseline_R2 <- results$M1_climate$metrics$R2_cv
  comparison_df$delta_R2 <- comparison_df$R2_cv - baseline_R2
} else {
  comparison_df$delta_R2 <- NA
}

print(comparison_df %>% 
        select(model_name, n_predictors, n_obs, R2_cv, RMSE_cv, delta_R2) %>%
        mutate(across(where(is.numeric), ~round(., 3))),
      row.names = FALSE)

write.csv(comparison_df, file.path(OUTPUT_DIR, "13_model_comparison.csv"),
          row.names = FALSE)

# =============================================================================
# VARIABLE IMPORTANCE
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  VARIABLE IMPORTANCE (Full Model)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

if ("M7_full" %in% names(results)) {
  
  importance_df <- data.frame(
    variable = names(results$M7_full$importance),
    importance = results$M7_full$importance
  ) %>%
    mutate(
      var_group = case_when(
        variable %in% CLIMATE_VARS ~ "Climate",
        variable %in% EDAPHIC_VARS ~ "Edaphic",
        variable %in% LANDUSE_VARS ~ "LandUse",
        variable %in% BIOLOGICAL_VARS ~ "Biological",
        TRUE ~ "Other"
      )
    ) %>%
    arrange(desc(importance))
  
  cat("Top 20 predictors (permutation importance):\n")
  print(head(importance_df, 20), row.names = FALSE)
  
  # Summary by group
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
  
  # Save
  write.csv(importance_df, file.path(OUTPUT_DIR, "13_variable_importance.csv"),
            row.names = FALSE)
}

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  SAVING OUTPUTS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# 1. Save all models and results
models_output <- file.path(OUTPUT_DIR, "13_rf_models.rds")
saveRDS(results, models_output)
cat("✓ Models saved:", models_output, "\n")

# 2. Save model comparison table
comparison_output <- file.path(OUTPUT_DIR, "13_model_comparison.csv")
write.csv(comparison_df, comparison_output, row.names = FALSE)
cat("✓ Model comparison saved:", comparison_output, "\n")

# 3. Save variable importance (for full model)
if ("M7_full" %in% names(results)) {
  importance_output <- file.path(OUTPUT_DIR, "13_variable_importance.csv")
  write.csv(importance_df, importance_output, row.names = FALSE)
  cat("✓ Variable importance saved:", importance_output, "\n")
}

# 4. Save variable groups configuration (for prediction and visualization)
var_groups_output <- file.path(OUTPUT_DIR, "13_var_groups.rds")
saveRDS(VAR_GROUPS, var_groups_output)
cat("✓ Variable groups saved:", var_groups_output, "\n")

# 5. Save model configuration for prediction step
model_config_for_prediction <- lapply(names(results), function(m) {
  list(
    model_id = m,
    model_name = results[[m]]$config$name,
    predictors = results[[m]]$predictors,
    raster_mapping = data.frame(
      predictor = results[[m]]$predictors,
      raster_file = NA
    )
  )
})
names(model_config_for_prediction) <- names(results)

config_output <- file.path(OUTPUT_DIR, "13_model_config.rds")
saveRDS(model_config_for_prediction, config_output)
cat("✓ Model config saved:", config_output, "\n")

# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================

cat("\nGenerating diagnostic plots...\n")

# Observed vs Predicted for full model
if ("M7_full" %in% names(results)) {
  
  cv_preds <- results$M7_full$predictions_cv
  
  # Back-transform from log scale for plotting
  cv_preds$observed_MRT <- exp(cv_preds$observed)
  cv_preds$predicted_MRT <- exp(cv_preds$predicted)
  
  p1 <- ggplot(cv_preds, aes(x = observed_MRT, y = predicted_MRT)) +
    geom_hex(bins = 50) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    scale_fill_viridis_c(trans = "log10") +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "M7 Full Model: Observed vs Predicted MRT",
      subtitle = paste0("R² (log-scale) = ", round(results$M7_full$metrics$R2_cv, 3),
                        ", RMSE (log) = ", round(results$M7_full$metrics$RMSE_cv, 3)),
      x = "Observed MRT (years)",
      y = "Predicted MRT (years)"
    ) +
    theme_bw()
  
  ggsave(file.path(PLOT_DIR, "13_obs_vs_pred_full.png"), p1, 
         width = 8, height = 7, dpi = 150)
  cat("✓ Saved:", file.path(PLOT_DIR, "13_obs_vs_pred_full.png"), "\n")
  
  # Model comparison barplot
  p2 <- ggplot(comparison_df, aes(x = reorder(model_name, R2_cv), y = R2_cv)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = round(R2_cv, 3)), hjust = -0.1, size = 3) +
    coord_flip() +
    ylim(0, max(comparison_df$R2_cv) * 1.15) +
    labs(
      title = "Model Comparison: Cross-validated R² (log-scale)",
      x = "",
      y = "R² (spatial CV)"
    ) +
    theme_bw()
  
  ggsave(file.path(PLOT_DIR, "13_model_comparison.png"), p2,
         width = 8, height = 5, dpi = 150)
  cat("✓ Saved:", file.path(PLOT_DIR, "13_model_comparison.png"), "\n")
  
  # Variable importance plot with new groups
  group_colors <- c(
    "Climate" = "#666666",
    "Edaphic" = "#E41A1C",
    "LandUse" = "#4DAF4A",
    "Biological" = "#377EB8",
    "Other" = "#999999"
  )
  
  p3 <- ggplot(head(importance_df, 20), 
               aes(x = reorder(variable, importance), y = importance, fill = var_group)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = group_colors) +
    labs(
      title = "Variable Importance (Full Model)",
      subtitle = "Permutation importance - top 20 predictors",
      x = "",
      y = "Importance",
      fill = "Variable Group"
    ) +
    theme_bw()
  
  ggsave(file.path(PLOT_DIR, "13_variable_importance.png"), p3,
         width = 10, height = 8, dpi = 150)
  cat("✓ Saved:", file.path(PLOT_DIR, "13_variable_importance.png"), "\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  STEP 13 COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Next step: Run 14_predict_global_MRT.R to generate global maps\n")
cat("Required: Spatialized predictor layers in ./Global_MRT_code/spatialized_layers/\n")




