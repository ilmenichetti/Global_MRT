################################################################################
# 14b_predict_global_MRT_xgboost.R
#
# Generate global MRT prediction maps from fitted XGBoost models
# Creates 6 maps (one per model) at 0.1° resolution
#
# NOTE: Models are trained on log(MRT), predictions are back-transformed
#       using exp() to return MRT in years.
#
# Input:  ./Global_MRT_code/outputs/13b_xgb_models.rds
#         ./Global_MRT_code/spatialized_layers/*
# Output: ./Global_MRT_code/outputs/MRT_predictions_xgb/*.tif
#
# Author: Lorenzo
# Date: 2026-01-08
################################################################################

library(terra)
library(xgboost)
library(dplyr)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
RASTER_DIR   <- file.path(PIPELINE_DIR, "spatialized_layers")
PRED_DIR     <- file.path(OUTPUT_DIR, "MRT_predictions_xgb")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots")

MODEL_FILE  <- file.path(OUTPUT_DIR, "13b_xgb_models.rds")

# Create output directory
dir.create(PRED_DIR, recursive = TRUE, showWarnings = FALSE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  GLOBAL MRT PREDICTION - XGBOOST (Step 14b)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD MODELS
# =============================================================================

cat("Loading XGBoost models...\n")
models <- readRDS(MODEL_FILE)
cat("  Loaded", length(models), "models:\n")
for (m in names(models)) {
  cat("    -", models[[m]]$config$name, "(", length(models[[m]]$predictors), "predictors)\n")
}
cat("\n")

# =============================================================================
# RASTER MAPPING
# =============================================================================

raster_mapping <- list(
  # Climate variables (from ERA5)
  temperature_2m_mean = "climate/temperature_2m_mean_0.1deg.tif",
  temperature_seasonality = "climate/temperature_seasonality_0.1deg.tif",
  soil_temperature_0_20cm = "climate/soil_temperature_0_20cm_0.1deg.tif",
  soil_moisture_0_20cm = "climate/soil_moisture_0_20cm_0.1deg.tif",
  snow_cover_mean = "climate/snow_cover_mean_0.1deg.tif",
  potential_evaporation = "climate/potential_evaporation_0.1deg.tif",
  lai_high_veg = "climate/lai_high_veg_0.1deg.tif",
  lai_low_veg = "climate/lai_low_veg_0.1deg.tif",
  
  # Physical variables - soil (from SoilGrids)
  sg_clay = "soilgrids/clay_0_20cm_global.tif",
  sg_sand = "soilgrids/sand_0_20cm_global.tif",
  sg_bdod = "soilgrids/bdod_0_20cm_global.tif",
  sg_cfvo = "soilgrids/cfvo_0_20cm_global.tif",
  
  # Physical variables - topography
  terrain_elev_mean = "topography/elev_mean_0p1deg.tif",
  terrain_slope_mean = "topography/slope_mean_0p1deg.tif",
  terrain_northness = "topography/northness_0p1deg.tif",
  terrain_eastness = "topography/eastness_0p1deg.tif",
  terrain_ruggedness = "topography/elev_ruggedness_0p1deg.tif",
  
  # Physical variables - landcover fractions (ESA WorldCover)
  lc_trees = "landcover/landcover_fractions_0p1deg.tif",
  lc_grassland = "landcover/landcover_fractions_0p1deg.tif",
  lc_shrubs = "landcover/landcover_fractions_0p1deg.tif",
  lc_cropland = "landcover/landcover_fractions_0p1deg.tif",
  lc_bare = "landcover/landcover_fractions_0p1deg.tif",
  lc_wetland = "landcover/landcover_fractions_0p1deg.tif",
  
  # Physical variables - disturbance
  hansen_any_loss = "disturbances/Hansen_forest_loss/Hansen_lossyear_global.tif",
  burn_count_before_obs = "disturbances/MODIS_burned/Burned_2020_global.tif",
  
  # Chemical variables (from SoilGrids)
  sg_phh2o = "soilgrids/phh2o_0_20cm_global.tif",
  sg_cec = "soilgrids/cec_0_20cm_global.tif",
  sg_nitrogen = "soilgrids/nitrogen_0_20cm_global.tif",
  
  # Chemical variables - soil classification
  soilclass_func_code = "soilclass/functional_group_0p1deg.tif",
  
  # Biological variables - microbial
  fungal_proportion = "microbial/fungal_proportion_0.1deg.tif",
  AM_roots_colonized = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  EcM_roots_colonized = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  EcM_AM_root_ratio = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  AM_richness = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_richness = "microbial/mycorrhiza_spun_0.1deg.tif",
  AM_endemism = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_endemism = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_AM_richness_ratio = "microbial/mycorrhiza_spun_0.1deg.tif"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

load_predictor_raster <- function(var_name, template, raster_dir, mapping) {
  
  if (!var_name %in% names(mapping) || is.null(mapping[[var_name]])) {
    return(NULL)
  }
  
  raster_path <- file.path(raster_dir, mapping[[var_name]])
  
  if (!file.exists(raster_path)) {
    warning("File not found: ", raster_path)
    return(NULL)
  }
  
  r <- rast(raster_path)
  
  # Handle multi-layer rasters
  if (nlyr(r) > 1) {
    layer_idx <- grep(var_name, names(r), ignore.case = TRUE)
    if (length(layer_idx) > 0) {
      r <- r[[layer_idx[1]]]
    } else {
      r <- r[[1]]
    }
  }
  
  # Resample if needed
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    r <- resample(r, template, method = "bilinear")
  }
  
  names(r) <- var_name
  return(r)
}

build_predictor_stack <- function(predictors, template, raster_dir, mapping) {
  
  cat("  Loading", length(predictors), "predictor layers...\n")
  
  layers <- list()
  missing <- c()
  
  for (var in predictors) {
    r <- load_predictor_raster(var, template, raster_dir, mapping)
    if (!is.null(r)) {
      layers[[var]] <- r
    } else {
      missing <- c(missing, var)
    }
  }
  
  if (length(missing) > 0) {
    cat("  WARNING:", length(missing), "predictors missing:", paste(missing, collapse = ", "), "\n")
  }
  
  if (length(layers) == 0) {
    stop("No predictor layers loaded!")
  }
  
  stack <- rast(layers)
  cat("  Loaded", nlyr(stack), "layers\n")
  
  return(stack)
}

# XGBoost prediction function
predict_xgb_from_stack <- function(model, stack, predictors, chunk_size = 1000000) {
  
  cat("  Predicting (log-scale, will back-transform)...\n")
  
  n_cells <- ncell(stack)
  n_chunks <- ceiling(n_cells / chunk_size)
  
  out <- rast(stack[[1]])
  names(out) <- "MRT_predicted"
  
  cat("    Processing", n_chunks, "chunks: ")
  
  for (i in 1:n_chunks) {
    cat(i, "")
    
    start_cell <- (i - 1) * chunk_size + 1
    end_cell <- min(i * chunk_size, n_cells)
    cells <- start_cell:end_cell
    
    vals <- values(stack)[cells, , drop = FALSE]
    colnames(vals) <- names(stack)
    
    complete <- complete.cases(vals)
    
    if (sum(complete) > 0) {
      # Ensure column order matches model expectation
      pred_matrix <- vals[complete, predictors, drop = FALSE]
      dtest <- xgb.DMatrix(data = as.matrix(pred_matrix))
      
      # Predict log_MRT
      log_predictions <- predict(model, dtest)
      
      # Back-transform
      predictions <- exp(log_predictions)
      
      out_vals <- rep(NA, length(cells))
      out_vals[complete] <- predictions
      out[cells] <- out_vals
    }
  }
  cat("done\n")
  
  return(out)
}

# =============================================================================
# CREATE TEMPLATE
# =============================================================================

cat("Creating global template at 0.1° resolution...\n")
template <- rast(resolution = 0.1, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
crs(template) <- "EPSG:4326"
cat("  Template:", ncol(template), "x", nrow(template), "cells\n\n")

# =============================================================================
# GENERATE PREDICTIONS
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  GENERATING PREDICTIONS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

prediction_rasters <- list()

for (model_id in names(models)) {
  
  cat("Processing", models[[model_id]]$config$name, "...\n")
  
  predictors <- models[[model_id]]$predictors
  model <- models[[model_id]]$model
  
  # Build predictor stack
  stack <- build_predictor_stack(predictors, template, RASTER_DIR, raster_mapping)
  
  # Check we have all required predictors
  available_predictors <- names(stack)
  missing_in_stack <- setdiff(predictors, available_predictors)
  
  if (length(missing_in_stack) > 0) {
    cat("  WARNING: Missing predictors in stack:", paste(missing_in_stack, collapse = ", "), "\n")
    # Use only available predictors
    predictors <- intersect(predictors, available_predictors)
  }
  
  # Predict
  pred_rast <- predict_xgb_from_stack(model, stack, predictors)
  
  # Save individual raster
  out_name <- gsub("_", "_", model_id)
  out_file <- file.path(PRED_DIR, paste0("MRT_", out_name, ".tif"))
  
  writeRaster(pred_rast, out_file, overwrite = TRUE)
  cat("  ✓ Saved:", out_file, "\n\n")
  
  prediction_rasters[[model_id]] <- pred_rast
}

# =============================================================================
# CREATE COMBINED OUTPUTS
# =============================================================================

cat("Creating combined outputs...\n")

# Stack all predictions
all_preds <- rast(prediction_rasters)
names(all_preds) <- names(prediction_rasters)

writeRaster(all_preds, file.path(PRED_DIR, "MRT_all_models_xgb.tif"), overwrite = TRUE)
cat("✓ Saved: MRT_all_models_xgb.tif\n")

# Difference map (M6 - M1)
if ("M6_full" %in% names(prediction_rasters) && "M1_climate" %in% names(prediction_rasters)) {
  diff_rast <- prediction_rasters$M6_full - prediction_rasters$M1_climate
  names(diff_rast) <- "MRT_diff_M6_M1"
  writeRaster(diff_rast, file.path(PRED_DIR, "MRT_difference_M6_M1_xgb.tif"), overwrite = TRUE)
  cat("✓ Saved: MRT_difference_M6_M1_xgb.tif\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  STEP 14b (XGBoost Prediction) COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Outputs saved to:", PRED_DIR, "\n")
cat("\nTo visualize, update PRED_DIR in 15_visualize_MRT_maps.R to:\n")
cat("  PRED_DIR <- file.path(PIPELINE_DIR, 'outputs/MRT_predictions_xgb')\n")