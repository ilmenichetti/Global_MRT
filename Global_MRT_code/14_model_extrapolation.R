################################################################################
# 14_predict_global_MRT.R
#
# Generate global MRT prediction maps from fitted Random Forest models
# Creates maps for each model configuration at 0.1° resolution
#
# Uses 4-group variable structure:
#   CLIMATE   - Temperature, moisture controls
#   EDAPHIC   - Intrinsic soil properties (texture, chemistry, terrain)
#   LANDUSE   - Land cover fractions & disturbance history
#   BIOLOGICAL - Soil ecology (mycorrhizal, fungal)
#
# NOTE: Models are trained on log(MRT), predictions are back-transformed
#       using exp() to return MRT in years.
#
# Input:  ./Global_MRT_code/outputs/13_rf_models.rds
#         ./Global_MRT_code/outputs/13_var_groups.rds
#         ./Global_MRT_code/spatialized_layers/*
# Output: ./Global_MRT_code/outputs/MRT_predictions/*.tif
#
# Author: Lorenzo
# Date: 2026-01-12
################################################################################

library(terra)
library(dplyr)
library(ranger)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
RASTER_DIR   <- file.path(PIPELINE_DIR, "spatialized_layers")
PRED_DIR     <- file.path(OUTPUT_DIR, "MRT_predictions")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots")

# Create output directory
dir.create(PRED_DIR, recursive = TRUE, showWarnings = FALSE)

# Target resolution
TARGET_RES <- 0.1  # degrees

cat("═══════════════════════════════════════════════════════════════\n")
cat("  GLOBAL MRT PREDICTION (Step 14)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD FITTED MODELS AND VARIABLE GROUPS
# =============================================================================

cat("Loading fitted models...\n")
models <- readRDS(file.path(OUTPUT_DIR, "13_rf_models.rds"))
cat("  Loaded", length(models), "models:\n")
for (m in names(models)) {
  cat("    -", models[[m]]$config$name, "(", length(models[[m]]$predictors), "predictors)\n")
}
cat("\n")

# Load variable groups
cat("Loading variable groups...\n")
VAR_GROUPS <- readRDS(file.path(OUTPUT_DIR, "13_var_groups.rds"))
cat("  Climate:", length(VAR_GROUPS$CLIMATE), "vars\n")
cat("  Edaphic:", length(VAR_GROUPS$EDAPHIC), "vars\n")
cat("  LandUse:", length(VAR_GROUPS$LANDUSE), "vars\n")
cat("  Biological:", length(VAR_GROUPS$BIOLOGICAL), "vars\n\n")

# =============================================================================
# DEFINE RASTER LAYER MAPPING
# =============================================================================

raster_mapping <- list(
  # ─────────────────────────────────────────────────────────────────────────────
  # CLIMATE variables (from ERA5)
  # ─────────────────────────────────────────────────────────────────────────────
  temperature_seasonality = "climate/temperature_seasonality_0.1deg.tif",
  soil_temperature_0_20cm = "climate/soil_temperature_0_20cm_0.1deg.tif",
  soil_moisture_0_20cm = "climate/soil_moisture_0_20cm_0.1deg.tif",
  snow_cover_mean = "climate/snow_cover_mean_0.1deg.tif",
  potential_evaporation = "climate/potential_evaporation_0.1deg.tif",
  koppen_value = "climate/koppen_detailed_0.1deg.tif",
  
  # ─────────────────────────────────────────────────────────────────────────────
  # EDAPHIC variables - Intrinsic soil properties
  # ─────────────────────────────────────────────────────────────────────────────
  # Texture (SoilGrids)
  sg_clay = "soilgrids/clay_0_20cm_global.tif",
  sg_sand = "soilgrids/sand_0_20cm_global.tif",
  
  # Bulk density & structure (SoilGrids)
  sg_bdod = "soilgrids/bdod_0_20cm_global.tif",
  sg_cfvo = "soilgrids/cfvo_0_20cm_global.tif",
  
  # Chemistry (SoilGrids)
  sg_phh2o = "soilgrids/phh2o_0_20cm_global.tif",
  sg_cec = "soilgrids/cec_0_20cm_global.tif",
  sg_nitrogen = "soilgrids/nitrogen_0_20cm_global.tif",
  
  # Soil classification
  soilclass_func_code = "soilclass/functional_group_0p1deg.tif",
  
  # Topography
  terrain_elev_mean = "topography/elev_mean_0p1deg.tif",
  terrain_slope_mean = "topography/slope_mean_0p1deg.tif",
  terrain_northness = "topography/northness_0p1deg.tif",
  terrain_eastness = "topography/eastness_0p1deg.tif",
  terrain_ruggedness = "topography/elev_ruggedness_0p1deg.tif",
  
  # ─────────────────────────────────────────────────────────────────────────────
  # LANDUSE variables - Land cover & disturbance
  # ─────────────────────────────────────────────────────────────────────────────
  # Land cover fractions (ESA WorldCover)
  lc_trees = "landcover/landcover_fractions_0p1deg.tif",
  lc_grassland = "landcover/landcover_fractions_0p1deg.tif",
  lc_shrubs = "landcover/landcover_fractions_0p1deg.tif",
  lc_cropland = "landcover/landcover_fractions_0p1deg.tif",
  lc_bare = "landcover/landcover_fractions_0p1deg.tif",
  lc_wetland = "landcover/landcover_fractions_0p1deg.tif",
  
  # Disturbance history
  hansen_any_loss = "disturbances/Hansen_forest_loss/Hansen_lossyear_global.tif",
  burn_count_before_obs = "disturbances/MODIS_burned/Burned_2020_global.tif",
  
  # ─────────────────────────────────────────────────────────────────────────────
  # BIOLOGICAL variables - Soil ecology
  # ─────────────────────────────────────────────────────────────────────────────
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

# Load and align a raster to template
load_predictor_raster <- function(var_name, template, raster_dir, mapping) {
  
  if (!var_name %in% names(mapping) || is.null(mapping[[var_name]])) {
    return(NULL)
  }
  
  fpath <- file.path(raster_dir, mapping[[var_name]])
  
  if (!file.exists(fpath)) {
    warning("File not found: ", fpath)
    return(NULL)
  }
  
  r <- rast(fpath)
  
  # Handle multi-layer rasters (e.g., mycorrhiza files)
  if (nlyr(r) > 1) {
    layer_idx <- grep(var_name, names(r), ignore.case = TRUE)
    if (length(layer_idx) > 0) {
      r <- r[[layer_idx[1]]]
    } else {
      r <- r[[1]]
    }
  }
  
  # Resample to template if needed
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    r <- resample(r, template, method = "bilinear")
  }
  
  names(r) <- var_name
  return(r)
}

# Build predictor stack for a model
build_predictor_stack <- function(predictors, template, raster_dir, mapping) {
  
  cat("  Loading", length(predictors), "predictor layers...\n")
  
  rast_list <- list()
  missing <- c()
  
  for (var in predictors) {
    r <- load_predictor_raster(var, template, raster_dir, mapping)
    
    if (!is.null(r)) {
      rast_list[[var]] <- r
      cat("    ✓", var, "\n")
    } else {
      missing <- c(missing, var)
      cat("    ✗", var, "(missing)\n")
    }
  }
  
  if (length(missing) > 0) {
    cat("  WARNING:", length(missing), "predictors missing\n")
  }
  
  if (length(rast_list) == 0) {
    return(NULL)
  }
  
  stack <- rast(rast_list)
  return(stack)
}

# Predict from raster stack using ranger model
predict_from_stack <- function(model, stack, chunk_size = 1000000) {
  
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
    complete <- complete.cases(vals)
    
    if (sum(complete) > 0) {
      pred_df <- as.data.frame(vals[complete, , drop = FALSE])
      log_predictions <- predict(model, data = pred_df)$predictions
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
# CREATE GLOBAL TEMPLATE
# =============================================================================

cat("Creating global template at", TARGET_RES, "degree resolution...\n")

template <- rast(
  xmin = -180, xmax = 180,
  ymin = -90, ymax = 90,
  res = TARGET_RES,
  crs = "EPSG:4326"
)

cat("  Dimensions:", nrow(template), "x", ncol(template), "\n")
cat("  Cells:", format(ncell(template), big.mark = ","), "\n\n")

# =============================================================================
# GENERATE PREDICTIONS FOR EACH MODEL
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  GENERATING PREDICTIONS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

prediction_rasters <- list()

for (model_id in names(models)) {
  
  model_info <- models[[model_id]]
  cat("─── ", model_info$config$name, " ───\n\n")
  
  predictors <- model_info$predictors
  
  pred_stack <- build_predictor_stack(predictors, template, RASTER_DIR, raster_mapping)
  
  if (is.null(pred_stack)) {
    cat("  SKIPPED - no predictor layers available\n\n")
    next
  }
  
  available <- names(pred_stack)
  missing <- setdiff(predictors, available)
  
  if (length(missing) > 0) {
    cat("  WARNING: Missing predictors:", paste(missing, collapse = ", "), "\n")
  }
  
  tryCatch({
    mrt_pred <- predict_from_stack(model_info$model, pred_stack)
    
    out_file <- file.path(PRED_DIR, paste0("MRT_", model_id, ".tif"))
    writeRaster(mrt_pred, out_file, overwrite = TRUE)
    cat("  ✓ Saved:", out_file, "\n\n")
    
    prediction_rasters[[model_id]] <- mrt_pred
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n\n")
  })
}

# =============================================================================
# CREATE GROUP-SPECIFIC PREDICTION LAYERS
# =============================================================================
#
# For RGB composite visualization, we need to create prediction layers
# that show the contribution of each variable group beyond climate
#
# Strategy: Calculate MRT predictions from:
#   - Climate-only model (baseline)
#   - Climate + each group (to get group contribution)
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CALCULATING GROUP CONTRIBUTIONS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Check which models exist for group contribution calculation
model_pairs <- list(
  EDAPHIC = c("M1_climate", "M2_climate_edaphic"),
  LANDUSE = c("M1_climate", "M3_climate_landuse"),
  BIOLOGICAL = c("M1_climate", "M4_climate_biological")
)

group_contributions <- list()

for (group_name in names(model_pairs)) {
  
  base_model <- model_pairs[[group_name]][1]
  full_model <- model_pairs[[group_name]][2]
  
  cat("Calculating", group_name, "contribution...\n")
  
  if (base_model %in% names(prediction_rasters) && 
      full_model %in% names(prediction_rasters)) {
    
    # Contribution = |prediction with group - prediction without|
    # Using ratio to make it scale-independent
    contrib <- prediction_rasters[[full_model]] / prediction_rasters[[base_model]]
    
    # Convert to log-ratio for better visualization
    # Positive = group increases MRT, Negative = decreases
    contrib_log <- log(contrib)
    
    # Also calculate absolute difference
    contrib_diff <- prediction_rasters[[full_model]] - prediction_rasters[[base_model]]
    
    names(contrib_log) <- paste0(group_name, "_contribution_logratio")
    names(contrib_diff) <- paste0(group_name, "_contribution_diff")
    
    group_contributions[[paste0(group_name, "_logratio")]] <- contrib_log
    group_contributions[[paste0(group_name, "_diff")]] <- contrib_diff
    
    # Save
    writeRaster(contrib_log, 
                file.path(PRED_DIR, paste0("contribution_", group_name, "_logratio.tif")),
                overwrite = TRUE)
    writeRaster(contrib_diff,
                file.path(PRED_DIR, paste0("contribution_", group_name, "_diff.tif")),
                overwrite = TRUE)
    
    cat("  ✓ Saved contribution rasters for", group_name, "\n")
    
  } else {
    cat("  SKIPPED - required models not available\n")
  }
}

# =============================================================================
# CALCULATE CLIMATE CONTRIBUTION (R² equivalent)
# =============================================================================
#
# Climate contribution = R² of M1 (climate-only) model
# This will be shown as grayscale intensity
# =============================================================================

if ("M1_climate" %in% names(models)) {
  
  cat("\nClimate model R² (CV):", round(models$M1_climate$metrics$R2_cv, 3), "\n")
  cat("This will be used as the grayscale intensity baseline.\n\n")
  
  # Store this for visualization step
  climate_R2 <- models$M1_climate$metrics$R2_cv
  saveRDS(climate_R2, file.path(OUTPUT_DIR, "14_climate_R2.rds"))
}

# =============================================================================
# CREATE COMPARISON MAPS
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CREATING COMPARISON OUTPUTS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

if (length(prediction_rasters) > 1) {
  
  # Stack all predictions
  all_preds <- rast(prediction_rasters)
  names(all_preds) <- names(prediction_rasters)
  
  combined_file <- file.path(PRED_DIR, "MRT_all_models.tif")
  writeRaster(all_preds, combined_file, overwrite = TRUE)
  cat("✓ Saved combined stack:", combined_file, "\n")
  
  # ─────────────────────────────────────────────────────────────────────────────
  # Difference maps: Each model minus climate baseline (M1)
  # ─────────────────────────────────────────────────────────────────────────────
  
  if ("M1_climate" %in% names(prediction_rasters)) {
    
    baseline <- prediction_rasters$M1_climate
    
    # M2 - M1: Edaphic contribution
    if ("M2_climate_edaphic" %in% names(prediction_rasters)) {
      diff_rast <- prediction_rasters$M2_climate_edaphic - baseline
      names(diff_rast) <- "MRT_diff_edaphic"
      writeRaster(diff_rast, file.path(PRED_DIR, "MRT_difference_M2_M1_edaphic.tif"), overwrite = TRUE)
      cat("✓ Saved difference map (M2 - M1): Edaphic contribution\n")
    }
    
    # M3 - M1: LandUse contribution
    if ("M3_climate_landuse" %in% names(prediction_rasters)) {
      diff_rast <- prediction_rasters$M3_climate_landuse - baseline
      names(diff_rast) <- "MRT_diff_landuse"
      writeRaster(diff_rast, file.path(PRED_DIR, "MRT_difference_M3_M1_landuse.tif"), overwrite = TRUE)
      cat("✓ Saved difference map (M3 - M1): LandUse contribution\n")
    }
    
    # M4 - M1: Biological contribution
    if ("M4_climate_biological" %in% names(prediction_rasters)) {
      diff_rast <- prediction_rasters$M4_climate_biological - baseline
      names(diff_rast) <- "MRT_diff_biological"
      writeRaster(diff_rast, file.path(PRED_DIR, "MRT_difference_M4_M1_biological.tif"), overwrite = TRUE)
      cat("✓ Saved difference map (M4 - M1): Biological contribution\n")
    }
    
    # M5 - M1: Edaphic + LandUse combined
    if ("M5_climate_edaphic_landuse" %in% names(prediction_rasters)) {
      diff_rast <- prediction_rasters$M5_climate_edaphic_landuse - baseline
      names(diff_rast) <- "MRT_diff_edaphic_landuse"
      writeRaster(diff_rast, file.path(PRED_DIR, "MRT_difference_M5_M1_edaphic_landuse.tif"), overwrite = TRUE)
      cat("✓ Saved difference map (M5 - M1): Edaphic + LandUse contribution\n")
    }
    
    # M6 - M1: Edaphic + Biological combined
    if ("M6_climate_edaphic_bio" %in% names(prediction_rasters)) {
      diff_rast <- prediction_rasters$M6_climate_edaphic_bio - baseline
      names(diff_rast) <- "MRT_diff_edaphic_bio"
      writeRaster(diff_rast, file.path(PRED_DIR, "MRT_difference_M6_M1_edaphic_bio.tif"), overwrite = TRUE)
      cat("✓ Saved difference map (M6 - M1): Edaphic + Biological contribution\n")
    }
    
    # M7 - M1: Full model (all groups)
    if ("M7_full" %in% names(prediction_rasters)) {
      diff_rast <- prediction_rasters$M7_full - baseline
      names(diff_rast) <- "MRT_diff_full"
      writeRaster(diff_rast, file.path(PRED_DIR, "MRT_difference_M7_M1_full.tif"), overwrite = TRUE)
      cat("✓ Saved difference map (M7 - M1): Full model contribution\n")
    }
    
    # ─────────────────────────────────────────────────────────────────────────────
    # Stack all differences for easy comparison
    # ─────────────────────────────────────────────────────────────────────────────
    
    diff_files <- list.files(PRED_DIR, pattern = "MRT_difference_M[2-7]_M1", full.names = TRUE)
    if (length(diff_files) > 1) {
      diff_stack <- rast(diff_files)
      writeRaster(diff_stack, file.path(PRED_DIR, "MRT_all_differences_vs_climate.tif"), overwrite = TRUE)
      cat("✓ Saved stacked differences:", file.path(PRED_DIR, "MRT_all_differences_vs_climate.tif"), "\n")
    }
  }
  
  # CV across models
  if (length(prediction_rasters) >= 3) {
    mean_pred <- mean(all_preds, na.rm = TRUE)
    sd_pred <- stdev(all_preds, na.rm = TRUE)
    cv_pred <- sd_pred / mean_pred * 100
    names(cv_pred) <- "MRT_CV_across_models"
    
    cv_file <- file.path(PRED_DIR, "MRT_CV_across_models.tif")
    writeRaster(cv_pred, cv_file, overwrite = TRUE)
    cat("✓ Saved CV across models:", cv_file, "\n")
  }
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  PREDICTION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (model_id in names(prediction_rasters)) {
  r <- prediction_rasters[[model_id]]
  vals <- values(r, na.rm = TRUE)
  
  cat(models[[model_id]]$config$name, ":\n")
  cat("  Mean MRT:", round(mean(vals), 1), "years\n")
  cat("  Median MRT:", round(median(vals), 1), "years\n")
  cat("  Range:", round(min(vals), 1), "-", round(max(vals), 1), "years\n")
  cat("  Valid cells:", format(length(vals), big.mark = ","), 
      "(", round(100 * length(vals) / ncell(r), 1), "%)\n\n")
}

# =============================================================================
# QUICK VISUALIZATION
# =============================================================================

cat("Generating preview plots...\n")

# Plot full model prediction
if ("M7_full" %in% names(prediction_rasters)) {
  
  png(file.path(PLOT_DIR, "14_MRT_prediction_M7_full.png"), 
      width = 1200, height = 600, res = 100)
  
  plot(prediction_rasters$M7_full,
       main = "Predicted Mean Residence Time (Full Model)",
       col = hcl.colors(50, "YlOrRd", rev = TRUE),
       range = c(0, 200))
  
  dev.off()
  cat("✓ Saved:", file.path(PLOT_DIR, "14_MRT_prediction_M7_full.png"), "\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  STEP 14 COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Output files in:", PRED_DIR, "\n")
list.files(PRED_DIR, pattern = "\\.tif$")