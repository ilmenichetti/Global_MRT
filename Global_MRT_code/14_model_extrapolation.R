################################################################################
# 14_predict_global_MRT.R
#
# Generate global MRT prediction maps from fitted Random Forest models.
# Creates maps for each model configuration at 0.1 degree resolution.
#
# 2026-04 FIX: SoilGrids rasters on disk are in native integer encoding
# (pH x10, bulk density x100, etc.). This script now applies the
# d_factor division at raster load, matching the physical units used at
# training time (see script 12b). Rasters on disk are left untouched.
#
# Uses 4-group variable structure:
#   CLIMATE    - Temperature, moisture controls
#   EDAPHIC    - Intrinsic soil properties (texture, chemistry, terrain)
#   LANDUSE    - Land cover fractions and disturbance history
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
# Date:   2026-01-12  (fix 2026-04)
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

dir.create(PRED_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_RES <- 0.1  # degrees

# --- SoilGrids v2.0 d_factor scale factors ----------------------------------
# See https://www.isric.org/explore/soilgrids/faq-soilgrids
# Variables not listed here are read at native units (scale = 1).
SG_SCALE <- c(
  sg_clay     = 10,    # g/kg     -> %
  sg_sand     = 10,    # g/kg     -> %
  sg_cfvo     = 10,    # cm3/dm3  -> %
  sg_cec      = 10,    # mmol(c)/kg -> cmol(c)/kg
  sg_phh2o    = 10,    # pH*10    -> pH
  sg_bdod     = 100,   # cg/cm3   -> g/cm3
  sg_nitrogen = 100    # cg/kg    -> g/kg
)

# Plausible physical-unit maxima for a post-scaling sanity check
SG_EXPECTED_MAX <- c(
  sg_clay = 100, sg_sand = 100, sg_cfvo = 100, sg_cec = 250,
  sg_phh2o = 12, sg_bdod = 3,   sg_nitrogen = 500
)

cat("\u2550\u2550\u2550 GLOBAL MRT PREDICTION (Step 14) \u2550\u2550\u2550\n\n")

# =============================================================================
# LOAD FITTED MODELS AND VARIABLE GROUPS
# =============================================================================

cat("Loading fitted models...\n")
models <- readRDS(file.path(OUTPUT_DIR, "13_rf_models.rds"))
cat("  Loaded", length(models), "models:\n")
for (m in names(models)) {
  cat("    -", models[[m]]$config$name,
      "(", length(models[[m]]$predictors), "predictors)\n")
}
cat("\n")

cat("Loading variable groups...\n")
VAR_GROUPS <- readRDS(file.path(OUTPUT_DIR, "13_var_groups.rds"))
cat("  Climate:",    length(VAR_GROUPS$CLIMATE),    "vars\n")
cat("  Edaphic:",    length(VAR_GROUPS$EDAPHIC),    "vars\n")
cat("  LandUse:",    length(VAR_GROUPS$LANDUSE),    "vars\n")
cat("  Biological:", length(VAR_GROUPS$BIOLOGICAL), "vars\n\n")

# =============================================================================
# DEFINE RASTER LAYER MAPPING
# =============================================================================

raster_mapping <- list(
  # CLIMATE (ERA5)
  temperature_seasonality = "climate/temperature_seasonality_0.1deg.tif",
  soil_temperature_0_20cm = "climate/soil_temperature_0_20cm_0.1deg.tif",
  soil_moisture_0_20cm    = "climate/soil_moisture_0_20cm_0.1deg.tif",
  snow_cover_mean         = "climate/snow_cover_mean_0.1deg.tif",
  potential_evaporation   = "climate/potential_evaporation_0.1deg.tif",
  koppen_value            = "climate/koppen_detailed_0.1deg.tif",
  
  # EDAPHIC -- SoilGrids (native integer encoding; scale factor applied at load)
  sg_clay     = "soilgrids/clay_0_20cm_global.tif",
  sg_sand     = "soilgrids/sand_0_20cm_global.tif",
  sg_bdod     = "soilgrids/bdod_0_20cm_global.tif",
  sg_cfvo     = "soilgrids/cfvo_0_20cm_global.tif",
  sg_phh2o    = "soilgrids/phh2o_0_20cm_global.tif",
  sg_cec      = "soilgrids/cec_0_20cm_global.tif",
  sg_nitrogen = "soilgrids/nitrogen_0_20cm_global.tif",
  
  # EDAPHIC -- Soil classification and topography
  soilclass_func_code = "soilclass/functional_group_0p1deg.tif",
  terrain_elev_mean   = "topography/elev_mean_0p1deg.tif",
  terrain_slope_mean  = "topography/slope_mean_0p1deg.tif",
  terrain_northness   = "topography/northness_0p1deg.tif",
  terrain_eastness    = "topography/eastness_0p1deg.tif",
  terrain_ruggedness  = "topography/elev_ruggedness_0p1deg.tif",
  
  # LANDUSE
  lc_trees     = "landcover/landcover_fractions_0p1deg.tif",
  lc_grassland = "landcover/landcover_fractions_0p1deg.tif",
  lc_shrubs    = "landcover/landcover_fractions_0p1deg.tif",
  lc_cropland  = "landcover/landcover_fractions_0p1deg.tif",
  lc_bare      = "landcover/landcover_fractions_0p1deg.tif",
  lc_wetland   = "landcover/landcover_fractions_0p1deg.tif",
  
  hansen_any_loss       = "disturbances/Hansen_forest_loss/Hansen_lossyear_global.tif",
  burn_count_before_obs = "disturbances/MODIS_burned/Burned_2020_global.tif",
  
  # BIOLOGICAL
  fungal_proportion     = "microbial/fungal_proportion_0.1deg.tif",
  AM_roots_colonized    = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  EcM_roots_colonized   = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  EcM_AM_root_ratio     = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  AM_richness           = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_richness          = "microbial/mycorrhiza_spun_0.1deg.tif",
  AM_endemism           = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_endemism          = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_AM_richness_ratio = "microbial/mycorrhiza_spun_0.1deg.tif"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load and align a raster to template, applying SoilGrids d_factor if needed.
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
  
  # Multi-layer raster: select layer whose name matches var_name
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
  
  # --- SoilGrids d_factor scaling --------------------------------------------
  # Apply after resampling so any terra-internal type promotions have settled.
  if (var_name %in% names(SG_SCALE)) {
    r <- r / SG_SCALE[[var_name]]
    
    # Sanity check: warn if post-scaling values look suspicious
    if (var_name %in% names(SG_EXPECTED_MAX)) {
      mx <- as.numeric(global(r, "max", na.rm = TRUE))
      if (is.finite(mx) && mx > SG_EXPECTED_MAX[[var_name]]) {
        warning(sprintf(
          "%s post-scaling max (%.2f) exceeds expected (%.2f); ",
          var_name, mx, SG_EXPECTED_MAX[[var_name]]),
          "check that raster is in native SoilGrids integer encoding.")
      }
    }
  }
  
  names(r) <- var_name
  return(r)
}

#' Build predictor stack for a model.
build_predictor_stack <- function(predictors, template, raster_dir, mapping) {
  
  cat("  Loading", length(predictors), "predictor layers...\n")
  
  rast_list <- list()
  missing   <- c()
  
  for (var in predictors) {
    r <- load_predictor_raster(var, template, raster_dir, mapping)
    
    if (!is.null(r)) {
      rast_list[[var]] <- r
      scale_note <- if (var %in% names(SG_SCALE))
        sprintf(" (/%g)", SG_SCALE[[var]]) else ""
      cat("    \u2713", var, scale_note, "\n")
    } else {
      missing <- c(missing, var)
      cat("    \u2717", var, "(missing)\n")
    }
  }
  
  if (length(missing) > 0) {
    cat("  WARNING:", length(missing), "predictors missing\n")
  }
  
  if (length(rast_list) == 0) return(NULL)
  
  rast(rast_list)
}

#' Predict from raster stack using ranger model.
predict_from_stack <- function(model, stack, chunk_size = 1e6) {
  
  cat("  Predicting (log-scale, will back-transform)...\n")
  
  n_cells  <- ncell(stack)
  n_chunks <- ceiling(n_cells / chunk_size)
  
  out <- rast(stack[[1]])
  names(out) <- "MRT_predicted"
  
  cat("    Processing", n_chunks, "chunks: ")
  
  for (i in 1:n_chunks) {
    cat(i, "")
    
    start_cell <- (i - 1) * chunk_size + 1
    end_cell   <- min(i * chunk_size, n_cells)
    cells      <- start_cell:end_cell
    
    vals     <- values(stack)[cells, , drop = FALSE]
    complete <- complete.cases(vals)
    
    if (sum(complete) > 0) {
      pred_df         <- as.data.frame(vals[complete, , drop = FALSE])
      log_predictions <- predict(model, data = pred_df)$predictions
      predictions     <- exp(log_predictions)
      
      out_vals           <- rep(NA, length(cells))
      out_vals[complete] <- predictions
      out[cells]         <- out_vals
    }
  }
  cat("done\n")
  
  out
}

# =============================================================================
# CREATE GLOBAL TEMPLATE
# =============================================================================

cat("Creating global template at", TARGET_RES, "degree resolution...\n")

template <- rast(
  xmin = -180, xmax = 180,
  ymin = -90,  ymax = 90,
  res  = TARGET_RES,
  crs  = "EPSG:4326"
)

cat("  Dimensions:", nrow(template), "x", ncol(template), "\n")
cat("  Cells:",      format(ncell(template), big.mark = ","), "\n\n")

# =============================================================================
# GENERATE PREDICTIONS FOR EACH MODEL
# =============================================================================

cat("\u2550\u2550\u2550 GENERATING PREDICTIONS \u2550\u2550\u2550\n\n")

prediction_rasters <- list()

for (model_id in names(models)) {
  
  model_info <- models[[model_id]]
  cat("--- ", model_info$config$name, " ---\n\n")
  
  predictors <- model_info$predictors
  pred_stack <- build_predictor_stack(predictors, template,
                                      RASTER_DIR, raster_mapping)
  
  if (is.null(pred_stack)) {
    cat("  SKIPPED - no predictor layers available\n\n")
    next
  }
  
  available <- names(pred_stack)
  missing   <- setdiff(predictors, available)
  if (length(missing) > 0) {
    cat("  WARNING: Missing predictors:",
        paste(missing, collapse = ", "), "\n")
  }
  
  tryCatch({
    mrt_pred <- predict_from_stack(model_info$model, pred_stack)
    
    out_file <- file.path(PRED_DIR, paste0("MRT_", model_id, ".tif"))
    writeRaster(mrt_pred, out_file, overwrite = TRUE)
    cat("  \u2713 Saved:", out_file, "\n\n")
    
    prediction_rasters[[model_id]] <- mrt_pred
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n\n")
  })
}

# =============================================================================
# GROUP CONTRIBUTION RASTERS
# =============================================================================

cat("\u2550\u2550\u2550 CALCULATING GROUP CONTRIBUTIONS \u2550\u2550\u2550\n\n")

model_pairs <- list(
  EDAPHIC    = c("M1_climate", "M2_climate_edaphic"),
  LANDUSE    = c("M1_climate", "M3_climate_landuse"),
  BIOLOGICAL = c("M1_climate", "M4_climate_biological")
)

group_contributions <- list()

for (group_name in names(model_pairs)) {
  
  base_model <- model_pairs[[group_name]][1]
  full_model <- model_pairs[[group_name]][2]
  
  cat("Calculating", group_name, "contribution...\n")
  
  if (base_model %in% names(prediction_rasters) &&
      full_model %in% names(prediction_rasters)) {
    
    contrib      <- prediction_rasters[[full_model]] /
      prediction_rasters[[base_model]]
    contrib_log  <- log(contrib)
    contrib_diff <- prediction_rasters[[full_model]] -
      prediction_rasters[[base_model]]
    
    names(contrib_log)  <- paste0(group_name, "_contribution_logratio")
    names(contrib_diff) <- paste0(group_name, "_contribution_diff")
    
    group_contributions[[paste0(group_name, "_logratio")]] <- contrib_log
    group_contributions[[paste0(group_name, "_diff")]]     <- contrib_diff
    
    writeRaster(contrib_log,
                file.path(PRED_DIR,
                          paste0("contribution_", group_name, "_logratio.tif")),
                overwrite = TRUE)
    writeRaster(contrib_diff,
                file.path(PRED_DIR,
                          paste0("contribution_", group_name, "_diff.tif")),
                overwrite = TRUE)
    
    cat("  \u2713 Saved contribution rasters for", group_name, "\n")
    
  } else {
    cat("  SKIPPED - required models not available\n")
  }
}

# =============================================================================
# CLIMATE R2 REFERENCE
# =============================================================================

if ("M1_climate" %in% names(models)) {
  cat("\nClimate model R2 (CV):",
      round(models$M1_climate$metrics$R2_cv, 3), "\n")
  climate_R2 <- models$M1_climate$metrics$R2_cv
  saveRDS(climate_R2, file.path(OUTPUT_DIR, "14_climate_R2.rds"))
}

# =============================================================================
# COMPARISON MAPS (M2-M7 minus M1 climate baseline)
# =============================================================================

cat("\u2550\u2550\u2550 CREATING COMPARISON OUTPUTS \u2550\u2550\u2550\n\n")

if (length(prediction_rasters) > 1) {
  
  all_preds <- rast(prediction_rasters)
  names(all_preds) <- names(prediction_rasters)
  writeRaster(all_preds, file.path(PRED_DIR, "MRT_all_models.tif"),
              overwrite = TRUE)
  cat("\u2713 Saved combined stack\n")
  
  if ("M1_climate" %in% names(prediction_rasters)) {
    
    baseline <- prediction_rasters$M1_climate
    
    diff_specs <- list(
      c("M2_climate_edaphic",           "edaphic"),
      c("M3_climate_landuse",           "landuse"),
      c("M4_climate_biological",        "biological"),
      c("M5_climate_edaphic_landuse",   "edaphic_landuse"),
      c("M6_climate_edaphic_bio",       "edaphic_bio"),
      c("M7_full",                      "full")
    )
    
    for (spec in diff_specs) {
      model_id <- spec[1]
      label    <- spec[2]
      if (model_id %in% names(prediction_rasters)) {
        diff_rast <- prediction_rasters[[model_id]] - baseline
        names(diff_rast) <- paste0("MRT_diff_", label)
        out_file <- file.path(PRED_DIR,
                              sprintf("MRT_difference_%s_M1_%s.tif",
                                      sub("_.*$", "", model_id), label))
        writeRaster(diff_rast, out_file, overwrite = TRUE)
        cat("\u2713", basename(out_file), "\n")
      }
    }
    
    diff_files <- list.files(PRED_DIR,
                             pattern = "MRT_difference_M[2-7]_M1",
                             full.names = TRUE)
    if (length(diff_files) > 1) {
      diff_stack <- rast(diff_files)
      writeRaster(diff_stack,
                  file.path(PRED_DIR,
                            "MRT_all_differences_vs_climate.tif"),
                  overwrite = TRUE)
      cat("\u2713 Saved stacked differences\n")
    }
  }
  
  if (length(prediction_rasters) >= 3) {
    mean_pred <- mean(all_preds,   na.rm = TRUE)
    sd_pred   <- stdev(all_preds,  na.rm = TRUE)
    cv_pred   <- sd_pred / mean_pred * 100
    names(cv_pred) <- "MRT_CV_across_models"
    writeRaster(cv_pred,
                file.path(PRED_DIR, "MRT_CV_across_models.tif"),
                overwrite = TRUE)
    cat("\u2713 Saved CV across models\n")
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\u2550\u2550\u2550 PREDICTION SUMMARY \u2550\u2550\u2550\n\n")

for (model_id in names(prediction_rasters)) {
  r <- prediction_rasters[[model_id]]
  vals <- values(r, na.rm = TRUE)
  cat(models[[model_id]]$config$name, ":\n")
  cat("  Mean MRT:",   round(mean(vals),   1), "years\n")
  cat("  Median MRT:", round(median(vals), 1), "years\n")
  cat("  Range:",      round(min(vals),    1), "-",
      round(max(vals),    1), "years\n\n")
}

# Quick preview
if ("M7_full" %in% names(prediction_rasters)) {
  png(file.path(PLOT_DIR, "14_MRT_prediction_M7_full.png"),
      width = 1200, height = 600, res = 100)
  plot(prediction_rasters$M7_full,
       main = "Predicted Mean Residence Time (Full Model)",
       col = hcl.colors(50, "YlOrRd", rev = TRUE),
       range = c(0, 200))
  dev.off()
  cat("\u2713 Preview plot saved\n")
}

cat("\n\u2550\u2550\u2550 STEP 14 COMPLETE \u2550\u2550\u2550\n\n")