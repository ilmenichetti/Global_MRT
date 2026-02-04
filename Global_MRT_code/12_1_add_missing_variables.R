################################################################################
# 12b_add_missing_predictors.R
#
# THIS IS A FIX!!!!! 
# Not all predictors were filled in the table, so they are added  from spatialized GeoTIFF layers
# Run this AFTER step 12 (MRT calculation) and BEFORE step 13 (model fitting)
#
# Extracts from:
#   - climate/ folder (ERA5 variables)
#   - topography/ folder (if missing)
#   - soilgrids/ folder (if missing)
#   - landcover/ folder (fractions)
#   - microbial/ folder (if missing)
#   - soilclass/ folder (functional groups)
#
# Input:  ./Global_MRT_code/outputs/12_with_MRT.rds
# Output: ./Global_MRT_code/outputs/12b_model_ready.rds
#
################################################################################

library(terra)
library(dplyr)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
RASTER_DIR   <- file.path(PIPELINE_DIR, "spatialized_layers")

INPUT_FILE  <- file.path(OUTPUT_DIR, "12_with_MRT.rds")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "12b_model_ready.rds")

cat("═══════════════════════════════════════════════════════════════\n")
cat("  ADD MISSING PREDICTORS (Step 12b)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
soil_data <- readRDS(INPUT_FILE)
n_start <- ncol(soil_data)
cat("  Loaded", format(nrow(soil_data), big.mark = ","), "observations\n")
cat("  Starting variables:", n_start, "\n\n")

# Create spatial points (once, reuse for all extractions)
cat("Creating spatial points...\n")
pts <- vect(soil_data,
            geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
            crs = "EPSG:4326")
cat("  Done\n\n")

# =============================================================================
# DEFINE VARIABLE TO RASTER MAPPING
# =============================================================================
#
# Format: variable_name = "subfolder/filename.tif"
# For multi-band rasters, specify band name after |
# e.g., "landcover/fractions.tif|lc_trees"
# =============================================================================

var_raster_map <- list(
  
  # ─────────────────────────────────────────────────────────────────────────
  # CLIMATE (ERA5)
  # ─────────────────────────────────────────────────────────────────────────
  temperature_2m_mean = "climate/temperature_2m_mean_0.1deg.tif",
  temperature_seasonality = "climate/temperature_seasonality_0.1deg.tif",
  soil_temperature_0_20cm = "climate/soil_temperature_0_20cm_0.1deg.tif",
  soil_moisture_0_20cm = "climate/soil_moisture_0_20cm_0.1deg.tif",
  snow_cover_mean = "climate/snow_cover_mean_0.1deg.tif",
  potential_evaporation = "climate/potential_evaporation_0.1deg.tif",
  lai_high_veg = "climate/lai_high_veg_0.1deg.tif",
  lai_low_veg = "climate/lai_low_veg_0.1deg.tif",
  
  # ─────────────────────────────────────────────────────────────────────────
  # TOPOGRAPHY
  # ─────────────────────────────────────────────────────────────────────────
  terrain_elev_mean = "topography/elev_mean_0p1deg.tif",
  terrain_slope_mean = "topography/slope_mean_0p1deg.tif",
  terrain_northness = "topography/northness_0p1deg.tif",
  terrain_eastness = "topography/eastness_0p1deg.tif",
  terrain_ruggedness = "topography/elev_ruggedness_0p1deg.tif",
  
  # ─────────────────────────────────────────────────────────────────────────
  # SOILGRIDS
  # ─────────────────────────────────────────────────────────────────────────
  sg_clay = "soilgrids/clay_0_20cm_global.tif",
  sg_sand = "soilgrids/sand_0_20cm_global.tif",
  sg_bdod = "soilgrids/bdod_0_20cm_global.tif",
  sg_cfvo = "soilgrids/cfvo_0_20cm_global.tif",
  sg_cec = "soilgrids/cec_0_20cm_global.tif",
  sg_phh2o = "soilgrids/phh2o_0_20cm_global.tif",
  sg_nitrogen = "soilgrids/nitrogen_0_20cm_global.tif",
  
  # ─────────────────────────────────────────────────────────────────────────
  # SOIL CLASSIFICATION
  # ─────────────────────────────────────────────────────────────────────────
  soilclass_func_code = "soilclass/functional_group_0p1deg.tif",
  
  # ─────────────────────────────────────────────────────────────────────────
  # MICROBIAL
  # ─────────────────────────────────────────────────────────────────────────
  fungal_proportion = "microbial/fungal_proportion_0.1deg.tif"
  
  # Note: AM/EcM variables are in multi-band rasters, handled separately below
)

# =============================================================================
# HELPER FUNCTION: Extract single-band raster
# =============================================================================

extract_from_raster <- function(var_name, raster_path, points, data) {
  
  # Skip if already present
  if (var_name %in% names(data)) {
    n_valid <- sum(!is.na(data[[var_name]]))
    if (n_valid > 0) {
      cat(sprintf("  %-30s: already present (%d valid)\n", var_name, n_valid))
      return(data)
    }
  }
  
  full_path <- file.path(RASTER_DIR, raster_path)
  
  if (!file.exists(full_path)) {
    cat(sprintf("  %-30s: FILE NOT FOUND\n", var_name))
    return(data)
  }
  
  cat(sprintf("  %-30s: extracting...", var_name))
  
  tryCatch({
    r <- rast(full_path)
    
    # If multi-layer, use first layer
    if (nlyr(r) > 1) {
      r <- r[[1]]
    }
    
    vals <- terra::extract(r, points, ID = FALSE)[[1]]
    data[[var_name]] <- vals
    cat(sprintf(" done (%d valid)\n", sum(!is.na(vals))))
    
  }, error = function(e) {
    cat(sprintf(" ERROR: %s\n", e$message))
  })
  
  return(data)
}

# =============================================================================
# EXTRACT SINGLE-BAND RASTERS
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  EXTRACTING SINGLE-BAND RASTERS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (var_name in names(var_raster_map)) {
  soil_data <- extract_from_raster(var_name, var_raster_map[[var_name]], pts, soil_data)
}

# =============================================================================
# EXTRACT MULTI-BAND RASTERS: LANDCOVER FRACTIONS
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  EXTRACTING LANDCOVER FRACTIONS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

lc_fractions_file <- file.path(RASTER_DIR, "landcover/landcover_fractions_0p1deg.tif")

if (file.exists(lc_fractions_file)) {
  
  cat("Loading landcover fractions raster...\n")
  lc_rast <- rast(lc_fractions_file)
  cat("  Bands:", nlyr(lc_rast), "\n")
  cat("  Band names:", paste(names(lc_rast), collapse = ", "), "\n\n")
  
  # Extract all bands at once
  cat("Extracting all bands...\n")
  lc_vals <- terra::extract(lc_rast, pts, ID = FALSE)
  
  # Map band names to expected variable names
  # Adjust this mapping based on your actual band names
  lc_mapping <- list(
    lc_trees = c("trees", "lc_trees", "Tree"),
    lc_grassland = c("grassland", "lc_grassland", "Grassland"),
    lc_shrubs = c("shrubs", "lc_shrubs", "Shrubland"),
    lc_cropland = c("cropland", "lc_cropland", "Cropland"),
    lc_bare = c("bare", "lc_bare", "Bare"),
    lc_wetland = c("wetland", "lc_wetland", "Wetland"),
    lc_built = c("built", "lc_built", "Built"),
    lc_water = c("water", "lc_water", "Water"),
    lc_snow = c("snow", "lc_snow", "Snow"),
    lc_mangroves = c("mangroves", "lc_mangroves", "Mangroves"),
    lc_moss = c("moss", "lc_moss", "Moss")
  )
  
  for (var_name in names(lc_mapping)) {
    
    # Skip if already present with valid data
    if (var_name %in% names(soil_data) && sum(!is.na(soil_data[[var_name]])) > 0) {
      cat(sprintf("  %-20s: already present\n", var_name))
      next
    }
    
    # Find matching band
    possible_names <- lc_mapping[[var_name]]
    band_idx <- NULL
    
    for (pn in possible_names) {
      idx <- grep(paste0("^", pn, "$"), names(lc_vals), ignore.case = TRUE)
      if (length(idx) > 0) {
        band_idx <- idx[1]
        break
      }
      # Also try partial match
      idx <- grep(pn, names(lc_vals), ignore.case = TRUE)
      if (length(idx) > 0) {
        band_idx <- idx[1]
        break
      }
    }
    
    if (!is.null(band_idx)) {
      soil_data[[var_name]] <- lc_vals[[band_idx]]
      cat(sprintf("  %-20s: extracted from band '%s' (%d valid)\n", 
                  var_name, names(lc_vals)[band_idx], sum(!is.na(lc_vals[[band_idx]]))))
    } else {
      cat(sprintf("  %-20s: no matching band found\n", var_name))
    }
  }
  
} else {
  cat("Landcover fractions file not found:", lc_fractions_file, "\n")
}

# =============================================================================
# EXTRACT MULTI-BAND RASTERS: MYCORRHIZA (BARCELO)
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  EXTRACTING MYCORRHIZA DATA (BARCELO)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

myc_barcelo_file <- file.path(RASTER_DIR, "microbial/mycorrhiza_barcelo_0.1deg.tif")

if (file.exists(myc_barcelo_file)) {
  
  cat("Loading mycorrhiza Barcelo raster...\n")
  myc_rast <- rast(myc_barcelo_file)
  cat("  Bands:", nlyr(myc_rast), "\n")
  cat("  Band names:", paste(names(myc_rast), collapse = ", "), "\n\n")
  
  myc_vals <- terra::extract(myc_rast, pts, ID = FALSE)
  
  # Expected variables and possible band names
  myc_mapping <- list(
    AM_roots_colonized = c("AM_roots", "AM_colonized", "AM"),
    EcM_roots_colonized = c("EcM_roots", "EcM_colonized", "EcM"),
    EcM_AM_root_ratio = c("EcM_AM_ratio", "ratio")
  )
  
  for (var_name in names(myc_mapping)) {
    
    if (var_name %in% names(soil_data) && sum(!is.na(soil_data[[var_name]])) > 0) {
      cat(sprintf("  %-25s: already present\n", var_name))
      next
    }
    
    # Try to find matching band
    for (pn in myc_mapping[[var_name]]) {
      idx <- grep(pn, names(myc_vals), ignore.case = TRUE)
      if (length(idx) > 0) {
        soil_data[[var_name]] <- myc_vals[[idx[1]]]
        cat(sprintf("  %-25s: extracted from band '%s'\n", var_name, names(myc_vals)[idx[1]]))
        break
      }
    }
  }
  
} else {
  cat("Barcelo mycorrhiza file not found\n")
}

# =============================================================================
# EXTRACT MULTI-BAND RASTERS: MYCORRHIZA (SPUN)
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  EXTRACTING MYCORRHIZA DATA (SPUN)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

myc_spun_file <- file.path(RASTER_DIR, "microbial/mycorrhiza_spun_0.1deg.tif")

if (file.exists(myc_spun_file)) {
  
  cat("Loading mycorrhiza SPUN raster...\n")
  myc_rast <- rast(myc_spun_file)
  cat("  Bands:", nlyr(myc_rast), "\n")
  cat("  Band names:", paste(names(myc_rast), collapse = ", "), "\n\n")
  
  myc_vals <- terra::extract(myc_rast, pts, ID = FALSE)
  
  myc_mapping <- list(
    AM_richness = c("AM_richness", "AM_Fungi_Richness"),
    EcM_richness = c("EcM_richness", "EcM_Fungi_Richness"),
    AM_endemism = c("AM_endemism", "AM_RWR", "AM_Fungi_RWR"),
    EcM_endemism = c("EcM_endemism", "EcM_RWR", "EcM_Fungi_RWR"),
    EcM_AM_richness_ratio = c("EcM_AM_richness_ratio", "richness_ratio")
  )
  
  for (var_name in names(myc_mapping)) {
    
    if (var_name %in% names(soil_data) && sum(!is.na(soil_data[[var_name]])) > 0) {
      cat(sprintf("  %-25s: already present\n", var_name))
      next
    }
    
    for (pn in myc_mapping[[var_name]]) {
      idx <- grep(pn, names(myc_vals), ignore.case = TRUE)
      if (length(idx) > 0) {
        soil_data[[var_name]] <- myc_vals[[idx[1]]]
        cat(sprintf("  %-25s: extracted from band '%s'\n", var_name, names(myc_vals)[idx[1]]))
        break
      }
    }
  }
  
} else {
  cat("SPUN mycorrhiza file not found\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

n_end <- ncol(soil_data)
n_added <- n_end - n_start

cat("Starting variables:", n_start, "\n")
cat("Final variables:", n_end, "\n")
cat("Variables added:", n_added, "\n\n")

# Check coverage of key predictor groups
check_coverage <- function(vars, group_name) {
  present <- vars[vars %in% names(soil_data)]
  cat(group_name, ":", length(present), "/", length(vars), "\n")
  for (v in vars) {
    if (v %in% names(soil_data)) {
      n_valid <- sum(!is.na(soil_data[[v]]))
      pct <- 100 * n_valid / nrow(soil_data)
      status <- ifelse(pct > 90, "✓", ifelse(pct > 50, "~", "!"))
      cat(sprintf("  %s %-28s: %5.1f%%\n", status, v, pct))
    } else {
      cat(sprintf("  ✗ %-28s: MISSING\n", v))
    }
  }
  cat("\n")
}

CLIMATE_VARS <- c("temperature_2m_mean", "temperature_seasonality", 
                  "soil_temperature_0_20cm", "soil_moisture_0_20cm",
                  "snow_cover_mean", "potential_evaporation", 
                  "lai_high_veg", "lai_low_veg")

PHYSICAL_VARS <- c("clay_total", "sand_total", "silt_total", "bd_filled",
                   "terrain_elev_mean", "terrain_slope_mean",
                   "lc_trees", "lc_grassland", "lc_cropland",
                   "hansen_any_loss", "burn_count_before_obs")

CHEMICAL_VARS <- c("ph_h2o", "sg_cec", "sg_nitrogen", "soilclass_func_code")

BIOLOGICAL_VARS <- c("fungal_proportion", "AM_roots_colonized", 
                     "EcM_roots_colonized", "AM_richness", "EcM_richness")

check_coverage(CLIMATE_VARS, "CLIMATE")
check_coverage(PHYSICAL_VARS, "PHYSICAL")
check_coverage(CHEMICAL_VARS, "CHEMICAL")
check_coverage(BIOLOGICAL_VARS, "BIOLOGICAL")

# =============================================================================
# SAVE OUTPUT
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  SAVING OUTPUT\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

saveRDS(soil_data, OUTPUT_FILE)
cat("✓ Saved:", OUTPUT_FILE, "\n")
cat("  Observations:", format(nrow(soil_data), big.mark = ","), "\n")
cat("  Variables:", ncol(soil_data), "\n")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  STEP 12b COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Next: source('13_model_fitting.R')\n")