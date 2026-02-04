# =============================================================================
# Step 3: Extract ERA5-Land climate data at point locations
# Uses pre-downloaded ERA5 NetCDF file (2001-2020 monthly means)
# Variables: t2m, stl1, stl2, swvl1, swvl2, snowc, pev, lai_hv, lai_lv
# =============================================================================

library(terra)
library(dplyr)

# --- Helper: Extract rasters at point locations ---
#' Extract annual climate rasters at observation point locations
#' @param annual_stats List of annual climate rasters
#' @param points_data Data frame with point locations
#' @return Data frame with extracted climate values
extract_climate_at_points <- function(annual_stats, points_data) {
  log_step("Extracting climate at point locations")
  
  # Create point geometry
  pts <- terra::vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  climate_extracted <- data.frame(
    dsiteid = points_data$dsiteid
  )
  
  # Extract each raster
  for (var_name in names(annual_stats)) {
    cat(sprintf("  Extracting: %s\n", var_name))
    
    tryCatch({
      vals <- terra::extract(annual_stats[[var_name]], pts, ID = FALSE)
      climate_extracted[[var_name]] <- vals[, 1]
    }, error = function(e) {
      warning(sprintf("    Failed: %s", e$message))
      climate_extracted[[var_name]] <- NA_real_
    })
  }
  
  message("  ✓ Extraction complete")
  return(climate_extracted)
}

# --- Helper: Rasterize point-based climate data to 0.5° grid ---
rasterize_climate_to_grid <- function(climate_extracted, points_data) {
  log_step("Rasterizing climate to 0.5° global grid")
  
  template <- create_global_grid(PARAMS$spatial_resolution)
  
  # Combine point data with climate
  points_with_climate <- cbind(points_data, climate_extracted)
  
  # Rasterize each variable
  var_names <- names(climate_extracted)[-1]  # Skip dsiteid
  
  for (var in var_names) {
    cat(sprintf("  Rasterizing: %s\n", var))
    
    output_file <- file.path(
      PATHS$raster_output_dir, "climate",
      sprintf("%s_0.5deg.tif", var)
    )
    
    tryCatch({
      aggregate_and_save(
        points_with_climate,
        lon_col = "longitude_decimal_degrees",
        lat_col = "latitude_decimal_degrees",
        value_col = var,
        output_filename = output_file,
        fun = "mean"
      )
    }, error = function(e) {
      warning(sprintf("    Failed to rasterize %s: %s", var, e$message))
    })
  }
  
  message("  ✓ Rasterization complete")
}

# =============================================================================
# Main execution
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║          Step 3: Climate Extraction (ERA5-Land)               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# --- Load input data ---
log_step("Loading standardized soil data")
data <- load_or_run(
  "02_5_standardized_0_20cm.rds",
  function() {
    stop("Run step 02_5 first to create standardized data")
  }
)
cat(sprintf("  Loaded %s sites\n", format(nrow(data), big.mark = ",")))

# --- Step 1: Check ERA5 data exists ---
era5_file <- "./era5_data/data_stream-moda.nc"

if (!file.exists(era5_file)) {
  stop("ERA5 file not found: ", era5_file, 
       "\nDownload from CDS: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means")
}
cat(sprintf("  ✓ Found ERA5 file: %s\n", era5_file))

# --- Step 2: Load and inspect ERA5 data ---
log_step("Loading ERA5 raster")
r <- rast(era5_file)
cat(sprintf("  Dimensions: %d x %d x %d\n", nrow(r), ncol(r), nlyr(r)))

# Extract variable names
var_names <- unique(gsub("_valid_time=.*", "", names(r)))
cat(sprintf("  Variables found: %s\n", paste(var_names, collapse = ", ")))

# --- Step 3: Extract monthly stacks and calculate annual stats ---
annual_stats <- load_or_run(
  "03_era5_annual_statistics.rds",
  function() {
    log_step("Extracting monthly stacks from ERA5")
    
    # Extract by variable prefix
    t2m_stack <- r[[grep("^t2m_", names(r))]]
    stl1_stack <- r[[grep("^stl1_", names(r))]]
    stl2_stack <- r[[grep("^stl2_", names(r))]]
    swvl1_stack <- r[[grep("^swvl1_", names(r))]]
    swvl2_stack <- r[[grep("^swvl2_", names(r))]]
    snowc_stack <- r[[grep("^snowc_", names(r))]]
    pev_stack <- r[[grep("^pev_", names(r))]]
    lai_hv_stack <- r[[grep("^lai_hv_", names(r))]]
    lai_lv_stack <- r[[grep("^lai_lv_", names(r))]]
    
    cat(sprintf("  t2m layers: %d\n", nlyr(t2m_stack)))
    cat(sprintf("  stl1 layers: %d\n", nlyr(stl1_stack)))
    cat(sprintf("  stl2 layers: %d\n", nlyr(stl2_stack)))
    cat(sprintf("  swvl1 layers: %d\n", nlyr(swvl1_stack)))
    cat(sprintf("  swvl2 layers: %d\n", nlyr(swvl2_stack)))
    cat(sprintf("  snowc layers: %d\n", nlyr(snowc_stack)))
    cat(sprintf("  pev layers: %d\n", nlyr(pev_stack)))
    cat(sprintf("  lai_hv layers: %d\n", nlyr(lai_hv_stack)))
    cat(sprintf("  lai_lv layers: %d\n", nlyr(lai_lv_stack)))
    
    # Calculate annual statistics
    log_step("Calculating annual statistics")
    
    stats <- list()
    
    # Temperature: mean and seasonality
    stats$temperature_2m_mean <- mean(t2m_stack, na.rm = TRUE) - 273.15  # K to C
    stats$temperature_seasonality <- terra::app(t2m_stack, sd)
    cat("  ✓ Temperature\n")
    
    # Soil temperature: depth-weighted 0-20cm mean
    # Layer 1: 0-7cm, Layer 2: 7-28cm
    stl1_mean <- mean(stl1_stack, na.rm = TRUE) - 273.15
    stl2_mean <- mean(stl2_stack, na.rm = TRUE) - 273.15
    stats$soil_temperature_0_20cm <- (stl1_mean * 0.35) + (stl2_mean * 0.65)
    cat("  ✓ Soil temperature (depth-weighted)\n")
    
    # Soil moisture: depth-weighted 0-20cm mean
    swvl1_mean <- mean(swvl1_stack, na.rm = TRUE)
    swvl2_mean <- mean(swvl2_stack, na.rm = TRUE)
    stats$soil_moisture_0_20cm <- (swvl1_mean * 0.35) + (swvl2_mean * 0.65)
    cat("  ✓ Soil moisture (depth-weighted)\n")
    
    # Snow cover
    if (nlyr(snowc_stack) > 0) {
      stats$snow_cover_mean <- mean(snowc_stack, na.rm = TRUE)
      cat("  ✓ Snow cover\n")
    }
    
    # Potential evaporation
    if (nlyr(pev_stack) > 0) {
      stats$potential_evaporation <- mean(pev_stack, na.rm = TRUE)
      cat("  ✓ Potential evaporation\n")
    }
    
    # LAI high vegetation
    if (nlyr(lai_hv_stack) > 0) {
      stats$lai_high_veg <- mean(lai_hv_stack, na.rm = TRUE)
      cat("  ✓ LAI high vegetation\n")
    }
    
    # LAI low vegetation
    if (nlyr(lai_lv_stack) > 0) {
      stats$lai_low_veg <- mean(lai_lv_stack, na.rm = TRUE)
      cat("  ✓ LAI low vegetation\n")
    }
    
    return(stats)
  }
)

cat(sprintf("  Annual statistics: %d variables\n", length(annual_stats)))

# --- Step 4: Extract at point locations ---
climate_extracted <- load_or_run(
  "03_climate_extracted.rds",
  function() {
    extract_climate_at_points(annual_stats, data)
  }
)

cat(sprintf("  Extracted %d climate variables at %s points\n", 
            ncol(climate_extracted) - 1, 
            format(nrow(climate_extracted), big.mark = ",")))

# --- Step 5: Merge with main dataset ---
data <- data %>%
  left_join(climate_extracted, by = "dsiteid")

cat(sprintf("  ✓ Merged: dataset now has %d columns\n", ncol(data)))

# --- Step 6: Save updated dataset ---
saveRDS(data, file.path(PATHS$output_dir, "03_with_climate.rds"))
cat(sprintf("  ✓ Saved: %s\n", file.path(PATHS$output_dir, "03_with_climate.rds")))

# --- Step 7: Summary ---
cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                    Climate Extraction Summary                  ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

climate_vars <- names(climate_extracted)[-1]
for (var in climate_vars) {
  vals <- data[[var]]
  n_valid <- sum(!is.na(vals))
  pct <- 100 * n_valid / nrow(data)
  if (n_valid > 0) {
    cat(sprintf("  %-30s: %d sites (%.1f%%), mean = %.2f\n",
                var, n_valid, pct, mean(vals, na.rm = TRUE)))
  } else {
    cat(sprintf("  %-30s: NO DATA\n", var))
  }
}

log_step("Climate extraction complete")



# --- Step 6b: Save climate rasters as GeoTIFFs ---
cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║           Saving Geotiff for model extrapolation               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

log_step("Saving climate rasters to GeoTIFF")

climate_raster_dir <- file.path(PATHS$raster_output_dir, "climate")
if (!dir.exists(climate_raster_dir)) {
  dir.create(climate_raster_dir, recursive = TRUE, showWarnings = FALSE)
}

for (var_name in names(annual_stats)) {
  output_file <- file.path(climate_raster_dir, sprintf("%s_0.1deg.tif", var_name))
  cat(sprintf("  Saving: %s\n", var_name))
  
  resample_and_save(
    annual_stats[[var_name]],
    output_file,
    resolution = 0.1,
    method = "bilinear"
  )
}
cat(sprintf("  ✓ Saved %d climate rasters to %s\n", length(annual_stats), climate_raster_dir))
