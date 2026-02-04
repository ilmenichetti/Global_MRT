# =============================================================================
# Utility functions for point extraction AND spatial rasterization
# =============================================================================

# --- Standard point extraction helper ---
extract_raster_at_points <- function(raster, coords_df, lon_col, lat_col) {
  pts <- terra::vect(
    coords_df,
    geom = c(lon_col, lat_col),
    crs = "EPSG:4326"
  )
  vals <- terra::extract(raster, pts, ID = FALSE)
  return(vals)
}



# --- SPATIAL WORKFLOW: Create/resample to standard global grid ---

#' Create a standardized 0.5° global raster grid
create_global_grid <- function(resolution = 0.5) {
  # Create empty raster at specified resolution
  r <- terra::rast(
    xmin = -180, xmax = 180,
    ymin = -90, ymax = 90,
    resolution = resolution,
    crs = "EPSG:4326"
  )
  return(r)
}

#' Resample any raster to standard resolution and save
resample_and_save <- function(input_raster, output_filename, 
                              resolution = PARAMS$spatial_resolution,
                              method = "bilinear") {
  
  # Create template grid
  template <- create_global_grid(resolution)
  
  # Resample input to template
  resampled <- terra::resample(
    input_raster,
    template,
    method = method
  )
  
  # Ensure output directory exists
  output_dir <- dirname(output_filename)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save as GeoTIFF (compressed)
  terra::writeRaster(
    resampled,
    output_filename,
    datatype = PARAMS$raster_datatype,
    filetype = "GTiff",
    overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", "ZLEVEL=9")  # high compression
  )
  
  message(sprintf("  Saved: %s", output_filename))
  return(resampled)
}

#' Aggregate point data to 0.5° grid
aggregate_points_to_grid <- function(points_sf, variable_col, 
                                     resolution = PARAMS$spatial_resolution,
                                     fun = "mean") {
  
  # Create template grid
  template <- create_global_grid(resolution)
  
  # Rasterize with aggregation function
  rasterized <- terra::rasterize(
    points_sf,
    template,
    values = variable_col,
    fun = fun
  )
  
  return(rasterized)
}

#' Aggregate and save point-based data
aggregate_and_save <- function(points_df, lon_col, lat_col, value_col,
                               output_filename, fun = "mean") {
  
  # Convert to sf
  pts_sf <- sf::st_as_sf(
    points_df,
    coords = c(lon_col, lat_col),
    crs = "EPSG:4326"
  )
  
  # Rasterize to grid
  rast_agg <- terra::rasterize(
    terra::vect(pts_sf),
    create_global_grid(PARAMS$spatial_resolution),
    values = points_df[[value_col]],
    fun = fun
  )
  
  # Save
  output_dir <- dirname(output_filename)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  terra::writeRaster(
    rast_agg,
    output_filename,
    datatype = PARAMS$raster_datatype,
    filetype = "GTiff",
    overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", "ZLEVEL=9")
  )
  
  message(sprintf("  Saved aggregated grid: %s", output_filename))
  return(rast_agg)
}

# --- Depth-weighted average for soil layers ---
calc_depth_weighted_avg <- function(values, top_depths, bot_depths, 
                                    target_top, target_bot) {
  overlaps <- (bot_depths > target_top) & (top_depths < target_bot)
  if (sum(overlaps & !is.na(values)) == 0) return(NA)
  
  weighted_sum <- 0
  total_weight <- 0
  
  for (i in which(overlaps)) {
    if (is.na(values[i])) next
    overlap_top <- max(top_depths[i], target_top)
    overlap_bot <- min(bot_depths[i], target_bot)
    thickness <- overlap_bot - overlap_top
    weighted_sum <- weighted_sum + values[i] * thickness
    total_weight <- total_weight + thickness
  }
  
  if (total_weight == 0) return(NA)
  return(weighted_sum / total_weight)
}

# --- Exponential depth standardization for SOC ---
calc_interval_avg_exp <- function(z1, z2, a, b) {
  ifelse(
    a <= 0,
    NA,
    (b / a) * (exp(-a * z1 / 100) - exp(-a * z2 / 100))
  )
}

# --- Caching wrapper: load if exists, else compute and save ---
load_or_run <- function(filename, compute_fn, subdir = NULL) {
  
  if (!is.null(subdir)) {
    filepath <- file.path(PATHS$output_dir, subdir, filename)
  } else {
    filepath <- file.path(PATHS$output_dir, filename)
  }
  
  if (file.exists(filepath)) {
    log_step(sprintf("Loading cached: %s", filename))
    return(readRDS(filepath))
  }
  
  log_step(sprintf("Computing: %s", filename))
  result <- compute_fn()
  
  # Create subdirectory if needed
  output_dir <- dirname(filepath)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  saveRDS(result, filepath)
  message(sprintf("  Saved: %s", filepath))
  return(result)
}

# --- Progress reporting for batch operations ---
batch_progress <- function(i, total, batch_size = 1000) {
  if (i %% batch_size == 0 || i == total) {
    cat(sprintf("\r  Progress: %d / %d (%.1f%%)",
                i, total, 100 * i / total))
    flush.console()
  }
}

# --- Logging utility ---
log_step <- function(msg) {
  # Print a timestamped log message to the console
  # Format: [HH:MM:SS] message
  # Used to track progress through pipeline steps
  cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

message("Utilities loaded.")