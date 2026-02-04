# =============================================================================
# Step 7: Extract Köppen-Geiger Climate Classification
# Source: Beck et al. (2023) - 1km resolution maps
# Period: 1991-2020
# Outputs: Point extraction + 0.1° spatialized rasters
# =============================================================================

library(terra)
library(dplyr)

# Source config and utilities
source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# --- Configuration ---
KOPPEN_FILE <- "./Global_MRT_code/koppen_data/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"
KOPPEN_OUTPUT_DIR <- "./Global_MRT_code/spatialized_layers/climate"
TARGET_RES <- 0.1  # degrees


# =============================================================================
# Köppen-Geiger Legend (Beck et al. 2023)
# Values 1-30 in raster → Köppen codes
# =============================================================================

KOPPEN_LEGEND <- data.frame(
  value = 1:30,
  code = c(
    "Af",  "Am",  "Aw",           # 1-3:   Tropical (A)
    "BWh", "BWk",                  # 4-5:   Arid desert (B)
    "BSh", "BSk",                  # 6-7:   Arid steppe (B)
    "Csa", "Csb", "Csc",          # 8-10:  Temperate dry summer (C)
    "Cwa", "Cwb", "Cwc",          # 11-13: Temperate dry winter (C)
    "Cfa", "Cfb", "Cfc",          # 14-16: Temperate no dry season (C)
    "Dsa", "Dsb", "Dsc", "Dsd",   # 17-20: Cold dry summer (D)
    "Dwa", "Dwb", "Dwc", "Dwd",   # 21-24: Cold dry winter (D)
    "Dfa", "Dfb", "Dfc", "Dfd",   # 25-28: Cold no dry season (D)
    "ET",  "EF"                    # 29-30: Polar (E)
  ),
  description = c(
    "Tropical rainforest",
    "Tropical monsoon",
    "Tropical savanna",
    "Hot desert",
    "Cold desert",
    "Hot steppe",
    "Cold steppe",
    "Mediterranean hot summer",
    "Mediterranean warm summer",
    "Mediterranean cold summer",
    "Humid subtropical dry winter",
    "Subtropical highland dry winter",
    "Subpolar oceanic dry winter",
    "Humid subtropical",
    "Oceanic",
    "Subpolar oceanic",
    "Continental hot dry summer",
    "Continental warm dry summer",
    "Continental cool dry summer",
    "Continental extreme dry summer",
    "Continental hot dry winter",
    "Continental warm dry winter",
    "Continental cool dry winter",
    "Continental extreme dry winter",
    "Hot-summer humid continental",
    "Warm-summer humid continental",
    "Subarctic",
    "Subarctic extreme",
    "Tundra",
    "Ice cap"
  ),
  # Main climate group (first letter)
  main_group = c(
    rep("A", 3),   # Tropical
    rep("B", 4),   # Arid
    rep("C", 9),   # Temperate
    rep("D", 12),  # Continental/Cold
    rep("E", 2)    # Polar
  ),
  main_group_num = c(
    rep(1, 3),     # A = 1
    rep(2, 4),     # B = 2
    rep(3, 9),     # C = 3
    rep(4, 12),    # D = 4
    rep(5, 2)      # E = 5
  ),
  main_description = c(
    rep("Tropical", 3),
    rep("Arid", 4),
    rep("Temperate", 9),
    rep("Continental", 12),
    rep("Polar", 2)
  ),
  stringsAsFactors = FALSE
)


# =============================================================================
# Helper Functions
# =============================================================================

#' Load Köppen raster and validate
load_koppen_raster <- function(filepath = KOPPEN_FILE) {
  log_step("Loading Köppen-Geiger raster")
  
  if (!file.exists(filepath)) {
    stop("Köppen raster not found: ", filepath)
  }
  
  r <- rast(filepath)
  cat(sprintf("  File: %s\n", basename(filepath)))
  cat(sprintf("  Resolution: %.6f° (~%.1f km at equator)\n", 
              res(r)[1], res(r)[1] * 111))
  cat(sprintf("  Extent: %.1f to %.1f lon, %.1f to %.1f lat\n",
              ext(r)[1], ext(r)[2], ext(r)[3], ext(r)[4]))
  cat(sprintf("  Values: %d to %d\n", 
              global(r, "min", na.rm = TRUE)[1,1],
              global(r, "max", na.rm = TRUE)[1,1]))
  
  return(r)
}


#' Create lookup vectors for fast value mapping
create_koppen_lookups <- function() {
  # Pre-allocate lookup vectors (index = raster value, value = code/group)
  # Maximum raster value is 30
  
  code_lookup <- rep(NA_character_, 30)
  code_lookup[KOPPEN_LEGEND$value] <- KOPPEN_LEGEND$code
  
  desc_lookup <- rep(NA_character_, 30)
  desc_lookup[KOPPEN_LEGEND$value] <- KOPPEN_LEGEND$description
  
  main_group_lookup <- rep(NA_character_, 30)
  main_group_lookup[KOPPEN_LEGEND$value] <- KOPPEN_LEGEND$main_group
  
  main_num_lookup <- rep(NA_integer_, 30)
  main_num_lookup[KOPPEN_LEGEND$value] <- KOPPEN_LEGEND$main_group_num
  
  main_desc_lookup <- rep(NA_character_, 30)
  main_desc_lookup[KOPPEN_LEGEND$value] <- KOPPEN_LEGEND$main_description
  
  list(
    code = code_lookup,
    description = desc_lookup,
    main_group = main_group_lookup,
    main_group_num = main_num_lookup,
    main_description = main_desc_lookup
  )
}


#' Extract Köppen values at points and map to codes
extract_koppen_at_points <- function(koppen_rast, points_data, lookups) {
  log_step("Extracting Köppen classification at points")
  
  # Create point geometry
  pts <- vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %s points...\n", format(nrow(points_data), big.mark = ",")))
  
  # Extract raw values
  raw_vals <- terra::extract(koppen_rast, pts, ID = FALSE)[, 1]
  
  # Map to codes using lookup vectors
  results <- data.frame(
    dsiteid = points_data$dsiteid,
    koppen_value = raw_vals,
    koppen_code = lookups$code[raw_vals],
    koppen_description = lookups$description[raw_vals],
    koppen_main_group = lookups$main_group[raw_vals],
    koppen_main_num = lookups$main_group_num[raw_vals],
    koppen_main_description = lookups$main_description[raw_vals],
    stringsAsFactors = FALSE
  )
  
  # Summary
  n_valid <- sum(!is.na(results$koppen_value))
  cat(sprintf("  Coverage: %s points (%.1f%%)\n", 
              format(n_valid, big.mark = ","),
              100 * n_valid / nrow(results)))
  
  # Distribution by main group
  cat("\n  Distribution by main climate group:\n")
  main_dist <- table(results$koppen_main_group, useNA = "ifany")
  for (grp in names(main_dist)) {
    grp_label <- if (is.na(grp)) "NA" else grp
    cat(sprintf("    %s: %s (%.1f%%)\n", 
                grp_label, 
                format(main_dist[grp], big.mark = ","),
                100 * main_dist[grp] / nrow(results)))
  }
  
  return(results)
}


#' Create spatialized Köppen rasters at target resolution
create_koppen_rasters <- function(koppen_rast, output_dir = KOPPEN_OUTPUT_DIR, 
                                  target_res = TARGET_RES) {
  log_step("Creating spatialized Köppen rasters")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create target grid at 0.1° resolution
  target_ext <- ext(-180, 180, -90, 90)
  target_rast <- rast(
    extent = target_ext,
    resolution = target_res,
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Target resolution: %.2f° (%d x %d pixels)\n", 
              target_res, ncol(target_rast), nrow(target_rast)))
  
  # --- Detailed classification (30 classes) ---
  cat("  Resampling detailed classification (modal aggregation)...\n")
  
  # Use modal (most common value) for categorical data
  koppen_detailed <- resample(koppen_rast, target_rast, method = "mode")
  
  detailed_file <- file.path(output_dir, "koppen_detailed_0.1deg.tif")
  writeRaster(koppen_detailed, detailed_file, overwrite = TRUE,
              datatype = "INT1U", gdal = c("COMPRESS=LZW"))
  cat(sprintf("  ✓ Saved: %s\n", basename(detailed_file)))
  
  # --- Main group classification (5 classes) ---
  cat("  Creating main group classification...\n")
  
  # Reclassify to main groups (1-5)
  # Build reclassification matrix: from, to, becomes
  rcl_matrix <- cbind(
    from = c(0.5, 3.5, 7.5, 16.5, 28.5),
    to   = c(3.5, 7.5, 16.5, 28.5, 30.5),
    becomes = c(1, 2, 3, 4, 5)
  )
  
  koppen_main <- classify(koppen_detailed, rcl_matrix)
  
  main_file <- file.path(output_dir, "koppen_main_0.1deg.tif")
  writeRaster(koppen_main, main_file, overwrite = TRUE,
              datatype = "INT1U", gdal = c("COMPRESS=LZW"))
  cat(sprintf("  ✓ Saved: %s\n", basename(main_file)))
  
  # File sizes
  cat(sprintf("\n  File sizes:\n"))
  cat(sprintf("    Detailed: %.1f MB\n", file.size(detailed_file) / 1e6))
  cat(sprintf("    Main:     %.1f MB\n", file.size(main_file) / 1e6))
  
  return(list(
    detailed = koppen_detailed,
    main = koppen_main,
    detailed_file = detailed_file,
    main_file = main_file
  ))
}


#' Create diagnostic plots
create_koppen_diagnostics <- function(data, koppen_rasters, output_dir = "./plots/07_koppen") {
  log_step("Creating diagnostic outputs")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- 1. Coverage table ---
  coverage_file <- file.path(output_dir, "koppen_coverage.csv")
  
  coverage_summary <- data %>%
    group_by(koppen_main_group, koppen_main_description) %>%
    summarise(
      n_sites = n(),
      pct = round(100 * n() / nrow(data), 2),
      .groups = "drop"
    ) %>%
    arrange(koppen_main_group)
  
  write.csv(coverage_summary, coverage_file, row.names = FALSE)
  cat(sprintf("  ✓ Coverage table: %s\n", basename(coverage_file)))
  
  # --- 2. Detailed distribution ---
  detailed_file <- file.path(output_dir, "koppen_detailed_distribution.csv")
  
  detailed_summary <- data %>%
    group_by(koppen_value, koppen_code, koppen_description) %>%
    summarise(n_sites = n(), .groups = "drop") %>%
    arrange(koppen_value)
  
  write.csv(detailed_summary, detailed_file, row.names = FALSE)
  cat(sprintf("  ✓ Detailed distribution: %s\n", basename(detailed_file)))
  
  # --- 3. Global map of main groups ---
  map_file <- file.path(output_dir, "koppen_main_map.png")
  
  png(map_file, width = 1600, height = 800, res = 150)
  
  # Define colors for 5 main groups
  main_colors <- c(
    "1" = "#0000FF",   # A - Tropical (blue)
    "2" = "#FF0000",   # B - Arid (red)
    "3" = "#00FF00",   # C - Temperate (green)
    "4" = "#00FFFF",   # D - Continental (cyan)
    "5" = "#FFFFFF"    # E - Polar (white)
  )
  
  plot(koppen_rasters$main, 
       main = "Köppen-Geiger Main Climate Groups (1991-2020)",
       col = main_colors,
       legend = FALSE,
       axes = TRUE)
  
  legend("bottomleft", 
         legend = c("A: Tropical", "B: Arid", "C: Temperate", 
                    "D: Continental", "E: Polar"),
         fill = main_colors,
         bg = "white",
         cex = 0.8)
  
  dev.off()
  cat(sprintf("  ✓ Main groups map: %s\n", basename(map_file)))
  
  # --- 4. Point distribution bar chart ---
  bar_file <- file.path(output_dir, "koppen_point_distribution.png")
  
  png(bar_file, width = 1200, height = 600, res = 150)
  
  par(mar = c(8, 4, 3, 1))
  
  # Count by detailed code
  code_counts <- table(data$koppen_code)
  code_counts <- sort(code_counts, decreasing = TRUE)
  
  # Assign colors based on main group
  bar_colors <- sapply(names(code_counts), function(code) {
    first_letter <- substr(code, 1, 1)
    switch(first_letter,
           "A" = "#0000FF",
           "B" = "#FF6600",
           "C" = "#00AA00",
           "D" = "#00AAAA",
           "E" = "#888888",
           "#000000")
  })
  
  barplot(code_counts,
          main = "Distribution of Soil Sites by Köppen-Geiger Class",
          ylab = "Number of sites",
          col = bar_colors,
          las = 2,
          cex.names = 0.7)
  
  dev.off()
  cat(sprintf("  ✓ Distribution bar chart: %s\n", basename(bar_file)))
}


# =============================================================================
# Main Execution
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║         Step 7: Köppen-Geiger Climate Classification          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# --- Step 1: Load input data ---
log_step("Loading input data")

input_file <- file.path(PATHS$output_dir, "06_with_soilgrids.rds")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, 
       "\nRun 06_extract_soilgrids.R first")
}

data <- readRDS(input_file)
cat(sprintf("  Loaded %s observations\n", format(nrow(data), big.mark = ",")))


# --- Step 2: Get unique site locations ---
log_step("Identifying unique locations")

sites <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees)

cat(sprintf("  %s unique site locations\n", format(nrow(sites), big.mark = ",")))


# --- Step 3: Load Köppen raster ---
koppen_rast <- load_koppen_raster(KOPPEN_FILE)


# --- Step 4: Create lookup tables ---
lookups <- create_koppen_lookups()


# --- Step 5: Extract at points ---
koppen_extracted <- load_or_run(
  "07_koppen_extracted.rds",
  function() {
    extract_koppen_at_points(koppen_rast, sites, lookups)
  }
)


# --- Step 6: Create spatialized rasters ---
koppen_rasters <- create_koppen_rasters(koppen_rast, KOPPEN_OUTPUT_DIR, TARGET_RES)


# --- Step 7: Merge with main dataset ---
log_step("Merging Köppen data with main dataset")

data <- data %>%
  left_join(koppen_extracted, by = "dsiteid")

cat(sprintf("  ✓ Added %d Köppen variables\n", 
            ncol(koppen_extracted) - 1))  # minus dsiteid


# --- Step 8: Save output ---
log_step("Saving output")

output_file <- file.path(PATHS$output_dir, "07_with_koppen.rds")
saveRDS(data, output_file)

cat(sprintf("  ✓ Saved: %s\n", basename(output_file)))
cat(sprintf("  Dataset: %s rows, %d columns\n", 
            format(nrow(data), big.mark = ","), ncol(data)))


# --- Step 9: Create diagnostics ---
create_koppen_diagnostics(data, koppen_rasters)


# --- Summary ---
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("                           SUMMARY                                  \n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Variables added:\n")
cat("  koppen_value          - Raw raster value (1-30)\n")
cat("  koppen_code           - Köppen code (e.g., 'Af', 'Cfb', 'Dfc')\n")
cat("  koppen_description    - Full description\n")
cat("  koppen_main_group     - Main group letter (A, B, C, D, E)\n")
cat("  koppen_main_num       - Main group numeric (1-5) for modeling\n")
cat("  koppen_main_description - Main group name\n")

cat("\nSpatialized rasters created:\n")
cat(sprintf("  %s\n", koppen_rasters$detailed_file))
cat(sprintf("  %s\n", koppen_rasters$main_file))

cat("\nNext step: Add DEM-derived covariates (elevation, slope, aspect)\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")







# =============================================================================
# Step 7b: Extract Terrain Metrics from Copernicus DEM
# Source: GEE-exported 0.1° rasters derived from Copernicus GLO-30
# =============================================================================

library(terra)
library(dplyr)

# Source config and utilities
source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# --- Configuration ---
TERRAIN_DIR <- "./Global_MRT_code/spatialized_layers/topography"

TERRAIN_FILES <- list(
  elev_mean  = "elev_mean_0p1deg.tif",
  slope_mean = "slope_mean_0p1deg.tif",
  northness  = "northness_0p1deg.tif",
  eastness   = "eastness_0p1deg.tif"
)


# =============================================================================
# Helper Functions
# =============================================================================

#' Load terrain rasters and validate
load_terrain_rasters <- function(terrain_dir = TERRAIN_DIR, files = TERRAIN_FILES) {
  log_step("Loading terrain rasters")
  
  rasters <- list()
  
  for (var in names(files)) {
    filepath <- file.path(terrain_dir, files[[var]])
    
    if (!file.exists(filepath)) {
      warning(sprintf("  ✗ %s: not found (%s)", var, files[[var]]))
      next
    }
    
    rasters[[var]] <- rast(filepath)
    
    # Get basic stats
    r_min <- global(rasters[[var]], "min", na.rm = TRUE)[1,1]
    r_max <- global(rasters[[var]], "max", na.rm = TRUE)[1,1]
    
    cat(sprintf("  ✓ %s: %.1f to %.1f\n", var, r_min, r_max))
  }
  
  if (length(rasters) == 0) {
    stop("No terrain files found in: ", terrain_dir)
  }
  
  # Check resolution consistency
  res_check <- sapply(rasters, function(r) res(r)[1])
  cat(sprintf("\n  Resolution: %.4f° (all layers)\n", res_check[1]))
  
  return(rasters)
}


#' Extract terrain values at point locations
extract_terrain_at_points <- function(rasters, points_data) {
  log_step("Extracting terrain metrics at points")
  
  # Create point geometry
  pts <- vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %s points...\n", format(nrow(points_data), big.mark = ",")))
  
  # Initialize results
  results <- data.frame(dsiteid = points_data$dsiteid)
  
  for (var in names(rasters)) {
    cat(sprintf("  Extracting: %s\n", var))
    
    vals <- terra::extract(rasters[[var]], pts, ID = FALSE)[, 1]
    
    # Store with terrain_ prefix
    col_name <- paste0("terrain_", var)
    results[[col_name]] <- vals
    
    # Coverage stats
    n_valid <- sum(!is.na(vals))
    cat(sprintf("    Coverage: %s (%.1f%%)\n", 
                format(n_valid, big.mark = ","),
                100 * n_valid / length(vals)))
  }
  
  return(results)
}


#' Create diagnostic outputs
create_terrain_diagnostics <- function(data, rasters, output_dir = "./plots/07b_terrain") {
  log_step("Creating diagnostic outputs")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- 1. Coverage summary ---
  coverage_file <- file.path(output_dir, "terrain_coverage.csv")
  
  terrain_cols <- grep("^terrain_", names(data), value = TRUE)
  
  coverage_summary <- data.frame(
    variable = terrain_cols,
    n_valid = sapply(terrain_cols, function(col) sum(!is.na(data[[col]]))),
    pct_valid = sapply(terrain_cols, function(col) 
      round(100 * sum(!is.na(data[[col]])) / nrow(data), 1)),
    mean = sapply(terrain_cols, function(col) round(mean(data[[col]], na.rm = TRUE), 2)),
    sd = sapply(terrain_cols, function(col) round(sd(data[[col]], na.rm = TRUE), 2)),
    min = sapply(terrain_cols, function(col) round(min(data[[col]], na.rm = TRUE), 2)),
    max = sapply(terrain_cols, function(col) round(max(data[[col]], na.rm = TRUE), 2))
  )
  
  write.csv(coverage_summary, coverage_file, row.names = FALSE)
  cat(sprintf("  ✓ Coverage summary: %s\n", basename(coverage_file)))
  
  # --- 2. Global elevation map ---
  elev_map_file <- file.path(output_dir, "terrain_elevation_map.png")
  
  png(elev_map_file, width = 1600, height = 800, res = 150)
  plot(rasters$elev_mean,
       main = "Mean Elevation (Copernicus DEM, 0.1°)",
       col = terrain.colors(100),
       axes = TRUE)
  dev.off()
  cat(sprintf("  ✓ Elevation map: %s\n", basename(elev_map_file)))
  
  # --- 3. Slope map ---
  slope_map_file <- file.path(output_dir, "terrain_slope_map.png")
  
  png(slope_map_file, width = 1600, height = 800, res = 150)
  plot(rasters$slope_mean,
       main = "Mean Slope (degrees, 0.1°)",
       col = hcl.colors(100, "YlOrRd"),
       axes = TRUE)
  dev.off()
  cat(sprintf("  ✓ Slope map: %s\n", basename(slope_map_file)))
  
  # --- 3b. Ruggedness map ---
  if ("ruggedness" %in% names(rasters)) {
    rugged_map_file <- file.path(output_dir, "terrain_ruggedness_map.png")
    
    png(rugged_map_file, width = 1600, height = 800, res = 150)
    plot(rasters$ruggedness,
         main = "Terrain Ruggedness (elevation SD, 0.1°)",
         col = hcl.colors(100, "Viridis"),
         axes = TRUE)
    dev.off()
    cat(sprintf("  ✓ Ruggedness map: %s\n", basename(rugged_map_file)))
  }
  
  # --- 4. Histograms of point values ---
  hist_file <- file.path(output_dir, "terrain_histograms.png")
  
  png(hist_file, width = 1200, height = 1000, res = 150)
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
  
  hist(data$terrain_elev_mean, breaks = 50, col = "steelblue",
       main = "Elevation at sites", xlab = "Elevation (m)")
  
  hist(data$terrain_slope_mean, breaks = 50, col = "tomato",
       main = "Slope at sites", xlab = "Slope (degrees)")
  
  hist(data$terrain_northness, breaks = 50, col = "forestgreen",
       main = "Northness at sites", xlab = "Northness (-1=S, +1=N)")
  
  hist(data$terrain_eastness, breaks = 50, col = "purple",
       main = "Eastness at sites", xlab = "Eastness (-1=W, +1=E)")
  
  hist(data$terrain_ruggedness, breaks = 50, col = "orange",
       main = "Ruggedness at sites", xlab = "Elevation SD (m)")
  
  dev.off()
  cat(sprintf("  ✓ Histograms: %s\n", basename(hist_file)))
  
  # --- 5. Elevation vs latitude scatter ---
  scatter_file <- file.path(output_dir, "terrain_elev_vs_lat.png")
  
  png(scatter_file, width = 1000, height = 600, res = 150)
  
  # Sample for plotting if too many points
  plot_data <- if (nrow(data) > 10000) data[sample(nrow(data), 10000), ] else data
  
  plot(plot_data$latitude_decimal_degrees, plot_data$terrain_elev_mean,
       pch = ".", col = rgb(0, 0, 0, 0.2),
       xlab = "Latitude", ylab = "Elevation (m)",
       main = "Site elevation vs latitude")
  
  dev.off()
  cat(sprintf("  ✓ Elevation vs latitude: %s\n", basename(scatter_file)))
}


# =============================================================================
# Main Execution
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║         Step 7b: Terrain Metrics (Copernicus DEM)              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# --- Step 1: Load input data ---
log_step("Loading input data")

input_file <- file.path(PATHS$output_dir, "07_with_koppen.rds")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, 
       "\nRun 07_extract_koppen.R first")
}

data <- readRDS(input_file)
cat(sprintf("  Loaded %s observations\n", format(nrow(data), big.mark = ",")))


# --- Step 2: Get unique site locations ---
log_step("Identifying unique locations")

sites <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees)

cat(sprintf("  %s unique site locations\n", format(nrow(sites), big.mark = ",")))


# --- Step 3: Load terrain rasters ---
terrain_rasters <- load_terrain_rasters(TERRAIN_DIR, TERRAIN_FILES)


# --- Step 3b: Compute ruggedness (SD of elevation in 3x3 window) ---
log_step("Computing terrain ruggedness")

ruggedness_file <- file.path(TERRAIN_DIR, "elev_ruggedness_0p1deg.tif")

if (file.exists(ruggedness_file)) {
  cat("  Loading cached ruggedness raster\n")
  terrain_rasters$ruggedness <- rast(ruggedness_file)
} else {
  cat("  Computing SD of elevation (3x3 focal window)...\n")
  
  # Focal SD in 3x3 window (~30km neighborhood at 0.1°)
  terrain_rasters$ruggedness <- focal(
    terrain_rasters$elev_mean,
    w = 3,
    fun = "sd",
    na.rm = TRUE
  )
  
  # Save for reuse
  writeRaster(terrain_rasters$ruggedness, ruggedness_file, 
              overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  cat(sprintf("  ✓ Saved: %s\n", basename(ruggedness_file)))
}

r_min <- global(terrain_rasters$ruggedness, "min", na.rm = TRUE)[1,1]
r_max <- global(terrain_rasters$ruggedness, "max", na.rm = TRUE)[1,1]
cat(sprintf("  Ruggedness range: %.1f to %.1f m\n", r_min, r_max))


# --- Step 4: Extract at points ---
terrain_extracted <- load_or_run(
  "07b_terrain_extracted.rds",
  function() {
    extract_terrain_at_points(terrain_rasters, sites)
  }
)


# --- Step 5: Merge with main dataset ---
log_step("Merging terrain data with main dataset")

data <- data %>%
  left_join(terrain_extracted, by = "dsiteid")

cat(sprintf("  ✓ Added %d terrain variables\n", 
            ncol(terrain_extracted) - 1))


# --- Step 6: Save output ---
log_step("Saving output")

output_file <- file.path(PATHS$output_dir, "07_with_covariates.rds")
saveRDS(data, output_file)

cat(sprintf("  ✓ Saved: %s\n", basename(output_file)))
cat(sprintf("  Dataset: %s rows, %d columns\n", 
            format(nrow(data), big.mark = ","), ncol(data)))


# --- Step 7: Create diagnostics ---
create_terrain_diagnostics(data, terrain_rasters)


# --- Summary ---
cat("\nCombined output (Köppen + Terrain):\n")
cat(sprintf("  %s\n", output_file))
