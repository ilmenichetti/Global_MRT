# =============================================================================
# Step 08: Extract Soil Classification (SoilGrids WRB)
# Source: SoilGrids WRB MostProbable (250m) via ISRIC COG/VRT
# Downloads, aggregates to 0.1°, and creates functional group raster
# =============================================================================

library(terra)
library(dplyr)

# Source config and utilities
source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# --- Configuration ---
SOILCLASS_DIR <- "./Global_MRT_code/spatialized_layers/soilclass"
WRB_FILE <- "soilgrids_wrb_0p1deg.tif"

# SoilGrids WRB source (250m resolution, COG via HTTPS)
SOILGRIDS_WRB_URL <- "/vsicurl/https://files.isric.org/soilgrids/latest/data/wrb/MostProbable.vrt"


# =============================================================================
# WRB Legend (SoilGrids numeric codes 0-29)
# =============================================================================

WRB_LEGEND <- data.frame(
  value = 0:29,
  wrb_code = c(
    "AC", "AB", "AL", "AN", "AR",  # 0-4
    "CL", "CM", "CH", "CR", "DU",  # 5-9
    "FR", "FL", "GL", "GY", "HS",  # 10-14
    "KS", "LP", "LX", "LV", "NT",  # 15-19
    "PH", "PL", "PT", "PZ", "RG",  # 20-24
    "SC", "SN", "ST", "UM", "VR"   # 25-29
  ),
  wrb_name = c(
    "Acrisols", "Albeluvisols", "Alisols", "Andosols", "Arenosols",
    "Calcisols", "Cambisols", "Chernozems", "Cryosols", "Durisols",
    "Ferralsols", "Fluvisols", "Gleysols", "Gypsisols", "Histosols",
    "Kastanozems", "Leptosols", "Lixisols", "Luvisols", "Nitisols",
    "Phaeozems", "Planosols", "Plinthosols", "Podzols", "Regosols",
    "Solonchaks", "Solonetz", "Stagnosols", "Umbrisols", "Vertisols"
  ),
  stringsAsFactors = FALSE
)


# =============================================================================
# Functional Groups - Mechanistically relevant for MRT
# =============================================================================

# Assign functional group to each WRB class
# Groups based on dominant controls on soil carbon dynamics
FUNCTIONAL_GROUPS <- data.frame(
  wrb_name = c(
    # Sandy/Coarse texture - fast drainage, low physical protection
    "Arenosols", "Podzols", "Regosols",
    
    # Temperate agricultural - moderate protection, high inputs
    "Luvisols", "Cambisols", "Phaeozems", "Chernozems", 
    "Kastanozems", "Albeluvisols",
    
    # High clay tropical - strong physical protection, weathered
    "Acrisols", "Ferralsols", "Lixisols", "Plinthosols",
    "Nitisols", "Alisols", "Vertisols",
    
    # Arid/chemical constraints - low inputs, mineral stabilization
    "Calcisols", "Gypsisols", "Solonchaks", "Solonetz", "Durisols",
    
    # Poorly drained/wetland - anaerobic, slow decomposition
    "Gleysols", "Histosols", "Fluvisols", "Stagnosols", "Planosols",
    
    # Volcanic - high allophane, strong organo-mineral bonds
    "Andosols",
    
    # Shallow/young - limited profile development
    "Leptosols", "Umbrisols",
    
    # Cold/permafrost - cryoturbation, frozen C
    "Cryosols"
  ),
  functional_group = c(
    rep("sandy_coarse", 3),
    rep("temperate_agricultural", 6),
    rep("high_clay_tropical", 7),
    rep("arid_chemical", 5),
    rep("poorly_drained", 5),
    rep("volcanic", 1),
    rep("shallow_young", 2),
    rep("cold_permafrost", 1)
  ),
  stringsAsFactors = FALSE
)

# Merge with WRB legend
WRB_LEGEND <- WRB_LEGEND %>%
  left_join(FUNCTIONAL_GROUPS, by = "wrb_name") %>%
  mutate(
    functional_group = ifelse(is.na(functional_group), "other", functional_group)
  )

# Numeric codes for functional groups (for raster and modeling)
FUNC_GROUP_CODES <- data.frame(
  functional_group = c(
    "sandy_coarse",
    "temperate_agricultural", 
    "high_clay_tropical",
    "arid_chemical",
    "poorly_drained",
    "volcanic",
    "shallow_young",
    "cold_permafrost",
    "other"
  ),
  func_code = 1:9,
  stringsAsFactors = FALSE
)

WRB_LEGEND <- WRB_LEGEND %>%
  left_join(FUNC_GROUP_CODES, by = "functional_group")


# =============================================================================
# Helper Functions
# =============================================================================

#' Load or create WRB raster at 0.1° resolution
#' Uses grid-sampling approach with progress monitoring
load_wrb_raster <- function(soilclass_dir = SOILCLASS_DIR, wrb_file = WRB_FILE) {
  log_step("Loading SoilGrids WRB raster")
  
  dir.create(soilclass_dir, recursive = TRUE, showWarnings = FALSE)
  filepath <- file.path(soilclass_dir, wrb_file)
  
  if (file.exists(filepath)) {
    cat("  Loading cached WRB raster\n")
    r <- rast(filepath)
  } else {
    cat("  Creating WRB raster via grid sampling (this will take several hours)...\n")
    cat(sprintf("  Source: %s\n", SOILGRIDS_WRB_URL))
    
    # Connect to remote VRT
    wrb_remote <- tryCatch({
      rast(SOILGRIDS_WRB_URL)
    }, error = function(e) {
      stop("Failed to connect to SoilGrids WRB: ", e$message,
           "\nCheck internet connection or try again later.")
    })
    
    cat(sprintf("  Native resolution: %.6f° (~%dm)\n", 
                res(wrb_remote)[1], round(res(wrb_remote)[1] * 111000)))
    
    # Grid sampling approach - extract at 0.1° grid cell centers
    target_res <- 0.1
    
    # Generate grid cell centers
    lons <- seq(-180 + target_res/2, 180 - target_res/2, by = target_res)
    lats <- seq(-60 + target_res/2, 85 - target_res/2, by = target_res)
    
    cat(sprintf("  Target grid: %d x %d = %s cells\n", 
                length(lons), length(lats), 
                format(length(lons) * length(lats), big.mark = ",")))
    cat("  Extracting WRB values by latitude bands...\n\n")
    
    # Process by latitude bands for progress monitoring
    n_lats <- length(lats)
    band_size <- 50  # Process 50 latitude rows at a time
    n_bands <- ceiling(n_lats / band_size)
    
    all_vals <- vector("list", n_bands)
    start_time <- Sys.time()
    
    for (b in seq_len(n_bands)) {
      lat_start <- (b - 1) * band_size + 1
      lat_end <- min(b * band_size, n_lats)
      band_lats <- lats[lat_start:lat_end]
      
      # Create grid for this band
      band_grid <- expand.grid(lon = lons, lat = band_lats)
      pts <- vect(band_grid, geom = c("lon", "lat"), crs = "EPSG:4326")
      
      # Extract values
      band_vals <- terra::extract(wrb_remote, pts, ID = FALSE)[, 1]
      all_vals[[b]] <- band_vals
      
      # Progress reporting
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      pct_done <- lat_end / n_lats
      est_total <- if (pct_done > 0) elapsed / pct_done else NA
      est_remain <- if (!is.na(est_total)) est_total - elapsed else NA
      
      cat(sprintf("\r  [%3.0f%%] Band %d/%d | Lat %.1f° to %.1f° | Elapsed: %.0f min | Remaining: ~%.0f min   ",
                  100 * pct_done, b, n_bands, 
                  min(band_lats), max(band_lats),
                  elapsed,
                  ifelse(is.na(est_remain), 0, est_remain)))
      flush.console()
    }
    cat("\n\n")
    
    # Combine all values
    all_vals <- unlist(all_vals)
    
    # Create raster from extracted values
    cat("  Building raster from extracted values...\n")
    r <- rast(
      xmin = -180, xmax = 180,
      ymin = -60, ymax = 85,
      resolution = target_res,
      crs = "EPSG:4326"
    )
    
    # expand.grid returns lon-major order, reshape to raster
    # Values go: all lons for lat1, all lons for lat2, etc. (south to north)
    val_matrix <- matrix(all_vals, nrow = length(lats), ncol = length(lons), byrow = TRUE)
    
    # Raster expects values from top (north) to bottom (south), so flip
    val_matrix <- val_matrix[nrow(val_matrix):1, ]
    
    # Fill raster
    values(r) <- as.vector(t(val_matrix))
    
    # Save cached version
    writeRaster(r, filepath, overwrite = TRUE,
                datatype = "INT1U", gdal = c("COMPRESS=LZW"))
    
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("  ✓ Saved: %s (%.1f minutes total)\n", basename(filepath), total_time))
  }
  
  r_min <- global(r, "min", na.rm = TRUE)[1,1]
  r_max <- global(r, "max", na.rm = TRUE)[1,1]
  
  cat(sprintf("  File: %s\n", wrb_file))
  cat(sprintf("  Resolution: %.4f°\n", res(r)[1]))
  cat(sprintf("  Value range: %d to %d\n", r_min, r_max))
  
  return(r)
}


#' Create functional group raster from WRB raster
create_functional_group_raster <- function(wrb_raster, output_dir = SOILCLASS_DIR) {
  log_step("Creating functional group raster")
  
  func_file <- file.path(output_dir, "functional_group_0p1deg.tif")
  
  if (file.exists(func_file)) {
    cat("  Loading cached functional group raster\n")
    return(rast(func_file))
  }
  
  cat("  Reclassifying WRB to functional groups...\n")
  
  # Build reclassification matrix: from, to, becomes
  # Since WRB values are exact integers, use from=to
  rcl_matrix <- cbind(
    from = WRB_LEGEND$value,
    to = WRB_LEGEND$value,
    becomes = WRB_LEGEND$func_code
  )
  
  func_raster <- classify(wrb_raster, rcl_matrix, right = NA)
  
  # Save
  writeRaster(func_raster, func_file, overwrite = TRUE,
              datatype = "INT1U", gdal = c("COMPRESS=LZW"))
  
  cat(sprintf("  ✓ Saved: %s\n", basename(func_file)))
  
  return(func_raster)
}


#' Extract soil classification at points
extract_soilclass_at_points <- function(wrb_raster, func_raster, points_data) {
  log_step("Extracting soil classification at points")
  
  pts <- vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %s points...\n", format(nrow(points_data), big.mark = ",")))
  
  # Extract WRB class (comes as factor/character with category labels)
  wrb_vals <- terra::extract(wrb_raster, pts, ID = FALSE)[, 1]
  wrb_names <- as.character(wrb_vals)
  
  # Build results - join by name instead of numeric code
  results <- data.frame(
    dsiteid = points_data$dsiteid,
    wrb_name = wrb_names,
    stringsAsFactors = FALSE
  )
  
  # Join WRB codes and functional groups by name
  results <- results %>%
    left_join(WRB_LEGEND[, c("value", "wrb_code", "wrb_name", "functional_group", "func_code")],
              by = "wrb_name")
  
  # Rename columns
  results <- results %>%
    rename(
      soilclass_wrb_value = value,
      soilclass_wrb_code = wrb_code,
      soilclass_wrb_name = wrb_name,
      soilclass_func_group = functional_group,
      soilclass_func_code = func_code
    )
  
  # Coverage stats
  n_valid <- sum(!is.na(results$soilclass_wrb_value))
  cat(sprintf("  WRB coverage: %s (%.1f%%)\n", 
              format(n_valid, big.mark = ","),
              100 * n_valid / nrow(results)))
  
  # Distribution by functional group
  cat("\n  Distribution by functional group:\n")
  func_dist <- table(results$soilclass_func_group, useNA = "ifany")
  func_dist <- sort(func_dist, decreasing = TRUE)
  for (grp in names(func_dist)) {
    grp_label <- if (is.na(grp)) "NA" else grp
    cat(sprintf("    %s: %s (%.1f%%)\n", 
                grp_label, 
                format(func_dist[grp], big.mark = ","),
                100 * func_dist[grp] / nrow(results)))
  }
  
  return(results)
}

#' Create diagnostic outputs
create_soilclass_diagnostics <- function(data, wrb_raster, func_raster, 
                                         output_dir = "./plots/08_soilclass") {
  log_step("Creating diagnostic outputs")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- 1. WRB distribution table ---
  wrb_dist_file <- file.path(output_dir, "wrb_distribution.csv")
  
  wrb_summary <- data %>%
    group_by(soilclass_wrb_value, soilclass_wrb_code, soilclass_wrb_name) %>%
    summarise(n_sites = n(), .groups = "drop") %>%
    arrange(desc(n_sites)) %>%
    mutate(pct = round(100 * n_sites / sum(n_sites), 1))
  
  write.csv(wrb_summary, wrb_dist_file, row.names = FALSE)
  cat(sprintf("  ✓ WRB distribution: %s\n", basename(wrb_dist_file)))
  
  # --- 2. Functional group distribution ---
  func_dist_file <- file.path(output_dir, "functional_group_distribution.csv")
  
  func_summary <- data %>%
    group_by(soilclass_func_code, soilclass_func_group) %>%
    summarise(n_sites = n(), .groups = "drop") %>%
    arrange(desc(n_sites)) %>%
    mutate(pct = round(100 * n_sites / sum(n_sites), 1))
  
  write.csv(func_summary, func_dist_file, row.names = FALSE)
  cat(sprintf("  ✓ Functional group distribution: %s\n", basename(func_dist_file)))
  
  # --- 3. WRB map ---
  wrb_map_file <- file.path(output_dir, "soilclass_wrb_map.png")
  
  png(wrb_map_file, width = 1600, height = 800, res = 150)
  plot(wrb_raster,
       main = "SoilGrids WRB Classification (0.1°)",
       col = hcl.colors(30, "Spectral"),
       axes = TRUE)
  dev.off()
  cat(sprintf("  ✓ WRB map: %s\n", basename(wrb_map_file)))
  
  # --- 4. Functional group map ---
  func_map_file <- file.path(output_dir, "soilclass_functional_map.png")
  
  # Define colors for functional groups
  func_colors <- c(
    "1" = "#FFD700",  # sandy_coarse - gold
    "2" = "#228B22",  # temperate_agricultural - forest green
    "3" = "#8B0000",  # high_clay_tropical - dark red
    "4" = "#DEB887",  # arid_chemical - burlywood
    "5" = "#4169E1",  # poorly_drained - royal blue
    "6" = "#FF4500",  # volcanic - orange red
    "7" = "#A0522D",  # shallow_young - sienna
    "8" = "#E0FFFF",  # cold_permafrost - light cyan
    "9" = "#808080"   # other - gray
  )
  
  png(func_map_file, width = 1600, height = 800, res = 150)
  plot(func_raster,
       main = "Soil Functional Groups (0.1°)",
       col = func_colors,
       axes = TRUE,
       legend = FALSE)
  
  legend("bottomleft",
         legend = FUNC_GROUP_CODES$functional_group,
         fill = func_colors,
         bg = "white",
         cex = 0.6,
         ncol = 2)
  dev.off()
  cat(sprintf("  ✓ Functional group map: %s\n", basename(func_map_file)))
  
  # --- 5. Bar chart of functional groups ---
  bar_file <- file.path(output_dir, "functional_group_barplot.png")
  
  png(bar_file, width = 1000, height = 600, res = 150)
  par(mar = c(8, 4, 3, 1))
  
  func_counts <- table(data$soilclass_func_group)
  func_counts <- sort(func_counts, decreasing = TRUE)
  
  bar_colors <- func_colors[as.character(FUNC_GROUP_CODES$func_code[
    match(names(func_counts), FUNC_GROUP_CODES$functional_group)])]
  
  barplot(func_counts,
          main = "Distribution of Soil Sites by Functional Group",
          ylab = "Number of sites",
          col = bar_colors,
          las = 2,
          cex.names = 0.8)
  dev.off()
  cat(sprintf("  ✓ Functional group bar chart: %s\n", basename(bar_file)))
}


# =============================================================================
# Main Execution
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║         Step 08: Soil Classification (SoilGrids WRB)           ║\n")
cat("║         Source: ISRIC SoilGrids COG (250m → 0.1°)              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# --- Step 1: Load input data ---
log_step("Loading input data")

input_file <- file.path(PATHS$output_dir, "07_with_covariates.rds")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, 
       "\nRun 07b_extract_terrain.R first")
}

data <- readRDS(input_file)
cat(sprintf("  Loaded %s observations\n", format(nrow(data), big.mark = ",")))


# --- Step 2: Get unique site locations ---
log_step("Identifying unique locations")

sites <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees)

cat(sprintf("  %s unique site locations\n", format(nrow(sites), big.mark = ",")))


# --- Step 3: Load WRB raster ---
wrb_raster <- load_wrb_raster(SOILCLASS_DIR, WRB_FILE)


# --- Step 4: Create functional group raster ---
func_raster <- create_functional_group_raster(wrb_raster, SOILCLASS_DIR)


# --- Step 5: Extract at points ---
soilclass_extracted <- load_or_run(
  "08_soilclass_extracted.rds",
  function() {
    extract_soilclass_at_points(wrb_raster, func_raster, sites)
  }
)


# Check raster values directly
print(wrb_raster)
global(wrb_raster, "range", na.rm = TRUE)

# Sample a few points manually
test_pts <- vect(
  data.frame(lon = c(0, 10, -50, 100), lat = c(50, 45, -20, 30)),
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)
terra::extract(wrb_raster, test_pts)

# Check if raster has any non-NA values
global(wrb_raster, "notNA")


# --- Step 6: Merge with main dataset ---
log_step("Merging soil classification with main dataset")

data <- data %>%
  left_join(soilclass_extracted, by = "dsiteid")

cat(sprintf("  ✓ Added %d soil classification variables\n", 
            ncol(soilclass_extracted) - 1))


# --- Step 7: Save output ---
log_step("Saving output")

output_file <- file.path(PATHS$output_dir, "08_with_soilclass.rds")
saveRDS(data, output_file)

cat(sprintf("  ✓ Saved: %s\n", basename(output_file)))
cat(sprintf("  Dataset: %s rows, %d columns\n", 
            format(nrow(data), big.mark = ","), ncol(data)))


# --- Step 8: Save legend files for reference ---
log_step("Saving legend files")

wrb_legend_file <- file.path(SOILCLASS_DIR, "wrb_legend.csv")
write.csv(WRB_LEGEND, wrb_legend_file, row.names = FALSE)
cat(sprintf("  ✓ WRB legend: %s\n", basename(wrb_legend_file)))

func_legend_file <- file.path(SOILCLASS_DIR, "functional_group_legend.csv")
write.csv(FUNC_GROUP_CODES, func_legend_file, row.names = FALSE)
cat(sprintf("  ✓ Functional group legend: %s\n", basename(func_legend_file)))


# --- Step 9: Create diagnostics ---
create_soilclass_diagnostics(data, wrb_raster, func_raster)


# --- Summary ---
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("                  WRB extraction summary                            \n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Variables added:\n")
cat("  soilclass_wrb_value   - WRB numeric code (0-29)\n")
cat("  soilclass_wrb_code    - WRB 2-letter code (e.g., CM, LV, PZ)\n")
cat("  soilclass_wrb_name    - WRB full name (e.g., Cambisols)\n")
cat("  soilclass_func_group  - Functional group name\n")
cat("  soilclass_func_code   - Functional group numeric (1-9) for modeling\n")

cat("\nFunctional groups:\n")
for (i in 1:nrow(FUNC_GROUP_CODES)) {
  cat(sprintf("  %d: %s\n", FUNC_GROUP_CODES$func_code[i], 
              FUNC_GROUP_CODES$functional_group[i]))
}

cat("\nSpatialized rasters:\n")
cat(sprintf("  %s/%s (downloaded from ISRIC, 0.1° modal)\n", SOILCLASS_DIR, WRB_FILE))
cat(sprintf("  %s/functional_group_0p1deg.tif (reclassified)\n", SOILCLASS_DIR))
cat("\n═══════════════════════════════════════════════════════════════════\n\n")
