# =============================================================================
# 09_extract_landcover.R
# Extract ESA WorldCover land cover fractions
# =============================================================================

library(terra)
library(dplyr)
library(geodata)
library(ggplot2)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

INPUT_FILE <- "./Global_MRT_code/outputs/08_with_soilclass.rds"
OUTPUT_DIR <- "./Global_MRT_code/outputs"
SPATIAL_DIR <- "./Global_MRT_code/spatialized_layers/landcover"
PLOT_DIR <- "./Global_MRT_code/plots/09_landcover"
TEMP_DIR <- "./Global_MRT_code/temp_data/landcover"

ESA_LC_CLASSES <- c("trees", "grassland", "shrubs", "cropland", "built",
                    "bare", "snow", "water", "wetland", "mangroves", "moss")

# Simplified groups for RF modeling
LC_GROUPS <- list(
  forest = "trees",
  grassland = "grassland", 
  shrubland = "shrubs",
  cropland = "cropland",
  wetland = c("wetland", "mangroves"),
  sparse_barren = c("bare", "snow", "moss"),
  other = c("built", "water")
)

for (d in c(OUTPUT_DIR, SPATIAL_DIR, PLOT_DIR, TEMP_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

load_or_run <- function(filename, fun, dir = OUTPUT_DIR) {
  path <- file.path(dir, filename)
  if (file.exists(path)) {
    log_step(paste("Loading cached:", filename))
    return(readRDS(path))
  }
  log_step(paste("Computing:", filename))
  result <- fun()
  saveRDS(result, path)
  cat(sprintf("  Saved: %s\n", path))
  return(result)
}

# -----------------------------------------------------------------------------
# Load landcover rasters (download if needed)
# -----------------------------------------------------------------------------

load_landcover_rasters <- function() {
  log_step("Loading ESA WorldCover layers")
  
  rasters <- list()
  for (lc in ESA_LC_CLASSES) {
    cat(sprintf("  Loading %s...\n", lc))
    rasters[[lc]] <- geodata::landcover(var = lc, path = TEMP_DIR)
  }
  
  # Stack all layers
  lc_stack <- rast(rasters)
  cat(sprintf("  Loaded %d layers at %.4f° resolution\n", 
              nlyr(lc_stack), res(lc_stack)[1]))
  
  return(lc_stack)
}

# -----------------------------------------------------------------------------
# Create spatialized 0.1° rasters
# -----------------------------------------------------------------------------

create_landcover_rasters <- function(lc_stack) {
  log_step("Creating 0.1° spatialized rasters")
  
  # Check for cached
  output_path <- file.path(SPATIAL_DIR, "landcover_fractions_0p1deg.tif")
  if (file.exists(output_path)) {
    log_step("Loading cached spatialized rasters")
    return(rast(output_path))
  }
  
  # Create target grid (0.1°, same extent as other layers)
  target <- rast(
    extent = ext(-180, 180, -60, 85),
    resolution = 0.1,
    crs = "EPSG:4326"
  )
  
  # Resample each layer (mean aggregation for fractions)
  cat("  Resampling to 0.1° (this may take a while)...\n")
  
  resampled <- list()
  for (i in seq_along(ESA_LC_CLASSES)) {
    lc <- ESA_LC_CLASSES[i]
    cat(sprintf("    [%d/%d] %s\n", i, length(ESA_LC_CLASSES), lc))
    resampled[[lc]] <- resample(lc_stack[[lc]], target, method = "average")
  }
  
  # Stack and save
  lc_0p1 <- rast(resampled)
  names(lc_0p1) <- ESA_LC_CLASSES
  
  writeRaster(lc_0p1, output_path, 
              overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=LZW"))
  cat(sprintf("  Saved: %s\n", output_path))
  
  return(lc_0p1)
}

# -----------------------------------------------------------------------------
# Create dominant class raster
# -----------------------------------------------------------------------------

create_dominant_raster <- function(lc_rasters) {
  log_step("Creating dominant land cover raster")
  
  output_path <- file.path(SPATIAL_DIR, "landcover_dominant_0p1deg.tif")
  if (file.exists(output_path)) {
    return(rast(output_path))
  }
  
  # Find dominant class (highest fraction) per cell
  dominant <- which.max(lc_rasters)
  names(dominant) <- "dominant_lc"
  
  writeRaster(dominant, output_path,
              overwrite = TRUE,
              datatype = "INT1U",
              gdal = c("COMPRESS=LZW"))
  cat(sprintf("  Saved: %s\n", output_path))
  
  return(dominant)
}

# -----------------------------------------------------------------------------
# Extract at points
# -----------------------------------------------------------------------------

extract_landcover_at_points <- function(lc_rasters, points_data) {
  log_step("Extracting land cover at points")
  
  pts <- vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %s points...\n", 
              format(nrow(points_data), big.mark = ",")))
  
  # Extract all fractions
  lc_vals <- terra::extract(lc_rasters, pts, ID = FALSE)
  
  # Build results
  results <- data.frame(dsiteid = points_data$dsiteid)
  
  # Add individual fractions with prefix
  for (lc in ESA_LC_CLASSES) {
    results[[paste0("lc_", lc)]] <- lc_vals[[lc]]
  }
  
  # Dominant class
  fraction_matrix <- as.matrix(lc_vals[, ESA_LC_CLASSES])
  
  get_dominant <- function(x) {
    if (all(is.na(x))) return(NA_character_)
    ESA_LC_CLASSES[which.max(x)]
  }
  
  results$lc_dominant <- apply(fraction_matrix, 1, get_dominant)
  results$lc_dominant_fraction <- apply(fraction_matrix, 1, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
  
  # Grouped fractions for RF
  for (grp_name in names(LC_GROUPS)) {
    classes <- LC_GROUPS[[grp_name]]
    if (length(classes) == 1) {
      results[[paste0("lc_grp_", grp_name)]] <- lc_vals[[classes]]
    } else {
      results[[paste0("lc_grp_", grp_name)]] <- rowSums(lc_vals[, classes], na.rm = TRUE)
    }
  }
  
  # Coverage stats
  n_valid <- sum(!is.na(results$lc_dominant))
  cat(sprintf("  Land cover coverage: %s (%.1f%%)\n",
              format(n_valid, big.mark = ","),
              100 * n_valid / nrow(results)))
  
  # Distribution
  cat("\n  Dominant land cover distribution:\n")
  lc_dist <- sort(table(results$lc_dominant, useNA = "ifany"), decreasing = TRUE)
  for (lc in names(lc_dist)) {
    lc_label <- if (is.na(lc)) "NA" else lc
    cat(sprintf("    %s: %s (%.1f%%)\n",
                lc_label,
                format(lc_dist[lc], big.mark = ","),
                100 * lc_dist[lc] / nrow(results)))
  }
  
  return(results)
}

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

create_landcover_diagnostics <- function(lc_rasters, dominant_raster, sites) {
  log_step("Creating land cover diagnostics")
  
  # --- Distribution CSV ---
  lc_dist <- data.frame(
    landcover = ESA_LC_CLASSES,
    n_dominant = sapply(ESA_LC_CLASSES, function(lc) sum(sites$lc_dominant == lc, na.rm = TRUE)),
    mean_fraction = sapply(ESA_LC_CLASSES, function(lc) mean(sites[[paste0("lc_", lc)]], na.rm = TRUE))
  ) %>%
    mutate(pct_dominant = 100 * n_dominant / sum(n_dominant, na.rm = TRUE)) %>%
    arrange(desc(n_dominant))
  
  write.csv(lc_dist, file.path(PLOT_DIR, "landcover_distribution.csv"), row.names = FALSE)
  
  # --- Dominant land cover map ---
  cat("  Creating dominant land cover map...\n")
  
  lc_colors <- c(
    "trees" = "#228B22",
    "grassland" = "#90EE90", 
    "shrubs" = "#DEB887",
    "cropland" = "#FFD700",
    "built" = "#DC143C",
    "bare" = "#D2B48C",
    "snow" = "#FFFFFF",
    "water" = "#4169E1",
    "wetland" = "#20B2AA",
    "mangroves" = "#006400",
    "moss" = "#9ACD32"
  )
  
  png(file.path(PLOT_DIR, "landcover_dominant_map.png"),
      width = 3600, height = 1800, res = 200)
  
  plot(dominant_raster, 
       col = lc_colors[ESA_LC_CLASSES],
       main = "ESA WorldCover - Dominant Land Cover (0.1°)",
       legend = FALSE,
       axes = TRUE)
  
  legend("bottomleft", 
         legend = ESA_LC_CLASSES,
         fill = lc_colors[ESA_LC_CLASSES],
         ncol = 3,
         cex = 0.8,
         bty = "n")
  
  dev.off()
  
  # --- Tree cover map ---
  cat("  Creating tree cover map...\n")
  
  png(file.path(PLOT_DIR, "landcover_trees_map.png"),
      width = 3600, height = 1800, res = 200)
  
  plot(lc_rasters[["trees"]], 
       col = colorRampPalette(c("#F5F5DC", "#228B22", "#006400"))(100),
       main = "Tree Cover Fraction",
       axes = TRUE)
  
  dev.off()
  
  # --- Cropland map ---
  cat("  Creating cropland map...\n")
  
  png(file.path(PLOT_DIR, "landcover_cropland_map.png"),
      width = 3600, height = 1800, res = 200)
  
  plot(lc_rasters[["cropland"]], 
       col = colorRampPalette(c("#F5F5DC", "#FFD700", "#FF8C00"))(100),
       main = "Cropland Fraction",
       axes = TRUE)
  
  dev.off()
  
  # --- Bar plot ---
  cat("  Creating distribution bar plot...\n")
  
  p <- ggplot(lc_dist, aes(x = reorder(landcover, -n_dominant), y = n_dominant, fill = landcover)) +
    geom_col() +
    scale_fill_manual(values = lc_colors, guide = "none") +
    labs(
      title = "Site Distribution by Dominant Land Cover",
      subtitle = sprintf("N = %s sites with land cover data", format(sum(lc_dist$n_dominant), big.mark = ",")),
      x = "Land Cover Class",
      y = "Number of Sites"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(PLOT_DIR, "landcover_barplot.png"), p, width = 10, height = 6, dpi = 200)
  
  cat(sprintf("  Diagnostics saved to: %s\n", PLOT_DIR))
}

# =============================================================================
# Main execution
# =============================================================================

log_step("Starting land cover extraction (Step 09)")

# Load input data
sites <- readRDS(INPUT_FILE)
cat(sprintf("  Loaded %s sites\n", format(nrow(sites), big.mark = ",")))

# Load/download landcover rasters (~5GB total download)
lc_native <- load_landcover_rasters()

# Create spatialized 0.1° rasters
lc_rasters <- create_landcover_rasters(lc_native)

# Create dominant class raster
dominant_raster <- create_dominant_raster(lc_rasters)

# Extract at points
lc_extracted <- load_or_run(
  "09_landcover_extracted.rds",
  function() extract_landcover_at_points(lc_rasters, sites)
)

# Merge with sites
sites <- sites %>%
  left_join(lc_extracted, by = "dsiteid")

# Save
log_step("Saving final dataset")
saveRDS(sites, file.path(OUTPUT_DIR, "09_with_landcover.rds"))
cat(sprintf("  Saved: %s observations, %s variables\n",
            format(nrow(sites), big.mark = ","), ncol(sites)))

# Diagnostics
create_landcover_diagnostics(lc_rasters, dominant_raster, sites)


# --- Summary ---
log_step("Land cover extraction complete!")
cat("\n  New variables added:\n")
cat("    - lc_trees, lc_grassland, ... (11 fraction columns)\n")
cat("    - lc_dominant (character: dominant class)\n")
cat("    - lc_dominant_fraction (0-1: fraction of dominant)\n")
cat("    - lc_grp_forest, lc_grp_cropland, ... (7 grouped fractions)\n")

