# =============================================================================
# Step 6: Extract SoilGrids data at point locations
# Variables: clay, sand, silt, bdod, cfvo, cec, phh2o, nitrogen (0-20cm)
# =============================================================================

library(terra)
library(dplyr)

# Source config and utilities
source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# --- Configuration ---
SOILGRIDS_DIR <- "./Global_MRT_code/spatialized_layers/soilgrids"

SOILGRIDS_VARS <- c(
  "clay",      # Clay content (%)
  "sand",      # Sand content (%)
  "silt",      # Silt content (%)
  "bdod",      # Bulk density (kg/dm³ × 100, i.e., cg/cm³)
  "cfvo",      # Coarse fragments (%)
  "cec",       # Cation exchange capacity (cmol(c)/kg × 10)
  "phh2o",     # pH in water (pH × 10)
  "nitrogen"   # Total nitrogen (g/kg × 100, i.e., cg/kg)
)

# Scale factors to convert to standard units
# SoilGrids stores values as integers, need to convert
SOILGRIDS_SCALE <- list(
  clay = 0.1,       # g/kg to %
  sand = 0.1,       # g/kg to %
  silt = 0.1,       # g/kg to %
  bdod = 0.01,      # cg/cm³ to kg/dm³
  cfvo = 0.1,       # cm³/dm³ to %
  cec = 0.1,        # mmol(c)/kg to cmol(c)/kg
  phh2o = 0.1,      # pH × 10 to pH
  nitrogen = 0.01   # cg/kg to g/kg
)

# =============================================================================
# Helper Functions
# =============================================================================

#' Load all SoilGrids rasters
#' @param sg_dir Directory containing the global mosaics
#' @param variables Vector of variable names
#' @return Named list of SpatRasters
load_soilgrids_rasters <- function(sg_dir = SOILGRIDS_DIR, variables = SOILGRIDS_VARS) {
  log_step("Loading SoilGrids rasters")
  
  rasters <- list()
  
  for (var in variables) {
    # Try different naming patterns
    patterns <- c(
      sprintf("%s_0_20cm_global.tif", var),
      sprintf("%s_0_20cm_0.1deg.tif", var),
      sprintf("%s_0_20cm.tif", var),
      sprintf("SG_%s_0_20cm_global.tif", var)
    )
    
    found <- FALSE
    for (pattern in patterns) {
      filepath <- file.path(sg_dir, pattern)
      if (file.exists(filepath)) {
        rasters[[var]] <- rast(filepath)
        cat(sprintf("  ✓ %s: %s\n", var, pattern))
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      warning(sprintf("  ✗ %s: not found", var))
    }
  }
  
  if (length(rasters) == 0) {
    files <- list.files(sg_dir, pattern = "\\.tif$")
    stop("No SoilGrids files found. Files in directory: ", paste(files, collapse = ", "))
  }
  
  cat(sprintf("\n  Loaded %d of %d variables\n", length(rasters), length(variables)))
  
  return(rasters)
}


#' Extract SoilGrids values at point locations
#' @param rasters Named list of SpatRasters
#' @param points_data Data frame with coordinates
#' @return Data frame with extracted values
extract_soilgrids_at_points <- function(rasters, points_data) {
  log_step("Extracting SoilGrids at point locations")
  
  # Create point geometry
  pts <- terra::vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %d points...\n", nrow(points_data)))
  
  # Initialize results
  results <- data.frame(dsiteid = points_data$dsiteid)
  
  for (var in names(rasters)) {
    cat(sprintf("  Extracting: %s\n", var))
    
    vals <- terra::extract(rasters[[var]], pts, ID = FALSE)[, 1]
    
    # Clean NoData
    vals[vals <= -9999] <- NA
    vals[vals < 0] <- NA  # SoilGrids shouldn't have negative values
    
    # Apply scale factor
    if (var %in% names(SOILGRIDS_SCALE)) {
      vals <- vals * SOILGRIDS_SCALE[[var]]
    }
    
    # Store with sg_ prefix to identify source
    col_name <- paste0("sg_", var)
    results[[col_name]] <- vals
    
    # Print coverage
    n_valid <- sum(!is.na(vals))
    cat(sprintf("    Coverage: %d (%.1f%%)\n", n_valid, 100 * n_valid / length(vals)))
  }
  
  return(results)
}


# =============================================================================
# Main Execution
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║           Step 6: SoilGrids Extraction (0-20cm)                ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# --- Step 1: Load input data ---
log_step("Loading input data")

input_file <- file.path(PATHS$output_dir, "05_with_disturbance.rds")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, 
       "\nRun 05_extract_disturbance.R first")
}

data <- readRDS(input_file)
cat(sprintf("  Loaded %d observations\n", nrow(data)))


# --- Step 2: Get unique site locations ---
log_step("Identifying unique locations")

sites <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees)

cat(sprintf("  %d unique site locations\n", nrow(sites)))


# --- Step 3: Load SoilGrids rasters ---
sg_rasters <- load_soilgrids_rasters(SOILGRIDS_DIR, SOILGRIDS_VARS)


# --- Step 4: Extract at points ---
sg_extracted <- load_or_run(
  "06_soilgrids_extracted.rds",
  function() {
    extract_soilgrids_at_points(sg_rasters, sites)
  }
)


# --- Step 5: Merge with main dataset ---
log_step("Merging SoilGrids data with main dataset")

data <- data %>%
  left_join(sg_extracted, by = "dsiteid")

cat(sprintf("  ✓ Added %d SoilGrids variables\n", ncol(sg_extracted) - 1))


# --- Step 6: Clean values ---
log_step("Cleaning SoilGrids values")

sg_cols <- grep("^sg_", names(data), value = TRUE)

for (col in sg_cols) {
  n_negative <- sum(data[[col]] < 0, na.rm = TRUE)
  if (n_negative > 0) {
    data[[col]][data[[col]] < 0] <- NA
    cat(sprintf("  Set %d negative values to NA in %s\n", n_negative, col))
  }
}

# Sanity checks
cat("\n  Sanity checks:\n")
if ("sg_clay" %in% names(data)) {
  cat(sprintf("  Clay range: %.1f - %.1f %%\n", 
              min(data$sg_clay, na.rm = TRUE), max(data$sg_clay, na.rm = TRUE)))
}
if ("sg_phh2o" %in% names(data)) {
  cat(sprintf("  pH range: %.1f - %.1f\n", 
              min(data$sg_phh2o, na.rm = TRUE), max(data$sg_phh2o, na.rm = TRUE)))
}
if ("sg_bdod" %in% names(data)) {
  cat(sprintf("  Bulk density range: %.2f - %.2f kg/dm³\n", 
              min(data$sg_bdod, na.rm = TRUE), max(data$sg_bdod, na.rm = TRUE)))
}


# --- Step 7: Save output ---
log_step("Saving output")

output_file <- file.path(PATHS$output_dir, "06_with_soilgrids.rds")
saveRDS(data, output_file)
cat(sprintf("  ✓ Saved: %s\n", output_file))


# =============================================================================
# Diagnostics
# =============================================================================

log_step("Generating diagnostics")

diag_dir <- "./Global_MRT_code/plots/06_extract_soilgrids"
if (!dir.exists(diag_dir)) {
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
}


# --- 1. Data coverage table ---
log_step("Creating data coverage table")

sg_vars <- grep("^sg_", names(data), value = TRUE)

coverage_df <- data.frame(
  variable = sg_vars,
  n_total = nrow(data),
  n_valid = sapply(sg_vars, function(v) sum(!is.na(data[[v]]))),
  pct_valid = sapply(sg_vars, function(v) round(100 * mean(!is.na(data[[v]])), 2)),
  mean = sapply(sg_vars, function(v) round(mean(data[[v]], na.rm = TRUE), 4)),
  sd = sapply(sg_vars, function(v) round(sd(data[[v]], na.rm = TRUE), 4)),
  min = sapply(sg_vars, function(v) round(min(data[[v]], na.rm = TRUE), 4)),
  max = sapply(sg_vars, function(v) round(max(data[[v]], na.rm = TRUE), 4))
)

write.csv(coverage_df, file.path(diag_dir, "soilgrids_data_coverage.csv"), row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "soilgrids_data_coverage.csv")))


# --- 2. World maps ---
log_step("Creating world maps")

sites_for_plot <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees,
           sg_clay, sg_phh2o, sg_cec, sg_nitrogen)

# Map: Clay content
png(file.path(diag_dir, "map_clay.png"), 
    width = 2400, height = 1200, res = 150)

par(mar = c(2, 2, 3, 1))
plot(c(-180, 180), c(-90, 90), type = "n", 
     xlab = "", ylab = "", asp = 1,
     main = "SoilGrids Clay Content (0-20cm)")

if (requireNamespace("maps", quietly = TRUE)) {
  maps::map("world", add = TRUE, col = "gray90", fill = TRUE, border = "gray70")
}

valid_sites <- sites_for_plot[!is.na(sites_for_plot$sg_clay), ]
if (nrow(valid_sites) > 0) {
  clay_colors <- colorRampPalette(c("sandybrown", "chocolate", "saddlebrown"))(100)
  clay_breaks <- seq(0, 60, length.out = 101)
  clay_idx <- cut(valid_sites$sg_clay, breaks = clay_breaks, labels = FALSE, include.lowest = TRUE)
  clay_idx[is.na(clay_idx)] <- 100
  
  points(valid_sites$longitude_decimal_degrees, 
         valid_sites$latitude_decimal_degrees,
         pch = 16, cex = 0.4, col = clay_colors[clay_idx])
}

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "map_clay.png")))


# Map: pH
png(file.path(diag_dir, "map_ph.png"), 
    width = 2400, height = 1200, res = 150)

par(mar = c(2, 2, 3, 1))
plot(c(-180, 180), c(-90, 90), type = "n", 
     xlab = "", ylab = "", asp = 1,
     main = "SoilGrids pH (0-20cm)")

if (requireNamespace("maps", quietly = TRUE)) {
  maps::map("world", add = TRUE, col = "gray90", fill = TRUE, border = "gray70")
}

valid_sites <- sites_for_plot[!is.na(sites_for_plot$sg_phh2o), ]
if (nrow(valid_sites) > 0) {
  ph_colors <- colorRampPalette(c("red", "yellow", "green", "blue"))(100)
  ph_breaks <- seq(4, 9, length.out = 101)
  ph_idx <- cut(valid_sites$sg_phh2o, breaks = ph_breaks, labels = FALSE, include.lowest = TRUE)
  ph_idx[is.na(ph_idx)] <- 50
  
  points(valid_sites$longitude_decimal_degrees, 
         valid_sites$latitude_decimal_degrees,
         pch = 16, cex = 0.4, col = ph_colors[ph_idx])
}

legend("bottomleft",
       legend = c("pH 4 (acidic)", "pH 6.5 (neutral)", "pH 9 (basic)"),
       pch = 16, col = c("red", "yellow", "blue"),
       bg = "white", cex = 0.9)

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "map_ph.png")))


# --- 3. Histograms ---
log_step("Creating histograms")

png(file.path(diag_dir, "hist_texture.png"), 
    width = 1000, height = 500, res = 150)

par(mfrow = c(1, 2), mar = c(5, 4, 3, 2))

# Clay
hist(data$sg_clay, breaks = 50, main = "Clay Content",
     xlab = "Clay (%)", col = "chocolate", border = "white")

# Sand  
hist(data$sg_sand, breaks = 50, main = "Sand Content",
     xlab = "Sand (%)", col = "sandybrown", border = "white")

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "hist_texture.png")))


png(file.path(diag_dir, "hist_chemistry.png"), 
    width = 1500, height = 500, res = 150)

par(mfrow = c(1, 3), mar = c(5, 4, 3, 2))

# pH
hist(data$sg_phh2o, breaks = 50, main = "pH",
     xlab = "pH", col = "steelblue", border = "white")

# CEC
hist(data$sg_cec, breaks = 50, main = "Cation Exchange Capacity",
     xlab = "CEC (cmol(c)/kg)", col = "darkgreen", border = "white")

# Nitrogen
hist(data$sg_nitrogen, breaks = 50, main = "Total Nitrogen",
     xlab = "N (g/kg)", col = "purple", border = "white")

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "hist_chemistry.png")))


# --- Summary ---
log_step("SoilGrids extraction complete")
