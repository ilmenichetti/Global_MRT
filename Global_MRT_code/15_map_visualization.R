################################################################################
# 15_visualize_MRT_maps.R
#
# Visualization of MRT prediction maps and difference maps
# Creates individual maps, multi-panel comparisons, and difference visualizations
#
# Outputs:
#   - Individual model prediction maps
#   - Multi-panel model comparison
#   - Individual difference maps (each group vs climate baseline)
#   - Multi-panel difference comparison
#
# Author: Lorenzo
# Date: 2026-01-12
################################################################################

library(terra)
library(sf)
library(rnaturalearth)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
PRED_DIR     <- file.path(PIPELINE_DIR, "outputs/MRT_predictions")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots")

# Create plot directory
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Target resolution for visualization (degrees)
VIS_RES <- 0.5

# MRT range for color scale (years)
MRT_MIN <- 0
MRT_MAX <- 200

# Difference range for color scale (years)
DIFF_RANGE <- c(-50, 50)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  VISUALIZE MRT PREDICTION MAPS (Step 15)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD HIGH-RESOLUTION WORLD BORDERS
# =============================================================================

cat("Loading country borders...\n")
sf_use_s2(FALSE)
world_borders <- ne_countries(scale = 50, returnclass = "sf")
world_borders <- st_geometry(world_borders)
cat("  ✓ Loaded Natural Earth 1:50m borders\n\n")

# =============================================================================
# FIND ALL PREDICTION MAPS
# =============================================================================

tif_files <- list.files(PRED_DIR, pattern = "^MRT_M[0-9].*\\.tif$", full.names = TRUE)
cat("Found", length(tif_files), "prediction maps:\n")
for (f in tif_files) cat("  -", basename(f), "\n")
cat("\n")

if (length(tif_files) == 0) {
  stop("No MRT prediction maps found in ", PRED_DIR)
}

# =============================================================================
# COLOR PALETTES
# =============================================================================

# For MRT predictions (sequential)
cols_mrt <- hcl.colors(100, "Spectral", rev = TRUE)

# For differences (diverging, centered at zero)
cols_diff <- hcl.colors(100, "Purple-Green")

# =============================================================================
# HELPER FUNCTION: Aggregate and trim raster
# =============================================================================

prep_raster <- function(r, agg_res = VIS_RES) {
  agg_factor <- round(agg_res / res(r)[1])
  if (agg_factor > 1) {
    r_agg <- aggregate(r, fact = agg_factor, fun = "mean", na.rm = TRUE)
  } else {
    r_agg <- r
  }
  return(r_agg)
}

# =============================================================================
# CREATE VISUALIZATION FOR EACH PREDICTION MAP
# =============================================================================

cat("Creating individual prediction maps...\n\n")

for (tif_file in tif_files) {
  
  model_name <- gsub("\\.tif$", "", basename(tif_file))
  cat("Processing", model_name, "...\n")
  
  r <- rast(tif_file)
  vals <- values(r, na.rm = TRUE)
  
  cat("  Value range:", round(min(vals), 1), "-", round(max(vals), 1), "years\n")
  
  r_agg <- prep_raster(r)
  r_vis <- clamp(r_agg, lower = MRT_MIN, upper = MRT_MAX)
  
  out_png <- file.path(PLOT_DIR, paste0(model_name, "_map.png"))
  
  png(out_png, width = 1600, height = 800, res = 120)
  
  plot(r_vis,
       main = paste0(model_name, " - Mean Residence Time (years)"),
       col = cols_mrt,
       range = c(MRT_MIN, MRT_MAX),
       mar = c(2, 2, 2, 4),
       axes = TRUE,
       plg = list(title = "MRT\n(years)", title.cex = 0.8))
  
  plot(world_borders, add = TRUE, col = NA, border = "grey40", lwd = 0.3)
  
  dev.off()
  
  cat("  ✓ Saved:", basename(out_png), "\n\n")
}

# =============================================================================
# CREATE MULTI-PANEL MODEL COMPARISON
# =============================================================================

cat("Creating multi-panel model comparison...\n")

# Create mask from M2 (to apply to M1 for consistent extent)
m2_file <- file.path(PRED_DIR, "MRT_M2_climate_edaphic.tif")
if (file.exists(m2_file)) {
  m2_rast <- rast(m2_file)
  m2_agg <- prep_raster(m2_rast)
  m1_mask <- !is.na(m2_agg)
  cat("  Created mask from M2\n")
} else {
  m1_mask <- NULL
}

png(file.path(PLOT_DIR, "MRT_all_models_comparison.png"), 
    width = 2400, height = 1200, res = 120)

n_models <- length(tif_files)
n_cols <- 4
n_rows <- ceiling(n_models / n_cols)

par(mfrow = c(n_rows, n_cols), mar = c(0.2, 0.2, 1.5, 0.2), oma = c(0, 0, 1, 0))

for (tif_file in tif_files) {
  
  model_name <- gsub("^MRT_", "", gsub("\\.tif$", "", basename(tif_file)))
  
  r <- rast(tif_file)
  r_agg <- prep_raster(r)
  
  # Apply mask to M1
  if (!is.null(m1_mask) && grepl("M1_climate$", model_name)) {
    r_agg <- mask(r_agg, m1_mask, maskvalues = 0)
  }
  
  r_vis <- clamp(r_agg, lower = MRT_MIN, upper = MRT_MAX)
  
  plot(r_vis,
       main = model_name,
       col = cols_mrt,
       range = c(MRT_MIN, MRT_MAX),
       axes = FALSE,
       legend = FALSE,
       mar = c(0.1, 0.1, 1, 0.1))
  
  plot(world_borders, add = TRUE, col = NA, border = "grey40", lwd = 0.3)
}

mtext("MRT Predictions - Model Comparison", outer = TRUE, cex = 1, font = 2, line = -0.5)

dev.off()

cat("✓ Saved: MRT_all_models_comparison.png\n\n")

# =============================================================================
# FIND ALL DIFFERENCE MAPS
# =============================================================================

diff_files <- list.files(PRED_DIR, pattern = "^MRT_difference_M[2-7]_M1.*\\.tif$", full.names = TRUE)
cat("Found", length(diff_files), "difference maps:\n")
for (f in diff_files) cat("  -", basename(f), "\n")
cat("\n")

# =============================================================================
# CREATE INDIVIDUAL DIFFERENCE MAP VISUALIZATIONS
# =============================================================================

if (length(diff_files) > 0) {
  
  cat("Creating individual difference maps...\n\n")
  
  # Labels for each difference type
  diff_labels <- list(
    "edaphic" = "Edaphic Contribution (M2 - M1)",
    "landuse" = "LandUse Contribution (M3 - M1)",
    "biological" = "Biological Contribution (M4 - M1)",
    "edaphic_landuse" = "Edaphic + LandUse (M5 - M1)",
    "edaphic_bio" = "Edaphic + Biological (M6 - M1)",
    "full" = "Full Model (M7 - M1)"
  )
  
  for (diff_file in diff_files) {
    
    file_name <- basename(diff_file)
    cat("Processing", file_name, "...\n")
    
    # Determine label
    label <- "Difference"
    for (key in names(diff_labels)) {
      if (grepl(key, file_name, ignore.case = TRUE)) {
        label <- diff_labels[[key]]
        break
      }
    }
    
    r <- rast(diff_file)
    vals <- values(r, na.rm = TRUE)
    
    cat("  Value range:", round(min(vals), 1), "to", round(max(vals), 1), "years\n")
    cat("  Mean:", round(mean(vals), 2), "years\n")
    
    r_agg <- prep_raster(r)
    r_vis <- clamp(r_agg, lower = DIFF_RANGE[1], upper = DIFF_RANGE[2])
    
    out_name <- gsub("\\.tif$", ".png", file_name)
    out_png <- file.path(PLOT_DIR, out_name)
    
    png(out_png, width = 1600, height = 800, res = 120)
    
    plot(r_vis,
         main = paste0("MRT Difference: ", label),
         col = cols_diff,
         range = DIFF_RANGE,
         mar = c(2, 2, 2, 4),
         axes = TRUE,
         plg = list(title = "Δ MRT\n(years)", title.cex = 0.8))
    
    plot(world_borders, add = TRUE, col = NA, border = "grey40", lwd = 0.3)
    
    # Add zero line note
    mtext("Positive = group increases MRT | Negative = group decreases MRT",
          side = 1, line = 0.5, cex = 0.7, col = "gray40")
    
    dev.off()
    
    cat("  ✓ Saved:", out_name, "\n\n")
  }
  
  # ===========================================================================
  # CREATE MULTI-PANEL DIFFERENCE COMPARISON
  # ===========================================================================
  
  cat("Creating multi-panel difference comparison...\n")
  
  # Order files logically (M2, M3, M4, M5, M6, M7)
  diff_order <- c("M2_M1", "M3_M1", "M4_M1", "M5_M1", "M6_M1", "M7_M1")
  diff_files_ordered <- character(0)
  for (pattern in diff_order) {
    match_file <- diff_files[grep(pattern, diff_files)]
    if (length(match_file) > 0) {
      diff_files_ordered <- c(diff_files_ordered, match_file[1])
    }
  }
  
  if (length(diff_files_ordered) > 0) {
    
    png(file.path(PLOT_DIR, "MRT_all_differences_comparison.png"), 
        width = 2400, height = 1200, res = 120)
    
    n_diffs <- length(diff_files_ordered)
    n_cols <- 3
    n_rows <- ceiling(n_diffs / n_cols)
    
    par(mfrow = c(n_rows, n_cols), mar = c(0.2, 0.2, 1.5, 0.2), oma = c(0, 0, 1, 0))
    
    short_labels <- c(
      "edaphic" = "Edaphic",
      "landuse" = "LandUse", 
      "biological" = "Biological",
      "edaphic_landuse" = "Edaphic+LandUse",
      "edaphic_bio" = "Edaphic+Bio",
      "full" = "Full Model"
    )
    
    for (diff_file in diff_files_ordered) {
      
      file_name <- basename(diff_file)
      
      # Get short label
      label <- "Difference"
      for (key in names(short_labels)) {
        if (grepl(key, file_name, ignore.case = TRUE)) {
          label <- short_labels[[key]]
          break
        }
      }
      
      r <- rast(diff_file)
      r_agg <- prep_raster(r)
      r_vis <- clamp(r_agg, lower = DIFF_RANGE[1], upper = DIFF_RANGE[2])
      
      plot(r_vis,
           main = label,
           col = cols_diff,
           range = DIFF_RANGE,
           axes = FALSE,
           legend = FALSE,
           mar = c(0.1, 0.1, 1, 0.1))
      
      plot(world_borders, add = TRUE, col = NA, border = "grey40", lwd = 0.3)
    }
    
    mtext("MRT Differences vs Climate Baseline (M1)", outer = TRUE, cex = 1, font = 2, line = -0.5)
    
    dev.off()
    
    cat("✓ Saved: MRT_all_differences_comparison.png\n\n")
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  VISUALIZATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Plots saved to:", PLOT_DIR, "\n")
cat("\nFiles created:\n")
plot_files <- list.files(PLOT_DIR, pattern = "MRT.*\\.png$")
for (f in plot_files) {
  size_kb <- round(file.size(file.path(PLOT_DIR, f)) / 1024)
  cat(sprintf("  %-50s %4d KB\n", f, size_kb))
}

cat("\nSettings used:\n")
cat("  MRT range: 0 -", MRT_MAX, "years\n")
cat("  Difference range:", DIFF_RANGE[1], "to", DIFF_RANGE[2], "years\n")
cat("  Visualization resolution:", VIS_RES, "°\n")