################################################################################
# 16_process_dominance_RGB.R
#
# Spatial magnitude of non-climate prediction adjustments.
#
# For each non-climate process group, compute the absolute prediction shift
# relative to the climate-only baseline:
#   |M2 - M1|  Edaphic
#   |M3 - M1|  LandUse
#   |M4 - M1|  Biological
#
# Three output figures (pick whichever works best for Fig. 5):
#
#   1) MRT_process_dominance_combined.png
#      Two panels: climate-dominance map (A) + CMY composite (B).
#
#   2) MRT_process_contribution_CMY.png
#      Single panel: CMY subtractive composite.
#      White = no modulation / climate dominates
#      Cyan / Magenta / Yellow = Edaphic / LandUse / Biological
#      Mixed hues = co-modulation; black = all three strong
#
#   3) MRT_process_contribution_RGB.png
#      Single panel: RGB additive composite.
#      Black = no modulation / climate dominates
#      Red / Green / Blue = Edaphic / LandUse / Biological
#      Mixed hues = co-modulation; white = all three strong
#
# IMPORTANT - what these figures show and do NOT show:
#   These maps display the MAGNITUDE of per-pixel prediction ADJUSTMENT that
#   each process group introduces relative to the climate baseline. They do
#   NOT quantify statistical "dominance" or variance explained - for that
#   see permutation importance (Fig. 4, Table 1). A large |Mi - M1| means
#   adding that process group shifts the prediction substantially at that
#   pixel; it does not mean that shift is more informative than others.
#
# Input:  ./Global_MRT_code/outputs/MRT_predictions/*.tif
# Output: ./Global_MRT_code/plots/MRT_process_dominance_combined.png
#         ./Global_MRT_code/plots/MRT_process_contribution_CMY.png
#         ./Global_MRT_code/plots/MRT_process_contribution_RGB.png
#         ./Global_MRT_code/outputs/MRT_predictions/process_dominance_RGB.tif
#
# Author: Lorenzo
# Date:   2026-01-12  (rewrite 2026-04)
################################################################################

library(terra)
library(dplyr)
library(sf)
library(rnaturalearth)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
PRED_DIR     <- file.path(PIPELINE_DIR, "outputs/MRT_predictions")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots")

VIS_RES         <- 0.5    # aggregation resolution (degrees) for display
LUMINANCE_BLEND <- 0.6    # strength of luminance modulation

cat("\u2550\u2550\u2550 PROCESS-MAGNITUDE COMPOSITES \u2550\u2550\u2550\n\n")

# =============================================================================
# LOAD PREDICTION MAPS
# =============================================================================

cat("Loading prediction maps...\n")

m1_file <- file.path(PRED_DIR, "MRT_M1_climate.tif")
m2_file <- file.path(PRED_DIR, "MRT_M2_climate_edaphic.tif")
m3_file <- file.path(PRED_DIR, "MRT_M3_climate_landuse.tif")
m4_file <- file.path(PRED_DIR, "MRT_M4_climate_biological.tif")
m7_file <- file.path(PRED_DIR, "MRT_M7_full.tif")

required_files <- c(m1_file, m2_file, m3_file, m4_file)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing required files:\n  ", paste(missing, collapse = "\n  "))
}

m1 <- rast(m1_file)
m2 <- rast(m2_file)
m3 <- rast(m3_file)
m4 <- rast(m4_file)

cat("  \u2713 Loaded M1, M2, M3, M4\n\n")

# =============================================================================
# PROCESS-GROUP PREDICTION ADJUSTMENTS
# =============================================================================

cat("Calculating prediction adjustments (|Mi - M1|)...\n")

edaphic_mag    <- abs(m2 - m1)
landuse_mag    <- abs(m3 - m1)
biological_mag <- abs(m4 - m1)

summarise_mag <- function(r, lbl) {
  m  <- as.numeric(global(r, "mean", na.rm = TRUE))
  rg <- as.numeric(unlist(global(r, "range", na.rm = TRUE)))
  cat(sprintf("  %-10s  mean = %5.2f  range = [%.1f, %.1f]\n",
              lbl, m, rg[1], rg[2]))
}
summarise_mag(edaphic_mag,    "Edaphic")
summarise_mag(landuse_mag,    "LandUse")
summarise_mag(biological_mag, "Biological")
cat("\n")

# =============================================================================
# NORMALISE EACH CHANNEL TO 0-255
# =============================================================================

cat("Normalising channels...\n")

all_vals <- c(values(edaphic_mag),
              values(landuse_mag),
              values(biological_mag))
max_val  <- as.numeric(quantile(all_vals, 0.98, na.rm = TRUE))
cat("  Common 98th-percentile max:", round(max_val, 2), "yr\n")

to_255 <- function(r, m) clamp(r / m, 0, 1) * 255

red   <- to_255(edaphic_mag,    max_val)
green <- to_255(landuse_mag,    max_val)
blue  <- to_255(biological_mag, max_val)

# =============================================================================
# COMMON LAND MASK
# =============================================================================

if (file.exists(m7_file)) {
  land_mask <- !is.na(rast(m7_file))
  cat("  Using M7 as common land mask\n")
} else {
  land_mask <- !is.na(m4)
  cat("  WARNING: M7 not found, falling back to M4 mask\n")
}

red   <- mask(red,   land_mask, maskvalues = 0)
green <- mask(green, land_mask, maskvalues = 0)
blue  <- mask(blue,  land_mask, maskvalues = 0)

# =============================================================================
# AGGREGATE FOR DISPLAY
# =============================================================================

agg_factor <- round(VIS_RES / res(red)[1])

if (agg_factor > 1) {
  cat("Aggregating by factor", agg_factor, "...\n")
  red   <- aggregate(red,   fact = agg_factor, fun = "mean", na.rm = TRUE)
  green <- aggregate(green, fact = agg_factor, fun = "mean", na.rm = TRUE)
  blue  <- aggregate(blue,  fact = agg_factor, fun = "mean", na.rm = TRUE)
  land_mask <- aggregate(land_mask, fact = agg_factor, fun = "max",
                         na.rm = TRUE)
}

# =============================================================================
# SAVE MULTI-BAND GEOTIFF (raw magnitudes)
# =============================================================================

rgb_stack <- c(red, green, blue)
names(rgb_stack) <- c("edaphic_mag", "landuse_mag", "biological_mag")
rgb_output <- file.path(PRED_DIR, "process_dominance_RGB.tif")
writeRaster(rgb_stack, rgb_output, overwrite = TRUE)
cat("\u2713 Saved GeoTIFF:", rgb_output, "\n\n")

# =============================================================================
# CLIMATE-DOMINANCE AUXILIARY LAYER
# =============================================================================

total_nonclim <- abs(m2 - m1) + abs(m3 - m1) + abs(m4 - m1)
climate_dom   <- max_val - total_nonclim
climate_dom   <- clamp(climate_dom, 0, max_val)
climate_dom   <- mask(climate_dom, !is.na(rast(m7_file)), maskvalues = 0)
if (agg_factor > 1) {
  climate_dom <- aggregate(climate_dom, fact = agg_factor, fun = "mean",
                           na.rm = TRUE)
}

# =============================================================================
# DATA EXTENT
# =============================================================================

trimmed     <- trim(red)
data_extent <- ext(trimmed)
data_extent <- ext(
  max(-180, xmin(data_extent) - 2),
  min( 180, xmax(data_extent) + 2),
  max( -90, ymin(data_extent) - 2),
  min(  90, ymax(data_extent) + 2)
)

# =============================================================================
# BUILD CMY AND RGB PLOTTING STACKS
# =============================================================================
#
# RGB (additive): raw normalised magnitudes.
#   Black = no contribution; full white = all three strong.
#   Luminance modulation pulls toward BLACK where climate dominates.
#
# CMY (subtractive): inverted magnitudes (255 - raw).
#   White = no contribution; black = all three strong.
#   Luminance modulation pulls toward WHITE where climate dominates.
# =============================================================================

climate_factor <- clamp(crop(climate_dom, data_extent) / max_val, 0, 1)

# --- RGB (additive) ---------------------------------------------------------

red_rgb   <- crop(red,   data_extent)
green_rgb <- crop(green, data_extent)
blue_rgb  <- crop(blue,  data_extent)

red_rgb   <- red_rgb   * (1 - climate_factor * LUMINANCE_BLEND)
green_rgb <- green_rgb * (1 - climate_factor * LUMINANCE_BLEND)
blue_rgb  <- blue_rgb  * (1 - climate_factor * LUMINANCE_BLEND)

rgb_plot_stack <- c(red_rgb, green_rgb, blue_rgb)

# --- CMY (subtractive) ------------------------------------------------------

red_cmy   <- 255 - crop(red,   data_extent)
green_cmy <- 255 - crop(green, data_extent)
blue_cmy  <- 255 - crop(blue,  data_extent)

red_cmy   <- red_cmy   + (255 - red_cmy)   * climate_factor * LUMINANCE_BLEND
green_cmy <- green_cmy + (255 - green_cmy) * climate_factor * LUMINANCE_BLEND
blue_cmy  <- blue_cmy  + (255 - blue_cmy)  * climate_factor * LUMINANCE_BLEND

cmy_plot_stack <- c(red_cmy, green_cmy, blue_cmy)

# =============================================================================
# COUNTRY BORDERS
# =============================================================================

cat("Loading Natural Earth borders...\n")
world <- st_geometry(ne_countries(scale = 50, returnclass = "sf"))

sf::sf_use_s2(FALSE)
bbox <- st_bbox(c(xmin = xmin(data_extent),
                  xmax = xmax(data_extent),
                  ymin = ymin(data_extent),
                  ymax = ymax(data_extent)),
                crs = st_crs(world))
world_crop <- st_crop(world, bbox)
sf::sf_use_s2(TRUE)

# =============================================================================
# FIGURE 1 -- TWO-PANEL COMBINED
# =============================================================================

cat("\nWriting two-panel combined figure...\n")

climate_plot <- mask(crop(climate_dom, data_extent),
                     !is.na(crop(cmy_plot_stack[[1]], data_extent)))

png_combined <- file.path(PLOT_DIR, "MRT_process_dominance_combined.png")
png(png_combined, width = 1200, height = 1300, res = 300, bg = "white")
layout(matrix(c(1, 2, 3), nrow = 3, byrow = TRUE), heights = c(4, 4, 1.2))

par(mar = c(0, 0, 1.2, 0), bg = "white")
climate_cols <- colorRampPalette(c("#5E3C99", "#E66101"))(100)
plot(climate_plot,
     col = climate_cols, legend = FALSE, axes = FALSE,
     mar = c(0, 0, 0.8, 0), background = "white",
     xlim = c(xmin(data_extent), xmax(data_extent)),
     ylim = c(ymin(data_extent), ymax(data_extent)))
plot(world_crop, add = TRUE, col = NA, border = "grey30", lwd = 0.3)
title("A) Climate explains most of the prediction",
      cex.main = 0.95, line = 0.3)

par(mar = c(0, 0, 1.2, 0))
plotRGB(cmy_plot_stack, r = 1, g = 2, b = 3, stretch = "lin",
        axes = FALSE, mar = c(0, 0, 0.8, 0), bgalpha = 0,
        xlim = c(xmin(data_extent), xmax(data_extent)),
        ylim = c(ymin(data_extent), ymax(data_extent)))
plot(world_crop, add = TRUE, col = NA, border = "grey30", lwd = 0.3)
title("B) Magnitude of non-climate prediction adjustment",
      cex.main = 0.95, line = 0.3)

par(mar = c(0, 0, 0, 0))
plot.new()
legend(x = 0.00, y = 0.95,
       legend = c("Climate alone explains most",
                  "Non-climate processes modulate"),
       fill = c("#5E3C99", "#E66101"), border = "gray50",
       title = "Climate contribution", title.font = 2,
       cex = 0.62, bty = "n")
legend(x = 0.35, y = 0.95,
       legend = c("Edaphic", "LandUse", "Biological",
                  "Edaph+Land", "Land+Bio", "Edaph+Bio",
                  "All three", "None (climate only)"),
       fill = c("cyan", "magenta", "yellow",
                "blue", "red", "green", "black", "white"),
       border = "gray50",
       title = "Non-climate channel (magnitude)", title.font = 2,
       ncol = 4, cex = 0.62, bty = "n")
dev.off()
cat("\u2713 Saved:", basename(png_combined), "\n")

# =============================================================================
# FIGURE 2 -- SINGLE-PANEL CMY
# =============================================================================

cat("Writing single-panel CMY figure...\n")

png_cmy <- file.path(PLOT_DIR, "MRT_process_contribution_CMY.png")
png(png_cmy, width = 1600, height = 900, res = 200, bg = "white")
layout(matrix(c(1, 2), nrow = 2), heights = c(5, 1))

par(mar = c(0.5, 0.5, 1.8, 0.5))
plotRGB(cmy_plot_stack, r = 1, g = 2, b = 3, stretch = "lin",
        axes = FALSE, mar = c(0.2, 0.2, 1.5, 0.2), bgalpha = 0,
        xlim = c(xmin(data_extent), xmax(data_extent)),
        ylim = c(ymin(data_extent), ymax(data_extent)))
plot(world_crop, add = TRUE, col = NA, border = "grey30", lwd = 0.3)
title("Magnitude of non-climate prediction adjustment (CMY)",
      cex.main = 1.15, line = 0.4)

par(mar = c(0.2, 1, 0.2, 1))
plot.new()
legend("center",
       legend = c("Edaphic", "LandUse", "Biological",
                  "Edaph+Land", "Land+Bio", "Edaph+Bio",
                  "All three", "None (climate only)"),
       fill = c("cyan", "magenta", "yellow",
                "blue", "red", "green", "black", "white"),
       border = "gray50", ncol = 4, cex = 0.85, bty = "n",
       title = "Channel intensity = |Mi - M1| per group")
dev.off()
cat("\u2713 Saved:", basename(png_cmy), "\n")

# =============================================================================
# FIGURE 3 -- SINGLE-PANEL RGB
# =============================================================================

cat("Writing single-panel RGB figure...\n")

png_rgb <- file.path(PLOT_DIR, "MRT_process_contribution_RGB.png")
png(png_rgb, width = 1600, height = 900, res = 200, bg = "white")
layout(matrix(c(1, 2), nrow = 2), heights = c(5, 1))

par(mar = c(0.5, 0.5, 1.8, 0.5), bg = "white")
plotRGB(rgb_plot_stack, r = 1, g = 2, b = 3, stretch = "lin",
        axes = FALSE, mar = c(0.2, 0.2, 1.5, 0.2), bgalpha = 0,
        colNA = "black",
        xlim = c(xmin(data_extent), xmax(data_extent)),
        ylim = c(ymin(data_extent), ymax(data_extent)))
plot(world_crop, add = TRUE, col = NA, border = "grey80", lwd = 0.3)
title("Magnitude of non-climate prediction adjustment (RGB)",
      cex.main = 1.15, line = 0.4)

par(mar = c(0.2, 1, 0.2, 1))
plot.new()
legend("center",
       legend = c("Edaphic", "LandUse", "Biological",
                  "Edaph+Land (yellow)", "Land+Bio (magenta)",
                  "Edaph+Bio (cyan)",
                  "All three (white)", "None (black, climate only)"),
       fill = c("red", "green", "blue",
                "yellow", "magenta", "cyan",
                "white", "black"),
       border = "gray50", ncol = 4, cex = 0.85, bty = "n",
       title = "Channel intensity = |Mi - M1| per group")
dev.off()
cat("\u2713 Saved:", basename(png_rgb), "\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\u2550\u2550\u2550 SUMMARY (for manuscript) \u2550\u2550\u2550\n\n")
cat("Per-pixel prediction adjustment magnitude |Mi - M1|:\n")
summarise_mag(edaphic_mag,    "Edaphic")
summarise_mag(landuse_mag,    "LandUse")
summarise_mag(biological_mag, "Biological")
cat("\nNote: these figures show MAGNITUDE of prediction shift, not\n")
cat("statistical dominance. For variance-explained rankings use the\n")
cat("permutation importance output from script 13.\n\n")

cat("\u2550\u2550\u2550 STEP 16 COMPLETE \u2550\u2550\u2550\n\n")