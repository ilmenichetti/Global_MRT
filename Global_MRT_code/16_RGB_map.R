################################################################################
# 16_process_dominance_RGB.R
#
# Create CMY composite map showing where different processes dominate MRT
#
# CMY Channels (beyond-climate contributions, subtractive color mixing):
#   Cyan    = Edaphic (intrinsic soil properties: texture, chemistry, terrain)
#   Magenta = LandUse (land cover fractions, disturbance history)
#   Yellow  = Biological (soil ecology: mycorrhiza, fungi)
#
# Color interpretation (subtractive/ink logic):
#   Cyan     = Edaphic dominates (soil texture, pH, topography)
#   Magenta  = LandUse dominates (vegetation cover, fire, forest loss)
#   Yellow   = Biological dominates (mycorrhiza, fungal communities)
#   Blue     = Edaphic + LandUse
#   Red      = LandUse + Biological
#   Green    = Edaphic + Biological
#   Black    = All three contribute equally
#   White    = Climate explains most (non-climate processes add little)
#
# Luminance modulation:
#   Colors are washed out (lighter) where climate dominates
#   Colors stay saturated (darker) where non-climate processes are strong
#
# Two-panel output:
#   Panel A: Climate contribution (colored: blue-red)
#   Panel B: CMY composite (Edaphic/LandUse/Biology) with luminance modulation
#
# Input:  ./Global_MRT_code/outputs/MRT_predictions/*.tif
# Output: ./Global_MRT_code/plots/MRT_process_dominance_combined.png
#         ./Global_MRT_code/plots/MRT_dominant_process.png
#         ./Global_MRT_code/outputs/MRT_predictions/process_dominance_RGB.tif
#
# Author: Lorenzo
# Date: 2026-01-12
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

# Visualization resolution
VIS_RES <- 0.5  # degrees

# Luminance blend strength: how much to wash out climate-dominated areas
# 0.0 = no luminance effect, 1.0 = full wash to white where climate dominates
LUMINANCE_BLEND <- 0.9

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PROCESS DOMINANCE RGB COMPOSITE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD PREDICTION MAPS
# =============================================================================

cat("Loading prediction maps...\n")

m1_file <- file.path(PRED_DIR, "MRT_M1_climate.tif")
m2_file <- file.path(PRED_DIR, "MRT_M2_climate_edaphic.tif")
m3_file <- file.path(PRED_DIR, "MRT_M3_climate_landuse.tif")
m4_file <- file.path(PRED_DIR, "MRT_M4_climate_biological.tif")

# Check files exist
required_files <- c(m1_file, m2_file, m3_file, m4_file)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing required files:\n  ", paste(missing, collapse = "\n  "))
}

m1 <- rast(m1_file)  # Climate only (baseline)
m2 <- rast(m2_file)  # Climate + Edaphic
m3 <- rast(m3_file)  # Climate + LandUse
m4 <- rast(m4_file)  # Climate + Biological

cat("  ✓ Loaded M1 (climate), M2 (edaphic), M3 (landuse), M4 (biological)\n\n")

# =============================================================================
# CALCULATE PROCESS CONTRIBUTIONS
# =============================================================================
#
# Each channel represents the CHANGE in predicted MRT when adding that 
# process group to climate-only model.
#
# Positive values = adding that process INCREASES predicted MRT
# Negative values = adding that process DECREASES predicted MRT
#
# For RGB visualization, we use ABSOLUTE contribution (magnitude matters,
# not direction), then normalize to 0-255 range.
# =============================================================================

cat("Calculating process contributions...\n")

# Raw differences (can be positive or negative)
edaphic_contrib <- m2 - m1    # What edaphic vars add
landuse_contrib <- m3 - m1    # What landuse vars add
biological_contrib <- m4 - m1 # What biological vars add

cat("  Edaphic contribution range:", 
    round(min(values(edaphic_contrib), na.rm = TRUE), 1), "to",
    round(max(values(edaphic_contrib), na.rm = TRUE), 1), "years\n")
cat("  LandUse contribution range:",
    round(min(values(landuse_contrib), na.rm = TRUE), 1), "to",
    round(max(values(landuse_contrib), na.rm = TRUE), 1), "years\n")
cat("  Biological contribution range:",
    round(min(values(biological_contrib), na.rm = TRUE), 1), "to",
    round(max(values(biological_contrib), na.rm = TRUE), 1), "years\n\n")

# =============================================================================
# OPTION 1: ABSOLUTE CONTRIBUTION (magnitude only)
# =============================================================================
# Use absolute values - we care about HOW MUCH each process matters,
# not whether it increases or decreases MRT

cat("Using absolute contributions (magnitude)...\n")

red_raw <- abs(edaphic_contrib)
green_raw <- abs(landuse_contrib)
blue_raw <- abs(biological_contrib)

# =============================================================================
# OPTION 2: SIGNED CONTRIBUTION (uncomment to use)
# =============================================================================
# # Positive = process increases MRT (slower turnover)
# # Negative = process decreases MRT (faster turnover)
# # Center at zero, stretch to show both directions
#
# red_raw <- edaphic_contrib
# green_raw <- landuse_contrib
# blue_raw <- biological_contrib

# =============================================================================
# NORMALIZE TO 0-255 FOR RGB
# =============================================================================

cat("Normalizing to RGB scale...\n")

# Find common scale (so colors are comparable)
# Use 98th percentile to avoid extreme outliers dominating
all_vals <- c(values(red_raw), values(green_raw), values(blue_raw))
max_val <- quantile(all_vals, 0.98, na.rm = TRUE)
cat("  Normalizing to max value:", round(max_val, 2), "\n")

normalize_to_255 <- function(r, max_val) {
  r_norm <- r / max_val
  r_norm <- clamp(r_norm, lower = 0, upper = 1)
  r_norm <- r_norm * 255
  return(r_norm)
}

red <- normalize_to_255(red_raw, max_val)
green <- normalize_to_255(green_raw, max_val)
blue <- normalize_to_255(blue_raw, max_val)

# =============================================================================
# CREATE COMMON MASK
# =============================================================================

# Use M4 (biological) as mask since it has most limited coverage
mask_layer <- !is.na(m4)
red <- mask(red, mask_layer, maskvalues = 0)
green <- mask(green, mask_layer, maskvalues = 0)
blue <- mask(blue, mask_layer, maskvalues = 0)

# =============================================================================
# AGGREGATE FOR VISUALIZATION
# =============================================================================

agg_factor <- round(VIS_RES / res(red)[1])
if (agg_factor > 1) {
  cat("Aggregating by factor", agg_factor, "...\n")
  red_agg <- aggregate(red, fact = agg_factor, fun = "mean", na.rm = TRUE)
  green_agg <- aggregate(green, fact = agg_factor, fun = "mean", na.rm = TRUE)
  blue_agg <- aggregate(blue, fact = agg_factor, fun = "mean", na.rm = TRUE)
} else {
  red_agg <- red
  green_agg <- green
  blue_agg <- blue
}

# =============================================================================
# CREATE RGB COMPOSITE
# =============================================================================

cat("Creating CMY composite (inverted RGB)...\n")

# Stack as RGB
rgb_stack <- c(red_agg, green_agg, blue_agg)
names(rgb_stack) <- c("red", "green", "blue")

# =============================================================================
# SAVE AS GEOTIFF
# =============================================================================

rgb_output <- file.path(PRED_DIR, "process_dominance_RGB.tif")
writeRaster(rgb_stack, rgb_output, overwrite = TRUE)
cat("✓ Saved GeoTIFF:", rgb_output, "\n")

# =============================================================================
# CALCULATE DATA EXTENT FOR CROPPING
# =============================================================================

cat("Calculating data extent...\n")

# Trim the actual data raster (not a boolean) to get true data bounds
trimmed_raster <- trim(red_agg)
data_extent <- ext(trimmed_raster)

cat("  Raw data extent: x [", round(xmin(data_extent), 1), ",", round(xmax(data_extent), 1), "]",
    " y [", round(ymin(data_extent), 1), ",", round(ymax(data_extent), 1), "]\n")

# Add small buffer (2 degrees) but don't exceed global bounds
data_extent <- ext(
  max(-180, xmin(data_extent) - 2),
  min(180, xmax(data_extent) + 2),
  max(-90, ymin(data_extent) - 2),
  min(90, ymax(data_extent) + 2)
)

cat("  Buffered extent: x [", round(xmin(data_extent), 1), ",", round(xmax(data_extent), 1), "]",
    " y [", round(ymin(data_extent), 1), ",", round(ymax(data_extent), 1), "]\n")

# =============================================================================
# LOAD HIGH-RESOLUTION WORLD BORDERS
# =============================================================================

cat("Loading high-resolution country borders...\n")

# Get 1:50m scale Natural Earth countries (higher detail than maps package)
world_borders <- ne_countries(scale = 50, returnclass = "sf")
world_borders <- st_geometry(world_borders)

cat("  ✓ Loaded Natural Earth 1:50m borders\n")

# =============================================================================
# PREPARE CLIMATE CONTRIBUTION LAYER (GRAYSCALE)
# =============================================================================
#
# Show how much variance climate alone explains
# Using predicted MRT from M1 as proxy for climate influence
# Higher MRT in climate-only model = climate predicts slower turnover
# 
# Alternative: Use total non-climate contribution (sum of RGB) as darkness
# Dark = climate dominates, Bright = non-climate processes matter
# =============================================================================

cat("Preparing climate contribution layer...\n")

# Option A: Climate model prediction (scaled)
# climate_layer <- m1
# climate_layer <- clamp(climate_layer, lower = 0, upper = 200)

# Option B: Inverse of total non-climate contribution
# Higher total contribution = brighter in RGB = less climate dominance
# So grayscale shows: bright = climate matters more, dark = climate matters less
total_nonclimate <- red_raw + green_raw + blue_raw
climate_dominance <- max_val - total_nonclimate  # Invert so high = climate dominant
climate_dominance <- clamp(climate_dominance, lower = 0, upper = max_val)

# Normalize to 0-255 grayscale
climate_gray <- (climate_dominance / max_val) * 255
climate_gray <- mask(climate_gray, mask_layer, maskvalues = 0)

# Aggregate
if (agg_factor > 1) {
  climate_gray_agg <- aggregate(climate_gray, fact = agg_factor, fun = "mean", na.rm = TRUE)
} else {
  climate_gray_agg <- climate_gray
}

# =============================================================================
# CREATE TWO-PANEL FIGURE: RGB + CLIMATE GRAYSCALE
# =============================================================================

cat("Creating two-panel figure...\n")

# ─────────────────────────────────────────────────────────────────────────────
# Get land mask from M7 (full model - requires ALL datasets)
# ─────────────────────────────────────────────────────────────────────────────

m7_file <- file.path(PRED_DIR, "MRT_M7_full.tif")
if (file.exists(m7_file)) {
  m7 <- rast(m7_file)
  land_mask <- !is.na(m7)
  cat("  Using M7 (full model) as common extent mask\n")
} else {
  # Fallback to M4 if M7 not available
  land_mask <- !is.na(m4)
  cat("  WARNING: M7 not found, using M4 as mask\n")
}

# Aggregate mask
if (agg_factor > 1) {
  land_mask_agg <- aggregate(land_mask, fact = agg_factor, fun = "max", na.rm = TRUE)
} else {
  land_mask_agg <- land_mask
}

# ─────────────────────────────────────────────────────────────────────────────
# Prepare RGB layers (set ocean to NA, which will plot as white)
# ─────────────────────────────────────────────────────────────────────────────

red_plot <- mask(red_agg, land_mask_agg, maskvalues = 0)
green_plot <- mask(green_agg, land_mask_agg, maskvalues = 0)
blue_plot <- mask(blue_agg, land_mask_agg, maskvalues = 0)

# Crop to data extent
red_plot <- crop(red_plot, data_extent)
green_plot <- crop(green_plot, data_extent)
blue_plot <- crop(blue_plot, data_extent)

# Invert to CMY logic: absence = white, all three = black
# This uses subtractive color mixing (like print/ink)
red_plot <- 255 - red_plot
green_plot <- 255 - green_plot
blue_plot <- 255 - blue_plot

# ─────────────────────────────────────────────────────────────────────────────
# Add luminance modulation based on climate dominance
# Areas where climate dominates → lighter (blend toward white)
# Areas where other processes dominate → keep saturated colors
# ─────────────────────────────────────────────────────────────────────────────

# Use the already-aggregated climate gray layer (same resolution as RGB)
climate_factor <- crop(climate_gray_agg, data_extent)

# Normalize to 0-1 (higher = more climate dominant)
climate_factor <- climate_factor / 255
climate_factor <- clamp(climate_factor, lower = 0, upper = 1)

# Blend each channel toward 255 (white) based on climate dominance
red_plot <- red_plot + (255 - red_plot) * climate_factor * LUMINANCE_BLEND
green_plot <- green_plot + (255 - green_plot) * climate_factor * LUMINANCE_BLEND
blue_plot <- blue_plot + (255 - blue_plot) * climate_factor * LUMINANCE_BLEND

cat("  ✓ Applied luminance modulation (blend strength:", LUMINANCE_BLEND, ")\n")

# Stack (keep NA as NA - will appear white in plot)
rgb_stack_plot <- c(red_plot, green_plot, blue_plot)

# ─────────────────────────────────────────────────────────────────────────────
# Prepare climate grayscale layer
# ─────────────────────────────────────────────────────────────────────────────

climate_plot <- mask(climate_gray_agg, land_mask_agg, maskvalues = 0)

# Crop to data extent
climate_plot <- crop(climate_plot, data_extent)

# Crop world borders to data extent for proper zoom
# Disable S2 spherical geometry to avoid edge crossing errors
sf_use_s2(FALSE)

crop_bbox <- st_bbox(c(xmin = xmin(data_extent), 
                       xmax = xmax(data_extent),
                       ymin = ymin(data_extent), 
                       ymax = ymax(data_extent)),
                     crs = st_crs(world_borders))

world_borders_crop <- st_crop(world_borders, crop_bbox)

sf_use_s2(TRUE)  # Re-enable S2

cat("  ✓ Cropped borders to data extent\n")

# ─────────────────────────────────────────────────────────────────────────────
# Create combined figure (Climate dominance + RGB)
# ─────────────────────────────────────────────────────────────────────────────

png_output <- file.path(PLOT_DIR, "MRT_process_dominance_combined.png")
png(png_output, width = 1200, height = 1300, res = 300, bg = "white")

# Layout: 2 maps stacked + 1 legend row at bottom
layout(matrix(c(1, 2, 3), nrow = 3, byrow = TRUE), heights = c(4, 4, 1.2))

# ─────────────────────────────────────────────────────────────────────────────
# Panel A: Climate Dominance (top) - colored palette
# ─────────────────────────────────────────────────────────────────────────────

par(mar = c(0, 0, 1.2, 0), bg = "white")

# Colored palette: pastel blue (low climate) to red (high climate)
#climate_cols <- colorRampPalette(c("#4575B4", "#D73027"))(100)
climate_cols <- colorRampPalette(c("#5E3C99", "#E66101"))(100)

plot(climate_plot,
     col = climate_cols,
     legend = FALSE,
     axes = FALSE,
     mar = c(0, 0, 0.8, 0),
     background = "white",
     xlim = c(xmin(data_extent), xmax(data_extent)),
     ylim = c(ymin(data_extent), ymax(data_extent)))

# Add high-resolution country borders (cropped)
plot(world_borders_crop, add = TRUE, col = NA, border = "grey30", lwd = 0.3)

title("A) Climate Dominance", cex.main = 1, line = 0.3)

# ─────────────────────────────────────────────────────────────────────────────
# Panel B: RGB Composite (bottom)
# ─────────────────────────────────────────────────────────────────────────────

par(mar = c(0, 0, 1.2, 0))

plotRGB(rgb_stack_plot, 
        r = 1, g = 2, b = 3,
        stretch = "lin",
        axes = FALSE,
        mar = c(0, 0, 0.8, 0),
        bgalpha = 0,
        xlim = c(xmin(data_extent), xmax(data_extent)),
        ylim = c(ymin(data_extent), ymax(data_extent)))

# Add high-resolution country borders (cropped)
plot(world_borders_crop, add = TRUE, col = NA, border = "grey30", lwd = 0.3)

title("B) Process Contribution", cex.main = 1, line = 0.3)

# ─────────────────────────────────────────────────────────────────────────────
# Panel C: Combined Legend (bottom)
# ─────────────────────────────────────────────────────────────────────────────

par(mar = c(0, 0, 0, 0))
plot.new()

# Left side: Climate dominance legend
legend(x = 0.0, y = 0.95,
       legend = c("High (climate dominates)", 
                  "Low (other processes)"),
       fill = c("#5E3C99", "#E66101"),
       border = "gray50",
       title = "Climate Dominance",
       title.font = 2,
       ncol = 1,
       cex = 0.65,
       bty = "n")



# Right side: CMY legend (inverted from RGB)
legend(x = 0.35, y = 0.95, 
       legend = c("Edaphic", "LandUse", "Biological",
                  "Edaph+Land", "Land+Bio", "Edaph+Bio",
                  "All three", "Climate only"),
       fill = c("cyan", "magenta", "yellow", 
                "blue", "red", "green",
                "black", "white"),
       border = "gray50",
       title = "Process Contribution",
       title.font = 2,
       ncol = 4,
       cex = 0.65,
       bty = "n")

dev.off()

cat("✓ Saved:", png_output, "\n")

# =============================================================================
# ALSO SAVE PROCESS CONTRIBUTION PLOT (standalone)
# =============================================================================

png(file.path(PLOT_DIR, "MRT_process_contribution_RGB.png"), 
    width = 1600, height = 1000, res = 150, bg = "white")

layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1))

par(mar = c(0.5, 0.5, 2, 0.5))
plotRGB(rgb_stack_plot, r = 1, g = 2, b = 3, stretch = "lin",
        axes = FALSE, mar = c(0.2, 0.2, 1.5, 0.2), bgalpha = 0,
        xlim = c(xmin(data_extent), xmax(data_extent)),
        ylim = c(ymin(data_extent), ymax(data_extent)))
plot(world_borders_crop, add = TRUE, col = NA, border = "grey30", lwd = 0.3)
title("MRT Process Contribution", cex.main = 1.2)

par(mar = c(0.5, 2, 0.5, 2))
plot.new()
legend("center", 
       legend = c("Edaphic", "LandUse", "Biological",
                  "Edaph + Land", "Land + Bio", "Edaph + Bio",
                  "All three", "Climate only"),
       fill = c("cyan", "magenta", "yellow", "blue", "red", "green", "black", "white"),
       border = "gray50", ncol = 4, cex = 0.9, bty = "n")

dev.off()
cat("✓ Saved: MRT_process_contribution_RGB.png\n")

# =============================================================================
# CREATE DOMINANT PROCESS MAP (SEPARATE FILE)
# =============================================================================

cat("Creating dominant process map...\n")

# ─────────────────────────────────────────────────────────────────────────────
# Prepare dominant process map
# ─────────────────────────────────────────────────────────────────────────────

# Which process has highest absolute contribution at each pixel?
dominant <- red_raw
values(dominant) <- NA

edaph_vals_raw <- values(red_raw)
land_vals_raw <- values(green_raw)
bio_vals_raw <- values(blue_raw)

# Find maximum
max_contrib <- pmax(edaph_vals_raw, land_vals_raw, bio_vals_raw, na.rm = TRUE)

# Assign category: 1=Edaphic, 2=LandUse, 3=Biological
dom_vals <- rep(NA, length(edaph_vals_raw))
dom_vals[edaph_vals_raw == max_contrib & !is.na(max_contrib)] <- 1
dom_vals[land_vals_raw == max_contrib & !is.na(max_contrib)] <- 2
dom_vals[bio_vals_raw == max_contrib & !is.na(max_contrib)] <- 3

values(dominant) <- dom_vals

# Aggregate
if (agg_factor > 1) {
  dominant_agg <- aggregate(dominant, fact = agg_factor, fun = "modal", na.rm = TRUE)
} else {
  dominant_agg <- dominant
}

# Apply common mask
dominant_agg <- mask(dominant_agg, land_mask_agg, maskvalues = 0)

# Crop to data extent
dominant_agg <- crop(dominant_agg, data_extent)

# ─────────────────────────────────────────────────────────────────────────────
# Save dominant process plot
# ─────────────────────────────────────────────────────────────────────────────

# Colors matching the RGB scheme
cols_cat <- c("#E41A1C",   # Edaphic - red
              "#4DAF4A",   # LandUse - green  
              "#377EB8")   # Biological - blue

png(file.path(PLOT_DIR, "MRT_dominant_process.png"), 
    width = 1600, height = 1000, res = 150, bg = "white")

layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1))

par(mar = c(0.5, 0.5, 2, 0.5))

# Plot data directly (NA will be white)
plot(dominant_agg, col = cols_cat, legend = FALSE, axes = FALSE,
     mar = c(0.2, 0.2, 1.5, 0.2), background = "white",
     xlim = c(xmin(data_extent), xmax(data_extent)),
     ylim = c(ymin(data_extent), ymax(data_extent)))
plot(world_borders_crop, add = TRUE, col = NA, border = "grey30", lwd = 0.3)
title("Dominant Non-Climate Process for MRT", cex.main = 1.2)

par(mar = c(0.5, 2, 0.5, 2))
plot.new()
legend("center", legend = c("Edaphic", "LandUse", "Biological"),
       fill = cols_cat, border = "gray50", ncol = 3, cex = 1.1, bty = "n")

dev.off()
cat("✓ Saved: MRT_dominant_process.png\n")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# What fraction of land is dominated by each process?
dom_table <- table(dom_vals, useNA = "no")
dom_pct <- round(100 * dom_table / sum(dom_table), 1)

cat("Dominant process (% of land area):\n")
cat("  Edaphic:    ", dom_pct[1], "%\n")
cat("  LandUse:    ", dom_pct[2], "%\n")
cat("  Biological: ", dom_pct[3], "%\n")

# Mean absolute contribution by process
cat("\nMean absolute contribution (years MRT change):\n")
cat("  Edaphic:    ", round(mean(edaph_vals_raw, na.rm = TRUE), 2), "\n")
cat("  LandUse:    ", round(mean(land_vals_raw, na.rm = TRUE), 2), "\n")
cat("  Biological: ", round(mean(bio_vals_raw, na.rm = TRUE), 2), "\n")

# Climate dominance stats
cat("\nClimate dominance (grayscale intensity):\n")
clim_dom_vals <- values(climate_dominance)
cat("  Mean:", round(mean(clim_dom_vals, na.rm = TRUE), 2), "\n")
cat("  Median:", round(median(clim_dom_vals, na.rm = TRUE), 2), "\n")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Outputs:\n")
cat("  ", file.path(PLOT_DIR, "MRT_process_dominance_combined.png"), "\n")
cat("  ", file.path(PLOT_DIR, "MRT_process_contribution_RGB.png"), "\n")
cat("  ", file.path(PLOT_DIR, "MRT_dominant_process.png"), "\n")
cat("  ", rgb_output, "\n")