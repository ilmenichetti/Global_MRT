# =============================================================================
# Step 08b: Peatland + open-water exclusion layers
# =============================================================================
# Produces two independent products used later in the pipeline:
#   (A) a per-observation peat flag  -> feeds the fit filter in 13c
#   (B) a 0.1deg peat/water fraction + combined exclusion mask -> feeds the maps
#
# Sources:
#   - Global Peatland Map 2.0 (Greifswald Mire Centre 2022; underlying dataset of
#     the UNEP Global Peatlands Assessment). CC BY-NC-SA 4.0. GeoTIFF, WGS84,
#     ~0.01deg (~1.1 km), values 1 = peat-dominated, 2 = peat in soil mosaic,
#     NoData = 255. Kept local (Datasets/ is gitignored); cite in the manuscript.
#   - ESA land-cover water fraction: band 8 ("water") of the existing
#     landcover_fractions_0p1deg.tif (already on the 0.1deg grid; no download).
#
# Decisions (see manuscript/decisions/2026-07-08_peat_exclusion_plan.md):
#   D1  exclude BOTH peat tiers (1 + 2)          -> MRT_PEAT_TIERS=1,2
#   D2  mask a 0.1deg cell at peat fraction >=0.5 -> MRT_PEAT_THRESH=0.5
#   D3  drop an observation whose GPM value is peat (tier per D1)
#   D5  keep the modal-Histosol flag alongside as a secondary cross-check
#   D6  exclude open water too: water fraction >=0.5 (same cutoff as D2)
# =============================================================================

library(terra)
library(dplyr)

source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# --- Configuration -----------------------------------------------------------
GPM_TIF <- "./Datasets/peat/GPM2022_GMC/GLOpeat_GPA22WGS_2cl_1x1km/peatGPA22WGS_2cl.tif"
LC_FRAC <- file.path(PATHS$raster_output_dir, "landcover", "landcover_fractions_0p1deg.tif")
LC_WATER_BAND <- 8L   # band 8 = "water" (verified via band descriptions)

PRED_DIR    <- file.path(PATHS$output_dir, "MRT_predictions")
TEMPLATE    <- file.path(PRED_DIR, "MRT_M7_full.tif")   # 0.1deg land grid (NA over non-land)
INPUT_TABLE <- file.path(PATHS$output_dir, "12b_model_ready.rds")

PLOT_DIR <- "./Global_MRT_code/plots/step_08b_peat"

# Env-configurable robustness knobs (defaults = the locked D1/D2/D6 choices)
PEAT_TIERS <- as.integer(strsplit(Sys.getenv("MRT_PEAT_TIERS", "1,2"), ",")[[1]])
THRESH     <- as.numeric(Sys.getenv("MRT_PEAT_THRESH", "0.5"))
stopifnot(all(PEAT_TIERS %in% c(1L, 2L)), THRESH > 0, THRESH <= 1)

cat("\n")
cat("+--------------------------------------------------------------+\n")
cat("|  Step 08b: Peat (GPM 2.0) + open-water exclusion layers       |\n")
cat(sprintf("|  Tiers excluded: %-8s  Mask threshold: %-5.2f            |\n",
            paste(PEAT_TIERS, collapse = "+"), THRESH))
cat("+--------------------------------------------------------------+\n\n")

# =============================================================================
# Part A -- per-observation peat flag (feeds the fit)
# =============================================================================
log_step("Part A: extracting per-observation peat flag from GPM 2.0")

if (!file.exists(GPM_TIF)) stop("GPM raster not found: ", GPM_TIF)
if (!file.exists(INPUT_TABLE)) stop("Model-ready table not found: ", INPUT_TABLE)

soil_data <- readRDS(INPUT_TABLE)
cat(sprintf("  Observations: %s (unique dsiteid: %s)\n",
            format(nrow(soil_data), big.mark = ","),
            format(length(unique(soil_data$dsiteid)), big.mark = ",")))

gpm <- rast(GPM_TIF)   # already WGS84 -> no reprojection needed
pts <- vect(soil_data,
            geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
            crs  = "EPSG:4326")

gpm_val   <- terra::extract(gpm, pts, ID = FALSE)[, 1]   # 1, 2, or NA (NoData 255)
peat_tier <- ifelse(is.na(gpm_val), 0L, as.integer(gpm_val))   # 0 = none
peat_flag_gpm <- peat_tier %in% PEAT_TIERS

# D5: keep the modal-Histosol flag alongside for cross-check
is_histosol <- !is.na(soil_data$soilclass_wrb_name) &
               soil_data$soilclass_wrb_name == "Histosols"

peat_flags <- data.frame(
  dsiteid        = soil_data$dsiteid,
  peat_gpm_value = gpm_val,
  peat_tier      = peat_tier,
  peat_flag_gpm  = peat_flag_gpm,
  is_histosol    = is_histosol,
  stringsAsFactors = FALSE
)

# Coverage + agreement diagnostics
n_pd  <- sum(peat_tier == 1L)
n_mos <- sum(peat_tier == 2L)
cat(sprintf("  GPM at obs: peat-dominated %s (%.2f%%), mosaic %s (%.2f%%)\n",
            format(n_pd, big.mark = ","),  100 * n_pd  / nrow(peat_flags),
            format(n_mos, big.mark = ","), 100 * n_mos / nrow(peat_flags)))
cat(sprintf("  Flagged peat (tiers %s): %s obs (%.2f%%) -> would be dropped from the fit\n",
            paste(PEAT_TIERS, collapse = "+"),
            format(sum(peat_flag_gpm), big.mark = ","),
            100 * mean(peat_flag_gpm)))
cat("  Cross-check GPM-peat vs modal-Histosol (D5):\n")
print(table(GPM_peat = peat_flag_gpm, Histosol = is_histosol))

flag_file <- file.path(PATHS$output_dir, "08b_peat_flag.rds")
saveRDS(peat_flags, flag_file)
write.csv(peat_flags, sub("\\.rds$", ".csv", flag_file), row.names = FALSE)
cat(sprintf("  -> Saved: %s\n", basename(flag_file)))

# =============================================================================
# Part B -- 0.1deg peat/water fraction rasters + combined exclusion mask (maps)
# =============================================================================
log_step("Part B: building 0.1deg peat/water fractions and the exclusion mask")

if (!file.exists(TEMPLATE)) stop("Template prediction raster not found: ", TEMPLATE)
template <- rast(TEMPLATE)

# --- peat fraction: presence (NoData/other -> 0), aggregate to ~0.1deg, snap to grid
# GPM NoData (255) covers BOTH ocean and non-peat land, so it must map to 0 (not
# peat), otherwise a mean over observed-only pixels would read ~1 everywhere peat
# occurs. Force NA -> 0 first, then mark the excluded tiers as present.
gpm0     <- ifel(is.na(gpm), 0L, gpm)
peat_bin <- subst(gpm0, from = PEAT_TIERS, to = rep(1L, length(PEAT_TIERS)), others = 0L)
cat("  Aggregating GPM presence to ~0.1deg (mean = peat fraction)...\n")
peat_agg  <- aggregate(peat_bin, fact = 10, fun = "mean", na.rm = TRUE)  # ~0.0998deg
peat_frac <- resample(peat_agg, template, method = "average")
names(peat_frac) <- "peat_fraction"

# --- water fraction: ESA land-cover band 8, snapped to the template grid
if (!file.exists(LC_FRAC)) stop("Landcover fraction stack not found: ", LC_FRAC)
water_raw  <- rast(LC_FRAC, lyrs = LC_WATER_BAND)
water_frac <- if (compareGeom(water_raw, template, stopOnError = FALSE)) {
  water_raw
} else {
  resample(water_raw, template, method = "average")
}
names(water_frac) <- "water_fraction"

# --- combined exclusion mask (1 = exclude), defined over the land template only
exclusion <- (peat_frac >= THRESH) | (water_frac >= THRESH)
exclusion <- ifel(is.na(template), NA, exclusion)   # keep only land cells
names(exclusion) <- "exclusion_mask"

writeRaster(peat_frac,  file.path(PRED_DIR, "peat_frac_0p1deg.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(water_frac, file.path(PRED_DIR, "water_frac_0p1deg.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(exclusion,  file.path(PRED_DIR, "exclusion_mask_0p1deg.tif"),
            overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW"))
cat(sprintf("  -> Saved: peat_frac / water_frac / exclusion_mask (0.1deg) in %s\n", PRED_DIR))

# =============================================================================
# SANITY GATE -- do not proceed until this looks right
# =============================================================================
log_step("Sanity gate: area check + regional overlays")

# (1) total peat area vs GPM's stated ~5.8 M km^2 (both tiers)
cell_km2 <- cellSize(template, unit = "km")
peat_km2  <- global(peat_frac  * cell_km2, "sum", na.rm = TRUE)[1, 1]
water_km2 <- global(water_frac * cell_km2, "sum", na.rm = TRUE)[1, 1]
excl_cells <- global(exclusion, "sum", na.rm = TRUE)[1, 1]
cat(sprintf("  Peat area (fraction-weighted): %.2f M km^2  (GPM reference ~5.8 M km^2)\n",
            peat_km2 / 1e6))
cat(sprintf("  Water area (fraction-weighted): %.2f M km^2\n", water_km2 / 1e6))
cat(sprintf("  Land cells masked (peat OR water >= %.2f): %s\n",
            THRESH, format(excl_cells, big.mark = ",")))

# (2) regional overlays against known peat regions
regions <- list(
  "Congo Cuvette Centrale" = ext(14, 26, -5, 4),
  "West Siberian Lowlands" = ext(60, 90, 55, 70),
  "Hudson Bay Lowlands"    = ext(-96, -74, 50, 60),
  "Fennoscandia"           = ext(4, 34, 58, 70),
  "SE Asian coastal peat"  = ext(95, 120, -5, 8)
)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
overlay_png <- file.path(PLOT_DIR, "08b_peat_sanity_overlay.png")
png(overlay_png, width = 1800, height = 1200, res = 140)
op <- par(mfrow = c(2, 3), mar = c(2, 2, 3, 4))
for (nm in names(regions)) {
  plot(crop(peat_frac, regions[[nm]]),
       main = nm, col = hcl.colors(20, "YlGnBu", rev = TRUE),
       range = c(0, 1), axes = TRUE)
}
# 6th panel: global exclusion mask
plot(exclusion, main = "Global exclusion mask (1 = peat/water)",
     col = c("grey90", "#B03030"), axes = TRUE, legend = FALSE)
par(op); dev.off()
cat(sprintf("  -> Overlay for review: %s\n", overlay_png))

cat("\n")
cat("===================================================================\n")
cat("  Step 08b complete. REVIEW the overlay + area check before Phase 2.\n")
cat("  Fit flag:  outputs/08b_peat_flag.rds  (join to 13c by dsiteid)\n")
cat("  Map masks: outputs/MRT_predictions/{peat_frac,water_frac,exclusion_mask}_0p1deg.tif\n")
cat("===================================================================\n\n")
