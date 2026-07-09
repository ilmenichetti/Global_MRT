# =============================================================================
# Step 14b: apply the peat + open-water exclusion mask to the prediction rasters
# =============================================================================
# Runs after 14 (extrapolation). Sets peat/water cells (exclusion_mask == 1) to NA
# in every MRT_* prediction raster, so all downstream maps (15/15b/17/20) and the
# manuscript number script (19) reflect the mineral-soil scope with no per-script
# path changes. Idempotent: re-masking already-NA cells is a no-op.
#
# Mask source: outputs/MRT_predictions/exclusion_mask_0p1deg.tif (built by 08b).
# =============================================================================

library(terra)

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
PRED_DIR     <- file.path(OUTPUT_DIR, "MRT_predictions")
MASK_FILE    <- file.path(PRED_DIR, "exclusion_mask_0p1deg.tif")

if (!file.exists(MASK_FILE)) stop("Exclusion mask not found: ", MASK_FILE,
                                  "\nRun 08b_extract_peat.R first.")
excl <- rast(MASK_FILE)   # 1 = exclude (peat OR water), NA over non-land

# Mask every prediction raster except the frac/mask layers themselves.
skip  <- c("peat_frac_0p1deg", "water_frac_0p1deg", "exclusion_mask_0p1deg")
files <- list.files(PRED_DIR, pattern = "\\.tif$", full.names = TRUE)
files <- files[!tools::file_path_sans_ext(basename(files)) %in% skip]

cat(sprintf("14b: masking %d prediction rasters (peat/water -> NA)...\n", length(files)))
n_excl <- global(excl == 1, "sum", na.rm = TRUE)[1, 1]
cat(sprintf("     exclusion mask flags %s land cells\n", format(n_excl, big.mark = ",")))

for (f in files) {
  r <- rast(f)
  m <- excl
  if (!compareGeom(m, r, stopOnError = FALSE, messages = FALSE))
    m <- resample(excl, r, method = "near")
  r2  <- mask(r, m, maskvalues = 1, updatevalue = NA)   # excluded cells -> NA
  tmp <- tempfile(fileext = ".tif")
  writeRaster(r2, tmp, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  file.copy(tmp, f, overwrite = TRUE); file.remove(tmp)
  cat("  masked:", basename(f), "\n")
}

# Sanity: report the new global tau summary over the masked (mineral-soil) domain
if (file.exists(file.path(PRED_DIR, "MRT_M7_full.tif"))) {
  m7  <- rast(file.path(PRED_DIR, "MRT_M7_full.tif"))
  v   <- values(m7, mat = FALSE); v <- v[!is.na(v)]
  cat(sprintf("\nMasked global tau (M7_full): mean %.1f  median %.1f  range %.1f-%.1f yr  (n=%s cells)\n",
              mean(v), median(v), min(v), max(v), format(length(v), big.mark = ",")))
}
cat("14b done: all MRT_* prediction rasters now exclude peat + open water.\n")
