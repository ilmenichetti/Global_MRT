################################################################################
# 15b_zonality_map.R
#
# The zonality map (framing A): where does NON-CLIMATE context push soil-carbon
# turnover longer or shorter than climate alone would predict? Two layers:
#   (a) ABIOTIC modulation  = log(M5 / M1) = log(tau_climate+edaphic+landuse /
#       tau_climate) -- near-global (soil-bearing land); the workhorse.
#   (b) BIOLOGICAL enrichment = biological-group log-ratio contribution -- the
#       partial (~29%) overlay where microbial layers exist.
#
# Rendered at MESO scale (2 deg aggregation): the residual is ~2/3 unstructured
# noise with a ~16 km range (see 13f), so native-pixel texture is not
# interpretable; the organised signal lives at coarse scale. Coverage is honest:
# ice sheets / Antarctica are blank (no soil), and the boreal high latitudes are
# sparse (SoilGrids), which we state rather than hide.
#
# Input:  outputs/MRT_predictions/zonality_modulation_M5_M1_logratio.tif
#         outputs/MRT_predictions/contribution_BIOLOGICAL_logratio.tif
# Output: plots/step_15b_zonality/zonality_map_abiotic.png
#         plots/step_15b_zonality/zonality_map_dual.png
#
# Author: Lorenzo   Date: 2026-06-26
################################################################################

suppressMessages({
  library(terra)
  library(rnaturalearth)
  library(sf)
})

OUTPUT_DIR <- "./Global_MRT_code/outputs/MRT_predictions"
PLOT_DIR   <- "./Global_MRT_code/plots/step_15b_zonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

AGG_DEG  <- 2        # meso aggregation (2 deg); 0.1 deg native -> fact 20
CLAMP    <- 0.6      # symmetric colour limit on the log-ratio (~2/98th pct)

# Diverging palette: blue = non-climate SHORTENS turnover (faster), red = LENGTHENS
pal <- colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(100)

borders <- st_geometry(ne_countries(scale = 50, returnclass = "sf"))

meso <- function(f) {
  r <- rast(file.path(OUTPUT_DIR, f))
  aggregate(r, fact = round(AGG_DEG / 0.1), fun = "mean", na.rm = TRUE)
}

plot_modulation <- function(r, main, sub) {
  rc <- clamp(r, -CLAMP, CLAMP, values = TRUE)
  cov_pct <- round(100 * global(!is.na(r), "sum")[1, 1] /
                     (ncell(r) * 0.292))   # ~29.2% of a 0.1->2deg global grid is land
  plot(rc, col = pal, range = c(-CLAMP, CLAMP),
       main = main, mar = c(2.2, 2.2, 2.6, 3.5),
       plg = list(title = "log(tau / climate-tau)", title.cex = 0.8))
  plot(borders, add = TRUE, col = NA, border = "grey45", lwd = 0.3)
  mtext(sub, side = 1, line = 0.6, cex = 0.7, col = "grey30")
}

# ---- (a) Abiotic zonality map (near-global) ---------------------------------
abio <- meso("zonality_modulation_M5_M1_logratio.tif")
png(file.path(PLOT_DIR, "zonality_map_abiotic.png"), width = 2200, height = 1200, res = 200)
plot_modulation(abio,
  main = "Zonality of non-climate (abiotic) control on soil-carbon turnover",
  sub = paste0("log(tau_M5 / tau_M1) at ", AGG_DEG,
               "°; red = abiotic context lengthens turnover, blue = shortens.",
               " Blank = ice/Antarctica (no soil); boreal high-lat sparse (SoilGrids)."))
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_map_abiotic.png"), "\n")

# ---- (b) Dual-layer: abiotic + biological enrichment overlay ----------------
bio_f <- "contribution_BIOLOGICAL_logratio.tif"
if (file.exists(file.path(OUTPUT_DIR, bio_f))) {
  bio <- meso(bio_f)
  png(file.path(PLOT_DIR, "zonality_map_dual.png"), width = 2200, height = 2100, res = 200)
  par(mfrow = c(2, 1))
  plot_modulation(abio,
    main = "(a) Abiotic non-climate modulation  (M5 - M1; near-global)",
    sub = paste0("log(tau_M5 / tau_M1) at ", AGG_DEG, "°"))
  plot_modulation(bio,
    main = "(b) Biological enrichment  (microbial layers; ~29% of land)",
    sub = paste0("biological log-ratio contribution at ", AGG_DEG,
                 "°; partial coverage where microbial data exist"))
  dev.off()
  cat("OK  ", file.path(PLOT_DIR, "zonality_map_dual.png"), "\n")
} else {
  cat("  (biology contribution raster not found; skipped dual panel)\n")
}

cat("\nDone. Meso scale =", AGG_DEG, "deg; colour clamp = +/-", CLAMP, "log-ratio.\n")
