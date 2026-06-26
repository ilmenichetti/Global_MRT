################################################################################
# 15b_zonality_map.R
#
# The zonality map (framing A): where does NON-CLIMATE context push soil-carbon
# turnover longer or shorter than the simpler model would predict? Built as a
# NESTED SUCCESSION of log-ratios of RF model predictions (additive:
# log(M7/M1) = log(M5/M1) + log(M7/M5)):
#   (a) ABIOTIC modulation    = log(M5 / M1)
#       = log( tau_{climate+edaphic+landuse} / tau_{climate} )
#       -- what edaphic + land-use add on top of CLIMATE; near-global.
#   (b) BIOLOGICAL modulation = log(M7 / M5)
#       = log( tau_{full} / tau_{climate+edaphic+landuse} )
#       -- what biology adds on top of the ABIOTIC model; ~29% footprint.
#
# Each panel uses its OWN diverging palette and a colorbar titled with the exact
# quantity, so the two are never conflated. Both quantities are signed (~50%
# +/-), hence diverging (not sequential) palettes. Rendered at MESO scale (2 deg)
# because the residual is ~2/3 noise with a ~16 km range (13f): native-pixel
# texture is not interpretable.
#
# Input:  outputs/MRT_predictions/zonality_modulation_M5_M1_logratio.tif  (log M5/M1)
#         outputs/MRT_predictions/MRT_M5_climate_edaphic_landuse.tif
#         outputs/MRT_predictions/MRT_M7_full.tif
# Output: plots/step_15b_zonality/zonality_map_abiotic.png   (headline, near-global)
#         plots/step_15b_zonality/zonality_map_dual.png      (a abiotic + b biological)
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

AGG_DEG <- 2          # meso aggregation (2 deg); native 0.1 deg -> fact 20
CLAMP_A <- 0.6        # colour limit for log(M5/M1)  (~2/98th pct)
CLAMP_B <- 0.25       # colour limit for log(M7/M5)  (smaller-magnitude signal)

# Two DIVERGING palettes (both quantities are signed): each panel its own family
# so abiotic and biological modulation are never visually conflated.
pal_abiotic <- hcl.colors(100, "Blue-Red 3")    # neg = blue (shorter), pos = red (longer)
pal_biology <- hcl.colors(100, "Purple-Green")  # neg = purple (shorter), pos = green (longer)

borders <- st_geometry(ne_countries(scale = 50, returnclass = "sf"))
meso_agg <- function(r) aggregate(r, fact = round(AGG_DEG / 0.1), fun = "mean", na.rm = TRUE)

# One method for every panel: log( tau[climate+domain] / tau[baseline] ) from the
# nested RF predictions. The colorbar title names the exact quantity.
plot_modulation <- function(r, main, sub, pal, clampval, legtitle) {
  rc <- clamp(r, -clampval, clampval, values = TRUE)
  plot(rc, col = pal, range = c(-clampval, clampval),
       main = main, mar = c(2.4, 2.2, 2.6, 4.8),
       plg = list(title = legtitle, title.cex = 0.8))
  plot(borders, add = TRUE, col = NA, border = "grey45", lwd = 0.3)
  mtext(sub, side = 1, line = 0.5, cex = 0.72, col = "grey30")
}

# ---- Build the two layers (nested succession) -------------------------------
abio <- meso_agg(rast(file.path(OUTPUT_DIR, "zonality_modulation_M5_M1_logratio.tif")))
m5   <- rast(file.path(OUTPUT_DIR, "MRT_M5_climate_edaphic_landuse.tif"))
m7   <- rast(file.path(OUTPUT_DIR, "MRT_M7_full.tif"))
bio  <- meso_agg(log(m7 / m5))

# ---- (1) Headline: abiotic modulation, near-global --------------------------
png(file.path(PLOT_DIR, "zonality_map_abiotic.png"), width = 2200, height = 1200, res = 200)
plot_modulation(abio,
  main = "Non-climate (abiotic) modulation of soil-carbon turnover",
  sub  = "log( turnover from climate+edaphic+land-use  /  turnover from climate alone ).  Red = longer turnover, blue = shorter.",
  pal = pal_abiotic, clampval = CLAMP_A, legtitle = "log(M5 / M1)")
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_map_abiotic.png"), "\n")

# ---- (2) Dual layer: abiotic (over climate) + biological (over abiotic) ------
png(file.path(PLOT_DIR, "zonality_map_dual.png"), width = 2200, height = 2100, res = 200)
par(mfrow = c(2, 1))
plot_modulation(abio,
  main = "(a) Abiotic modulation: edaphic + land-use added to climate",
  sub  = "log( turnover from climate+edaphic+land-use  /  turnover from climate ).  Red = longer, blue = shorter.",
  pal = pal_abiotic, clampval = CLAMP_A, legtitle = "log(M5 / M1)")
plot_modulation(bio,
  main = "(b) Biological modulation: biology added on top of the abiotic model",
  sub  = "log( turnover from full model  /  turnover from climate+edaphic+land-use ).  Green = longer, purple = shorter.",
  pal = pal_biology, clampval = CLAMP_B, legtitle = "log(M7 / M5)")
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_map_dual.png"), "\n")

cat(sprintf("\nDone. Meso = %g deg; clamps: abiotic +/-%.2f, biological +/-%.2f.\n",
            AGG_DEG, CLAMP_A, CLAMP_B))
