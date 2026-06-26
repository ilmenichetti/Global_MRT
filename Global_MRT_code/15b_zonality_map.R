################################################################################
# 15b_zonality_map.R
#
# The zonality map (framing A): where does NON-CLIMATE context push soil-carbon
# turnover longer or shorter than the simpler model would predict? Maps are the
# log-ratio of nested RF model predictions, rendered at MESO scale (2 deg)
# because the residual is ~2/3 noise with a ~16 km range (13f) -- native-pixel
# texture is not interpretable.
#
# MAIN (manuscript) = SYMMETRIC formulation -- each domain relative to climate,
#   no abiotic-before-biology ordering (consistent with the order-independent
#   Shapley attribution and the "beyond-climate" framing):
#     (a) abiotic   = log(M5 / M1) = log(tau_{clim+edaphic+landuse} / tau_{clim})
#     (b) biology   = log(M4 / M1) = log(tau_{clim+biology}        / tau_{clim})
# APPENDIX = NESTED formulation (additive but order-dependent):
#     (a) abiotic   = log(M5 / M1)
#     (b) biology   = log(M7 / M5)  -- biology on top of abiotic (its minimum)
#
# Each dual figure carries a 3rd panel: the abiotic-vs-biology scatter with the
# spatial correlation (native + meso). Interpretation: the SYMMETRIC pair
# co-varies (r~0.4-0.5; both tap the shared beyond-climate signal) while the
# NESTED pair is ~orthogonal (biology conditional on abiotic = its unique axis).
# Together this is the "correlated, not redundant" result, spatially. The
# correlations are computed here (not by hand) and saved for repeatability.
#
# Input:  outputs/MRT_predictions/{zonality_modulation_M5_M1_logratio,
#         MRT_M1_climate, MRT_M4_climate_biological,
#         MRT_M5_climate_edaphic_landuse, MRT_M7_full}.tif
# Output: plots/step_15b_zonality/zonality_map_abiotic.png          (headline, near-global)
#         plots/step_15b_zonality/zonality_map_dual_symmetric.png   (MAIN; a,b,c)
#         plots/step_15b_zonality/zonality_map_dual.png             (APPENDIX nested; a,b,c)
#         outputs/15b_modulation_correlations.csv
#
# Author: Lorenzo   Date: 2026-06-26
################################################################################

suppressMessages({
  library(terra)
  library(rnaturalearth)
  library(sf)
})

PRED_DIR     <- "./Global_MRT_code/outputs/MRT_predictions"
OUTPUTS_ROOT <- "./Global_MRT_code/outputs"
PLOT_DIR     <- "./Global_MRT_code/plots/step_15b_zonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

AGG_DEG <- 2          # meso aggregation (2 deg); native 0.1 deg -> fact 20
CLAMP_A <- 0.6        # colour limit for log(M5/M1) and log(M4/M1)
CLAMP_B <- 0.25       # colour limit for the smaller-magnitude log(M7/M5)

pal_abiotic <- hcl.colors(100, "Blue-Red 3")    # neg = blue (shorter), pos = red (longer)
pal_biology <- hcl.colors(100, "Purple-Green")  # neg = purple (shorter), pos = green (longer)

borders  <- st_geometry(ne_countries(scale = 50, returnclass = "sf"))
meso_agg <- function(r) aggregate(r, fact = round(AGG_DEG / 0.1), fun = "mean", na.rm = TRUE)

plot_modulation <- function(r, main, sub, pal, clampval, legtitle) {
  rc <- clamp(r, -clampval, clampval, values = TRUE)
  plot(rc, col = pal, range = c(-clampval, clampval),
       main = main, mar = c(2.4, 2.2, 2.6, 4.8),
       plg = list(title = legtitle, title.cex = 0.8))
  plot(borders, add = TRUE, col = NA, border = "grey45", lwd = 0.3)
  mtext(sub, side = 1, line = 0.5, cex = 0.72, col = "grey30")
}

# Köppen main zones (for colouring the scatter): 1=A Tropical, 2=B Arid,
# 3=C Temperate, 4=D Continental, 5=E Polar.
KOP_COL <- c("1" = "#2c7fb8", "2" = "#d95f02", "3" = "#1b9e77",
             "4" = "#7570b3", "5" = "#737373")
KOP_LAB <- c("1" = "A Tropical", "2" = "B Arid", "3" = "C Temperate",
             "4" = "D Continental", "5" = "E Polar")

# (c) abiotic-vs-biology scatter at meso scale, full-width (matches the map
# panels), points coloured by Köppen zone, annotated with native + meso r.
scatter_panel <- function(x_meso, y_meso, kop_meso, r_nat, r_meso, xlab, ylab) {
  df <- data.frame(x = values(x_meso)[, 1], y = values(y_meso)[, 1],
                   k = values(kop_meso)[, 1])
  df <- df[!is.na(df$x) & !is.na(df$y), ]
  pcol <- KOP_COL[as.character(df$k)]; pcol[is.na(pcol)] <- "#cccccc"
  op <- par(mar = c(4.2, 4.4, 2.6, 2)); on.exit(par(op))   # no pty="s": fill width
  plot(df$x, df$y, pch = 16, col = adjustcolor(pcol, 0.4), cex = 0.5,
       xlab = xlab, ylab = ylab,
       main = "(c) Do the two axes co-vary?  (points by Köppen zone)")
  abline(h = 0, v = 0, col = "grey75"); abline(lm(df$y ~ df$x), col = "black", lwd = 2)
  legend("topleft", bty = "n", cex = 0.95,
         legend = c(sprintf("r (2° meso) = %.2f", r_meso),
                    sprintf("r (native)  = %.2f", r_nat)))
  legend("bottomright", bty = "n", cex = 0.85, pch = 16, horiz = FALSE,
         col = KOP_COL, legend = KOP_LAB, title = "Köppen main")
}

corpair <- function(x, y) { v <- na.omit(cbind(values(x), values(y))); cor(v[, 1], v[, 2]) }

# ---- Build layers (native + meso) -------------------------------------------
abio_nat  <- rast(file.path(PRED_DIR, "zonality_modulation_M5_M1_logratio.tif"))  # log M5/M1
m1 <- rast(file.path(PRED_DIR, "MRT_M1_climate.tif"))
m4 <- rast(file.path(PRED_DIR, "MRT_M4_climate_biological.tif"))
m5 <- rast(file.path(PRED_DIR, "MRT_M5_climate_edaphic_landuse.tif"))
m7 <- rast(file.path(PRED_DIR, "MRT_M7_full.tif"))
bnest_nat <- log(m7 / m5)   # nested:    biology on top of abiotic
bclim_nat <- log(m4 / m1)   # symmetric: biology relative to climate

abio  <- meso_agg(abio_nat)
bnest <- meso_agg(bnest_nat)
bclim <- meso_agg(bclim_nat)

# Köppen main zone per meso cell (majority class), on the same 2° grid
kop <- aggregate(rast("./Global_MRT_code/spatialized_layers/climate/koppen_main_0.1deg.tif"),
                 fact = round(AGG_DEG / 0.1), fun = "modal", na.rm = TRUE)

# ---- Correlations (computed in-code, saved for repeatability) ----------------
cors <- data.frame(
  formulation = c("symmetric (a=M5/M1, b=M4/M1)", "nested (a=M5/M1, b=M7/M5)"),
  r_native    = c(corpair(abio_nat, bclim_nat), corpair(abio_nat, bnest_nat)),
  r_meso_2deg = c(corpair(abio, bclim),         corpair(abio, bnest))
)
write.csv(cors, file.path(OUTPUTS_ROOT, "15b_modulation_correlations.csv"), row.names = FALSE)
cat("Abiotic-vs-biology modulation correlations:\n"); print(cors, row.names = FALSE)

# ---- (1) Headline: abiotic modulation, near-global --------------------------
png(file.path(PLOT_DIR, "zonality_map_abiotic.png"), width = 2200, height = 1200, res = 200)
plot_modulation(abio,
  main = "Non-climate (abiotic) modulation of soil-carbon turnover",
  sub  = "log( turnover from climate+edaphic+land-use  /  turnover from climate alone ).  Red = longer turnover, blue = shorter.",
  pal = pal_abiotic, clampval = CLAMP_A, legtitle = "log(M5 / M1)")
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_map_abiotic.png"), "\n")

# ---- (2) MAIN: symmetric dual (each domain over climate) + scatter ----------
png(file.path(PLOT_DIR, "zonality_map_dual_symmetric.png"), width = 2000, height = 2700, res = 200)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(1, 1, 1.05))
plot_modulation(abio,
  main = "(a) Abiotic effect relative to climate (edaphic + land-use)",
  sub  = "log( turnover from climate+edaphic+land-use  /  turnover from climate ).  Red = longer, blue = shorter.",
  pal = pal_abiotic, clampval = CLAMP_A, legtitle = "log(M5 / M1)")
plot_modulation(bclim,
  main = "(b) Biological effect relative to climate",
  sub  = "log( turnover from climate+biology  /  turnover from climate ).  Green = longer, purple = shorter.",
  pal = pal_biology, clampval = CLAMP_A, legtitle = "log(M4 / M1)")
scatter_panel(abio, bclim, kop, cors$r_native[1], cors$r_meso_2deg[1],
  xlab = "Abiotic modulation  log(M5/M1)", ylab = "Biological modulation  log(M4/M1)")
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_map_dual_symmetric.png"), "\n")

# ---- (3) APPENDIX: nested dual (biology on top of abiotic) + scatter ---------
png(file.path(PLOT_DIR, "zonality_map_dual.png"), width = 2000, height = 2700, res = 200)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(1, 1, 1.05))
plot_modulation(abio,
  main = "(a) Abiotic modulation: edaphic + land-use added to climate",
  sub  = "log( turnover from climate+edaphic+land-use  /  turnover from climate ).  Red = longer, blue = shorter.",
  pal = pal_abiotic, clampval = CLAMP_A, legtitle = "log(M5 / M1)")
plot_modulation(bnest,
  main = "(b) Biological modulation: biology added on top of the abiotic model",
  sub  = "log( turnover from full model  /  turnover from climate+edaphic+land-use ).  Green = longer, purple = shorter.",
  pal = pal_biology, clampval = CLAMP_B, legtitle = "log(M7 / M5)")
scatter_panel(abio, bnest, kop, cors$r_native[2], cors$r_meso_2deg[2],
  xlab = "Abiotic modulation  log(M5/M1)", ylab = "Biological modulation  log(M7/M5)")
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_map_dual.png"), "\n")

cat(sprintf("\nDone. Meso = %g deg. MAIN = symmetric; APPENDIX = nested.\n", AGG_DEG))
