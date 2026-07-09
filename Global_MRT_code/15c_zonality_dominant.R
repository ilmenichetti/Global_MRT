################################################################################
# 15c_zonality_dominant.R  --  PROTOTYPE, not wired into the manuscript
#
# Discrete "dominant modulator" zonality map: at each meso (2 deg) block, which
# non-climate domain (Edaphic / LandUse / Biological) most modulates tau relative
# to climate-only, and in which direction (longer / shorter). Three definitions
# of "dominant" are compared, because raw magnitude over-credits LandUse (whose
# large per-pixel shift is mostly the land-cover allocation-coefficient imprint,
# yet carries little explained variance):
#   A  raw          dominant = max |log ratio|            (magnitude)
#   C  standardised dominant = max |log ratio| / SD_domain (per-domain z-score)
#   D  Shapley-wtd  dominant = max |log ratio| x Shapley%  (attribution-consistent)
# The climate-dominated mask (raw |log ratio| < THRESH) and the direction (sign of
# the dominant domain's raw log ratio) are identical across methods; only the
# ranking of "which domain dominates" changes.
#
# Standalone review figure. Not referenced by manuscript.tex; output stays
# gitignored (no un-ignore line), so it is not committed until we adopt it.
#
# Input : outputs/MRT_predictions/contribution_{EDAPHIC,LANDUSE,BIOLOGICAL}_logratio.tif
#         outputs/13c_decomposition.rds   (Shapley shares, for method D)
# Output: plots/step_15b_zonality/zonality_dominant_comparison.png
#
# Author: Lorenzo   Date: 2026-07-07
################################################################################

suppressMessages({ library(terra); library(rnaturalearth); library(sf) })

PRED_DIR <- "./Global_MRT_code/outputs/MRT_predictions"
PLOT_DIR <- "./Global_MRT_code/plots/step_15b_zonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

AGG_DEG <- 2
THRESH  <- 0.05     # |log ratio| below this at 2 deg = climate already explains tau

meso_agg <- function(r) aggregate(r, fact = round(AGG_DEG / 0.1), fun = "mean", na.rm = TRUE)
borders  <- st_geometry(ne_countries(scale = 50, returnclass = "sf"))

# --- per-domain log-ratio of tau vs climate-only, aggregated to meso ------------
ed <- rast(file.path(PRED_DIR, "contribution_EDAPHIC_logratio.tif"))
lu <- rast(file.path(PRED_DIR, "contribution_LANDUSE_logratio.tif"))
bi <- rast(file.path(PRED_DIR, "contribution_BIOLOGICAL_logratio.tif"))
s  <- meso_agg(c(ed, lu, bi)); names(s) <- c("Edaphic", "LandUse", "Biological")

# --- shared pieces: climate mask + no-data + per-domain SD + Shapley weights -----
absS   <- abs(s)
maxRaw <- max(absS, na.rm = TRUE)         # NA where all three domains missing
nodata <- is.na(maxRaw)
clim   <- maxRaw < THRESH                  # climate-dominated

sds <- global(s, "sd", na.rm = TRUE)[, 1]  # per-domain spatial spread (method C)
dec <- readRDS("./Global_MRT_code/outputs/13c_decomposition.rds")
shp <- setNames(dec$shapley_df$shapley_pct, dec$shapley_df$group)
w   <- shp[c("EDAPHIC", "LANDUSE", "BIOLOGICAL")]   # Shapley weights (method D)

# ranking-magnitude stacks (all >= 0), one per method
mag_raw <- absS
mag_std <- rast(lapply(1:3, function(i) absS[[i]] / sds[i]))
mag_shp <- rast(lapply(1:3, function(i) absS[[i]] * w[i]))

# --- classify: dominant domain (by a given magnitude) x direction ---------------
classify <- function(mag) {
  mf   <- ifel(is.na(mag), -1, mag)         # so which.max ignores missing domains
  dom  <- which.max(mf)                      # 1=Edaphic, 2=LandUse, 3=Biological
  domv <- selectRange(s, dom)                # signed raw log ratio of that domain
  code <- (dom - 1) * 2 + 1 + (domv < 0)     # 1/2 E, 3/4 L, 5/6 B (long/short)
  code <- ifel(clim, 0L, code)               # 0 = climate-dominated
  code <- ifel(nodata, NA, code)
  as.int(code)
}

LAB <- c("None (climate alone)",
         "Edaphic  (longer τ)",    "Edaphic  (shorter τ)",
         "LandUse  (longer τ)",    "LandUse  (shorter τ)",
         "Biological  (longer τ)", "Biological  (shorter τ)")
PAL <- setNames(c("#BDBDBD", "#A50F15", "#FB9A99", "#238B45",
                  "#A1D99B", "#08519C", "#9ECAE1"), LAB)

set_levels <- function(code) { levels(code) <- data.frame(value = 0:6, class = LAB); code }
codes <- list(
  "A  Raw magnitude  (max |log ratio|)"                 = set_levels(classify(mag_raw)),
  "C  Standardised  (max |log ratio| / SD)"             = set_levels(classify(mag_std)),
  "D  Shapley-weighted  (max |log ratio| x Shapley %)"  = set_levels(classify(mag_shp)))

# --- plot: 3 stacked panels, one shared legend ---------------------------------
png(file.path(PLOT_DIR, "zonality_dominant_comparison.png"),
    width = 2000, height = 2650, res = 190)
layout(matrix(c(1, 2, 3, 4), ncol = 1), heights = c(1, 1, 1, 0.22))
for (nm in names(codes)) {
  cd  <- codes[[nm]]
  pres <- levels(cd)[[1]]$class %in% unique(values(cd, na.rm = TRUE))  # not used; keep all
  par(mar = c(1.6, 2.0, 2.4, 0.5))
  plot(cd, col = unname(PAL), type = "classes", legend = FALSE,
       axes = TRUE, main = nm, cex.main = 1.05)
  plot(borders, add = TRUE, col = NA, border = "grey40", lwd = 0.3)
}
# shared legend strip
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = LAB, fill = unname(PAL), border = "grey60",
       ncol = 4, bty = "n", cex = 1.0, title = "Dominant non-climate modulator of τ")
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_dominant_comparison.png"), "\n\n")

# --- 4-class shares (method D, domain only): console diagnostic only ------------
# The standalone clean 4-class map was superseded by the combined 3-panel Fig 3
# (panel c uses the direction-shaded 7-class version); code4 is kept for the shares.
dom_D <- which.max(ifel(is.na(mag_shp), -1, mag_shp))   # 1=Edaphic,2=LandUse,3=Biological
code4 <- as.int(ifel(nodata, NA, ifel(clim, 0L, dom_D)))
ft4 <- freq(code4); ft4$pct <- round(100 * ft4$count / sum(ft4$count), 1)
cat("\n4-class (D) shares:\n"); print(ft4[, c("value", "pct")], row.names = FALSE)

# ============================================================================
# COMBINED 3-panel MANUSCRIPT figure (Fig 3), single image so the panels align.
#   (a) abiotic modulation over climate    = log(M5/M1), meso   [Blue-Red 3]
#   (b) biological modulation over climate = log(M4/M1), meso   [Purple-Green]
#   (c) discrete dominant BEYOND-climate modulator (method D), shaded by
#       direction (darker = longer tau, lighter = shorter), grey = climate alone.
# Every panel uses the SAME right margin (MAR) so the map boxes line up exactly.
#   -> plots/step_15b_zonality/zonality_three_panel.png
# ============================================================================
CLAMP       <- 0.6
pal_abiotic <- hcl.colors(100, "Blue-Red 3")     # red = longer tau, blue = shorter
pal_biology <- hcl.colors(100, "Purple-Green")   # green = longer, purple = shorter
MAR         <- c(3.4, 3.4, 4.2, 13)              # identical on all panels -> aligned
CEX_MAIN <- 2.4; CEX_AXIS <- 1.9; CEX_LEG <- 1.8

abio  <- meso_agg(rast(file.path(PRED_DIR, "zonality_modulation_M5_M1_logratio.tif")))
m1    <- rast(file.path(PRED_DIR, "MRT_M1_climate.tif"))
m4    <- rast(file.path(PRED_DIR, "MRT_M4_climate_biological.tif"))
bclim <- meso_agg(log(m4 / m1))

plot_mod <- function(r, main, pal, legtitle) {
  rc <- clamp(r, -CLAMP, CLAMP, values = TRUE)
  plot(rc, col = pal, range = c(-CLAMP, CLAMP), main = main, mar = MAR,
       cex.main = CEX_MAIN, pax = list(cex.axis = CEX_AXIS),
       plg = list(title = legtitle, title.cex = 1.7, cex = CEX_LEG))
  plot(borders, add = TRUE, col = NA, border = "grey45", lwd = 0.3)
}

code7 <- set_levels(classify(mag_shp))   # method D, 7-class (domain x direction)

png(file.path(PLOT_DIR, "zonality_three_panel.png"), width = 1900, height = 2650, res = 170)
par(mfrow = c(3, 1))
plot_mod(abio,  "(a) Abiotic effect relative to climate (edaphic + land-use)",
         pal_abiotic, "Abiotic\n(log τ ratio)")
plot_mod(bclim, "(b) Biological effect relative to climate",
         pal_biology, "Biological\n(log τ ratio)")
plot(code7, col = unname(PAL), type = "classes", mar = MAR,
     cex.main = CEX_MAIN, pax = list(cex.axis = CEX_AXIS),
     plg = list(cex = CEX_LEG, title = "Beyond-climate\nmodulator", title.cex = 1.3),
     main = "(c) Dominant beyond-climate modulator")
plot(borders, add = TRUE, col = NA, border = "grey45", lwd = 0.3)
dev.off()
cat("OK  ", file.path(PLOT_DIR, "zonality_three_panel.png"), "\n")

# --- class shares per method (console) -----------------------------------------
for (nm in names(codes)) {
  ft <- freq(codes[[nm]]); ft$pct <- round(100 * ft$count / sum(ft$count), 1)
  cat("===", nm, "===\n"); print(ft[, c("value", "pct")], row.names = FALSE)
}
