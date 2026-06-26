################################################################################
# 13f_scale_noisefloor.R
#
# Scale spectrum & noise floor of the beyond-climate residual (appendix).
# House-style reconstruction of the (now-deleted) scratchpad `scale_noisefloor.R`
# so the appendix figure/numbers are regenerable.
#
# Three diagnostics on the beyond-climate residual (observed log(MRT) - climate-
# only spatial-block-CV prediction; same recipe/seed as 13e):
#   (a) Empirical variogram -> nugget/sill ratio (unstructured-noise fraction)
#       and the spatial range. Expectation: ~2/3 nugget, near-flat beyond ~15 km.
#   (b) Multiscale organisation: observed between-block variance / permutation
#       null (obs/null ratio) at 1-20 deg -> coarse structure grows with scale,
#       but the absolute coarse signal is a small share of total variance.
#   (c) Block-mean covariate skill vs scale: linear R^2 of block-mean covariates
#       on the block-mean residual at 1-20 deg -> peaks at a meso scale (~3-4 deg),
#       far above the native-point skill. Motivates building maps at meso scale.
#
# Input:  outputs/12b_model_ready.rds, outputs/13_var_groups.rds
# Output: outputs/13f_scale_noisefloor.csv, outputs/13f_scale_noisefloor.rds
#         plots/appendix/scale_noisefloor.png
#
# Env override: MRT_13F_NTREES (default 500).
# Author: Lorenzo   Date: 2026-06-26
################################################################################

library(dplyr)
library(ranger)
library(gstat)
library(sp)
library(ggplot2)
library(patchwork)

OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR   <- "./Global_MRT_code/plots/appendix"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

RF_NTREES   <- as.integer(Sys.getenv("MRT_13F_NTREES", "500"))
RF_THREADS  <- 3
CV_FOLDS    <- 10
BLOCK_CV    <- 2
MIN_BLOCK_N <- 5
NPERM       <- 199
VG_NSAMP    <- 20000       # subsample for the variogram (full n too large)
set.seed(42)

BIO_MAIN <- c("fungal_proportion", "AM_richness", "EcM_richness",
              "AM_endemism", "EcM_endemism", "EcM_AM_richness_ratio")

cat("===============================================================\n")
cat("  SCALE SPECTRUM & NOISE FLOOR (Step 13f) -", RF_NTREES, "trees\n")
cat("===============================================================\n\n")

# ---- Data prep + residual (same as 13e/13g) ---------------------------------
d <- readRDS(file.path(OUTPUT_DIR, "12b_model_ready.rds"))
g <- readRDS(file.path(OUTPUT_DIR, "13_var_groups.rds"))
CLIM <- g$CLIMATE
COVS <- c(g$CLIMATE, g$EDAPHIC, g$LANDUSE, BIO_MAIN)

d <- d %>% filter(MRT_QC == "valid", !is.na(MRT_years), MRT_years > 0, MRT_years < Inf)
q <- quantile(d$MRT_years, c(.01, .99), na.rm = TRUE)
d <- d %>% filter(MRT_years > q[1] & MRT_years < q[2])
d$log_MRT <- log(d$MRT_years)

need <- c("log_MRT", COVS, "longitude_decimal_degrees", "latitude_decimal_degrees")
cv <- d[complete.cases(d[, need]), ]
cv$lon <- cv$longitude_decimal_degrees; cv$lat <- cv$latitude_decimal_degrees
cat("MAIN common-row n:", format(nrow(cv), big.mark = ","), "\n")

cv$grid <- paste(floor(cv$lon / BLOCK_CV), floor(cv$lat / BLOCK_CV))
ug <- unique(cv$grid); set.seed(42)
cv$fold <- data.frame(grid = ug, f = sample(CV_FOLDS, length(ug), TRUE))$f[match(cv$grid, ug)]
fclim <- as.formula(paste("log_MRT ~", paste(CLIM, collapse = " + ")))
cv$clim_pred <- NA_real_
cat("Fitting climate-only model for OOS residuals...\n")
for (k in sort(unique(cv$fold))) {
  te <- cv$fold == k
  m <- ranger(fclim, cv[!te, ], num.trees = RF_NTREES, num.threads = RF_THREADS, seed = 42)
  cv$clim_pred[te] <- predict(m, cv[te, ])$predictions
}
cv$resid <- cv$log_MRT - cv$clim_pred
resid_var <- var(cv$resid)
cat("Climate-only block-CV R2 =",
    round(1 - sum(cv$resid^2) / sum((cv$log_MRT - mean(cv$log_MRT))^2), 4),
    "| residual variance =", round(resid_var, 4), "\n\n")

# =============================================================================
# (a) VARIOGRAM -> nugget/sill ratio + range
# =============================================================================
cat("Computing empirical variogram (n =", VG_NSAMP, "subsample)...\n")
set.seed(42)
samp <- cv[sample(nrow(cv), min(VG_NSAMP, nrow(cv))), c("lon", "lat", "resid")]
coordinates(samp) <- ~lon + lat
proj4string(samp) <- CRS("+proj=longlat +datum=WGS84")   # great-circle distances (km)
# Fine short-range binning: the structured part lives < ~15 km, so coarse bins
# cannot resolve the nugget. Use 4 km bins out to 150 km.
vg <- variogram(resid ~ 1, samp, cutoff = 150, width = 4)

# Model-free nugget: extrapolate the first few short-lag bins to distance 0
# (robust to fit.variogram non-convergence). Sill ~ total residual variance.
lo <- head(vg[order(vg$dist), ], 4)
nugget   <- max(0, unname(coef(lm(gamma ~ dist, data = lo))[1]))
sill     <- resid_var
nug_sill <- nugget / sill

# Fitted model (for the plotted line + range estimate); fallback if no convergence
vfit <- tryCatch(
  fit.variogram(vg, vgm(psill = sill - nugget, model = "Exp",
                        range = 20, nugget = nugget)),
  warning = function(w) NULL, error = function(e) NULL)
rng <- if (!is.null(vfit)) vfit$range[vfit$model != "Nug"][1] else NA
psill <- sill - nugget
cat(sprintf("  nugget = %.3f, sill = %.3f, nugget/sill = %.2f, range ≈ %s km\n\n",
            nugget, sill, nug_sill, ifelse(is.na(rng), "n/a", round(rng))))

# =============================================================================
# (b) MULTISCALE ORGANISATION (obs/null) and (c) BLOCK-MEAN COVARIATE R^2
# =============================================================================
block_mean_var <- function(r, blk) {
  bm <- tapply(r, blk, function(x) if (length(x) >= MIN_BLOCK_N) mean(x) else NA_real_)
  var(bm, na.rm = TRUE)
}
org_ratio <- function(deg) {
  blk <- paste(floor(cv$lon / deg), floor(cv$lat / deg))
  obs  <- block_mean_var(cv$resid, blk)
  null <- mean(replicate(NPERM, block_mean_var(sample(cv$resid), blk)))
  obs / null
}
# Out-of-sample (5-fold-CV) R^2 of block-mean covariates on the block-mean
# residual. CV is essential here: at coarse scales the block count falls toward
# the predictor count, so an in-sample R^2 overfits and rises spuriously. The
# CV R^2 instead reveals the scale at which covariates genuinely track the
# block-mean (zonal) residual.
block_cov_R2 <- function(deg, kf = 5) {
  df <- cv[, COVS]; df$r <- cv$resid
  df$blk <- paste(floor(cv$lon / deg), floor(cv$lat / deg))
  agg <- df %>% group_by(blk) %>% filter(n() >= MIN_BLOCK_N) %>%
    summarise(across(all_of(c("r", COVS)), mean), .groups = "drop")
  n <- nrow(agg)
  if (n < 40) return(list(R2 = NA_real_, n_blocks = n))   # too few blocks
  set.seed(42); fold <- sample(rep_len(seq_len(kf), n))
  pred <- rep(NA_real_, n)
  for (k in seq_len(kf)) {
    te <- fold == k
    m  <- lm(reformulate(COVS, response = "r"), data = agg[!te, ])
    pred[te] <- predict(m, agg[te, ])
  }
  list(R2 = 1 - sum((agg$r - pred)^2) / sum((agg$r - mean(agg$r))^2),
       n_blocks = n)
}

scales <- c(1, 2, 3, 4, 5, 7, 10, 15, 20)
cat("Scanning scales (organisation + covariate skill)...\n")
spectrum <- do.call(rbind, lapply(scales, function(s) {
  bc <- block_cov_R2(s)
  data.frame(scale_deg = s, org_obs_null = org_ratio(s),
             block_cov_R2 = bc$R2, n_blocks = bc$n_blocks)
}))
print(transform(spectrum, org_obs_null = round(org_obs_null, 2),
                block_cov_R2 = round(block_cov_R2, 3)), row.names = FALSE)

write.csv(spectrum, file.path(OUTPUT_DIR, "13f_scale_noisefloor.csv"), row.names = FALSE)
saveRDS(list(
  config = list(rf_ntrees = RF_NTREES, n_common = nrow(cv), nperm = NPERM,
                vg_nsamp = VG_NSAMP, seed = 42, date = Sys.Date()),
  variogram = list(empirical = vg, fit = vfit, nugget = nugget, psill = psill,
                   range_km = rng, nugget_sill_ratio = nug_sill,
                   resid_var = resid_var),
  spectrum = spectrum),
  file.path(OUTPUT_DIR, "13f_scale_noisefloor.rds"))

# =============================================================================
# FIGURE (3 panels)
# =============================================================================
pa <- ggplot(vg, aes(dist, gamma)) +
  geom_point(size = 1.6) +
  { if (!is.null(vfit)) geom_line(data = variogramLine(vfit, maxdist = max(vg$dist)),
                                  aes(dist, gamma), colour = "#B2182B") } +
  geom_hline(yintercept = resid_var, linetype = "dotted", colour = "grey50") +
  labs(tag = "a", title = "Residual variogram",
       subtitle = sprintf("nugget/sill = %.2f, range ≈ %.0f km (dotted = total variance)",
                          nug_sill, rng),
       x = "Distance (km)", y = "Semivariance") +
  theme_bw(base_size = 11) + theme(plot.tag = element_text(face = "bold"))

pb <- ggplot(spectrum, aes(scale_deg, org_obs_null)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_line() + geom_point() +
  labs(tag = "b", title = "Coarse organisation grows with scale",
       subtitle = "Between-block variance / spatial-randomness null",
       x = "Block size (degrees)", y = "Observed / null ratio") +
  theme_bw(base_size = 11) + theme(plot.tag = element_text(face = "bold"))

peak_scale <- spectrum$scale_deg[which.max(spectrum$block_cov_R2)]
pc <- ggplot(spectrum, aes(scale_deg, block_cov_R2)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line() + geom_point() +
  geom_vline(xintercept = peak_scale, linetype = "dashed", colour = "#377EB8") +
  labs(tag = "c", title = "Covariate skill on block means vs scale",
       subtitle = sprintf("Out-of-sample (5-fold-CV) R²; peaks near %g° (coarse scales block-limited)",
                          peak_scale),
       x = "Block size (degrees)", y = expression("Block-mean CV " * R^2)) +
  theme_bw(base_size = 11) + theme(plot.tag = element_text(face = "bold"))

fig <- (pa | pb | pc) +
  plot_annotation(title = "Scale spectrum & noise floor of the beyond-climate residual",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(PLOT_DIR, "scale_noisefloor.png"), fig, width = 14, height = 4.6, dpi = 200)
cat("\nOK  ", file.path(PLOT_DIR, "scale_noisefloor.png"), "\n")
cat("OK  outputs/13f_scale_noisefloor.csv + .rds\n")
cat(sprintf("\nVerdict: ~%.0f%% of the beyond-climate residual is unstructured noise (nugget),\n",
            100 * nug_sill))
cat("structured only at short range (~15 km). Coarse organisation (obs/null) grows\n")
cat("with scale, but the OUT-OF-SAMPLE linear skill of block-mean covariates on the\n")
cat(sprintf("block-mean residual is weak (max CV R2 ~ %.2f, near %g deg) and turns negative\n",
            max(spectrum$block_cov_R2, na.rm = TRUE), peak_scale))
cat("at coarse scales (overfitting; few blocks). NOTE: this lm proxy is not the RF\n")
cat("mapping skill; the meso-scale mapping recommendation must be re-assessed with an\n")
cat("RF-based, scale-resolved skill metric before it is asserted.\n")
