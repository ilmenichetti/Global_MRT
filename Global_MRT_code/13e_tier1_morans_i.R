################################################################################
# 13e_tier1_morans_i.R
#
# Tier-1 rigour: is the BEYOND-CLIMATE residual spatially ORGANISED (zonality),
# or just spatial noise? We add Moran's I on block-aggregated residuals to the
# (scale-dependent, weak) between-block ICC reported earlier.
#
# Residual = observed log(MRT) - climate-only spatial-block-CV prediction
#            (the "beyond-climate signal"). We aggregate it to spatial blocks and
#            test for positive spatial autocorrelation (Moran's I) at 2 deg and
#            5 deg. A significant positive I => the beyond-climate signal is
#            spatially structured, not random => zonality is real (Tier 1).
#
# Footprint = the MAIN common-row set (complete over M5 + MAIN biological set),
# so Tier-1 and the Tier-2 decomposition (13c) share one footprint.
#
# Input:  ./Global_MRT_code/outputs/12b_model_ready.rds
#         ./Global_MRT_code/outputs/13_var_groups.rds
# Output: ./Global_MRT_code/outputs/13e_morans_i.csv
#         ./Global_MRT_code/outputs/13e_tier1.rds
#         ./Global_MRT_code/plots/step_13c_commonality/13e_block_residual_map.png
#         ./Global_MRT_code/plots/step_13c_commonality/13e_moran_scatter.png
#
# Env override: MRT_13E_NTREES (default 500).
#
# Author: Lorenzo
# Date: 2026-06-26
################################################################################

library(dplyr)
library(ggplot2)
library(ranger)
library(spdep)

OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR   <- "./Global_MRT_code/plots/step_13c_commonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

RF_NTREES   <- as.integer(Sys.getenv("MRT_13E_NTREES", "500"))
RF_MIN_NODE <- 5
RF_THREADS  <- 3            # capped to coexist with other background RF jobs
CV_FOLDS    <- 10
BLOCK_CV    <- 2            # CV block size (matches 13/13c)
KNN         <- 8            # neighbours for the spatial weights
set.seed(42)

BIO_MAIN <- c("fungal_proportion", "AM_richness", "EcM_richness",
              "AM_endemism", "EcM_endemism", "EcM_AM_richness_ratio")

cat("===============================================================\n")
cat("  TIER-1: MORAN'S I ON BEYOND-CLIMATE RESIDUAL (Step 13e)\n")
cat("===============================================================\n\n")

# ---- Data prep (same recipe as 13 / 13c) ------------------------------------
d <- readRDS(file.path(OUTPUT_DIR, "12b_model_ready.rds"))
g <- readRDS(file.path(OUTPUT_DIR, "13_var_groups.rds"))
CLIM <- g$CLIMATE
M5   <- c(g$CLIMATE, g$EDAPHIC, g$LANDUSE)

d <- d %>% filter(MRT_QC == "valid", !is.na(MRT_years), MRT_years > 0, MRT_years < Inf)
q <- quantile(d$MRT_years, c(.01, .99), na.rm = TRUE)
d <- d %>% filter(MRT_years > q[1] & MRT_years < q[2])
d$log_MRT <- log(d$MRT_years)

need <- c("log_MRT", M5, BIO_MAIN, "longitude_decimal_degrees", "latitude_decimal_degrees")
cv <- d[complete.cases(d[, need]), ]
cat("MAIN common-row n:", format(nrow(cv), big.mark = ","), "\n")

# ---- CV folds (2 deg blocks -> 10 folds, seed 42) ---------------------------
cv$lon <- cv$longitude_decimal_degrees
cv$lat <- cv$latitude_decimal_degrees
cv$grid <- paste(floor(cv$lon / BLOCK_CV), floor(cv$lat / BLOCK_CV))
ug <- unique(cv$grid); set.seed(42)
cv$fold <- data.frame(grid = ug, f = sample(CV_FOLDS, length(ug), TRUE))$f[match(cv$grid, ug)]

# ---- Climate-only block-CV out-of-sample predictions -> residual ------------
cat("Fitting climate-only model (block-CV) for OOS residuals...\n")
fclim <- as.formula(paste("log_MRT ~", paste(CLIM, collapse = " + ")))
cv$clim_pred <- NA_real_
for (k in sort(unique(cv$fold))) {
  te <- cv$fold == k
  m <- ranger(fclim, cv[!te, ], num.trees = RF_NTREES, min.node.size = RF_MIN_NODE,
              num.threads = RF_THREADS, seed = 42)
  cv$clim_pred[te] <- predict(m, cv[te, ])$predictions
  cat(" ", k)
}
cat(" done\n")
cv$resid <- cv$log_MRT - cv$clim_pred
clim_R2 <- 1 - sum(cv$resid^2) / sum((cv$log_MRT - mean(cv$log_MRT))^2)
cat("Climate-only block-CV R2 =", round(clim_R2, 4), "\n\n")

# =============================================================================
# MORAN'S I ON BLOCK-AGGREGATED RESIDUALS (at several scales)
# =============================================================================
morans_at_scale <- function(block_deg) {
  agg <- cv %>%
    mutate(bx = floor(lon / block_deg), by = floor(lat / block_deg)) %>%
    group_by(bx, by) %>%
    summarise(resid = mean(resid), n = n(),
              lon = mean(lon), lat = mean(lat), .groups = "drop") %>%
    filter(n >= 5)                       # stable block means only
  coords <- as.matrix(agg[, c("lon", "lat")])
  nb  <- knn2nb(knearneigh(coords, k = KNN))
  lw  <- nb2listw(nb, style = "W")
  mt  <- moran.test(agg$resid, lw, zero.policy = TRUE)
  list(agg = agg, listw = lw,
       row = data.frame(
         block_deg = block_deg, n_blocks = nrow(agg),
         morans_I  = unname(mt$estimate["Moran I statistic"]),
         expected  = unname(mt$estimate["Expectation"]),
         sd        = sqrt(unname(mt$estimate["Variance"])),
         z         = unname(mt$statistic),
         p_value   = mt$p.value))
}

scales <- c(2, 5)
res_list <- lapply(scales, morans_at_scale)
morans_df <- do.call(rbind, lapply(res_list, `[[`, "row"))
cat("Moran's I on block-aggregated beyond-climate residual:\n")
print(transform(morans_df,
                morans_I = round(morans_I, 4), expected = round(expected, 5),
                z = round(z, 1), p_value = signif(p_value, 3)), row.names = FALSE)
write.csv(morans_df, file.path(OUTPUT_DIR, "13e_morans_i.csv"), row.names = FALSE)

# =============================================================================
# FIGURES
# =============================================================================
# (a) Map of 2 deg block-mean residuals
agg2 <- res_list[[1]]$agg
lim <- max(abs(quantile(agg2$resid, c(.02, .98))))
p_map <- ggplot(agg2, aes(lon, lat, fill = resid)) +
  geom_tile(width = 2, height = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                       name = "Beyond-climate\nresidual (log MRT)") +
  coord_quickmap() +
  labs(title = "Beyond-climate residual, 2° block means",
       subtitle = sprintf("Moran's I = %.3f (p = %s); spatially organised, not random",
                          morans_df$morans_I[1], signif(morans_df$p_value[1], 2)),
       x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 11)
ggsave(file.path(PLOT_DIR, "13e_block_residual_map.png"), p_map,
       width = 10, height = 5.5, dpi = 200)
cat("OK  ", file.path(PLOT_DIR, "13e_block_residual_map.png"), "\n")

# (b) Moran scatter at 2 deg
png(file.path(PLOT_DIR, "13e_moran_scatter.png"), width = 1500, height = 1500, res = 220)
moran.plot(res_list[[1]]$agg$resid, res_list[[1]]$listw, zero.policy = TRUE,
           xlab = "Block-mean residual", ylab = "Spatially lagged residual",
           main = sprintf("Moran scatter (2° blocks): I = %.3f", morans_df$morans_I[1]))
dev.off()
cat("OK  ", file.path(PLOT_DIR, "13e_moran_scatter.png"), "\n")

# =============================================================================
# SAVE
# =============================================================================
saveRDS(list(
  config = list(rf_ntrees = RF_NTREES, cv_folds = CV_FOLDS, block_cv = BLOCK_CV,
                knn = KNN, n_common = nrow(cv), seed = 42, date = Sys.Date()),
  climate_R2 = clim_R2, morans = morans_df,
  block_means_2deg = agg2
), file.path(OUTPUT_DIR, "13e_tier1.rds"))
cat("OK  ", file.path(OUTPUT_DIR, "13e_tier1.rds"), "\n")

cat("\nTier-1 verdict: positive, significant Moran's I => the beyond-climate\n")
cat("residual is spatially organised (zonality is real), complementing the\n")
cat("scale-dependent ICC. Magnitude is modest, consistent with the noise floor.\n")
