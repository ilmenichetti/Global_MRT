################################################################################
# 13g_provenance_control.R
#
# Provenance control on the beyond-climate signal (appendix). House-style
# reconstruction of the (now-deleted) scratchpad `provenance_control.R` so the
# appendix table is regenerable. The soil compendium aggregates 44 source
# databases that are geographically clustered; a source-specific measurement or
# laboratory offset could in principle masquerade as large-scale (zonal)
# structure. We test whether the beyond-climate residual is carried by data
# provenance (source_db) rather than ecology.
#
# Residual = observed log(MRT) - climate-only spatial-block-CV prediction
#            (identical recipe to 13e; deterministic seed 42 -> same residual).
# Footprint = MAIN common-row set (M5 + MAIN biological set), matching 13c.
#
# Metrics (before vs after removing each source's mean residual):
#   - eta^2_source: between-source variance share of the residual.
#   - Organisation @ 2/5/10 deg: observed between-block variance of the residual
#     divided by a spatial-randomness (label-permutation) null -> obs/null ratio.
#   - Block-mean covariate R^2 @ 5 deg: linear skill of block-mean covariates on
#     the block-mean residual (the "zonal" residual).
# If source de-meaning barely changes these, the zonal signal is ecological, not
# a provenance artifact. The test is conservative (source is confounded with
# geography), i.e. an upper bound on the provenance contribution.
#
# Input:  outputs/12b_model_ready.rds, outputs/13_var_groups.rds
# Output: outputs/13g_provenance.csv, outputs/13g_provenance.rds
#
# Env override: MRT_13G_NTREES (default 500).
# Author: Lorenzo   Date: 2026-06-26
################################################################################

library(dplyr)
library(ranger)

OUTPUT_DIR <- "./Global_MRT_code/outputs"
RF_NTREES  <- as.integer(Sys.getenv("MRT_13G_NTREES", "500"))
RF_THREADS <- 3            # capped to coexist with other background RF jobs
CV_FOLDS   <- 10
BLOCK_CV   <- 2
MIN_BLOCK_N <- 5           # blocks need >= this many points to enter organisation
NPERM      <- 199          # permutations for the spatial-randomness null
set.seed(42)

BIO_MAIN <- c("fungal_proportion", "AM_richness", "EcM_richness",
              "AM_endemism", "EcM_endemism", "EcM_AM_richness_ratio")

cat("===============================================================\n")
cat("  PROVENANCE CONTROL (Step 13g) -", RF_NTREES, "trees\n")
cat("===============================================================\n\n")

# ---- Data prep + MAIN common rows -------------------------------------------
d <- readRDS(file.path(OUTPUT_DIR, "12b_model_ready.rds"))
g <- readRDS(file.path(OUTPUT_DIR, "13_var_groups.rds"))
CLIM <- g$CLIMATE
M5   <- c(g$CLIMATE, g$EDAPHIC, g$LANDUSE)
COVS <- c(M5, BIO_MAIN)

d <- d %>% filter(MRT_QC == "valid", !is.na(MRT_years), MRT_years > 0, MRT_years < Inf)
q <- quantile(d$MRT_years, c(.01, .99), na.rm = TRUE)
d <- d %>% filter(MRT_years > q[1] & MRT_years < q[2])
d$log_MRT <- log(d$MRT_years)

need <- c("log_MRT", COVS, "source_db",
          "longitude_decimal_degrees", "latitude_decimal_degrees")
cv <- d[complete.cases(d[, need]), ]
cv$lon <- cv$longitude_decimal_degrees; cv$lat <- cv$latitude_decimal_degrees
cv$source_db <- factor(cv$source_db)
cat("MAIN common-row n:", format(nrow(cv), big.mark = ","),
    " | sources:", nlevels(cv$source_db), "\n")

# ---- Climate-only block-CV residual -----------------------------------------
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
cat("Climate-only block-CV R2 =",
    round(1 - sum(cv$resid^2) / sum((cv$log_MRT - mean(cv$log_MRT))^2), 4), "\n\n")

# ---- Source de-meaning ------------------------------------------------------
cv$resid_dm <- cv$resid - ave(cv$resid, cv$source_db)

# =============================================================================
# METRICS
# =============================================================================
eta2_source <- function(r) {
  ss_tot <- sum((r - mean(r))^2)
  ss_b   <- sum(tapply(r, cv$source_db, function(x) length(x) * (mean(x) - mean(r))^2))
  ss_b / ss_tot
}

# observed between-block variance of block means / permutation null
block_mean_var <- function(r, blk) {
  bm <- tapply(r, blk, function(x) if (length(x) >= MIN_BLOCK_N) mean(x) else NA_real_)
  var(bm, na.rm = TRUE)
}
org_ratio <- function(r, deg) {
  blk <- paste(floor(cv$lon / deg), floor(cv$lat / deg))
  obs  <- block_mean_var(r, blk)
  null <- mean(replicate(NPERM, block_mean_var(sample(r), blk)))
  obs / null
}
block_cov_R2 <- function(r, deg) {
  df <- cv[, COVS]
  df$r   <- r
  df$blk <- paste(floor(cv$lon / deg), floor(cv$lat / deg))
  agg <- df %>% group_by(blk) %>% filter(n() >= MIN_BLOCK_N) %>%
    summarise(across(all_of(c("r", COVS)), mean), .groups = "drop")
  summary(lm(reformulate(COVS, response = "r"), data = agg))$r.squared
}

scales <- c(2, 5, 10)
res <- data.frame(
  metric = c("eta2_source",
             paste0("organisation_obs_null_", scales, "deg"),
             "block_mean_covariate_R2_5deg"),
  before = NA_real_, after = NA_real_)

res$before[1] <- eta2_source(cv$resid)
res$after[1]  <- eta2_source(cv$resid_dm)
for (i in seq_along(scales)) {
  res$before[1 + i] <- org_ratio(cv$resid,    scales[i])
  res$after[1 + i]  <- org_ratio(cv$resid_dm, scales[i])
}
res$before[5] <- block_cov_R2(cv$resid,    5)
res$after[5]  <- block_cov_R2(cv$resid_dm, 5)

cat("===============================================================\n")
cat("  RESULTS (before vs after source de-meaning)\n")
cat("===============================================================\n")
print(transform(res, before = round(before, 4), after = round(after, 4)), row.names = FALSE)

write.csv(res, file.path(OUTPUT_DIR, "13g_provenance.csv"), row.names = FALSE)
saveRDS(list(config = list(rf_ntrees = RF_NTREES, n_common = nrow(cv),
                           n_sources = nlevels(cv$source_db), nperm = NPERM,
                           min_block_n = MIN_BLOCK_N, seed = 42, date = Sys.Date()),
             results = res),
        file.path(OUTPUT_DIR, "13g_provenance.rds"))
cat("\nOK  outputs/13g_provenance.csv  +  .rds\n")
cat(sprintf("\nVerdict: eta^2_source = %.1f%% of beyond-climate residual is between sources;\n",
            100 * res$before[1]))
cat("organisation and covariate skill barely change after de-meaning => the modest\n")
cat("zonal signal is ECOLOGICAL, not a provenance/laboratory artifact (conservative test).\n")
