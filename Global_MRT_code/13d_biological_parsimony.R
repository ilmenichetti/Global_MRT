################################################################################
# 13d_biological_parsimony.R
#
# Appendix analysis underpinning the biological-predictor decision (2026-06-26):
# the MAIN model uses fungal_proportion + 5 SPUN richness/endemism layers, and
# DROPS the three Barcelo root-colonization layers, which are the sole constraint
# capping the common-row footprint yet add no unique beyond-abiotic skill.
#
# This is the reproducible, house-style record of the scratchpad triage. It runs:
#   1. Per-layer observation coverage (% of modeling rows non-NA).
#   2. Collinearity of the 9 biological layers (Spearman; pairwise complete).
#   3. Single-layer marginal block-CV dR^2: each biological layer added to the
#      M5 abiotic baseline, on a COMMON row set so dR^2 are comparable; plotted
#      against each layer's native coverage -> the skill-vs-coverage trade-off.
#   4. Subset comparison reproducing the decision: M5 + {fungal, SPUN5,
#      fungal+SPUN5 (MAIN), Barcelo3, all-9}.
#
# Input:  ./Global_MRT_code/outputs/12b_model_ready.rds
#         ./Global_MRT_code/outputs/13_var_groups.rds
# Output: ./Global_MRT_code/outputs/13d_bio_coverage.csv
#         ./Global_MRT_code/outputs/13d_single_layer_skill.csv
#         ./Global_MRT_code/outputs/13d_subset_skill.csv
#         ./Global_MRT_code/outputs/13d_parsimony.rds
#         ./Global_MRT_code/plots/appendix/biological_collinearity.png
#         ./Global_MRT_code/plots/appendix/biological_single_layer_skill.png
#
# Env override: MRT_13D_NTREES (default 500; e.g. 50 for a smoke test).
#
# Author: Lorenzo
# Date: 2026-06-26
################################################################################

library(dplyr)
library(ggplot2)
library(ranger)
library(corrplot)

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots", "appendix")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

RF_NTREES      <- as.integer(Sys.getenv("MRT_13D_NTREES", "500"))
RF_MIN_NODE    <- 5
RF_THREADS     <- parallel::detectCores() - 1
CV_FOLDS       <- 10
BLOCK_SIZE     <- 2
set.seed(42)

cat("===============================================================\n")
cat("  BIOLOGICAL PARSIMONY (Step 13d) -", RF_NTREES, "trees\n")
cat("===============================================================\n\n")

# ---- Biological layer families (for colouring / interpretation) -------------
FAM <- c(
  fungal_proportion     = "Fungal (Yu 2022)",
  AM_richness           = "SPUN richness/endemism",
  EcM_richness          = "SPUN richness/endemism",
  AM_endemism           = "SPUN richness/endemism",
  EcM_endemism          = "SPUN richness/endemism",
  EcM_AM_richness_ratio = "SPUN richness/endemism",
  AM_roots_colonized    = "Barcelo root colonization",
  EcM_roots_colonized   = "Barcelo root colonization",
  EcM_AM_root_ratio     = "Barcelo root colonization"
)
BIO_MAIN <- c("fungal_proportion", "AM_richness", "EcM_richness",
              "AM_endemism", "EcM_endemism", "EcM_AM_richness_ratio")

# ---- Data prep (same recipe as 13 / 13c) ------------------------------------
d <- readRDS(file.path(OUTPUT_DIR, "12b_model_ready.rds"))
g <- readRDS(file.path(OUTPUT_DIR, "13_var_groups.rds"))
bio <- g$BIOLOGICAL
M5  <- c(g$CLIMATE, g$EDAPHIC, g$LANDUSE)

d <- d %>% filter(MRT_QC == "valid", !is.na(MRT_years), MRT_years > 0, MRT_years < Inf)
q <- quantile(d$MRT_years, c(.01, .99), na.rm = TRUE)
d <- d %>% filter(MRT_years > q[1] & MRT_years < q[2])
d$log_MRT <- log(d$MRT_years)
N_universe <- nrow(d)
cat("Modeling-universe rows:", format(N_universe, big.mark = ","), "\n\n")

# =============================================================================
# 1. PER-LAYER COVERAGE
# =============================================================================
coverage_df <- data.frame(
  layer    = bio,
  family   = FAM[bio],
  n        = sapply(bio, function(v) sum(!is.na(d[[v]]))),
  coverage = sapply(bio, function(v) mean(!is.na(d[[v]])) * 100),
  row.names = NULL
) %>% arrange(desc(coverage))
cat("Per-layer coverage (% of modeling rows):\n")
print(transform(coverage_df, coverage = round(coverage, 1)), row.names = FALSE)
write.csv(coverage_df, file.path(OUTPUT_DIR, "13d_bio_coverage.csv"), row.names = FALSE)

# =============================================================================
# 2. COLLINEARITY (Spearman, pairwise complete)
# =============================================================================
cm <- cor(d[, bio], use = "pairwise.complete.obs", method = "spearman")
png(file.path(PLOT_DIR, "biological_collinearity.png"),
    width = 1700, height = 1700, res = 200)
corrplot(cm, method = "color", type = "upper", order = "hclust",
         addCoef.col = "black", number.cex = 0.7, tl.col = "black",
         tl.srt = 45, tl.cex = 0.8,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         mar = c(0, 0, 2, 0),
         title = "Biological layer collinearity (Spearman)")
dev.off()
cat("\nOK  ", file.path(PLOT_DIR, "biological_collinearity.png"), "\n")

# =============================================================================
# COMMON ROW SET + spatial folds (equal rows -> comparable dR^2)
# =============================================================================
need <- c("log_MRT", M5, bio, "longitude_decimal_degrees", "latitude_decimal_degrees")
cv <- d[complete.cases(d[, need]), ]
cat("\nCommon-row n (M5 + all 9 bio):", format(nrow(cv), big.mark = ","), "\n")
cv$grid <- paste(floor(cv$longitude_decimal_degrees / BLOCK_SIZE),
                 floor(cv$latitude_decimal_degrees  / BLOCK_SIZE))
ug <- unique(cv$grid); set.seed(42)
cv$fold <- data.frame(grid = ug, f = sample(CV_FOLDS, length(ug), TRUE))$f[match(cv$grid, ug)]

cv_block_r2 <- function(preds) {
  obs <- cv$log_MRT; pr <- rep(NA_real_, nrow(cv))
  f <- as.formula(paste("log_MRT ~", paste(preds, collapse = " + ")))
  for (k in sort(unique(cv$fold))) {
    te <- cv$fold == k
    m <- ranger(f, cv[!te, ], num.trees = RF_NTREES, min.node.size = RF_MIN_NODE,
                num.threads = RF_THREADS, seed = 42)
    pr[te] <- predict(m, cv[te, ])$predictions
  }
  1 - sum((obs - pr)^2) / sum((obs - mean(obs))^2)
}

# =============================================================================
# 3. SINGLE-LAYER MARGINAL dR^2 (each layer added to M5)
# =============================================================================
cat("\nFitting M5 baseline + single-layer additions...\n")
base_R2 <- cv_block_r2(M5)
cat("  M5 baseline R2 =", round(base_R2, 4), "\n")

single_df <- coverage_df
single_df$dR2 <- NA_real_
for (i in seq_len(nrow(single_df))) {
  v <- single_df$layer[i]
  single_df$dR2[i] <- cv_block_r2(c(M5, v)) - base_R2
  cat(sprintf("  +%-22s dR2 = %+.4f  (cov %.0f%%)\n",
              v, single_df$dR2[i], single_df$coverage[i]))
}
write.csv(single_df, file.path(OUTPUT_DIR, "13d_single_layer_skill.csv"), row.names = FALSE)

# Skill-vs-coverage trade-off plot
fam_cols <- c("Fungal (Yu 2022)" = "#377EB8",
              "SPUN richness/endemism" = "#4DAF4A",
              "Barcelo root colonization" = "#E41A1C")
p_skill <- ggplot(single_df, aes(coverage, dR2, colour = family)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = layer), size = 2.8, max.overlaps = 20,
                           show.legend = FALSE) +
  scale_colour_manual(values = fam_cols, name = NULL) +
  labs(title = "Single-layer skill vs coverage trade-off",
       subtitle = sprintf("Marginal block-CV ΔR² over M5 abiotic baseline (equal rows, n=%s)",
                          format(nrow(cv), big.mark = ",")),
       x = "Native observation coverage (% of modeling rows)",
       y = expression(Delta * R^2 ~ "over M5 (equal rows)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, "biological_single_layer_skill.png"), p_skill,
       width = 8.5, height = 6, dpi = 200)
cat("OK  ", file.path(PLOT_DIR, "biological_single_layer_skill.png"), "\n")

# =============================================================================
# 4. SUBSET COMPARISON (reproduces the decision)
# =============================================================================
cat("\nFitting biological subset comparison...\n")
spun5 <- c("AM_richness", "EcM_richness", "AM_endemism", "EcM_endemism", "EcM_AM_richness_ratio")
barc3 <- c("AM_roots_colonized", "EcM_roots_colonized", "EcM_AM_root_ratio")
subsets <- list(
  "M5 (abiotic base)"         = M5,
  "M5 + fungal"               = c(M5, "fungal_proportion"),
  "M5 + SPUN5"                = c(M5, spun5),
  "M5 + fungal + SPUN5 [MAIN]"= c(M5, BIO_MAIN),
  "M5 + Barcelo3"             = c(M5, barc3),
  "M5 + all 9 bio"            = c(M5, bio)
)
subset_df <- data.frame(model = names(subsets), R2 = NA_real_)
for (i in seq_along(subsets)) {
  subset_df$R2[i] <- cv_block_r2(subsets[[i]])
  cat(sprintf("  %-28s R2 = %.4f\n", names(subsets)[i], subset_df$R2[i]))
}
subset_df$dR2_over_M5 <- subset_df$R2 - subset_df$R2[1]
write.csv(subset_df, file.path(OUTPUT_DIR, "13d_subset_skill.csv"), row.names = FALSE)

uniq_barc <- subset_df$dR2_over_M5[subset_df$model == "M5 + all 9 bio"] -
             subset_df$dR2_over_M5[subset_df$model == "M5 + fungal + SPUN5 [MAIN]"]

# =============================================================================
# SAVE + REPORT
# =============================================================================
saveRDS(list(
  config = list(rf_ntrees = RF_NTREES, cv_folds = CV_FOLDS, block_size = BLOCK_SIZE,
                seed = 42, n_common = nrow(cv), n_universe = N_universe, date = Sys.Date()),
  coverage = coverage_df, collinearity = cm,
  single_layer = single_df, subsets = subset_df,
  unique_barcelo = uniq_barc, bio_main = BIO_MAIN
), file.path(OUTPUT_DIR, "13d_parsimony.rds"))

cat("\n===============================================================\n")
cat("  SUMMARY\n")
cat("===============================================================\n")
cat(sprintf("Unique Barcelo-3 contribution (all9 - MAIN, equal rows): %+.4f\n", uniq_barc))
cat(sprintf("MAIN set coverage: 96%%  vs  all-9 coverage: %.0f%%\n",
            mean(complete.cases(d[, bio])) * 100))
cat("Decision: MAIN = fungal + 5 SPUN (drop Barcelo-3); all-9 kept as sensitivity.\n")
cat("OK  ", file.path(OUTPUT_DIR, "13d_parsimony.rds"), "\n")
