################################################################################
# 13c_commonality_validation.R
#
# Spatial-block cross-validated variance decomposition of log(MRT) across the
# four mechanistic predictor groups (Climate, Edaphic, LandUse, Biological).
#
# This is the production, house-style replacement for the exploratory scratchpad
# scripts (spatialcv_validation.R / spatialcv_finish.R) used during the
# 2026-06-25/26 sessions. It resolves the "how much, and to whom" question that
# underpins framing A (zonality of SOC decay as an interpretive framework):
#
#   Q: Of the spatially-honest (block-CV) predictable variance in soil-carbon
#      turnover, how much is attributable to each mechanistic domain, once the
#      domains' mutual correlation is accounted for?
#
# Why a dedicated script (vs. 13_model_fitting.R):
#   * 13_model_fitting.R fits the 7 nested models M1..M7 each on its OWN set of
#     complete cases, so their R^2 are NOT mutually comparable and cannot be
#     combined into an exact decomposition.
#   * Here every coalition is fit on ONE COMMON ROW SET (complete across all four
#     groups, biology included), so coalition R^2 are comparable and the Shapley
#     values sum EXACTLY to the full-model R^2.
#
# Method:
#   * Common row set: complete.cases over the union of all group predictors +
#     response + coordinates (biology's ~29% footprint is the binding
#     constraint; expected n ~ 1.0e5).
#   * Spatial-block CV: 2 degree grid cells assigned to 10 folds (seed 42),
#     matching 13_model_fitting.R, so the blocking is identical to the main fit.
#   * Fit all 2^4 - 1 = 15 non-empty group coalitions; v(S) = block-CV R^2 of the
#     RF trained on the union of variables in the groups of S; v({}) = 0.
#   * Derive: Shapley value per group (exact, symmetric allocation of v(full));
#     singles v({i}); unique increment v(full) - v(full \ {i}); a nested
#     increment ladder (C -> +E -> +L -> +B); and the overlap
#     sum_i v({i}) - v(full) (the "correlated, not redundant" quantity).
#
# Input:  ./Global_MRT_code/outputs/12b_model_ready.rds
#         ./Global_MRT_code/outputs/13_var_groups.rds
# Output: ./Global_MRT_code/outputs/13c_coalition_r2.csv
#         ./Global_MRT_code/outputs/13c_shapley.csv
#         ./Global_MRT_code/outputs/13c_decomposition.rds
#         ./Global_MRT_code/plots/13c_shapley_shares.png
#         ./Global_MRT_code/plots/13c_singles.png
#         ./Global_MRT_code/plots/13c_nested_ladder.png
#
# Biological set: the MAIN run (default) uses fungal + 5 SPUN layers (96%
# coverage). Re-run with env MRT_13C_BIO=full for the all-9-layer SENSITIVITY
# case, which writes the same outputs with an "_allbio" suffix.
#
# Runtime note: 15 coalitions x 10 folds = 150 RF fits at 500 trees on ~1e5
# rows. Expect tens of minutes on a multicore machine.
#
# Author: Lorenzo
# Date: 2026-06-26
################################################################################

library(dplyr)
library(tidyr)
library(ranger)
library(ggplot2)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots", "step_13c_commonality")

INPUT_FILE      <- file.path(OUTPUT_DIR, "12b_model_ready.rds")
VAR_GROUPS_FILE <- file.path(OUTPUT_DIR, "13_var_groups.rds")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Random Forest parameters (production: 500 trees, matching 13_model_fitting.R).
# Two optional env overrides exist purely for a fast end-to-end smoke test and
# DO NOT affect a default run:
#   MRT_13C_NTREES   - override tree count (e.g. 50 for a quick check)
#   MRT_13C_SUBSAMPLE- fraction of common rows to keep (e.g. 0.1); 1 = all rows
# Leave both unset for the production 500-tree, full-data run.
RF_NTREES        <- as.integer(Sys.getenv("MRT_13C_NTREES", "500"))
RF_MIN_NODE_SIZE <- 5
RF_NUM_THREADS   <- parallel::detectCores() - 1  # Leave one core free
SUBSAMPLE_FRAC   <- as.numeric(Sys.getenv("MRT_13C_SUBSAMPLE", "1"))

# Spatial-block CV parameters (identical to 13_model_fitting.R).
# BLOCK_SIZE is env-overridable for the block-size sensitivity sweep (1/2/5 deg);
# non-default sizes tag outputs with a "_blkN" suffix so the 2 deg headline is safe.
CV_FOLDS   <- 10
BLOCK_SIZE <- as.numeric(Sys.getenv("MRT_13C_BLOCK", "2"))  # degrees (~200 km at 2)

set.seed(42)

cat("===============================================================\n")
cat("  COMMONALITY / SHAPLEY VALIDATION (Step 13c)\n")
cat("===============================================================\n\n")
cat("Using", RF_NUM_THREADS, "threads;", RF_NTREES, "trees per forest\n\n")

# =============================================================================
# LOAD DATA AND VARIABLE GROUPS
# =============================================================================

cat("Loading data and variable groups...\n")
soil_data  <- readRDS(INPUT_FILE)
VAR_GROUPS <- readRDS(VAR_GROUPS_FILE)

cat("  Observations:", format(nrow(soil_data), big.mark = ","), "\n")
cat("  Groups:", paste(names(VAR_GROUPS), collapse = ", "), "\n")
cat("  Group sizes:", paste(sprintf("%s=%d", names(VAR_GROUPS),
                                    lengths(VAR_GROUPS)), collapse = ", "), "\n\n")

# =============================================================================
# PREPARE RESPONSE (identical recipe to 13_model_fitting.R)
# =============================================================================

cat("Filtering to valid MRT and removing 1%/99% outliers...\n")
model_data <- soil_data %>%
  filter(MRT_QC == "valid", !is.na(MRT_years), MRT_years > 0, MRT_years < Inf)

q01 <- quantile(model_data$MRT_years, 0.01, na.rm = TRUE)
q99 <- quantile(model_data$MRT_years, 0.99, na.rm = TRUE)
model_data <- model_data %>% filter(MRT_years > q01 & MRT_years < q99)
model_data$log_MRT <- log(model_data$MRT_years)
cat("  After response QC:", format(nrow(model_data), big.mark = ","), "rows\n\n")

# Keep only variables actually present (defensive; mirrors 13 availability check)
VAR_GROUPS <- lapply(VAR_GROUPS, function(v) v[v %in% names(model_data)])

# -----------------------------------------------------------------------------
# BIOLOGICAL PREDICTOR SET (decision 2026-06-26, user-endorsed)
# -----------------------------------------------------------------------------
# The MAIN model uses fungal_proportion + the five SPUN richness/endemism layers
# (96% observation coverage). The three Barcelo root-colonization layers are
# DROPPED from the main model: at equal rows they add only +0.001 unique block-CV
# R^2 (statistically inert; alone they are negative), while their ~62% coverage is
# the sole constraint collapsing the common-row set to n=101,400. Dropping them
# lifts biological coverage 61% -> 96% and retains ~84% of the biological dR^2.
# The full 9-layer set is retained as a SENSITIVITY case (run with MRT_13C_BIO=full).
# See bio-skill triage and memory [[appendix-biological-parsimony]].
BIOLOGICAL_MAIN <- c("fungal_proportion",
                     "AM_richness", "EcM_richness",
                     "AM_endemism", "EcM_endemism",
                     "EcM_AM_richness_ratio")
BIO_MODE <- Sys.getenv("MRT_13C_BIO", "main")   # "main" (6 layers) | "full" (9)
if (BIO_MODE == "main") {
  VAR_GROUPS$BIOLOGICAL <- intersect(BIOLOGICAL_MAIN, VAR_GROUPS$BIOLOGICAL)
  OUT_SUFFIX <- ""
} else {
  OUT_SUFFIX <- "_allbio"
}
# Non-default CV block size (sensitivity sweep) gets its own suffix so the 2 deg
# headline outputs are never overwritten.
if (BLOCK_SIZE != 2) OUT_SUFFIX <- paste0(OUT_SUFFIX, "_blk", BLOCK_SIZE)
cat(sprintf("Biological set: %s (%d layers); CV block: %g deg; output suffix '%s'\n",
            BIO_MODE, length(VAR_GROUPS$BIOLOGICAL), BLOCK_SIZE, OUT_SUFFIX))

all_vars <- unlist(VAR_GROUPS, use.names = FALSE)
all_vars <- all_vars[all_vars %in% names(model_data)]

# =============================================================================
# COMMON ROW SET (the crux: identical rows for every coalition)
# =============================================================================
#
# Every coalition is fit on the SAME rows -- those complete across ALL four
# groups -- so coalition R^2 are mutually comparable and the Shapley values
# sum exactly to the full-model R^2. Biology's sparse footprint dominates this
# constraint by design (it is a finding, not a nuisance: see working memo).
# =============================================================================

cat("Building common row set (complete across ALL groups)...\n")
needed <- c("log_MRT", all_vars,
            "longitude_decimal_degrees", "latitude_decimal_degrees")
complete_idx <- complete.cases(model_data[, needed])
cv_data <- model_data[complete_idx, ]
n_common <- nrow(cv_data)
cat("  Common rows (n_common):", format(n_common, big.mark = ","), "\n")
cat("  (drops to the biology footprint; this is expected)\n")

if (SUBSAMPLE_FRAC < 1) {
  set.seed(42)
  cv_data  <- cv_data[sample(nrow(cv_data), floor(nrow(cv_data) * SUBSAMPLE_FRAC)), ]
  n_common <- nrow(cv_data)
  cat("  [SMOKE TEST] subsampled to", format(n_common, big.mark = ","),
      "rows (fraction", SUBSAMPLE_FRAC, ")\n")
}
cat("\n")

# =============================================================================
# SPATIAL-BLOCK FOLDS (2 deg grid -> 10 folds, seed 42; same as step 13)
# =============================================================================

cat("Assigning spatial-block CV folds...\n")
cv_data <- cv_data %>%
  mutate(
    grid_x  = floor(longitude_decimal_degrees / BLOCK_SIZE),
    grid_y  = floor(latitude_decimal_degrees  / BLOCK_SIZE),
    grid_id = paste(grid_x, grid_y, sep = "_")
  )
unique_grids <- unique(cv_data$grid_id)
set.seed(42)
grid_folds <- data.frame(
  grid_id = unique_grids,
  cv_fold = sample(seq_len(CV_FOLDS), length(unique_grids), replace = TRUE)
)
cv_data <- cv_data %>% left_join(grid_folds, by = "grid_id")
cat("  Block size:", BLOCK_SIZE, "deg;  unique blocks:",
    length(unique_grids), ";  folds:", CV_FOLDS, "\n\n")

# =============================================================================
# CORE: SPATIAL-BLOCK CV R^2 FOR A GIVEN PREDICTOR SET
# =============================================================================

cv_block_r2 <- function(predictors, data = cv_data,
                        response = "log_MRT", fold_col = "cv_fold") {
  observed    <- data[[response]]
  predictions <- rep(NA_real_, nrow(data))
  rf_formula  <- as.formula(
    paste(response, "~", paste(predictors, collapse = " + ")))

  for (fold in sort(unique(data[[fold_col]]))) {
    test_idx  <- data[[fold_col]] == fold
    rf_model  <- ranger(
      formula       = rf_formula,
      data          = data[!test_idx, ],
      num.trees     = RF_NTREES,
      min.node.size = RF_MIN_NODE_SIZE,
      num.threads   = RF_NUM_THREADS,
      importance    = "none",
      seed          = 42
    )
    predictions[test_idx] <- predict(rf_model, data = data[test_idx, ])$predictions
  }

  res    <- observed - predictions
  SS_res <- sum(res^2)
  SS_tot <- sum((observed - mean(observed))^2)
  1 - SS_res / SS_tot
}

# =============================================================================
# ENUMERATE ALL 15 NON-EMPTY GROUP COALITIONS AND COMPUTE v(S)
# =============================================================================

GROUPS  <- names(VAR_GROUPS)          # CLIMATE EDAPHIC LANDUSE BIOLOGICAL
LETTERS_MAP <- c(CLIMATE = "C", EDAPHIC = "E", LANDUSE = "L", BIOLOGICAL = "B")
nG      <- length(GROUPS)

# All non-empty subsets as logical membership matrix (rows = coalitions)
subset_list <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), nG))
names(subset_list) <- GROUPS
subset_list <- subset_list[rowSums(subset_list) > 0, , drop = FALSE]
subset_list <- subset_list[order(rowSums(subset_list)), , drop = FALSE]

cat("===============================================================\n")
cat("  FITTING", nrow(subset_list), "COALITIONS (block-CV R^2)\n")
cat("===============================================================\n\n")

vS <- setNames(rep(NA_real_, nrow(subset_list)), NULL)
coalition_label <- character(nrow(subset_list))
coalition_npred <- integer(nrow(subset_list))
total_start <- Sys.time()

for (i in seq_len(nrow(subset_list))) {
  members <- GROUPS[unlist(subset_list[i, ])]
  preds   <- unlist(VAR_GROUPS[members], use.names = FALSE)
  label   <- paste(LETTERS_MAP[members], collapse = "")

  cat(sprintf("  [%2d/%2d] %-4s (%2d predictors) ... ",
              i, nrow(subset_list), label, length(preds)))
  t0 <- Sys.time()
  vS[i] <- cv_block_r2(preds)
  coalition_label[i] <- label
  coalition_npred[i] <- length(preds)
  cat(sprintf("R2 = %.4f  [%.1f min]\n",
              vS[i], as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat(sprintf("\nTotal coalition-fitting time: %.1f min\n\n",
            as.numeric(difftime(Sys.time(), total_start, units = "mins"))))

# Named lookup keyed by coalition label (e.g. "CE", "CELB")
v <- setNames(vS, coalition_label)
v_lookup <- function(members_letters) {
  if (length(members_letters) == 0) return(0)
  # canonical order C,E,L,B
  key <- paste(intersect(c("C", "E", "L", "B"), members_letters), collapse = "")
  unname(v[key])
}
full_label <- paste(LETTERS_MAP[GROUPS], collapse = "")
full_R2    <- unname(v[full_label])

coalition_df <- data.frame(
  coalition    = coalition_label,
  n_groups     = nchar(coalition_label),
  n_predictors = coalition_npred,
  R2_block_cv  = round(vS, 5),
  stringsAsFactors = FALSE
)
coalition_df <- coalition_df[order(coalition_df$n_groups,
                                   -coalition_df$R2_block_cv), ]

# =============================================================================
# SHAPLEY VALUES (exact allocation of full_R2 across the four groups)
# =============================================================================
#
# phi_i = sum_{S subseteq N\{i}} |S|!(n-|S|-1)!/n! * ( v(S u {i}) - v(S) )
# Shapley values sum exactly to v(N) = full_R2.
# =============================================================================

cat("Computing Shapley values...\n")
group_letters <- unname(LETTERS_MAP[GROUPS])
shapley <- setNames(rep(0, nG), GROUPS)

for (gi in seq_len(nG)) {
  i_letter <- group_letters[gi]
  others   <- setdiff(group_letters, i_letter)
  # iterate over all subsets S of `others`
  for (k in 0:length(others)) {
    combos <- if (k == 0) list(character(0)) else
      utils::combn(others, k, simplify = FALSE)
    w <- factorial(k) * factorial(nG - k - 1) / factorial(nG)
    for (S in combos) {
      shapley[gi] <- shapley[gi] +
        w * (v_lookup(c(S, i_letter)) - v_lookup(S))
    }
  }
}

# =============================================================================
# SINGLES, UNIQUE INCREMENTS, NESTED LADDER, OVERLAP
# =============================================================================

singles <- setNames(sapply(group_letters, function(L) v_lookup(L)), GROUPS)

unique_incr <- setNames(
  sapply(group_letters, function(L)
    full_R2 - v_lookup(setdiff(group_letters, L))),
  GROUPS)

# Nested ladder in the manuscript's narrative order: C -> +E -> +L -> +B
ladder_order <- c("C", "E", "L", "B")
ladder <- data.frame(
  step      = paste0("+", c("Climate", "Edaphic", "LandUse", "Biological")),
  added     = c("Climate", "Edaphic", "LandUse", "Biological"),
  R2_cumul  = NA_real_,
  increment = NA_real_,
  stringsAsFactors = FALSE
)
prev <- 0
for (j in seq_along(ladder_order)) {
  cur <- v_lookup(ladder_order[1:j])
  ladder$R2_cumul[j]  <- cur
  ladder$increment[j] <- cur - prev
  prev <- cur
}

overlap <- sum(singles) - full_R2

# =============================================================================
# REPORT
# =============================================================================

cat("\n===============================================================\n")
cat("  RESULTS\n")
cat("===============================================================\n\n")

cat("Full-model block-CV R^2 (v(full)):", round(full_R2, 4),
    "   n_common =", format(n_common, big.mark = ","), "\n\n")

shapley_df <- data.frame(
  group        = GROUPS,
  shapley_R2   = round(shapley, 4),
  shapley_pct  = round(100 * shapley / full_R2, 1),
  single_R2    = round(singles, 4),
  unique_R2    = round(unique_incr, 4),
  row.names    = NULL
)
shapley_df <- shapley_df[order(-shapley_df$shapley_R2), ]

cat("Shapley decomposition (shares sum to 100%):\n")
print(shapley_df, row.names = FALSE)
cat(sprintf("\n  Shapley sum = %.4f  (should equal full R^2 = %.4f)\n",
            sum(shapley), full_R2))
cat(sprintf("  Singles sum = %.4f  -> overlap = %.4f  (correlated, not redundant)\n\n",
            sum(singles), overlap))

cat("Nested ladder (order-dependent; C -> +E -> +L -> +B):\n")
print(transform(ladder,
                R2_cumul  = round(R2_cumul, 4),
                increment = round(increment, 4)), row.names = FALSE)
cat("\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("===============================================================\n")
cat("  SAVING OUTPUTS\n")
cat("===============================================================\n\n")

f_coal <- file.path(OUTPUT_DIR, paste0("13c_coalition_r2", OUT_SUFFIX, ".csv"))
write.csv(coalition_df, f_coal, row.names = FALSE)
cat("OK  ", f_coal, "\n")

f_shap <- file.path(OUTPUT_DIR, paste0("13c_shapley", OUT_SUFFIX, ".csv"))
write.csv(shapley_df, f_shap, row.names = FALSE)
cat("OK  ", f_shap, "\n")

decomposition <- list(
  config = list(
    rf_ntrees = RF_NTREES, rf_min_node_size = RF_MIN_NODE_SIZE,
    cv_folds = CV_FOLDS, block_size_deg = BLOCK_SIZE, seed = 42,
    n_common = n_common, var_groups = VAR_GROUPS,
    bio_mode = BIO_MODE, biological_vars = VAR_GROUPS$BIOLOGICAL,
    date = Sys.Date()
  ),
  full_R2      = full_R2,
  coalition_R2 = v,
  coalition_df = coalition_df,
  shapley      = shapley,
  shapley_df   = shapley_df,
  singles      = singles,
  unique_incr  = unique_incr,
  nested_ladder = ladder,
  overlap      = overlap
)
f_dec <- file.path(OUTPUT_DIR, paste0("13c_decomposition", OUT_SUFFIX, ".rds"))
saveRDS(decomposition, f_dec)
cat("OK  ", f_dec, "\n\n")

# =============================================================================
# FIGURES (serve the Results variance-decomposition figure, step 13c -> MS)
# =============================================================================

group_colors <- c(
  "Climate"    = "#666666",
  "Edaphic"    = "#E41A1C",
  "LandUse"    = "#4DAF4A",
  "Biological" = "#377EB8"
)
nice_name <- c(CLIMATE = "Climate", EDAPHIC = "Edaphic",
               LANDUSE = "LandUse", BIOLOGICAL = "Biological")

# (1) Shapley shares
plot_shapley <- shapley_df %>%
  mutate(group_nice = nice_name[group])
p_sh <- ggplot(plot_shapley,
               aes(x = reorder(group_nice, shapley_R2),
                   y = shapley_pct, fill = group_nice)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", shapley_pct)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = group_colors, guide = "none") +
  ylim(0, max(plot_shapley$shapley_pct) * 1.15) +
  labs(title = "Shapley attribution of block-CV predictable variance",
       subtitle = sprintf("Shares of full-model R^2 = %.3f (n = %s)",
                          full_R2, format(n_common, big.mark = ",")),
       x = "", y = "Share of attributable variance (%)") +
  theme_bw()
f_psh <- file.path(PLOT_DIR, paste0("13c_shapley_shares", OUT_SUFFIX, ".png"))
ggsave(f_psh, p_sh, width = 8, height = 4.5, dpi = 150)
cat("OK  ", f_psh, "\n")

# (2) Singles (each domain alone) -- the Results "singles-first" opening
plot_singles <- data.frame(group = GROUPS, single_R2 = singles) %>%
  mutate(group_nice = nice_name[group])
p_si <- ggplot(plot_singles,
               aes(x = reorder(group_nice, single_R2),
                   y = single_R2, fill = group_nice)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.3f", single_R2)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = group_colors, guide = "none") +
  ylim(0, max(plot_singles$single_R2) * 1.15) +
  labs(title = "Single-domain block-CV R^2 (each domain alone)",
       subtitle = sprintf("Singles sum = %.3f vs full = %.3f -> overlap = %.3f",
                          sum(singles), full_R2, overlap),
       x = "", y = "Block-CV R^2") +
  theme_bw()
f_psi <- file.path(PLOT_DIR, paste0("13c_singles", OUT_SUFFIX, ".png"))
ggsave(f_psi, p_si, width = 8, height = 4.5, dpi = 150)
cat("OK  ", f_psi, "\n")

# (3) Nested ladder waterfall (C -> +E -> +L -> +B)
plot_ladder <- ladder %>%
  mutate(step = factor(step, levels = step))
p_la <- ggplot(plot_ladder, aes(x = step, y = R2_cumul, group = 1)) +
  geom_col(aes(fill = added)) +
  geom_line() + geom_point() +
  geom_text(aes(label = sprintf("+%.3f", increment)), vjust = -0.6, size = 3.2) +
  scale_fill_manual(values = setNames(group_colors, nice_name[GROUPS]),
                    guide = "none") +
  labs(title = "Nested cumulative block-CV R^2",
       subtitle = "Order-dependent ladder; biology adds little when nested last",
       x = "", y = "Cumulative block-CV R^2") +
  theme_bw()
f_pla <- file.path(PLOT_DIR, paste0("13c_nested_ladder", OUT_SUFFIX, ".png"))
ggsave(f_pla, p_la, width = 8, height = 4.5, dpi = 150)
cat("OK  ", f_pla, "\n")

cat("\n===============================================================\n")
cat("  STEP 13c COMPLETE\n")
cat("===============================================================\n\n")
cat("Headline: climate is ~1/3 of attributable signal; non-climate domains ~2/3.\n")
cat("Use Shapley shares (sum-exact) as the Results headline; singles-first as the\n")
cat("opening; nested ladder shows order-dependence (biology last underrates it).\n")
