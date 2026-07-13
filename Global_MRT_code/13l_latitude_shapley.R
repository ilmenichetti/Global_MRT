################################################################################
# 13l_latitude_shapley.R  --  EXPLORATORY (not yet wired into the manuscript)
#
# "Explanatory power by latitude": the band-wise counterpart to the global
# Shapley of step 13c. Within each latitude band we REFIT all 15 group coalitions
# with local spatial-block cross-validation and compute a band-specific Shapley
# attribution of the log(tau) variance across the four process classes
# (Climate / Edaphic / Land use / Biological). This tests the hypothesis that the
# ATTRIBUTION of tau (not just its modulation, Fig. 2) is zonally structured:
# climate's grip loosening at the extremes, edaphic rising in the tropics,
# land use / biology rising toward the boreal.
#
# Why local refit (not band-restricted evaluation of the global model): refitting
# inside the band lets each domain express the LOCAL tau relationship, which is
# exactly the "which class controls tau here" question. Because this is variance
# EXPLAINED (cross-validated R^2), it also strips the land-cover allocation-
# coefficient imprint that inflates land use in the per-pixel modulation maps.
#
# INTERPRETATION CAVEAT: shares are WITHIN-band (composition of that band's own
# full-model R^2). Conditioning on latitude removes the poleward gradient, so
# these measure sub-band spatial structure, not the cross-latitude level; a class
# that sets a band's mean but is uniform within it earns little local credit.
# Read alongside n and the band's full-model R^2 (both reported).
#
# Prep, CV and Shapley MUST MATCH 13c (production MAIN: below input, peat dropped
# via GPM, 6-layer biology, 2 deg blocks, 500 trees, seed 42).
#
# Output: outputs/13l_latitude_shapley.csv   (fine 15 deg bands + 4 biome bands)
#         outputs/13l_latitude_shapley.rds
#         plots/step_13c_commonality/13l_latitude_shapley.png
#
# Author: Lorenzo   Date: 2026-07-10
################################################################################

suppressMessages({ library(dplyr); library(ranger) })

# ---- CONFIG (identical to 13c) -------------------------------------------------
PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots", "step_13c_commonality")
INPUT_FILE      <- file.path(OUTPUT_DIR, "12b_model_ready.rds")
VAR_GROUPS_FILE <- file.path(OUTPUT_DIR, "13_var_groups.rds")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

RF_NTREES        <- as.integer(Sys.getenv("MRT_13C_NTREES", "500"))
RF_MIN_NODE_SIZE <- 5
RF_NUM_THREADS   <- as.integer(Sys.getenv("MRT_RF_THREADS",
                              as.character(parallel::detectCores() - 1)))
CV_FOLDS   <- 10
BLOCK_SIZE <- 2
set.seed(42)

# band-analysis reliability guards (env-overridable only to keep smoke tests exercising the path)
MIN_N      <- as.integer(Sys.getenv("MRT_13L_MIN_N", "800"))
MIN_BLOCKS <- as.integer(Sys.getenv("MRT_13L_MIN_BLOCKS", "40"))
FINE_STEP  <- as.numeric(Sys.getenv("MRT_13L_STEP", "10"))   # fine contiguous band width (deg)
# smoke-test hooks (leave unset for the production run): MRT_13C_NTREES=20 +
# MRT_13L_SUBSAMPLE=0.05 verifies the whole path in a couple of minutes.
SUBSAMPLE_FRAC <- as.numeric(Sys.getenv("MRT_13L_SUBSAMPLE", "1"))

# ---- LOAD + PREP (mirrors 13c: peat drop GPM, QC, 6-layer biology) -------------
cat("Loading 12b + variable groups...\n")
soil_data  <- readRDS(INPUT_FILE)
VAR_GROUPS <- readRDS(VAR_GROUPS_FILE)

source(file.path(PIPELINE_DIR, "peat_filter.R"))
soil_data <- apply_peat_filter(soil_data, OUTPUT_DIR)      # MAIN default = drop, GPM

model_data <- soil_data %>%
  filter(MRT_QC == "valid", !is.na(MRT_years), MRT_years > 0, MRT_years < Inf)
q01 <- quantile(model_data$MRT_years, 0.01, na.rm = TRUE)
q99 <- quantile(model_data$MRT_years, 0.99, na.rm = TRUE)
model_data <- model_data %>% filter(MRT_years > q01 & MRT_years < q99)
model_data$log_MRT <- log(model_data$MRT_years)

VAR_GROUPS <- lapply(VAR_GROUPS, function(v) v[v %in% names(model_data)])
BIOLOGICAL_MAIN <- c("fungal_proportion", "AM_richness", "EcM_richness",
                     "AM_endemism", "EcM_endemism", "EcM_AM_richness_ratio")
VAR_GROUPS$BIOLOGICAL <- intersect(BIOLOGICAL_MAIN, VAR_GROUPS$BIOLOGICAL)
all_vars <- unlist(VAR_GROUPS, use.names = FALSE)
all_vars <- all_vars[all_vars %in% names(model_data)]

# common row set (complete across ALL groups), with coords for banding + folds
needed <- c("log_MRT", all_vars, "longitude_decimal_degrees", "latitude_decimal_degrees")
cv_data <- model_data[complete.cases(model_data[, needed]), ]
cv_data <- cv_data %>%
  mutate(grid_x  = floor(longitude_decimal_degrees / BLOCK_SIZE),
         grid_y  = floor(latitude_decimal_degrees  / BLOCK_SIZE),
         grid_id = paste(grid_x, grid_y, sep = "_"),
         lat     = latitude_decimal_degrees)
if (SUBSAMPLE_FRAC < 1) {
  set.seed(42); cv_data <- cv_data[sample(nrow(cv_data), floor(nrow(cv_data) * SUBSAMPLE_FRAC)), ]
  cat(sprintf("  [SMOKE] subsampled to %s rows (frac %.3g)\n",
              format(nrow(cv_data), big.mark = ","), SUBSAMPLE_FRAC))
}
cat(sprintf("  Common rows: %s;  latitude span %.0f..%.0f\n\n",
            format(nrow(cv_data), big.mark = ","), min(cv_data$lat), max(cv_data$lat)))

GROUPS      <- names(VAR_GROUPS)                                    # CLIMATE EDAPHIC LANDUSE BIOLOGICAL
LETTERS_MAP <- c(CLIMATE = "C", EDAPHIC = "E", LANDUSE = "L", BIOLOGICAL = "B")
nG          <- length(GROUPS)
subset_list <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), nG))
names(subset_list) <- GROUPS
subset_list <- subset_list[rowSums(subset_list) > 0, , drop = FALSE]
subset_list <- subset_list[order(rowSums(subset_list)), , drop = FALSE]

# ---- CORE: block-CV R^2 for a predictor set on a given data subset -------------
cv_block_r2 <- function(predictors, data, response = "log_MRT", fold_col = "cv_fold") {
  observed    <- data[[response]]
  predictions <- rep(NA_real_, nrow(data))
  rf_formula  <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  for (fold in sort(unique(data[[fold_col]]))) {
    test_idx <- data[[fold_col]] == fold
    if (all(test_idx) || !any(test_idx)) next
    rf_model <- ranger(formula = rf_formula, data = data[!test_idx, ],
                       num.trees = RF_NTREES, min.node.size = RF_MIN_NODE_SIZE,
                       num.threads = RF_NUM_THREADS, importance = "none", seed = 42)
    predictions[test_idx] <- predict(rf_model, data = data[test_idx, ])$predictions
  }
  ok <- !is.na(predictions)
  1 - sum((observed[ok] - predictions[ok])^2) / sum((observed[ok] - mean(observed[ok]))^2)
}

# ---- exact Shapley allocation of full_R2 across the 4 groups (from 13c) ---------
shapley_from_v <- function(v) {                       # v: named coalition-R^2 lookup (keys "C","CE",...)
  v_lookup <- function(L) if (length(L) == 0) 0 else
    unname(v[paste(intersect(c("C","E","L","B"), L), collapse = "")])
  gl <- unname(LETTERS_MAP[GROUPS]); phi <- setNames(rep(0, nG), GROUPS)
  for (gi in seq_len(nG)) {
    il <- gl[gi]; others <- setdiff(gl, il)
    for (k in 0:length(others)) {
      combos <- if (k == 0) list(character(0)) else utils::combn(others, k, simplify = FALSE)
      w <- factorial(k) * factorial(nG - k - 1) / factorial(nG)
      for (S in combos) phi[gi] <- phi[gi] + w * (v_lookup(c(S, il)) - v_lookup(S))
    }
  }
  phi
}

# ---- decompose one band (local refit of all 15 coalitions) ---------------------
decompose_band <- function(d, band_name, lat_center) {
  nb <- length(unique(d$grid_id))
  if (nrow(d) < MIN_N || nb < MIN_BLOCKS) {
    cat(sprintf("  %-14s n=%-6d blocks=%-4d  -> INSUFFICIENT (skipped)\n",
                band_name, nrow(d), nb)); return(NULL)
  }
  nf <- min(CV_FOLDS, nb)
  ug <- unique(d$grid_id); set.seed(42)
  gf <- data.frame(grid_id = ug, cv_fold = sample(seq_len(nf), length(ug), replace = TRUE))
  d  <- left_join(d, gf, by = "grid_id")

  vlab <- character(nrow(subset_list)); vval <- numeric(nrow(subset_list))
  for (i in seq_len(nrow(subset_list))) {
    members <- GROUPS[unlist(subset_list[i, ])]
    preds   <- unlist(VAR_GROUPS[members], use.names = FALSE)
    vlab[i] <- paste(LETTERS_MAP[members], collapse = "")
    vval[i] <- cv_block_r2(preds, data = d)
  }
  v <- setNames(vval, vlab)
  full_R2 <- unname(v[paste(LETTERS_MAP[GROUPS], collapse = "")])
  phi <- shapley_from_v(v)
  cat(sprintf("  %-14s n=%-6d blocks=%-4d full_R2=%.3f  | C %.0f E %.0f L %.0f B %.0f %%\n",
              band_name, nrow(d), nb, full_R2,
              100*phi["CLIMATE"]/full_R2, 100*phi["EDAPHIC"]/full_R2,
              100*phi["LANDUSE"]/full_R2, 100*phi["BIOLOGICAL"]/full_R2))
  data.frame(band = band_name, lat_center = lat_center, n = nrow(d), n_blocks = nb,
             full_R2 = full_R2,
             CLIMATE = phi["CLIMATE"], EDAPHIC = phi["EDAPHIC"],
             LANDUSE = phi["LANDUSE"], BIOLOGICAL = phi["BIOLOGICAL"],
             row.names = NULL)
}

# ---- run: fine 15 deg bands (heatmap) + 4 biome bands (robust table) -----------
cat(sprintf("=== FINE %g-degree bands ===\n", FINE_STEP))
fine_edges <- seq(-60, 70, by = FINE_STEP)
fine <- do.call(rbind, lapply(seq_len(length(fine_edges) - 1), function(k) {
  lo <- fine_edges[k]; hi <- fine_edges[k + 1]
  decompose_band(cv_data[cv_data$lat >= lo & cv_data$lat < hi, , drop = FALSE],
                 sprintf("[%d,%d)", lo, hi), (lo + hi) / 2)
}))
fine$scheme <- "fine"

cat("\n=== 4 biome bands ===\n")
biome <- list("S-temperate" = c(-60, -23.5), "Tropics" = c(-23.5, 23.5),
              "N-temperate" = c(23.5, 50),   "Boreal-Arctic" = c(50, 90))
biome_df <- do.call(rbind, lapply(names(biome), function(nm) {
  lo <- biome[[nm]][1]; hi <- biome[[nm]][2]
  decompose_band(cv_data[cv_data$lat >= lo & cv_data$lat < hi, , drop = FALSE],
                 nm, (lo + hi) / 2)
}))
biome_df$scheme <- "biome4"

res <- rbind(fine, biome_df)
# within-band shares (% of that band's full-model R^2) + beyond-climate renorm
for (g in GROUPS) res[[paste0("pct_", g)]] <- 100 * res[[g]] / res$full_R2
bc <- res$EDAPHIC + res$LANDUSE + res$BIOLOGICAL
for (g in c("EDAPHIC","LANDUSE","BIOLOGICAL")) res[[paste0("bcpct_", g)]] <- 100 * res[[g]] / bc

write.csv(res, file.path(OUTPUT_DIR, "13l_latitude_shapley.csv"), row.names = FALSE)
saveRDS(res, file.path(OUTPUT_DIR, "13l_latitude_shapley.rds"))
cat("\nSaved outputs/13l_latitude_shapley.{csv,rds}\n")

# ============================================================================
# FIGURE: latitude x process-class heatmap of within-band Shapley share.
# Each domain column is a white->domain-hue ramp keyed to its % share; a right
# strip reports the band's full-model R^2 and n so unreliable bands are visible.
# ============================================================================
DOM     <- c("CLIMATE","EDAPHIC","LANDUSE","BIOLOGICAL")
DOM_LAB <- c("Climate","Edaphic","Land use","Biological")
# established manuscript domain palette (matches 13c_results_figure.R / 13k / Fig. 2c hue families)
DOM_HUE <- c(CLIMATE="#666666", EDAPHIC="#E41A1C", LANDUSE="#4DAF4A", BIOLOGICAL="#377EB8")
SHARE_MAX <- 70
ramp <- lapply(DOM_HUE, function(h) colorRampPalette(c("white", h))(100))
cellcol <- function(pct, dom) {
  idx <- pmin(pmax(round(pct / SHARE_MAX * 99) + 1, 1), 100); ramp[[dom]][idx]
}

f <- res[res$scheme == "fine", ]; f <- f[order(f$lat_center), ]
hh <- FINE_STEP / 2

png(file.path(PLOT_DIR, "13l_latitude_shapley.png"), width = 1700, height = 1500, res = 190)
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1.3))
par(mar = c(5, 5, 4, 1))
plot(NA, xlim = c(0.5, 4.5), ylim = c(min(f$lat_center) - hh, max(f$lat_center) + hh),
     xaxt = "n", xlab = "", ylab = "Latitude (deg)", las = 1,
     main = "Within-band attribution of log(tau) variance\n(local Shapley share by latitude)")
for (r in seq_len(nrow(f))) for (j in seq_along(DOM)) {
  pct <- f[[paste0("pct_", DOM[j])]][r]; y <- f$lat_center[r]
  rect(j - 0.5, y - hh, j + 0.5, y + hh, col = cellcol(pct, DOM[j]), border = "grey88")
  text(j, y, sprintf("%.0f", pct), cex = 0.9,
       col = if (pct > SHARE_MAX * 0.55) "white" else "grey20")
}
axis(1, at = seq_along(DOM), labels = DOM_LAB, cex.axis = 1.0)
mtext("share of within-band explained variance (%)", side = 1, line = 3, cex = 0.95)

# right strip: full-model R^2 (bar) and n per band
par(mar = c(5, 0.5, 4, 4))
plot(NA, xlim = c(0, max(f$full_R2) * 1.15), ylim = c(min(f$lat_center) - hh, max(f$lat_center) + hh),
     yaxt = "n", ylab = "", xlab = expression("band full-model " * R^2), main = "predictability")
abline(v = pretty(c(0, max(f$full_R2))), col = "grey92")
rect(0, f$lat_center - hh * 0.7, f$full_R2, f$lat_center + hh * 0.7, col = "grey65", border = NA)
text(max(f$full_R2) * 1.13, f$lat_center, sprintf("n=%s", format(f$n, big.mark = ",")),
     adj = 1, cex = 0.72, col = "grey35", xpd = NA)
dev.off()
cat("OK  ", file.path(PLOT_DIR, "13l_latitude_shapley.png"), "\n")

# ---- companion: beyond-climate composition (Edaphic/Land-use/Biological) --------
# The climate-free mix, renormalised to 100% within each band, pairing with Fig. 2c.
# Negative Shapley (suppression) is clamped to 0 before renormalising.
BC_HUE <- c(EDAPHIC = "#E41A1C", LANDUSE = "#4DAF4A", BIOLOGICAL = "#377EB8")
BC_LAB <- c("Edaphic", "Land use", "Biological")
Ep <- pmax(f$EDAPHIC, 0); Lp <- pmax(f$LANDUSE, 0); Bp <- pmax(f$BIOLOGICAL, 0)
comp <- cbind(EDAPHIC = Ep, LANDUSE = Lp, BIOLOGICAL = Bp) / (Ep + Lp + Bp) * 100

png(file.path(PLOT_DIR, "13l_beyondclimate_composition.png"), width = 1500, height = 1400, res = 190)
par(mar = c(5, 5, 4, 8))
plot(NA, xlim = c(0, 100), ylim = c(min(f$lat_center) - hh, max(f$lat_center) + hh),
     xlab = "share of BEYOND-climate explained variance (%)", ylab = "Latitude (deg)", las = 1,
     main = "Beyond-climate composition by latitude\n(Edaphic / Land use / Biological, renormalised)")
for (r in seq_len(nrow(f))) {
  x0 <- 0; y <- f$lat_center[r]
  for (g in names(BC_HUE)) {
    w <- comp[r, g]; if (is.na(w)) next
    rect(x0, y - hh * 0.9, x0 + w, y + hh * 0.9, col = BC_HUE[g], border = "white")
    if (w >= 8) text(x0 + w/2, y, sprintf("%.0f", w), col = "white", cex = 0.8)
    x0 <- x0 + w
  }
}
legend(x = 103, y = mean(range(f$lat_center)), legend = BC_LAB, fill = BC_HUE, border = NA,
       bty = "n", xpd = NA, cex = 1.1, title = "Process class")
dev.off()
cat("OK  ", file.path(PLOT_DIR, "13l_beyondclimate_composition.png"), "\n")
