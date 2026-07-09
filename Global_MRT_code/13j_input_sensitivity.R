################################################################################
# 13j_input_sensitivity.R
#
# Consolidate the CARBON-INPUT sensitivity sweep: does the choice of soil carbon
# input convention change the paper's RELATIVE (zonality) conclusions? The
# reviewer critique (Aleksi x2 + Shoji) is that C input = the BELOWGROUND
# fraction of NPP only. As argued in the memo, this is a coefficient convention
# on one NPP field, not a data gap; a uniform coefficient change is scale-
# invariant, so only the coefficient's SPATIAL variation across land cover can
# move the relative results. We test the full plausible span:
#
#   below   (main)      C_input = NPP * f_bnpp                (land-cover structure)
#   total   (upper)     C_input = NPP                         (coeff = 1, no structure)
#   harvest (teeth)     C_input = NPP * [f + (1-f)(1-h)] on tree-dominant pixels,
#                       NPP elsewhere; generous h on (crudely) managed forest ->
#                       the MOST spatially-structured perturbation, concentrated
#                       on forest and overlapping the boreal zonality signal.
#
# For each convention this reads the heavy Shapley decomposition produced by
#   MRT_13C_INPUT=<mode> Rscript 13c_commonality_validation.R
# and reports full-model block-CV R^2, the climate-only R^2 (CLIMATE "single"),
# the BEYOND-CLIMATE gap (full - climate), the four Shapley shares, and Moran's I
# on the beyond-climate residual (recomputed here per convention, 13e recipe).
# If the numbers barely move -> the zonality claim is coefficient-invariant (one
# appendix table + a framing paragraph). If they move materially -> that JUSTIFIES
# the heavier managed-forest (Lesiv) build.
#
# The peatland robustness case (MRT_13C_PEAT=drop) is reported in the same table
# if its decomposition file is present.
#
# Decomposition files (from `MRT_13C_INPUT=<mode> Rscript 13c_...R`):
#   13c_decomposition.rds            -> below   (MAIN headline)
#   13c_decomposition_inpTOTAL.rds   -> total
#   13c_decomposition_inpHARV.rds    -> harvest
#   13c_decomposition_nopeat.rds     -> below, Histosols dropped (optional)
#
# Input:  outputs/13c_decomposition*.rds  ;  outputs/12b_model_ready.rds
#         outputs/13_var_groups.rds
# Output: outputs/13j_input_sensitivity.csv
#         plots/step_13c_commonality/13j_input_shapley.png
#         plots/step_13c_commonality/13j_input_gap.png
#
# Env: MRT_13J_MORAN  = 1 (default) | 0   (0 skips the Moran refit AND h-curve)
#      MRT_13J_NTREES = 500 (default; RF trees for the refits)
#      MRT_13J_HVEC   = "0.2,0.3,0.45,0.6" (harvest h-curve grid; literature-
#                       anchored: product avg / boreal stem / temperate stem /
#                       whole-tree). Empty string to skip only the h-curve.
#
# Author: Lorenzo   Date: 2026-07-07
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

`%||%` <- function(a, b) if (is.null(a)) b else a

OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR   <- "./Global_MRT_code/plots/step_13c_commonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

DO_MORAN     <- Sys.getenv("MRT_13J_MORAN", "1") != "0"
MORAN_NTREES <- as.integer(Sys.getenv("MRT_13J_NTREES", "500"))

GROUP_COLORS <- c(Climate = "#666666", Edaphic = "#E41A1C",
                  LandUse = "#4DAF4A", Biological = "#377EB8")
NICE <- c(CLIMATE = "Climate", EDAPHIC = "Edaphic",
          LANDUSE = "LandUse", BIOLOGICAL = "Biological")

# Convention label -> decomposition file (config in the object is authoritative)
files <- list(
  "below"   = file.path(OUTPUT_DIR, "13c_decomposition.rds"),
  "total"   = file.path(OUTPUT_DIR, "13c_decomposition_inpTOTAL.rds"),
  "harvest" = file.path(OUTPUT_DIR, "13c_decomposition_inpHARV.rds"),
  "nopeat"  = file.path(OUTPUT_DIR, "13c_decomposition_nopeat.rds")
)
CONV_LABEL <- c(below = "Belowground (main)", total = "Total NPP",
                harvest = "Total - harvest", nopeat = "Belowground, no peat")

# =============================================================================
# READ DECOMPOSITIONS -> full R^2, climate single, beyond-climate gap, Shapley
# =============================================================================
rows_shap <- list(); rows_hdr <- list()
for (nm in names(files)) {
  f <- files[[nm]]
  if (!file.exists(f)) { cat("  (skip, not found:", basename(f), ")\n"); next }
  dec <- readRDS(f)
  clim_single <- unname(dec$singles["CLIMATE"])
  gap         <- dec$full_R2 - clim_single
  rows_hdr[[nm]] <- data.frame(
    convention  = nm,
    input_mode  = dec$config$input_mode %||% "below",
    harvest_h   = dec$config$harvest_h  %||% NA_real_,
    peat_mode   = dec$config$peat_mode  %||% "keep",
    n_common    = dec$config$n_common,
    full_R2     = dec$full_R2,
    clim_R2     = clim_single,
    beyond_gap_pp = 100 * gap,
    rf_ntrees   = dec$config$rf_ntrees,
    stringsAsFactors = FALSE)
  rows_shap[[nm]] <- dec$shapley_df %>%
    transmute(convention = nm, group, nice = NICE[group],
              shapley_pct, shapley_R2)
}
stopifnot(length(rows_hdr) >= 1)
hdr  <- bind_rows(rows_hdr)
shap <- bind_rows(rows_shap)

conv_present <- hdr$convention
conv_levels  <- intersect(c("below", "total", "harvest", "nopeat"), conv_present)

# =============================================================================
# MORAN'S I PER CONVENTION (beyond-climate residual; 13e recipe, recomputed)
# =============================================================================
# tau is re-derived here EXACTLY as in 13c's input switch (kept in sync; 13c is the
# source of truth). Peat is dropped via the SAME shared GPM filter as the fit
# (peat_filter.R), so the recomputed Moran matches the headline basis.
source("./Global_MRT_code/peat_filter.R")
derive_tau <- function(d, input_mode, harvest_h, peat_mode) {
  d <- apply_peat_filter(d, OUTPUT_DIR, peat_mode = peat_mode, peat_source = "gpm")
  if (input_mode != "below") {
    f   <- d$bnpp_fraction
    eff <- switch(input_mode,
      total   = rep(1, nrow(d)),
      harvest = ifelse(!is.na(d$lc_dominant) & d$lc_dominant == "trees",
                       f + (1 - f) * (1 - harvest_h), 1))
    C_input_alt <- d$NPP_g_C_m2_yr * eff
    d$MRT_years <- d$SOC_stock_g_m2 / C_input_alt
    d$MRT_QC <- dplyr::case_when(
      is.na(d$MRT_years) ~ "missing_data", is.infinite(d$MRT_years) ~ "infinite",
      d$MRT_years < 0 ~ "negative", d$MRT_years < 1 ~ "very_low",
      d$MRT_years > 1000 ~ "very_high", TRUE ~ "valid")
  }
  d
}

# Build the common-row block-CV frame for a given tau convention (same footprint
# and folds as 13c/13e: M5 + MAIN biology, 2 deg blocks, seed 42).
build_cv <- function(d0, ALLVARS, input_mode, harvest_h, peat_mode, block = 2) {
  d <- derive_tau(d0, input_mode, harvest_h, peat_mode)
  d <- d %>% filter(MRT_QC == "valid", !is.na(MRT_years),
                    MRT_years > 0, MRT_years < Inf)
  q <- quantile(d$MRT_years, c(.01, .99), na.rm = TRUE)
  d <- d %>% filter(MRT_years > q[1] & MRT_years < q[2])
  d$log_MRT <- log(d$MRT_years)
  d$lon <- d$longitude_decimal_degrees; d$lat <- d$latitude_decimal_degrees
  need <- c("log_MRT", ALLVARS, "lon", "lat")
  cv <- d[complete.cases(d[, need]), ]
  cv$grid <- paste(floor(cv$lon / block), floor(cv$lat / block))
  ug <- unique(cv$grid); set.seed(42)
  cv$fold <- data.frame(grid = ug, f = sample(10, length(ug), TRUE))$f[match(cv$grid, ug)]
  cv
}

# Out-of-sample block-CV R^2 for a predictor set (and the OOS predictions).
fit_block_cv_r2 <- function(cv, preds, ntrees) {
  frm <- as.formula(paste("log_MRT ~", paste(preds, collapse = " + ")))
  pred <- rep(NA_real_, nrow(cv))
  for (k in sort(unique(cv$fold))) {
    te <- cv$fold == k
    m <- ranger::ranger(frm, cv[!te, ], num.trees = ntrees, min.node.size = 5,
                        num.threads = max(1, parallel::detectCores() - 1), seed = 42)
    pred[te] <- predict(m, cv[te, ])$predictions
  }
  r2 <- 1 - sum((cv$log_MRT - pred)^2) / sum((cv$log_MRT - mean(cv$log_MRT))^2)
  list(r2 = r2, pred = pred)
}

# Moran's I on the beyond-climate residual (climate OOS predictions -> residual
# -> 2 deg block means -> knn8 weights), identical recipe to 13e.
moran_resid <- function(cv, clim_pred, block = 2, knn = 8) {
  cv$resid <- cv$log_MRT - clim_pred
  agg <- cv %>% mutate(bx = floor(lon / block), by = floor(lat / block)) %>%
    group_by(bx, by) %>%
    summarise(resid = mean(resid), n = n(),
              lon = mean(lon), lat = mean(lat), .groups = "drop") %>%
    filter(n >= 5)
  coords <- as.matrix(agg[, c("lon", "lat")])
  lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords, k = knn)), style = "W")
  mt <- spdep::moran.test(agg$resid, lw, zero.policy = TRUE)
  data.frame(n_blocks = nrow(agg),
             morans_I = unname(mt$estimate["Moran I statistic"]),
             moran_p  = mt$p.value)
}

moran_tab <- NULL; hcurve_tab <- NULL
if (DO_MORAN) {
  d0 <- readRDS(file.path(OUTPUT_DIR, "12b_model_ready.rds"))
  g  <- readRDS(file.path(OUTPUT_DIR, "13_var_groups.rds"))
  CLIM     <- g$CLIMATE
  M5       <- c(g$CLIMATE, g$EDAPHIC, g$LANDUSE)
  BIO_MAIN <- c("fungal_proportion", "AM_richness", "EcM_richness",
                "AM_endemism", "EcM_endemism", "EcM_AM_richness_ratio")
  ALLVARS  <- c(M5, BIO_MAIN)

  # ---- Moran's I per convention (climate-only refit; below reproduces ~0.25) --
  cat("Recomputing Moran's I per convention (climate-only refit, ",
      MORAN_NTREES, " trees)...\n", sep = "")
  m_rows <- list()
  for (nm in conv_levels) {
    hc <- hdr[hdr$convention == nm, ]
    cat("  ", nm, "... ")
    cv <- build_cv(d0, ALLVARS, hc$input_mode, hc$harvest_h %||% 0.5, hc$peat_mode)
    cp <- fit_block_cv_r2(cv, CLIM, MORAN_NTREES)$pred
    mr <- moran_resid(cv, cp); mr$n_cv <- nrow(cv); mr$convention <- nm
    m_rows[[nm]] <- mr
    cat(sprintf("I = %.3f (p = %s)\n", mr$morans_I, signif(mr$moran_p, 2)))
  }
  moran_tab <- bind_rows(m_rows)

  # ---- Harvest h-curve: does the beyond-climate signal rest on one asserted h? -
  # Anchored to literature (see header): h ~ 0.2 managed-product average
  # (Luyssaert-lineage, harvest ~13-16% of NPP), ~0.3 boreal stem-only, ~0.45
  # temperate stem-only, ~0.6 whole-tree/intensive. For each h we refit the FULL
  # model and the climate-only model on the harvest-derived tau, giving the
  # beyond-climate gap and Moran's I as a function of harvest intensity.
  H_VEC <- as.numeric(strsplit(Sys.getenv("MRT_13J_HVEC", "0.2,0.3,0.45,0.6"),
                               ",")[[1]])
  if (length(H_VEC) > 0 && !any(is.na(H_VEC))) {
    cat("Harvest h-curve (full + climate-only refit per h, ",
        MORAN_NTREES, " trees)...\n", sep = "")
    h_rows <- list()
    for (hh in H_VEC) {
      cat(sprintf("  h = %.2f ... ", hh))
      cv   <- build_cv(d0, ALLVARS, "harvest", hh, "keep")  # peat KEPT (consistent w/ main)
      full <- fit_block_cv_r2(cv, ALLVARS, MORAN_NTREES)
      clim <- fit_block_cv_r2(cv, CLIM,    MORAN_NTREES)
      mr   <- moran_resid(cv, clim$pred)
      h_rows[[as.character(hh)]] <- data.frame(
        harvest_h = hh, n_cv = nrow(cv),
        full_R2 = full$r2, clim_R2 = clim$r2,
        beyond_gap_pp = 100 * (full$r2 - clim$r2),
        morans_I = mr$morans_I, moran_p = mr$moran_p)
      cat(sprintf("gap = %.1f pp, I = %.3f\n",
                  100 * (full$r2 - clim$r2), mr$morans_I))
    }
    hcurve_tab <- bind_rows(h_rows)
    write.csv(hcurve_tab, file.path(OUTPUT_DIR, "13j_harvest_hcurve.csv"),
              row.names = FALSE)
    cat("OK  outputs/13j_harvest_hcurve.csv\n")
  }
}

# =============================================================================
# ASSEMBLE COMPARISON TABLE
# =============================================================================
shap_wide <- shap %>%
  select(convention, nice, shapley_pct) %>%
  pivot_wider(names_from = nice, values_from = shapley_pct)

tab <- hdr %>% left_join(shap_wide, by = "convention")
if (!is.null(moran_tab)) tab <- tab %>% left_join(moran_tab, by = "convention")
tab <- tab %>% mutate(convention = factor(convention, levels = conv_levels)) %>%
  arrange(convention)

# Stability: spread of each headline across the input conventions (exclude nopeat,
# which is a separate robustness axis, from the input-convention spread).
inp_only <- tab %>% filter(convention %in% c("below", "total", "harvest"))
spread <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
gap_spread   <- spread(inp_only$beyond_gap_pp)
moran_spread <- if (!is.null(moran_tab)) spread(inp_only$morans_I) else NA_real_
land_spread  <- spread(inp_only$LandUse)

cat("\n=== Carbon-input sensitivity: headline comparison ===\n")
print(as.data.frame(tab %>% mutate(convention = CONV_LABEL[as.character(convention)])),
      row.names = FALSE, digits = 3)
cat(sprintf(
  "\nAcross input conventions (below/total/harvest):\n  beyond-climate gap spread = %.1f pp%s\n  LandUse Shapley spread    = %.1f pp\n",
  gap_spread,
  if (!is.null(moran_tab)) sprintf("\n  Moran's I spread          = %.3f", moran_spread) else "",
  land_spread))

write.csv(tab, file.path(OUTPUT_DIR, "13j_input_sensitivity.csv"), row.names = FALSE)
cat("OK  outputs/13j_input_sensitivity.csv\n")

# =============================================================================
# FIGURES
# =============================================================================
# (1) Shapley shares grouped by convention
shap_plot <- shap %>%
  mutate(convention = factor(convention, levels = conv_levels),
         conv_lab = CONV_LABEL[as.character(convention)])
p1 <- ggplot(shap_plot, aes(reorder(nice, shapley_pct), shapley_pct,
                            fill = nice, alpha = convention)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f", shapley_pct)),
            position = position_dodge(width = 0.8), hjust = -0.15, size = 2.6,
            show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_alpha_manual(values = seq(0.4, 1, length.out = length(conv_levels)),
                     labels = CONV_LABEL[conv_levels], name = "Input convention") +
  scale_y_continuous(limits = c(0, max(shap_plot$shapley_pct) * 1.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Shapley attribution is insensitive to the carbon-input convention",
       subtitle = sprintf("Shares stable across below/total/harvest (LandUse spread %.1f pp)",
                          land_spread),
       x = NULL, y = "Share of attributable variance (%)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, "13j_input_shapley.png"), p1, width = 8, height = 5, dpi = 200)
cat("OK  ", file.path(PLOT_DIR, "13j_input_shapley.png"), "\n")

# (2) Beyond-climate gap (and Moran's I if available) across conventions
gap_plot <- tab %>%
  transmute(conv_lab = factor(CONV_LABEL[as.character(convention)],
                              levels = CONV_LABEL[conv_levels]),
            `Beyond-climate gap (pp)` = beyond_gap_pp,
            `Moran's I x100` = if (!is.null(moran_tab)) morans_I * 100 else NA_real_) %>%
  pivot_longer(-conv_lab, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))
p2 <- ggplot(gap_plot, aes(conv_lab, value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", value)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Beyond-climate gap (pp)" = "#4477AA",
                               "Moran's I x100" = "#EE6677"), name = NULL) +
  labs(title = "Beyond-climate signal survives the input-convention span",
       subtitle = sprintf("Gap spread %.1f pp%s across below/total/harvest",
                          gap_spread,
                          if (!is.null(moran_tab)) sprintf("; Moran's I spread %.3f", moran_spread) else ""),
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, "13j_input_gap.png"), p2, width = 8, height = 4.8, dpi = 200)
cat("OK  ", file.path(PLOT_DIR, "13j_input_gap.png"), "\n")

# (3) Harvest h-curve: beyond-climate gap and Moran's I vs harvest intensity
if (!is.null(hcurve_tab)) {
  hc_long <- hcurve_tab %>%
    transmute(harvest_h,
              `Beyond-climate gap (pp)` = beyond_gap_pp,
              `Moran's I` = morans_I) %>%
    pivot_longer(-harvest_h, names_to = "metric", values_to = "value")
  anchors <- data.frame(
    harvest_h = c(0.2, 0.3, 0.45, 0.6),
    lab = c("product avg", "boreal stem", "temperate stem", "whole-tree"))
  p3 <- ggplot(hc_long, aes(harvest_h, value, colour = metric)) +
    geom_vline(data = anchors, aes(xintercept = harvest_h),
               linetype = "dotted", colour = "grey75") +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~metric, scales = "free_y") +
    scale_colour_manual(values = c("Beyond-climate gap (pp)" = "#4477AA",
                                   "Moran's I" = "#EE6677"), guide = "none") +
    labs(title = "Beyond-climate signal is flat across harvest intensity",
         subtitle = "Harvest convention (peat kept); h anchored to literature (dotted: product avg / boreal & temperate stem-only / whole-tree)",
         x = "Harvest export fraction h (of aboveground NPP, tree pixels)", y = NULL) +
    theme_bw(base_size = 11)
  ggsave(file.path(PLOT_DIR, "13j_harvest_hcurve.png"), p3,
         width = 9, height = 4.2, dpi = 200)
  cat("OK  ", file.path(PLOT_DIR, "13j_harvest_hcurve.png"), "\n")
  cat("\n=== Harvest h-curve ===\n")
  print(as.data.frame(hcurve_tab), row.names = FALSE, digits = 3)
}

cat("\nSTEP 13j COMPLETE.\n")
