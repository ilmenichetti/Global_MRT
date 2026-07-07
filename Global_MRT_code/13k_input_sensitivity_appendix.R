################################################################################
# 13k_input_sensitivity_appendix.R
#
# Presentation layer for the carbon-input sensitivity sweep (compute = 13j).
# Builds the appendix TABLE + a 4-panel FIGURE that carries the three-beat
# argument:
#   (A) the input change is a REAL, land-cover-structured shift in tau
#       (median tau by land cover under each convention) -- so the invariance
#       below is not because we nudged tau imperceptibly;
#   (B) the SIGNATURE: Shapley shares are flat across conventions except LandUse,
#       which is exactly the axis a land-cover-structured coefficient projects onto;
#   (C) the HEADLINE is invariant: beyond-climate gap barely moves across the four
#       conventions (Moran's I annotated);
#   (D) and it does not rest on one harvest intensity: gap is flat-to-rising across
#       the literature-anchored h range.
#
# One axis per panel (no dual-axis); Moran's I is annotation/table, never a second
# y-scale. Domain palette matches 13c/13h (manuscript Figs); convention palette is
# an Okabe-Ito triple (validated: CVD-safe, in-band).
#
# Input:  outputs/13j_input_sensitivity.csv
#         outputs/13j_harvest_hcurve.csv
#         outputs/12b_model_ready.rds   (for panel A: tau by land cover)
# Output: outputs/13k_input_sensitivity_table.tex   (booktabs, appendix-ready)
#         plots/step_13c_commonality/13k_input_sensitivity_appendix.png
#
# Author: Lorenzo   Date: 2026-07-07
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR   <- "./Global_MRT_code/plots/step_13c_commonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Established manuscript domain palette (matches 13c/13h). Climate = neutral grey
# by design (the reference domain); every domain line is also direct-labelled, so
# identity never rests on colour alone.
DOMAIN_COL <- c(Climate = "#666666", Edaphic = "#E41A1C",
                LandUse = "#4DAF4A", Biological = "#377EB8")
# Convention palette: Okabe-Ito triple, validator-passing (CVD-safe, in-band).
CONV_COL <- c(below = "#0072B2", total = "#009E73", harvest = "#E69F00")
# Only the three INPUT conventions are shown (belowground / total / harvest); the
# peat-drop case is a different robustness axis (a data filter, not an input
# definition) and is not shown here unless peat is filtered consistently
# downstream. It stays flag-not-mask on the maps + a one-line robustness note.
CONV_LAB <- c(below = "Belowground\n(main)", total = "Total NPP",
              harvest = "Total − harvest\n(h=0.5)")
CONV_ORDER <- c("below", "total", "harvest")

theme_app <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.7, colour = "grey30"),
        legend.position = "none")

# =============================================================================
# READ SWEEP RESULTS
# =============================================================================
tab_all <- read.csv(file.path(OUTPUT_DIR, "13j_input_sensitivity.csv"),
                    stringsAsFactors = FALSE)   # all conventions incl. nopeat
hc  <- read.csv(file.path(OUTPUT_DIR, "13j_harvest_hcurve.csv"),
                stringsAsFactors = FALSE)
tab <- tab_all
tab$convention <- factor(tab$convention, levels = CONV_ORDER)
tab <- tab %>% filter(!is.na(convention)) %>% arrange(convention)  # 3 input conventions
gap_main <- tab$beyond_gap_pp[tab$convention == "below"]

# =============================================================================
# PANEL A -- tau by land cover under each convention (the perturbation is REAL)
# =============================================================================
# Re-derive tau three ways per row (same recipe as 13c/13j) and take the median
# per land-cover class. Only vegetated classes with a defined allocation fraction.
d <- readRDS(file.path(OUTPUT_DIR, "12b_model_ready.rds"))
keep_lc <- c("trees", "shrubs", "grassland", "cropland",
             "wetland", "mangroves", "moss", "bare")
LC_NICE <- c(trees = "Tree cover", shrubs = "Shrubland", grassland = "Grassland",
             cropland = "Cropland", wetland = "Herb. wetland",
             mangroves = "Mangroves", moss = "Moss/lichen", bare = "Bare/sparse")

f   <- d$bnpp_fraction
eff_h <- ifelse(!is.na(d$lc_dominant) & d$lc_dominant == "trees",
                f + (1 - f) * (1 - 0.5), 1)
tau_df <- data.frame(
  lc      = d$lc_dominant,
  below   = d$SOC_stock_g_m2 / (d$NPP_g_C_m2_yr * f),
  total   = d$SOC_stock_g_m2 /  d$NPP_g_C_m2_yr,
  harvest = d$SOC_stock_g_m2 / (d$NPP_g_C_m2_yr * eff_h)
) %>%
  filter(lc %in% keep_lc) %>%
  pivot_longer(c(below, total, harvest), names_to = "convention",
               values_to = "tau") %>%
  filter(is.finite(tau), tau > 0, tau < 1000) %>%
  group_by(lc, convention) %>%
  summarise(tau_med = median(tau, na.rm = TRUE), .groups = "drop")

lc_ord <- tau_df %>% filter(convention == "below") %>% arrange(tau_med) %>% pull(lc)
tau_df <- tau_df %>%
  mutate(lc = factor(lc, levels = lc_ord),
         convention = factor(convention, levels = c("below", "total", "harvest")))
seg_df <- tau_df %>% group_by(lc) %>%
  summarise(lo = min(tau_med), hi = max(tau_med), .groups = "drop")

pA <- ggplot(tau_df, aes(tau_med, lc)) +
  geom_segment(data = seg_df, aes(x = lo, xend = hi, y = lc, yend = lc),
               inherit.aes = FALSE, colour = "grey75", linewidth = 0.9) +
  geom_point(aes(colour = convention), size = 2.6) +
  scale_colour_manual(values = CONV_COL,
                      labels = c(below = "Belowground (main)", total = "Total NPP",
                                 harvest = "Total − harvest"),
                      name = NULL) +
  scale_y_discrete(labels = LC_NICE) +
  scale_x_log10() +
  labs(title = "A  The input change is a real, land-cover-structured shift in τ",
       subtitle = "Median τ per land-cover class; harvest departs from total only on tree cover",
       x = "Median τ (years, log scale)", y = NULL) +
  theme_app + theme(legend.position = c(0.98, 0.03),
                    legend.justification = c(1, 0),
                    legend.background = element_rect(fill = "white", colour = "grey85"),
                    legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7.5))

# =============================================================================
# PANEL B -- Shapley shares across conventions (the SIGNATURE: only LandUse moves)
# =============================================================================
shp <- tab %>% filter(convention %in% c("below", "total", "harvest")) %>%
  select(convention, Climate, Edaphic, Biological, LandUse) %>%
  pivot_longer(-convention, names_to = "domain", values_to = "pct") %>%
  mutate(convention = factor(convention, levels = c("below", "total", "harvest")),
         cx = as.integer(convention),
         domain = factor(domain, levels = c("Climate", "Edaphic", "Biological", "LandUse")))
lab_r <- shp %>% filter(convention == "harvest")

pB <- ggplot(shp, aes(cx, pct, colour = domain, group = domain)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.3) +
  geom_text(data = lab_r, aes(label = sprintf("%s  %.0f%%", domain, pct)),
            hjust = 0, nudge_x = 0.08, size = 3, fontface = "bold") +
  scale_colour_manual(values = DOMAIN_COL, guide = "none") +
  scale_x_continuous(breaks = 1:3,
                     labels = c("Belowground\n(main)", "Total NPP", "Total −\nharvest"),
                     limits = c(0.9, 4.4)) +
  labs(title = "B  Only LandUse responds — the predicted signature",
       subtitle = "Shapley shares are flat except LandUse (where a land-cover input lands)",
       x = NULL, y = "Share of attributable variance (%)") +
  theme_app

# =============================================================================
# PANEL C -- beyond-climate gap across conventions (HEADLINE is invariant)
# =============================================================================
gap_df <- tab %>%
  mutate(clab = factor(as.character(convention), levels = CONV_ORDER))
pC <- ggplot(gap_df, aes(clab, beyond_gap_pp)) +
  geom_hline(yintercept = gap_main, linetype = "dashed", colour = "grey65") +
  geom_col(width = 0.62, fill = "#4477AA") +
  geom_text(aes(label = sprintf("%.1f", beyond_gap_pp)),
            vjust = -0.5, size = 3.1, fontface = "bold") +
  geom_text(aes(y = 0.5, label = sprintf("I=%.2f", morans_I)),
            colour = "white", size = 2.7) +
  scale_x_discrete(labels = function(x) gsub("\n", " ", CONV_LAB[x])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(title = "C  Beyond-climate gap is invariant across the three conventions",
       subtitle = "Gap = full − climate-only block-CV R² (pp); Moran's I of the residual in white",
       x = NULL, y = "Beyond-climate gap (pp)") +
  theme_app + theme(axis.text.x = element_text(size = 7.6))

# =============================================================================
# PANEL D -- gap vs harvest intensity h (does not rest on one h)
# =============================================================================
anchors <- data.frame(h = c(0.2, 0.3, 0.45, 0.6),
                      lab = c("product\navg", "boreal\nstem", "temperate\nstem", "whole\ntree"))
pD <- ggplot(hc, aes(harvest_h, beyond_gap_pp)) +
  geom_hline(yintercept = gap_main, linetype = "dashed", colour = "grey65") +
  geom_vline(data = anchors, aes(xintercept = h), linetype = "dotted", colour = "grey80") +
  geom_text(data = anchors, aes(x = h, y = max(hc$beyond_gap_pp) + 0.15, label = lab),
            inherit.aes = FALSE, size = 2.4, colour = "grey45", lineheight = 0.85, vjust = 1) +
  geom_line(linewidth = 1.1, colour = "#E69F00") +
  geom_point(size = 2.4, colour = "#E69F00") +
  annotate("text", x = min(hc$harvest_h), y = gap_main, vjust = 1.5, hjust = 0,
           label = sprintf("main gap %.1f pp", gap_main), size = 2.6, colour = "grey45") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(title = "D  Harvest convention: boreal signal does not rest on one h",
       subtitle = "Third convention (Total − harvest), peat kept; gap vs export fraction h; anchors = literature h",
       x = "Harvest export fraction h (of aboveground NPP, tree pixels)",
       y = "Beyond-climate gap (pp)") +
  theme_app

# =============================================================================
# COMBINE + SAVE
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Carbon-input sensitivity: the zonality is invariant across the plausible input range",
    theme = theme(plot.title = element_text(face = "bold", size = 12)))
f_png <- file.path(PLOT_DIR, "13k_input_sensitivity_appendix.png")
ggsave(f_png, fig, width = 11, height = 8.4, dpi = 200)
cat("OK  ", f_png, "\n")

# =============================================================================
# APPENDIX TABLE (booktabs LaTeX + console)
# =============================================================================
tt <- tab %>% transmute(
  Convention = c(below = "Belowground (main)", total = "Total NPP",
                 harvest = "Total $-$ harvest ($h{=}0.5$)")[as.character(convention)],
  n = formatC(n_common, big.mark = ",", format = "d"),
  R2full = sprintf("%.3f", full_R2),
  R2clim = sprintf("%.3f", clim_R2),
  Gap = sprintf("%.1f", beyond_gap_pp),
  Moran = sprintf("%.3f", morans_I),
  C = sprintf("%.1f", Climate), E = sprintf("%.1f", Edaphic),
  B = sprintf("%.1f", Biological), L = sprintf("%.1f", LandUse))

tex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\small",
  "\\caption{Carbon-input sensitivity of the beyond-climate signal. The soil carbon input coefficient is re-derived on the same NPP field under three conventions (belowground, total NPP, and total minus a generous forest-harvest export $h$). Block-CV $R^2$, the beyond-climate gap (full $-$ climate-only), Moran's $I$ of the beyond-climate residual, and the four Shapley shares are all but unchanged; only LandUse responds, as expected for a land-cover-structured coefficient.}",
  "\\label{tab:input_sensitivity}",
  "\\begin{tabular}{lrccccrrrr}",
  "\\toprule",
  " & & \\multicolumn{2}{c}{Block-CV $R^2$} & Gap & Moran's & \\multicolumn{4}{c}{Shapley share (\\%)} \\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){7-10}",
  "Convention & $n$ & full & clim. & (pp) & $I$ & C & E & B & L \\\\",
  "\\midrule",
  paste0(apply(tt, 1, function(r) paste(paste(r, collapse = " & "), "\\\\")), collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}")
f_tex <- file.path(OUTPUT_DIR, "13k_input_sensitivity_table.tex")
writeLines(tex, f_tex)
cat("OK  ", f_tex, "\n\n")

cat("=== Appendix table (preview) ===\n")
print(as.data.frame(tt), row.names = FALSE)

# =============================================================================
# PEAT-SENSITIVITY TABLE (separate SI table; peat is a data-filter axis, not an
# input definition, so it is not in the input-convention table above)
# =============================================================================
peat_rows <- tab_all %>% filter(convention %in% c("below", "nopeat"))
if (all(c("below", "nopeat") %in% peat_rows$convention)) {
  peat_rows <- peat_rows[match(c("below", "nopeat"), peat_rows$convention), ]
  n_drop <- peat_rows$n_common[peat_rows$convention == "below"] -
            peat_rows$n_common[peat_rows$convention == "nopeat"]
  pp <- peat_rows %>% transmute(
    Run = c(below = "Main (Histosols retained)",
            nopeat = "Histosols dropped")[convention],
    n = formatC(n_common, big.mark = ",", format = "d"),
    R2full = sprintf("%.3f", full_R2), R2clim = sprintf("%.3f", clim_R2),
    Gap = sprintf("%.1f", beyond_gap_pp), Moran = sprintf("%.3f", morans_I),
    C = sprintf("%.1f", Climate), E = sprintf("%.1f", Edaphic),
    B = sprintf("%.1f", Biological), L = sprintf("%.1f", LandUse))
  ptex <- c(
    "\\begin{table}[t]", "\\centering", "\\small",
    sprintf("\\caption{Peatland sensitivity of the beyond-climate signal (flag-not-mask decision). Re-running the full commonality decomposition with Histosol-classified observations dropped (%s common rows removed, $\\sim$1.3\\%% of the full dataset) leaves the beyond-climate gap and its spatial organisation essentially unchanged; the gap in fact decreases slightly, i.e. organic soils are not inflating the signal. Histosols are therefore retained in the fit and flagged (not interpreted) on the maps.}", formatC(n_drop, big.mark = ",", format = "d")),
    "\\label{tab:peat_sensitivity}",
    "\\begin{tabular}{lrccccrrrr}", "\\toprule",
    " & & \\multicolumn{2}{c}{Block-CV $R^2$} & Gap & Moran's & \\multicolumn{4}{c}{Shapley share (\\%)} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){7-10}",
    "Run & $n$ & full & clim. & (pp) & $I$ & C & E & B & L \\\\", "\\midrule",
    paste0(apply(pp, 1, function(r) paste(paste(r, collapse = " & "), "\\\\")), collapse = "\n"),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  f_ptex <- file.path(OUTPUT_DIR, "13k_peat_sensitivity_table.tex")
  writeLines(ptex, f_ptex)
  cat("\nOK  ", f_ptex, "\n")
  cat("=== Peat-sensitivity table (preview) ===\n")
  print(as.data.frame(pp), row.names = FALSE)
} else {
  cat("\n(peat table skipped: need both 'below' and 'nopeat' rows; run MRT_13C_PEAT=drop 13c)\n")
}

cat("\nSTEP 13k COMPLETE.\n")
