################################################################################
# 13c_results_figure.R
#
# Build the consolidated variance-decomposition figure for the Results section,
# from the outputs of 13c_commonality_validation.R. Narrative order (endorsed):
#   (a) SINGLES   - each mechanistic domain alone (the opener)
#   (b) BUILD-UP  - all 15 coalitions by size: the big jump is 1->2 domains,
#                   then saturation = "correlated, not redundant" (overlap)
#   (c) SHAPLEY   - symmetric attribution of the full-model predictable variance
#                   (the headline: climate ~1/3, non-climate ~2/3)
#
# Also emits a MAIN (fungal+SPUN, 6-layer) vs SENSITIVITY (all-9) comparison.
#
# Input:  ./Global_MRT_code/outputs/13c_decomposition.rds        (MAIN headline)
#         ./Global_MRT_code/outputs/13c_decomposition_allbio.rds  (sensitivity)
# Output: ./Global_MRT_code/plots/13c_decomposition.png
#         ./Global_MRT_code/plots/13c_decomposition_sensitivity.png
#
# Author: Lorenzo
# Date: 2026-06-26
################################################################################

library(dplyr)
library(ggplot2)
library(patchwork)

OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR   <- "./Global_MRT_code/plots/step_13c_commonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_COLORS <- c(Climate = "#666666", Edaphic = "#E41A1C",
                  LandUse = "#4DAF4A", Biological = "#377EB8")
NICE <- c(CLIMATE = "Climate", EDAPHIC = "Edaphic",
          LANDUSE = "LandUse", BIOLOGICAL = "Biological")

dec      <- readRDS(file.path(OUTPUT_DIR, "13c_decomposition.rds"))
sens_f   <- file.path(OUTPUT_DIR, "13c_decomposition_allbio.rds")
have_sens<- file.exists(sens_f)
if (have_sens) dec_sens <- readRDS(sens_f)

full_R2  <- dec$full_R2
n_common <- dec$config$n_common

# =============================================================================
# PANEL (a) SINGLES -- each domain alone
# =============================================================================

singles_df <- data.frame(group = names(dec$singles),
                         R2 = as.numeric(dec$singles)) %>%
  mutate(nice = NICE[group])

pa <- ggplot(singles_df, aes(reorder(nice, R2), R2, fill = nice)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", R2)), hjust = -0.15, size = 3.4) +
  coord_flip() +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, max(singles_df$R2) * 1.18),
                     expand = expansion(mult = c(0, 0))) +
  labs(tag = "a", title = "Each domain alone",
       subtitle = "Single-domain block-CV R²",
       x = NULL, y = expression("Block-CV " * R^2)) +
  theme_bw(base_size = 11) +
  theme(plot.tag = element_text(face = "bold"))

# =============================================================================
# PANEL (b) BUILD-UP -- all 15 coalitions by number of domains
# =============================================================================
# Each coalition is a point at x = n_groups; the full model (CELB) is marked.
# The jump from 1->2 then saturation shows overlap (correlated, not redundant).

coal <- dec$coalition_df %>%
  mutate(is_full = coalition == "CELB")

singles_sum <- sum(dec$singles)

pb <- ggplot(coal, aes(n_groups, R2_block_cv)) +
  # naive additive expectation (sum of singles) for reference -> shows overlap
  geom_hline(yintercept = singles_sum, linetype = "dotted", colour = "grey55") +
  annotate("text", x = 1, y = singles_sum, vjust = -0.6, hjust = 0,
           size = 3, colour = "grey40",
           label = sprintf("Σ singles = %.2f (if independent)", singles_sum)) +
  geom_line(stat = "summary", fun = max, colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = is_full, size = is_full)) +
  ggrepel::geom_text_repel(aes(label = coalition), size = 2.7,
                           max.overlaps = 20, segment.size = 0.2) +
  scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = "#D55E00"),
                      guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.2), guide = "none") +
  scale_x_continuous(breaks = 1:4) +
  labs(tag = "b", title = "Adding domains: big jump, then saturation",
       subtitle = "All 15 coalitions (C/E/L/B); orange = full model",
       x = "Number of domains in model", y = expression("Block-CV " * R^2)) +
  theme_bw(base_size = 11) +
  theme(plot.tag = element_text(face = "bold"))

# =============================================================================
# PANEL (c) SHAPLEY -- attribution of full-model predictable variance
# =============================================================================

shap_df <- dec$shapley_df %>%
  mutate(nice = NICE[group])

pc <- ggplot(shap_df, aes(reorder(nice, shapley_R2), shapley_pct, fill = nice)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%  (%.3f)", shapley_pct, shapley_R2)),
            hjust = -0.1, size = 3.4) +
  coord_flip() +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, max(shap_df$shapley_pct) * 1.32),
                     expand = expansion(mult = c(0, 0))) +
  labs(tag = "c",
       title = "Attribution of predictable variance (Shapley)",
       subtitle = sprintf("Shares of full-model R² = %.3f (n = %s); climate ≈ 1/3, non-climate ≈ 2/3",
                          full_R2, format(n_common, big.mark = ",")),
       x = NULL, y = "Share of attributable variance (%)") +
  theme_bw(base_size = 11) +
  theme(plot.tag = element_text(face = "bold"))

# =============================================================================
# ASSEMBLE: (a | b) over (c)
# =============================================================================

fig <- (pa | pb) / pc +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Spatially-honest variance decomposition of soil-carbon turnover",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(PLOT_DIR, "13c_decomposition.png"), fig,
       width = 11, height = 8.5, dpi = 200)
cat("OK  ", file.path(PLOT_DIR, "13c_decomposition.png"), "\n")

# =============================================================================
# SENSITIVITY: MAIN (6-layer) vs all-9 Shapley shares
# =============================================================================

if (have_sens) {
  cmp <- bind_rows(
    dec$shapley_df      %>% mutate(set = sprintf("MAIN: fungal+SPUN (n=%s)",
                                   format(dec$config$n_common, big.mark = ","))),
    dec_sens$shapley_df %>% mutate(set = sprintf("Sensitivity: all 9 (n=%s)",
                                   format(dec_sens$config$n_common, big.mark = ",")))
  ) %>% mutate(nice = NICE[group])

  ps <- ggplot(cmp, aes(reorder(nice, shapley_pct), shapley_pct, fill = nice,
                        alpha = set)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", shapley_pct)),
              position = position_dodge(width = 0.75), hjust = -0.12, size = 3,
              show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_alpha_manual(values = c(1, 0.5), name = NULL) +
    scale_y_continuous(limits = c(0, max(cmp$shapley_pct) * 1.25),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = "Shapley shares are insensitive to the biological-set choice",
         subtitle = sprintf("Full R²: MAIN %.3f vs all-9 %.3f — same two-tier story",
                            dec$full_R2, dec_sens$full_R2),
         x = NULL, y = "Share of attributable variance (%)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(PLOT_DIR, "13c_decomposition_sensitivity.png"), ps,
         width = 8.5, height = 5, dpi = 200)
  cat("OK  ", file.path(PLOT_DIR, "13c_decomposition_sensitivity.png"), "\n")
}

# =============================================================================
# CONSOLE SUMMARY (for discussion)
# =============================================================================

cat("\n--- HEADLINE (MAIN, fungal+SPUN, n =", format(n_common, big.mark = ","),
    ") ---\n")
cat("Full-model block-CV R^2:", round(full_R2, 4), "\n\n")
print(shap_df %>% transmute(group = nice,
                            shapley_pct, shapley_R2 = round(shapley_R2, 3),
                            single_R2 = round(single_R2, 3),
                            unique_R2 = round(unique_R2, 3)), row.names = FALSE)
cat(sprintf("\nSum of singles = %.3f vs full = %.3f -> overlap = %.3f\n",
            singles_sum, full_R2, dec$overlap))
cat("Climate share:", round(shap_df$shapley_pct[shap_df$group == "CLIMATE"], 1),
    "%  -> non-climate =",
    round(100 - shap_df$shapley_pct[shap_df$group == "CLIMATE"], 1), "%\n")
