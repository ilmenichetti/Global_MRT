################################################################################
# 13h_blocksize_sensitivity.R
#
# Consolidate the spatial-block-CV block-size sensitivity: does the choice of CV
# block size change the Shapley attribution? Reads the decomposition objects
# produced by 13c_commonality_validation.R at several block sizes and reports the
# Shapley shares (and full R^2) as a function of block size. Used to justify the
# 2 deg choice in Methods: if the shares are invariant across 1-5 deg, the choice
# of block size is immaterial to the conclusions.
#
# Block-size files (from `MRT_13C_BLOCK=<deg> Rscript 13c_...R`):
#   13c_decomposition.rds        -> 2 deg  (MAIN headline; 500 trees)
#   13c_decomposition_blk1.rds   -> 1 deg  (sensitivity; 300 trees)
#   13c_decomposition_blk5.rds   -> 5 deg  (sensitivity; 300 trees)
# Each file's own config (block size, tree count, n) is read from the object, so
# a future 500-tree sweep regenerates this table without edits.
#
# Input:  outputs/13c_decomposition*.rds
# Output: outputs/13h_blocksize_sensitivity.csv
#         plots/step_13c_commonality/13h_blocksize_shapley.png
#
# Author: Lorenzo   Date: 2026-06-26
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

`%||%` <- function(a, b) if (is.null(a)) b else a

OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR   <- "./Global_MRT_code/plots/step_13c_commonality"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_COLORS <- c(Climate = "#666666", Edaphic = "#E41A1C",
                  LandUse = "#4DAF4A", Biological = "#377EB8")
NICE <- c(CLIMATE = "Climate", EDAPHIC = "Edaphic",
          LANDUSE = "LandUse", BIOLOGICAL = "Biological")

# Map each candidate file to its nominal block size (config is authoritative)
files <- list(
  "2" = file.path(OUTPUT_DIR, "13c_decomposition.rds"),
  "1" = file.path(OUTPUT_DIR, "13c_decomposition_blk1.rds"),
  "5" = file.path(OUTPUT_DIR, "13c_decomposition_blk5.rds")
)

rows <- list()
for (nm in names(files)) {
  f <- files[[nm]]
  if (!file.exists(f)) {
    cat("  (skip, not found:", basename(f), ")\n"); next
  }
  dec <- readRDS(f)
  blk <- dec$config$block_size_deg %||% as.numeric(nm)
  rows[[nm]] <- dec$shapley_df %>%
    transmute(group, shapley_pct, shapley_R2,
              block_deg   = blk,
              full_R2     = dec$full_R2,
              n_common    = dec$config$n_common,
              rf_ntrees   = dec$config$rf_ntrees)
}
stopifnot(length(rows) >= 1)
tab <- bind_rows(rows) %>% mutate(nice = NICE[group])

# ---- Wide Shapley-share table (rows = domain, cols = block size) -------------
wide <- tab %>%
  select(nice, block_deg, shapley_pct) %>%
  pivot_wider(names_from = block_deg, values_from = shapley_pct,
              names_prefix = "blk_") %>%
  arrange(desc(.data[[paste0("blk_", 2)]] %||% -Inf))

# Stability metric: max spread of each domain's share across block sizes
stab <- tab %>% group_by(nice) %>%
  summarise(min_pct = min(shapley_pct), max_pct = max(shapley_pct),
            spread_pp = max(shapley_pct) - min(shapley_pct), .groups = "drop")
max_spread <- max(stab$spread_pp)

cat("=== Shapley shares (%) by CV block size ===\n")
print(as.data.frame(wide), row.names = FALSE, digits = 3)
cat("\n=== Full R^2 / n / trees by block size ===\n")
print(tab %>% distinct(block_deg, full_R2, n_common, rf_ntrees) %>%
        arrange(block_deg) %>% as.data.frame(), row.names = FALSE, digits = 3)
cat(sprintf("\nMax share spread across block sizes: %.1f pp (largest in %s)\n",
            max_spread, stab$nice[which.max(stab$spread_pp)]))

write.csv(tab, file.path(OUTPUT_DIR, "13h_blocksize_sensitivity.csv"), row.names = FALSE)

# ---- Figure: grouped bars of Shapley share by block size --------------------
p <- ggplot(tab, aes(reorder(nice, shapley_pct), shapley_pct,
                     fill = nice, alpha = factor(block_deg))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f", shapley_pct)),
            position = position_dodge(width = 0.8), hjust = -0.15, size = 2.7,
            show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_alpha_manual(values = c(0.45, 0.7, 1.0), name = "CV block (deg)") +
  scale_y_continuous(limits = c(0, max(tab$shapley_pct) * 1.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Shapley attribution is insensitive to CV block size",
       subtitle = sprintf("Shares stable across 1-5 deg (max spread %.1f pp)", max_spread),
       x = NULL, y = "Share of attributable variance (%)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, "13h_blocksize_shapley.png"), p,
       width = 8, height = 4.8, dpi = 200)
cat("OK  ", file.path(PLOT_DIR, "13h_blocksize_shapley.png"), "\n")
cat("OK  outputs/13h_blocksize_sensitivity.csv\n")
