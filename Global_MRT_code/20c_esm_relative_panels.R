################################################################################
# 20c_esm_relative_panels.R
#
# Two-panel relative-turnover figure (replaces the single biome-zone bar figure):
#   (a) relative tau by biome zone (discrete)         <- cmip6_biome_relative_tau.csv
#   (b) the SAME relative tau by continuous latitude  <- cmip6_zonal_means.csv
#       (each product normalised to its own global mean; profiles read more
#        naturally than the zonal bars and carry the same information).
# Reads existing step-20 outputs only (no CMIP6/RF recomputation).
# Output: plots/step_20/fig_esm_relative_panels.png
################################################################################
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

CMIP_DIR <- "./Global_MRT_code/outputs/cmip6"
PLOT_DIR <- "./Global_MRT_code/plots/step_20"

biome <- read.csv(file.path(CMIP_DIR, "cmip6_biome_relative_tau.csv"))
zm    <- read.csv(file.path(CMIP_DIR, "cmip6_zonal_means.csv"))

src   <- c("RF_M7", "CESM2", "UKESM1-0-LL", "IPSL-CM6A-LR",
           "MPI-ESM1-2-LR", "CanESM5", "MIROC-ES2L")
labs  <- c(RF_M7 = "RF (this study)", CESM2 = "CESM2",
           "UKESM1-0-LL" = "UKESM1-0-LL", "IPSL-CM6A-LR" = "IPSL-CM6A-LR",
           "MPI-ESM1-2-LR" = "MPI-ESM1-2-LR", CanESM5 = "CanESM5",
           "MIROC-ES2L" = "MIROC-ES2L")
pal   <- c("RF (this study)" = "#B2182B", "CESM2" = "#332288",
           "UKESM1-0-LL" = "#117733", "IPSL-CM6A-LR" = "#88CCEE",
           "MPI-ESM1-2-LR" = "#DDCC77", "CanESM5" = "#CC6677",
           "MIROC-ES2L" = "#AA4499")
plev  <- unname(labs[src])

# ---- panel (a): biome zones ----
ba <- biome %>% filter(source %in% src) %>%
  mutate(zone = factor(zone, levels = c("Tropical", "Temperate", "Boreal", "Polar")),
         prod = factor(labs[source], levels = plev),
         is_rf = source == "RF_M7")
pa <- ggplot(ba, aes(zone, tau_rel, colour = prod, group = prod)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
  geom_line(aes(linewidth = is_rf), alpha = 0.9) +
  geom_point(data = ~filter(.x, is_rf), size = 2.2, show.legend = FALSE) +
  scale_y_log10() +
  scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.4), guide = "none") +
  scale_colour_manual(values = pal, name = NULL) +
  labs(tag = "a", title = "By biome zone",
       x = NULL, y = "Relative turnover (product / own global mean)") +
  theme_bw(base_size = 10) +
  theme(plot.tag = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1))

# ---- panel (b): continuous latitude (same normalisation) ----
gm <- biome %>% distinct(source, global_mean)
zb <- zm %>% filter(source %in% src, !is.na(tau_mean)) %>%
  left_join(gm, by = "source") %>%
  mutate(tau_rel = tau_mean / global_mean,
         prod = factor(labs[source], levels = plev),
         is_rf = source == "RF_M7")
pb <- ggplot(zb, aes(lat, tau_rel, colour = prod, group = prod)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
  geom_line(aes(linewidth = is_rf), alpha = 0.9) +
  scale_y_log10() +
  scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.4), guide = "none") +
  scale_colour_manual(values = pal, name = NULL) +
  labs(tag = "b", title = "By latitude (continuous)",
       x = "Latitude (deg)", y = NULL) +
  theme_bw(base_size = 10) + theme(plot.tag = element_text(face = "bold"))

fig <- (pa | pb) + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Relative soil-carbon turnover: CMIP6 ESMs vs observation-constrained estimate",
    subtitle = "Each product normalised to its own global mean; RF (this study) in bold red",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 9, colour = "grey30")))
ggsave(file.path(PLOT_DIR, "fig_esm_relative_panels.png"), fig,
       width = 10, height = 4.8, dpi = 200)
cat("OK  ", file.path(PLOT_DIR, "fig_esm_relative_panels.png"), "\n")
