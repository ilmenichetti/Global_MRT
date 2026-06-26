################################################################################
# 20b_zonal_amplification.R
#
# Over-zonalization contrast: how much more steeply do ESMs amplify polar vs
# tropical soil-carbon turnover than the observation-constrained RF? This is the
# quantitative payoff of the (structural) ESM comparison and is scale-free
# (ratios within each product), so it is immune to the absolute depth/pool
# incomparability between products.
#
# Reads the already-computed step-20 outputs (no CMIP6 nc re-processing):
#   - cmip6_biome_relative_tau.csv : per-product mean tau by biome zone
#   - cmip6_tau_ensemble_mean.tif  : CMIP6 ensemble-mean tau
#   - MRT_M7_full.tif              : RF M7 tau
#
# Output: outputs/cmip6/cmip6_zonal_amplification.csv
#
# Author: Lorenzo   Date: 2026-06-26
################################################################################

suppressMessages({ library(dplyr); library(tidyr); library(terra) })

CMIP_DIR <- "./Global_MRT_code/outputs/cmip6"
PRED_DIR <- "./Global_MRT_code/outputs/MRT_predictions"

# --- (1) Biome-based polar/tropical amplification (per product) ---------------
biome <- read.csv(file.path(CMIP_DIR, "cmip6_biome_relative_tau.csv"))
w <- biome %>%
  select(source, zone, tau, tau_rel) %>%
  pivot_wider(names_from = zone, values_from = c(tau, tau_rel)) %>%
  mutate(
    polar_over_tropical = tau_Polar / tau_Tropical,   # scale-free zonal amplification
    polar_rel_to_mean   = tau_rel_Polar               # polar tau relative to product's own mean
  ) %>%
  mutate(kind = ifelse(source == "RF_M7", "observation (RF)", "ESM")) %>%
  select(source, kind, polar_over_tropical, polar_rel_to_mean,
         tau_Tropical, tau_Boreal, tau_Polar) %>%
  arrange(kind, desc(polar_over_tropical))

rf_pt  <- w$polar_over_tropical[w$source == "RF_M7"]
esm_pt <- w$polar_over_tropical[w$kind == "ESM"]

# --- (2) Raster-based high-latitude (>50 deg) / tropical ratio ----------------
band_ratio <- function(f) {
  r <- rast(f); lat <- init(r, "y")
  trop <- global(ifel(abs(lat) < 23.5, r, NA), "mean", na.rm = TRUE)[1, 1]
  high <- global(ifel(abs(lat) > 50,   r, NA), "mean", na.rm = TRUE)[1, 1]
  high / trop
}
rf_rast  <- band_ratio(file.path(PRED_DIR, "MRT_M7_full.tif"))
ens_rast <- band_ratio(file.path(CMIP_DIR, "cmip6_tau_ensemble_mean.tif"))

# --- (3) Summary --------------------------------------------------------------
summary_tbl <- data.frame(
  metric = c("RF polar/tropical (biome)",
             "ESM median polar/tropical (biome)",
             "ESM range polar/tropical (biome)",
             "over-zonalization factor (ESM median / RF)",
             "RF high-lat(>50)/tropical (raster)",
             "CMIP6 ensemble high-lat(>50)/tropical (raster)",
             "ensemble/RF (raster)"),
  value = c(sprintf("%.2f", rf_pt),
            sprintf("%.2f", median(esm_pt)),
            sprintf("%.2f - %.2f", min(esm_pt), max(esm_pt)),
            sprintf("%.1fx", median(esm_pt) / rf_pt),
            sprintf("%.2f", rf_rast),
            sprintf("%.2f", ens_rast),
            sprintf("%.1fx", ens_rast / rf_rast))
)

cat("=== Per-product zonal amplification ===\n")
print(w %>% mutate(across(where(is.numeric), ~round(., 2))), row.names = FALSE)
cat("\n=== Summary contrast ===\n")
print(summary_tbl, row.names = FALSE)

write.csv(w, file.path(CMIP_DIR, "cmip6_zonal_amplification.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(CMIP_DIR, "cmip6_zonal_amplification_summary.csv"), row.names = FALSE)
cat("\nOK  outputs/cmip6/cmip6_zonal_amplification{,_summary}.csv\n")
