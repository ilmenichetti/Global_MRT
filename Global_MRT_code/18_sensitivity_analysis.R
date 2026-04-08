################################################################################
# 18_sensitivity_analysis.R
#
# Accumulated Local Effects (ALE) analysis for M7 predictor variables.
#
# For each continuous predictor, ALE computes the marginal effect on MRT
# by evaluating prediction differences within narrow quantile-based bins.
# Within each bin, only pixels whose predictor value falls in that bin are
# used, and the variable is nudged to the bin's upper and lower edge.
# The local prediction difference is accumulated across bins to produce
# a response curve in native variable units.
#
# This corrects for predictor correlations: because each bin only uses
# pixels where the variable is naturally at that level, all other variables
# take realistic co-occurring values rather than arbitrary held-constant ones.
#
# Y-axis: accumulated local effect on MRT (years), centered at zero at the
# left edge of the distribution. Interpreted as: "relative to the global
# baseline, how much does MRT change as this variable increases?"
#
# Categorical predictors (koppen, soil class) are excluded.
#
# Input:  ./Global_MRT_code/outputs/13_rf_models.rds
#         ./Global_MRT_code/spatialized_layers/*
# Output: ./Global_MRT_code/outputs/18_ale_results.rds
#         ./Global_MRT_code/plots/step_18_sensitivity/*.png
#
# Author: Lorenzo
# Date:   2026-04-08
################################################################################

library(terra)
library(ranger)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR  <- "./Global_MRT_code"
OUTPUT_DIR    <- file.path(PIPELINE_DIR, "outputs")
RASTER_DIR    <- file.path(PIPELINE_DIR, "spatialized_layers")
PLOT_DIR      <- file.path(PIPELINE_DIR, "plots", "step_18_sensitivity")
CACHE_FILE    <- file.path(OUTPUT_DIR, "18_ale_results.rds")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Spatially stratified subsample size (5 degree lat/lon strata)
N_SAMPLE      <- 100000

# Number of quantile-based ALE bins per variable
# Each bin will contain approximately N_SAMPLE / N_BINS pixels
N_BINS        <- 40

# Variables to skip (categorical / integer-coded)
SKIP_VARS     <- c("koppen_value", "soilclass_func_code")

# Latitude bands for stratified summaries
LAT_BREAKS    <- c(-90, -23.5, 23.5, 66.5, 90)
LAT_LABELS    <- c("Southern Extratropics", "Tropics",
                   "Northern Temperate",    "Boreal/Arctic")

cat("═══════════════════════════════════════════════════════════════\n")
cat("  SOC MRT ALE ANALYSIS (Step 18)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD MODEL (M7 only)
# =============================================================================

cat("Loading M7 model...\n")
models <- readRDS(file.path(OUTPUT_DIR, "13_rf_models.rds"))

if (!"M7_full" %in% names(models)) {
  stop("M7_full not found in 13_rf_models.rds. Available: ",
       paste(names(models), collapse = ", "))
}

m7       <- models[["M7_full"]]$model
m7_vars  <- models[["M7_full"]]$predictors
cat("  M7 predictors:", length(m7_vars), "\n\n")

cont_vars <- setdiff(m7_vars, SKIP_VARS)
cat_vars  <- intersect(m7_vars, SKIP_VARS)
cat("  Continuous predictors to analyse:", length(cont_vars), "\n")
cat("  Categorical predictors (skipped): ", length(cat_vars),  "\n\n")

# =============================================================================
# RASTER MAPPING  (mirrored from step 14)
# =============================================================================

raster_mapping <- list(
  # Climate
  temperature_seasonality  = "climate/temperature_seasonality_0.1deg.tif",
  soil_temperature_0_20cm  = "climate/soil_temperature_0_20cm_0.1deg.tif",
  soil_moisture_0_20cm     = "climate/soil_moisture_0_20cm_0.1deg.tif",
  snow_cover_mean          = "climate/snow_cover_mean_0.1deg.tif",
  potential_evaporation    = "climate/potential_evaporation_0.1deg.tif",
  koppen_value             = "climate/koppen_detailed_0.1deg.tif",
  # Edaphic
  sg_clay                  = "soilgrids/clay_0_20cm_global.tif",
  sg_sand                  = "soilgrids/sand_0_20cm_global.tif",
  sg_bdod                  = "soilgrids/bdod_0_20cm_global.tif",
  sg_cfvo                  = "soilgrids/cfvo_0_20cm_global.tif",
  sg_phh2o                 = "soilgrids/phh2o_0_20cm_global.tif",
  sg_cec                   = "soilgrids/cec_0_20cm_global.tif",
  sg_nitrogen              = "soilgrids/nitrogen_0_20cm_global.tif",
  soilclass_func_code      = "soilclass/functional_group_0p1deg.tif",
  terrain_elev_mean        = "topography/elev_mean_0p1deg.tif",
  terrain_slope_mean       = "topography/slope_mean_0p1deg.tif",
  terrain_northness        = "topography/northness_0p1deg.tif",
  terrain_eastness         = "topography/eastness_0p1deg.tif",
  terrain_ruggedness       = "topography/elev_ruggedness_0p1deg.tif",
  # LandUse
  lc_trees                 = "landcover/landcover_fractions_0p1deg.tif",
  lc_grassland             = "landcover/landcover_fractions_0p1deg.tif",
  lc_shrubs                = "landcover/landcover_fractions_0p1deg.tif",
  lc_cropland              = "landcover/landcover_fractions_0p1deg.tif",
  lc_bare                  = "landcover/landcover_fractions_0p1deg.tif",
  lc_wetland               = "landcover/landcover_fractions_0p1deg.tif",
  hansen_any_loss          = "disturbances/Hansen_forest_loss/Hansen_lossyear_global.tif",
  burn_count_before_obs    = "disturbances/MODIS_burned/Burned_2020_global.tif",
  # Biological
  fungal_proportion        = "microbial/fungal_proportion_0.1deg.tif",
  AM_roots_colonized       = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  EcM_roots_colonized      = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  EcM_AM_root_ratio        = "microbial/mycorrhiza_barcelo_0.1deg.tif",
  AM_richness              = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_richness             = "microbial/mycorrhiza_spun_0.1deg.tif",
  AM_endemism              = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_endemism             = "microbial/mycorrhiza_spun_0.1deg.tif",
  EcM_AM_richness_ratio    = "microbial/mycorrhiza_spun_0.1deg.tif"
)

# =============================================================================
# LOAD RASTER STACK
# =============================================================================

cat("Building M7 predictor raster stack...\n")

template <- rast(
  xmin = -180, xmax = 180,
  ymin = -90,  ymax = 90,
  res  = 0.1,
  crs  = "EPSG:4326"
)

loaded_layers <- list()

for (var in m7_vars) {
  
  if (!var %in% names(raster_mapping)) {
    cat("  ✗", var, "(no mapping defined)\n"); next
  }
  
  fpath <- file.path(RASTER_DIR, raster_mapping[[var]])
  if (!file.exists(fpath)) {
    cat("  ✗", var, "(file not found)\n"); next
  }
  
  tryCatch({
    r <- rast(fpath)
    
    if (nlyr(r) > 1) {
      idx <- grep(var, names(r), ignore.case = TRUE)
      if (length(idx) > 0) {
        r <- r[[idx[1]]]
      } else {
        cat("  !", var, "— no layer name match, using layer 1\n")
        cat("    Available:", paste(names(r), collapse = ", "), "\n")
        r <- r[[1]]
      }
    }
    
    if (ncell(r) == 0 || nlyr(r) == 0) {
      cat("  ✗", var, "(empty after selection)\n"); return()
    }
    
    if (!compareGeom(r, template, stopOnError = FALSE)) {
      r <- resample(r, template, method = "bilinear")
    }
    
    names(r) <- var
    loaded_layers[[var]] <- r
    cat("  ✓", var, "\n")
    
  }, error = function(e) {
    cat("  ✗", var, "— ERROR:", conditionMessage(e), "\n")
  })
}

pred_stack <- rast(loaded_layers)
cat("\n  Stack ready:", nlyr(pred_stack), "layers,",
    format(ncell(pred_stack), big.mark = ","), "cells\n\n")

# =============================================================================
# SPATIALLY STRATIFIED SUBSAMPLE
# =============================================================================

cat("Drawing spatially stratified subsample (n =",
    format(N_SAMPLE, big.mark = ","), ")...\n")

xy       <- crds(pred_stack[[1]], na.rm = FALSE)
vals_all <- values(pred_stack)

complete_idx <- which(complete.cases(vals_all))
cat("  Complete pixels:", format(length(complete_idx), big.mark = ","), "\n")

lat_strat <- cut(xy[complete_idx, 2], breaks = seq(-90,  90,  by = 5))
lon_strat <- cut(xy[complete_idx, 1], breaks = seq(-180, 180, by = 5))
strat_id  <- paste(lat_strat, lon_strat)

set.seed(42)
n_strata    <- length(unique(strat_id))
per_stratum <- max(1, floor(N_SAMPLE / n_strata))

sampled_local <- unlist(tapply(
  seq_along(complete_idx),
  strat_id,
  function(i) sample(i, min(length(i), per_stratum))
))

sampled_global <- complete_idx[sampled_local]
cat("  Sampled pixels:", format(length(sampled_global), big.mark = ","), "\n\n")

df_base  <- as.data.frame(vals_all[sampled_global, , drop = FALSE])
lat_sub  <- xy[sampled_global, 2]
lat_band <- cut(lat_sub,
                breaks         = LAT_BREAKS,
                labels         = LAT_LABELS,
                include.lowest = TRUE)

df_base_m7 <- df_base[, m7_vars, drop = FALSE]

# =============================================================================
# ALE COMPUTATION FUNCTION
# =============================================================================
#
# For each variable:
#   1. Define N_BINS quantile-based bin edges from the variable's distribution
#   2. For each bin, identify pixels whose variable value falls within it
#   3. Predict MRT with variable set to upper bin edge, then lower bin edge
#   4. The bin's local effect = mean(pred_upper - pred_lower) for those pixels
#   5. Accumulate effects across bins (cumulative sum)
#   6. Centre the curve so the weighted mean effect = 0
#
# Returns a data frame with columns: bin_mid, ale_effect, n_pixels
# ale_effect is in MRT years (back-transformed from log scale)
# =============================================================================

compute_ale <- function(var, df, model, n_bins, lat_band_vec) {
  
  x     <- df[[var]]
  probs <- seq(0, 1, length.out = n_bins + 1)
  edges <- unique(quantile(x, probs = probs, na.rm = TRUE))
  
  if (length(edges) < 3) {
    cat("    ! Too few unique values for ALE binning, skipping\n")
    return(NULL)
  }
  
  n_bins_actual <- length(edges) - 1
  bin_results   <- vector("list", n_bins_actual)
  
  for (k in seq_len(n_bins_actual)) {
    
    lo  <- edges[k]
    hi  <- edges[k + 1]
    mid <- (lo + hi) / 2
    
    # Pixels naturally in this bin
    idx <- which(x >= lo & x < hi)
    if (k == n_bins_actual) idx <- which(x >= lo & x <= hi)  # include right edge on last bin
    
    if (length(idx) < 2) {
      bin_results[[k]] <- data.frame(
        bin_mid    = mid,
        bin_lo     = lo,
        bin_hi     = hi,
        local_effect = 0,
        n_pixels   = length(idx),
        lat_band   = "Global"
      )
      next
    }
    
    df_lo <- df[idx, , drop = FALSE]
    df_hi <- df[idx, , drop = FALSE]
    df_lo[[var]] <- lo
    df_hi[[var]] <- hi
    
    pred_lo <- exp(predict(model, data = df_lo)$predictions)
    pred_hi <- exp(predict(model, data = df_hi)$predictions)
    diff    <- pred_hi - pred_lo
    
    # Global bin row
    global_row <- data.frame(
      bin_mid      = mid,
      bin_lo       = lo,
      bin_hi       = hi,
      local_effect = mean(diff, na.rm = TRUE),
      n_pixels     = length(idx),
      lat_band     = "Global"
    )
    
    # Per latitude-band rows
    lb_rows <- lapply(levels(lat_band_vec), function(lb) {
      lb_idx <- which(lat_band_vec[idx] == lb)
      if (length(lb_idx) < 2) {
        data.frame(bin_mid = mid, bin_lo = lo, bin_hi = hi,
                   local_effect = NA_real_,
                   n_pixels = length(lb_idx), lat_band = lb)
      } else {
        data.frame(bin_mid = mid, bin_lo = lo, bin_hi = hi,
                   local_effect = mean(diff[lb_idx], na.rm = TRUE),
                   n_pixels = length(lb_idx), lat_band = lb)
      }
    })
    
    bin_results[[k]] <- rbind(global_row, do.call(rbind, lb_rows))
  }
  
  bin_df <- do.call(rbind, bin_results)
  
  # Accumulate and centre separately for each lat_band level
  bin_df <- bin_df %>%
    group_by(lat_band) %>%
    arrange(bin_mid, .by_group = TRUE) %>%
    mutate(
      ale_raw    = cumsum(replace(local_effect, is.na(local_effect), 0)),
      # Centre: subtract weighted mean (weighted by n_pixels in global band)
      ale_effect = ale_raw - weighted.mean(ale_raw, w = pmax(n_pixels, 1),
                                           na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(variable = var)
  
  return(bin_df)
}

# =============================================================================
# RUN ALE FOR ALL CONTINUOUS VARIABLES
# =============================================================================

cat("Computing ALE curves...\n")
cat("  Variables:", length(cont_vars), "  Bins:", N_BINS, "\n\n")

ale_list <- list()

for (var in cont_vars) {
  
  cat("  ALE:", var, "...")
  
  result <- tryCatch(
    compute_ale(var, df_base_m7, m7, N_BINS, lat_band),
    error = function(e) {
      cat(" ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(result)) {
    ale_list[[var]] <- result
    n_bins_used <- length(unique(result$bin_mid))
    cat(" ✓ (", n_bins_used, "bins )\n")
  }
}

cat("\n")

ale_df <- do.call(rbind, ale_list)
rownames(ale_df) <- NULL

# =============================================================================
# CACHE
# =============================================================================

cat("Saving results...\n")
saveRDS(ale_df, CACHE_FILE)
cat("  ✓ Saved:", CACHE_FILE, "\n\n")

# =============================================================================
# VARIABLE GROUP METADATA
# =============================================================================

group_map <- c(
  temperature_seasonality = "Climate",  soil_temperature_0_20cm = "Climate",
  soil_moisture_0_20cm    = "Climate",  snow_cover_mean         = "Climate",
  potential_evaporation   = "Climate",
  sg_clay                 = "Edaphic",  sg_sand                 = "Edaphic",
  sg_bdod                 = "Edaphic",  sg_cfvo                 = "Edaphic",
  sg_phh2o                = "Edaphic",  sg_cec                  = "Edaphic",
  sg_nitrogen             = "Edaphic",  terrain_elev_mean       = "Edaphic",
  terrain_slope_mean      = "Edaphic",  terrain_northness       = "Edaphic",
  terrain_eastness        = "Edaphic",  terrain_ruggedness      = "Edaphic",
  lc_trees                = "LandUse",  lc_grassland            = "LandUse",
  lc_shrubs               = "LandUse",  lc_cropland             = "LandUse",
  lc_bare                 = "LandUse",  lc_wetland              = "LandUse",
  hansen_any_loss         = "LandUse",  burn_count_before_obs   = "LandUse",
  fungal_proportion       = "Biological", AM_roots_colonized    = "Biological",
  EcM_roots_colonized     = "Biological", EcM_AM_root_ratio     = "Biological",
  AM_richness             = "Biological", EcM_richness          = "Biological",
  AM_endemism             = "Biological", EcM_endemism          = "Biological",
  EcM_AM_richness_ratio   = "Biological"
)

group_colours <- c(
  Climate    = "#2166AC",
  Edaphic    = "#B2182B",
  LandUse    = "#4DAC26",
  Biological = "#E08214"
)

# =============================================================================
# PLOTTING
# =============================================================================

cat("Generating plots...\n")

# Global-only data with group annotation and scaled x for spaghetti plots
plot_df <- ale_df %>%
  filter(lat_band == "Global") %>%
  mutate(
    group = group_map[variable],
    group = factor(group, levels = names(group_colours))
  ) %>%
  filter(!is.na(group)) %>%
  group_by(variable) %>%
  mutate(x_scaled = (bin_mid - min(bin_mid)) / (max(bin_mid) - min(bin_mid))) %>%
  ungroup()

# ── PLOT A: ALE effect in MRT years, one panel per variable ───────────────────

p_ale <- ggplot(plot_df,
                aes(x      = bin_mid,
                    y      = ale_effect,
                    colour = group,
                    fill   = group,
                    group  = variable)) +
  geom_hline(yintercept = 0, colour = "grey60",
             linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = ale_effect - abs(ale_effect) * 0.1,
                  ymax = ale_effect + abs(ale_effect) * 0.1),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.75, alpha = 0.85) +
  scale_colour_manual(values = group_colours, name = "Variable group") +
  scale_fill_manual(values   = group_colours, name = "Variable group") +
  facet_wrap(~ variable, scales = "free", ncol = 6) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    strip.text       = element_text(face = "bold", size = 8),
    legend.position  = "bottom",
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 7),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  ) +
  labs(
    title    = "ALE Response Curves — Effect on Global Mean MRT (M7)",
    subtitle = sprintf(
      "Accumulated Local Effects · %d quantile bins · n = %s pixels · y centred at weighted mean · dashed = zero",
      N_BINS, format(length(sampled_global), big.mark = ",")),
    x = "Variable value (native units)",
    y = "ALE effect on MRT (years)"
  )

ggsave(file.path(PLOT_DIR, "18A_ALE_curves_MRT.png"),
       p_ale, width = 20, height = 18, dpi = 300, bg = "white")
cat("  ✓ Saved: 18A_ALE_curves_MRT.png\n")

# ── PLOT B: Spaghetti — all variables overlaid, x as percentile rank ──────────
# Rescale x to [0,1] quantile rank within each variable for comparability

p_spaghetti <- ggplot(plot_df,
                      aes(x      = x_scaled,
                          y      = ale_effect,
                          colour = group,
                          group  = variable)) +
  geom_hline(yintercept = 0, colour = "grey60",
             linewidth = 0.4, linetype = "dashed") +
  geom_line(linewidth = 0.6, alpha = 0.6) +
  scale_colour_manual(values = group_colours, name = "Variable group") +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("p5", "p50", "p95")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    legend.position  = "right",
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.title       = element_text(size = 10)
  ) +
  labs(
    title    = "ALE Response Curves — All Variables Overlaid",
    subtitle = "x rescaled to p5-p95 range for cross-variable comparison · each line = one predictor",
    x        = "Variable percentile range",
    y        = "ALE effect on MRT (years)"
  )

ggsave(file.path(PLOT_DIR, "18B_ALE_spaghetti.png"),
       p_spaghetti, width = 12, height = 7, dpi = 300, bg = "white")
cat("  ✓ Saved: 18B_ALE_spaghetti.png\n")

# ── PLOT C: Latitude-band breakdown, top 10 most responsive variables ──────────

ale_range <- ale_df %>%
  filter(lat_band == "Global") %>%
  group_by(variable) %>%
  summarise(ale_range = diff(range(ale_effect, na.rm = TRUE)),
            .groups = "drop") %>%
  arrange(desc(ale_range))

top_vars <- ale_range$variable[1:min(10, nrow(ale_range))]

cat("\n  Top 10 variables by ALE range (years):\n")
print(ale_range[1:min(10, nrow(ale_range)), ])
cat("\n")

band_lty  <- c("Global"               = "solid",
               "Southern Extratropics" = "dashed",
               "Tropics"               = "dotdash",
               "Northern Temperate"    = "dotted",
               "Boreal/Arctic"         = "longdash")
band_alpha <- c("Global"               = 1.00,
                "Southern Extratropics" = 0.65,
                "Tropics"               = 0.65,
                "Northern Temperate"    = 0.65,
                "Boreal/Arctic"         = 0.65)

lat_df <- ale_df %>%
  filter(variable %in% top_vars) %>%
  mutate(
    group    = group_map[variable],
    group    = factor(group, levels = names(group_colours)),
    lat_band = factor(lat_band, levels = c("Global", LAT_LABELS))
  ) %>%
  group_by(variable, lat_band) %>%
  mutate(x_scaled = (bin_mid - min(bin_mid)) / (max(bin_mid) - min(bin_mid))) %>%
  ungroup()

p_latbands <- ggplot(lat_df,
                     aes(x        = x_scaled,
                         y        = ale_effect,
                         colour   = group,
                         linetype = lat_band,
                         alpha    = lat_band,
                         group    = interaction(variable, lat_band))) +
  geom_hline(yintercept = 0, colour = "grey60",
             linewidth = 0.4, linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = group_colours, name = "Group") +
  scale_linetype_manual(values = band_lty,   name = "Latitude band") +
  scale_alpha_manual(values = band_alpha,    name = "Latitude band") +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("p5", "p50", "p95")) +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    strip.text       = element_text(size = 8, face = "bold"),
    legend.position  = "bottom",
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9,  colour = "grey40"),
    axis.text        = element_text(size = 7),
    axis.title       = element_text(size = 9)
  ) +
  guides(colour   = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1),
         alpha    = guide_legend(nrow = 1)) +
  labs(
    title    = "ALE by Latitude Band — Top 10 Most Responsive Predictors",
    subtitle = "Solid = Global; dashed lines = latitude band means",
    x        = "Variable percentile range",
    y        = "ALE effect on MRT (years)"
  )

ggsave(file.path(PLOT_DIR, "18C_ALE_by_latband.png"),
       p_latbands, width = 18, height = 9, dpi = 300, bg = "white")
cat("  ✓ Saved: 18C_ALE_by_latband.png\n\n")

# ── PLOT D: Top 10 — ALE effect in MRT years, native x-axis ──────────────────

plot_df_top <- plot_df %>%
  filter(variable %in% top_vars) %>%
  mutate(variable = factor(variable, levels = top_vars))

p50_df_top <- plot_df_top %>%
  distinct(variable, group) %>%
  left_join(
    plot_df_top %>%
      group_by(variable) %>%
      slice(which.min(abs(x_scaled - 0.5))) %>%
      select(variable, bin_mid),
    by = "variable"
  ) %>%
  rename(x_ref = bin_mid)

p_ale_top <- ggplot(plot_df_top,
                    aes(x      = bin_mid,
                        y      = ale_effect,
                        colour = group,
                        fill   = group,
                        group  = variable)) +
  geom_hline(yintercept = 0, colour = "grey60",
             linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = ale_effect - abs(ale_effect) * 0.1,
                  ymax = ale_effect + abs(ale_effect) * 0.1),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  scale_colour_manual(values = group_colours, name = "Variable group") +
  scale_fill_manual(values   = group_colours, name = "Variable group") +
  facet_wrap(~ variable, scales = "free", ncol = 5) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 8),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  ) +
  labs(
    title    = "ALE Response Curves — Top 10 Variables, Global Mean MRT (M7)",
    subtitle = sprintf(
      "Accumulated Local Effects · %d quantile bins · n = %s pixels · dashed = zero",
      N_BINS, format(length(sampled_global), big.mark = ",")),
    x = "Variable value (native units)",
    y = "ALE effect on MRT (years)"
  )

ggsave(file.path(PLOT_DIR, "18D_ALE_top10_MRT.png"),
       p_ale_top, width = 16, height = 8, dpi = 300, bg = "white")
cat("  ✓ Saved: 18D_ALE_top10_MRT.png\n")

# ── PLOT E: Top 10 — spaghetti with native x-axis (one colour per variable) ───

# Assign a distinct colour per variable within its group
top_vars_df <- data.frame(variable = top_vars) %>%
  mutate(group = group_map[variable],
         group = factor(group, levels = names(group_colours)))

plot_df_top_scaled <- plot_df_top %>%
  group_by(variable) %>%
  mutate(x_scaled = (bin_mid - min(bin_mid)) / (max(bin_mid) - min(bin_mid))) %>%
  ungroup()

p_top_spaghetti <- ggplot(plot_df_top_scaled,
                          aes(x      = x_scaled,
                              y      = ale_effect,
                              colour = group,
                              group  = variable,
                              label  = variable)) +
  geom_hline(yintercept = 0, colour = "grey60",
             linewidth = 0.4, linetype = "dashed") +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  scale_colour_manual(values = group_colours, name = "Variable group") +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("p5", "p50", "p95")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    legend.position  = "right",
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.title       = element_text(size = 10)
  ) +
  labs(
    title    = "ALE Response Curves — Top 10 Variables Overlaid",
    subtitle = "x rescaled to p5-p95 range · colour = variable group",
    x        = "Variable percentile range",
    y        = "ALE effect on MRT (years)"
  )

ggsave(file.path(PLOT_DIR, "18E_ALE_top10_spaghetti.png"),
       p_top_spaghetti, width = 12, height = 7, dpi = 300, bg = "white")
cat("  ✓ Saved: 18E_ALE_top10_spaghetti.png\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  ALE RANKING (by total effect range in MRT years)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
print(ale_range, n = nrow(ale_range))

cat("\n✓ Step 18 complete.\n")
cat("  Results cached : ", CACHE_FILE, "\n")
cat("  Plots saved to : ", PLOT_DIR,   "\n\n")
