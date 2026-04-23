# =============================================================================
# Step 21: RF M7 tau vs Incubation-derived tau (Zhang et al. 2025)
# =============================================================================
#
# PURPOSE:
#   Compares RF M7 global MRT predictions against turnover times derived from
#   laboratory soil incubation experiments (Zhang et al. 2025, ESSD).
#   The incubation-derived tau provides a methodologically independent
#   estimate of decomposition kinetics, based on direct measurement of CO2
#   flux rather than equilibrium mass balance.
#
# APPROACH:
#   1. Filter incubation data to surface mineral soils (0-30 cm), no C inputs
#   2. Fit first-order decay model per experiment: Cumu_SOC ~ A*(1-exp(-k*day))
#   3. Temperature-normalise k to field MAT using Q10 = 2.0
#   4. Convert k_field (day-1) to tau_incubation (years)
#   5. Extract RF M7 tau at each site coordinate
#   6. Scatterplot: RF tau vs tau_incubation, coloured by ecosystem type
#
# KEY CAVEAT:
#   Incubation tau (decomposability under controlled conditions) is not
#   equivalent to field equilibrium tau. Incubations remove substrate inputs
#   and standardise temperature/moisture. The comparison therefore tests
#   spatial congruence in decomposition kinetics, not absolute equivalence.
#
# SOURCE:
#   Zhang S., Wang M., Zheng J., Luo Z. (2025). A global dataset of soil
#   organic carbon mineralization under various incubation conditions. ESSD.
#   https://doi.org/10.6084/m9.figshare.25808698
#
# OUTPUTS:
#   outputs/incubation/incubation_tau_sites.csv   — per-site tau estimates
#   plots/step_21/fig_scatterplot_rf_vs_incu.pdf  — main scatterplot
#   plots/step_21/fig_scatterplot_rf_vs_incu.png
# =============================================================================

source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(minpack.lm)   # for robust nonlinear least squares (nlsLM)

# =============================================================================
# 0. Configuration
# =============================================================================

ZHANG_DIR  <- "./Global_MRT_code/Zhang_et_al_2025"
OUTPUT_DIR <- "./Global_MRT_code/outputs/incubation"
PLOT_DIR   <- "./Global_MRT_code/plots/step_21"
RF_M7_PATH <- "./Global_MRT_code/outputs/MRT_predictions/MRT_M7_full.tif"

# Q10 for temperature normalisation (standard value; sensitivity tested below)
Q10 <- 2.0

# Depth categories to retain (surface/topsoil only)
SURFACE_DEPTHS <- c("0_5", "0_6", "0_7", "0_7.1", "0_7.5", "0_8", "0_8.4",
                    "0_10", "0_11", "0_12", "0_13", "0_15", "0_20", "0_28",
                    "0_30", "Top", "0_1.5")

# tau limits for display (years) — very short or very long tau unreliable
TAU_MIN <- 0.5
TAU_MAX <- 500

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,   recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load and filter incubation data
# =============================================================================

cat("\n=== Loading Zhang et al. 2025 incubation data ===\n\n")

data_raw <- read.csv(file.path(ZHANG_DIR, "data.csv"),
                     stringsAsFactors = FALSE)
cat(sprintf("Raw rows: %d | Unique sites: %d\n",
            nrow(data_raw),
            nrow(unique(data_raw[, c("Latitude", "Longitude")]))))

data_filt <- data_raw %>%
  filter(
    Soil_depth %in% SURFACE_DEPTHS,         # surface soils only
    C_input    %in% c("NONE", "NONE "),     # no substrate additions
    !is.na(Latitude),                        # must have coordinates
    !is.na(Longitude),
    !is.na(Cumu_SOC),                        # must have cumulative CO2
    !is.na(Measure.day),
    Measure.day > 0,
    Cumu_SOC   >= 0,
    !is.na(Incu_temp),
    !is.na(MAT)                              # need field MAT for normalisation
  ) %>%
  mutate(
    Exp_ID_full = paste(Reference_ID, Soil_ID, Profile_ID,
                        Exp_ID, Incu_temp, sep = "_")
  )

cat(sprintf("After filtering: %d rows | %d experiments | %d unique sites\n",
            nrow(data_filt),
            length(unique(data_filt$Exp_ID_full)),
            nrow(unique(data_filt[, c("Latitude", "Longitude")]))))

# =============================================================================
# 2. Fit first-order decay model per experiment
# =============================================================================
#
# Model: Cumu_SOC(t) = A * (1 - exp(-k * t))
#   A = potentially mineralisable C (mg CO2-C g-1 SOC)
#   k = first-order rate constant (day-1)
#
# We need at least 3 time points to fit 2 parameters reliably.
# =============================================================================

cat("\n=== Fitting first-order decay models ===\n\n")

fit_decay <- function(days, cumu, exp_id) {
  
  # Need at least 3 points
  if (sum(!is.na(cumu) & !is.na(days)) < 3) return(NULL)
  
  # Starting values: A ~ max observed, k ~ 0.01 day-1
  A_start <- max(cumu, na.rm = TRUE) * 1.2
  k_start <- 0.01
  
  tryCatch({
    fit <- minpack.lm::nlsLM(
      cumu ~ A * (1 - exp(-k * days)),
      start   = list(A = A_start, k = k_start),
      lower   = c(A = 0, k = 1e-6),
      upper   = c(A = Inf, k = 10),
      control = nls.lm.control(maxiter = 100)
    )
    coefs <- coef(fit)
    data.frame(
      Exp_ID_full = exp_id,
      A_fitted    = coefs["A"],
      k_day       = coefs["k"],    # day-1
      n_points    = sum(!is.na(cumu)),
      r2          = 1 - sum(residuals(fit)^2) /
        sum((cumu - mean(cumu, na.rm = TRUE))^2)
    )
  }, error = function(e) NULL)
}

# Run fitting per experiment
exp_list <- split(data_filt, data_filt$Exp_ID_full)

cat(sprintf("Fitting %d experiments...\n", length(exp_list)))

fit_results <- lapply(seq_along(exp_list), function(i) {
  if (i %% 500 == 0) cat(sprintf("  %d / %d\n", i, length(exp_list)))
  exp_data <- exp_list[[i]]
  fit_decay(exp_data$Measure.day, exp_data$Cumu_SOC,
            unique(exp_data$Exp_ID_full))
})

fit_df <- bind_rows(fit_results)
cat(sprintf("\nSuccessful fits: %d / %d (%.1f%%)\n",
            nrow(fit_df), length(exp_list),
            100 * nrow(fit_df) / length(exp_list)))

# =============================================================================
# 3. Merge fits with site metadata and temperature-normalise
# =============================================================================

cat("\n=== Temperature normalising k to field MAT ===\n\n")

# Site-level metadata (one row per experiment)
site_meta <- data_filt %>%
  group_by(Exp_ID_full) %>%
  summarise(
    Latitude   = first(Latitude),
    Longitude  = first(Longitude),
    MAT        = first(MAT),
    MAP        = first(MAP),
    Incu_temp  = first(Incu_temp),
    Eco_type   = first(Eco_type),
    Soil_depth = first(Soil_depth),
    SOC        = first(SOC),
    pH         = first(pH),
    Clay       = first(Clay),
    .groups    = "drop"
  )

# Join and normalise
results <- fit_df %>%
  left_join(site_meta, by = "Exp_ID_full") %>%
  filter(
    r2 > 0.7,              # only keep well-fitted experiments
    k_day > 0
  ) %>%
  mutate(
    # Temperature normalisation: k_field = k_lab * Q10^((MAT - Incu_temp)/10)
    k_field_day = k_day * Q10^((MAT - Incu_temp) / 10),
    k_field_yr  = k_field_day * 365,
    tau_incu    = 1 / k_field_yr,    # years
    
    # Q10 sensitivity: also compute with Q10 = 1.5 and 2.5
    tau_q10_15  = 1 / (k_day * 1.5^((MAT - Incu_temp) / 10) * 365),
    tau_q10_25  = 1 / (k_day * 2.5^((MAT - Incu_temp) / 10) * 365)
  ) %>%
  filter(
    tau_incu >= TAU_MIN,
    tau_incu <= TAU_MAX,
    !is.na(tau_incu)
  )

cat(sprintf("Sites after quality filtering: %d\n", nrow(results)))
cat(sprintf("tau_incu range: %.1f - %.1f yr (median: %.1f yr)\n",
            min(results$tau_incu),
            max(results$tau_incu),
            median(results$tau_incu)))

# =============================================================================
# 4. Extract RF M7 tau at incubation site coordinates
# =============================================================================

cat("\n=== Extracting RF M7 tau at incubation sites ===\n\n")

rf_m7 <- rast(RF_M7_PATH)

site_vect <- vect(results[, c("Longitude", "Latitude")],
                  geom = c("Longitude", "Latitude"),
                  crs  = "EPSG:4326")

rf_extracted <- terra::extract(rf_m7, site_vect)
results$tau_rf <- rf_extracted[, 2]

# Remove sites where RF has no prediction (ocean, ice, etc.)
results_final <- results %>%
  filter(!is.na(tau_rf), tau_rf > 0)

cat(sprintf("Sites with RF prediction: %d\n", nrow(results_final)))

# Save site-level results
write.csv(results_final,
          file.path(OUTPUT_DIR, "incubation_tau_sites.csv"),
          row.names = FALSE)

# =============================================================================
# 5. Scatterplot: RF M7 tau vs tau_incubation
# =============================================================================

cat("\n=== Plot: RF M7 vs incubation tau scatterplot ===\n\n")

# Correlation statistics
r_spearman <- cor(log10(results_final$tau_rf),
                  log10(results_final$tau_incu),
                  method = "spearman",
                  use    = "complete.obs")

r_pearson  <- cor(log10(results_final$tau_rf),
                  log10(results_final$tau_incu),
                  method = "pearson",
                  use    = "complete.obs")

cat(sprintf("Pearson r (log-log):   %.3f\n", r_pearson))
cat(sprintf("Spearman r (log-log):  %.3f\n", r_spearman))

# Ecosystem colour palette
eco_colors <- c(
  "Forest"    = "#2d6a4f",
  "Grassland" = "#95d5b2",
  "Cropland"  = "#f4a261",
  "Wetland"   = "#457b9d",
  "Shrubland" = "#a8dadc",
  "Tundra"    = "#8ecae6",
  "Desert"    = "#e9c46a"
)

# Clean up Eco_type labels
results_final <- results_final %>%
  mutate(Eco_label = case_when(
    grepl("forest|Forest|woodland|Woodland", Eco_type, ignore.case = TRUE) ~ "Forest",
    grepl("grass|Grass|meadow|Meadow|steppe|Steppe|savann",
          Eco_type, ignore.case = TRUE)                                     ~ "Grassland",
    grepl("crop|Crop|agri|Agri|farm|Farm|rice|Rice|wheat|maize",
          Eco_type, ignore.case = TRUE)                                     ~ "Cropland",
    grepl("wetland|Wetland|marsh|bog|peat|Peat",
          Eco_type, ignore.case = TRUE)                                     ~ "Wetland",
    grepl("shrub|Shrub|heath|Heath",
          Eco_type, ignore.case = TRUE)                                     ~ "Shrubland",
    grepl("tundra|Tundra|arctic|Arctic|alpine|Alpine",
          Eco_type, ignore.case = TRUE)                                     ~ "Tundra",
    grepl("desert|Desert|arid|Arid",
          Eco_type, ignore.case = TRUE)                                     ~ "Desert",
    TRUE                                                                    ~ "Other"
  ))

# 1:1 reference line limits
lim_min <- floor(log10(min(c(results_final$tau_rf,
                             results_final$tau_incu), na.rm = TRUE)))
lim_max <- ceiling(log10(max(c(results_final$tau_rf,
                               results_final$tau_incu), na.rm = TRUE)))

# Annotation label
ann_label <- sprintf("Pearson r = %.2f\nSpearman \u03c1 = %.2f\nn = %d",
                     r_pearson, r_spearman, nrow(results_final))

p_scatter <- ggplot(results_final,
                    aes(x = tau_rf, y = tau_incu, color = Eco_label)) +
  
  # 1:1 line
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey40", linewidth = 0.7) +
  
  # Q10 uncertainty ribbon (range from Q10=1.5 to Q10=2.5)
  geom_segment(aes(xend = tau_rf, y = tau_q10_15, yend = tau_q10_25),
               color = "grey80", linewidth = 0.3, alpha = 0.4) +
  
  # Points
  geom_point(alpha = 0.6, size = 1.8) +
  
  # Loess smoother
  geom_smooth(method  = "lm",
              formula = y ~ x,
              se      = TRUE,
              color   = "black",
              linewidth = 0.8,
              linetype  = "solid") +
  
  # Correlation annotation
  annotate("text",
           x     = 10^(lim_min + 0.2),
           y     = 10^(lim_max - 0.3),
           label = ann_label,
           hjust = 0, vjust = 1,
           size  = 3.2,
           color = "grey20") +
  
  scale_x_log10(
    limits = 10^c(lim_min, lim_max),
    breaks = 10^(lim_min:lim_max),
    labels = function(x) ifelse(x >= 1, as.character(round(x)), as.character(x))
  ) +
  scale_y_log10(
    limits = 10^c(lim_min, lim_max),
    breaks = 10^(lim_min:lim_max),
    labels = function(x) ifelse(x >= 1, as.character(round(x)), as.character(x))
  ) +
  scale_color_manual(
    values   = c(eco_colors, Other = "grey60"),
    name     = "Ecosystem",
    na.value = "grey60"
  ) +
  
  labs(
    x        = expression(paste("RF M7 predicted ", tau, " (years)")),
    y        = expression(paste("Incubation-derived ", tau, " (years)")),
    title    = "Comparison of RF M7 MRT predictions with incubation-derived decomposition times",
    subtitle = paste0("Incubation \u03c4 = 1/k, temperature-normalised to field MAT (Q10 = ",
                      Q10, "). Grey bars = Q10 uncertainty (1.5\u20132.5).\n",
                      "Dashed line = 1:1; solid line = linear regression on log-log axes.")
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.text       = element_text(size = 9),
    panel.grid.minor  = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 8, color = "grey40"),
    aspect.ratio      = 1
  )

ggsave(file.path(PLOT_DIR, "fig_scatterplot_rf_vs_incu.pdf"),
       p_scatter, width = 7, height = 7, dpi = 300)
ggsave(file.path(PLOT_DIR, "fig_scatterplot_rf_vs_incu.png"),
       p_scatter, width = 7, height = 7, dpi = 300)

cat("Saved: fig_scatterplot_rf_vs_incu\n")

# =============================================================================
# 6. Summary statistics by ecosystem type
# =============================================================================

summary_eco <- results_final %>%
  group_by(Eco_label) %>%
  summarise(
    n          = n(),
    median_rf  = round(median(tau_rf,   na.rm = TRUE), 1),
    median_inc = round(median(tau_incu, na.rm = TRUE), 1),
    r_spearman = round(cor(log10(tau_rf), log10(tau_incu),
                           method = "spearman",
                           use = "complete.obs"), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

cat("\nSummary by ecosystem:\n")
print(summary_eco)

write.csv(summary_eco,
          file.path(OUTPUT_DIR, "incubation_tau_summary_by_eco.csv"),
          row.names = FALSE)

# =============================================================================
# 7. Done
# =============================================================================

cat("\n\u2550\u2550\u2550 Step 21 complete \u2550\u2550\u2550\n")
cat(sprintf("Sites analysed : %d\n", nrow(results_final)))
cat(sprintf("Outputs        : %s\n", OUTPUT_DIR))
cat(sprintf("Plots          : %s\n", PLOT_DIR))