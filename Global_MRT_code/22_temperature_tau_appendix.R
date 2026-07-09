# =============================================================================
# Step 22: Temperature-tau scatter (appendix figure)
# =============================================================================
# Observed transit time (tau) against mean annual air temperature, the summary
# view commonly shown for this kind of study (cf. Varney et al. 2020, Fig. 3a).
# Points are the QC-passed profiles used for modelling, coloured by coarse Koppen
# zone; a single overall trend line summarises the temperature response, and we
# report the fitted slope and the implied apparent Q10 of bulk turnover.
#
# QC recipe and Koppen palette are kept identical to Step 13e for consistency.
#
# Input:  ./Global_MRT_code/outputs/12b_model_ready.rds
# Output: ./Global_MRT_code/plots/appendix/22_temperature_tau_scatter.png
#         ./Global_MRT_code/outputs/22_temperature_tau_fit.csv
# =============================================================================

library(dplyr)
library(ggplot2)

OUTPUT_DIR <- "./Global_MRT_code/outputs"
PLOT_DIR   <- "./Global_MRT_code/plots/appendix"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("===============================================================\n")
cat("  Step 22: TEMPERATURE-TAU SCATTER (appendix)\n")
cat("===============================================================\n\n")

# ---- Load and apply the Step 13e QC recipe ---------------------------------
d <- readRDS(file.path(OUTPUT_DIR, "12b_model_ready.rds"))
source("./Global_MRT_code/peat_filter.R")   # MAIN: drop peat (GPM 2.0); scoped to mineral soils
d <- apply_peat_filter(d, OUTPUT_DIR)
d <- d %>% filter(MRT_QC == "valid", !is.na(MRT_years), MRT_years > 0, MRT_years < Inf)
q <- quantile(d$MRT_years, c(.01, .99), na.rm = TRUE)
d <- d %>% filter(MRT_years > q[1] & MRT_years < q[2])

# temperature and zone; drop rows missing either
d <- d %>%
  filter(!is.na(temperature_2m_mean), !is.na(koppen_main_group)) %>%
  mutate(koppen = factor(koppen_main_group,
                         levels = c("A", "B", "C", "D", "E"),
                         labels = c("A Tropical", "B Arid", "C Temperate",
                                    "D Cold", "E Polar")))
cat("Profiles plotted:", nrow(d), "\n")

# ---- Overall temperature response (log tau vs MAT) -------------------------
fit <- lm(log(MRT_years) ~ temperature_2m_mean, data = d)
slope <- unname(coef(fit)[2])                 # d ln(tau) / dT  (per degC)
# turnover rate k = 1/tau, so d ln(k)/dT = -slope; apparent Q10 = exp(10 * dlnk/dT)
q10  <- exp(-10 * slope)
r2   <- summary(fit)$r.squared
cat(sprintf("ln(tau) ~ MAT slope = %.4f per degC  (R2 = %.3f)\n", slope, r2))
cat(sprintf("Implied apparent Q10 of bulk turnover = %.2f\n", q10))
write.csv(data.frame(slope_per_degC = slope, apparent_Q10 = q10, r2 = r2),
          file.path(OUTPUT_DIR, "22_temperature_tau_fit.csv"), row.names = FALSE)

# ---- Plot (Koppen palette identical to 13e) --------------------------------
pal <- c("A Tropical" = "#117733", "B Arid" = "#DDCC77", "C Temperate" = "#332288",
         "D Cold" = "#AA4499", "E Polar" = "#88CCEE")

p <- ggplot(d, aes(temperature_2m_mean, MRT_years)) +
  geom_point(aes(colour = koppen), alpha = 0.08, size = 0.4) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 1) +
  scale_y_log10() +
  scale_colour_manual(values = pal, name = "Koppen zone",
                      na.translate = FALSE,
                      guide = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
           label = sprintf("apparent~Q[10] == %.2f", q10), parse = TRUE) +
  labs(x = "Mean annual air temperature (°C)",
       y = expression("Transit time "*tau*" (yr, log scale)"),
       title = "Transit time versus temperature",
       subtitle = "QC-passed profiles; black line = overall loess trend") +
  theme_bw(base_size = 12) + theme(legend.position = "right")

ggsave(file.path(PLOT_DIR, "22_temperature_tau_scatter.png"),
       p, width = 7.6, height = 5.2, dpi = 200, bg = "white")
cat("OK  ", file.path(PLOT_DIR, "22_temperature_tau_scatter.png"), "\n")
