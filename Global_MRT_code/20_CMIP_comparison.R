# =============================================================================
# Step 20: CMIP6 Soil Carbon Turnover Time — Derivation and Comparison
# =============================================================================
# STREAMLINED VERSION (PNG only)
#
# Changes from previous version:
#   - CESM2 now uses cSoilAbove1m (0-1 m) for depth comparability with obs
#   - Removed duplicate heatmap block at end of script (kept log-ratio form)
#   - Area-weighted global means in summary table
#   - All library() calls moved to top
#   - All plots saved as PNG only
#   - Biome plot: single geom_col (rect drawn first as background)
#   - Added cmip6_dispersion_stats.csv (CV, 5/95 percentile ratio, log-SD)
#
# PURPOSE:
#   Derives tau = cSoil / rh from 6 structurally independent CMIP6 ESMs
#   (1980-2014) and compares the ensemble against the RF-predicted global
#   tau map (M7).
#
# CMIP6 MODELS (one per LSM family, following Varney et al. 2022):
#   CESM2          CLM5        (NCAR, USA)          -- uses cSoilAbove1m
#   UKESM1-0-LL    JULES       (MOHC, UK)
#   IPSL-CM6A-LR   ORCHIDEE    (IPSL, France)
#   MPI-ESM1-2-LR  JSBACH      (MPI-M, Germany)
#   CanESM5        CLASS-CTEM  (CCCma, Canada)
#   MIROC-ES2L     VISIT-e     (MIROC, Japan)
#
# DATA LAYOUT:
#   ./Global_MRT_code/cmip6/[MODEL_ID]/cSoil/         *.nc
#   ./Global_MRT_code/cmip6/[MODEL_ID]/cSoilAbove1m/  *.nc  (CESM2 only)
#   ./Global_MRT_code/cmip6/[MODEL_ID]/rh/            *.nc
#
# REFERENCE:
#   Varney et al. (2022) Biogeosciences 19:4671-4704
#   tau = cSoil / rh under quasi-equilibrium (Koven et al. 2017)
# =============================================================================

# --- Libraries (all at top) ---------------------------------------------------

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(sf)
library(rnaturalearth)
library(patchwork)

source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# =============================================================================
# 0. Configuration
# =============================================================================

CMIP6_DATA_DIR <- "./Global_MRT_code/cmip6"
OUTPUT_DIR     <- "./Global_MRT_code/outputs/cmip6"
PLOT_DIR       <- "./Global_MRT_code/plots/step_20"
RF_M7_PATH     <- "./Global_MRT_code/outputs/MRT_predictions/MRT_M7_full.tif"

YEAR_START <- 1980
YEAR_END   <- 2014

TARGET_RES <- 1.0
TARGET_EXT <- ext(-180, 180, -90, 90)
TARGET_CRS <- "EPSG:4326"

SEC_PER_YEAR <- 60 * 60 * 24 * 365.25

# FIX: CESM2 uses cSoilAbove1m for 0-1 m comparability with obs-based estimates
MODELS <- data.frame(
  model_id    = c("CESM2",        "UKESM1-0-LL", "IPSL-CM6A-LR",
                  "MPI-ESM1-2-LR", "CanESM5",     "MIROC-ES2L"),
  lsm         = c("CLM5",         "JULES",       "ORCHIDEE",
                  "JSBACH",       "CLASS-CTEM",  "VISIT-e"),
  institution = c("NCAR",         "MOHC",        "IPSL",
                  "MPI-M",        "CCCma",       "MIROC"),
  csoil_var   = c("cSoilAbove1m", "cSoil",       "cSoil",
                  "cSoil",        "cSoil",       "cSoil"),
  stringsAsFactors = FALSE
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,   recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Helpers
# =============================================================================

#' Load all NetCDF files for a model/variable, subset to reference period,
#' return the temporal mean.
load_cmip6_mean <- function(model_id, variable, year_start, year_end) {
  
  var_dir  <- file.path(CMIP6_DATA_DIR, model_id, variable)
  nc_files <- list.files(var_dir, pattern = "\\.nc$", full.names = TRUE)
  
  if (length(nc_files) == 0)
    stop(sprintf("No NetCDF files found for %s / %s in:\n  %s",
                 model_id, variable, var_dir))
  
  cat(sprintf("    Loading %s / %s (%d file(s))...\n",
              model_id, variable, length(nc_files)))
  
  r     <- rast(nc_files)
  times <- time(r)
  years <- as.integer(format(times, "%Y"))
  idx   <- which(years >= year_start & years <= year_end)
  
  if (length(idx) == 0)
    stop(sprintf("No time steps in %d-%d for %s / %s",
                 year_start, year_end, model_id, variable))
  
  r_mean        <- mean(r[[idx]], na.rm = TRUE)
  names(r_mean) <- variable
  
  if (is.na(crs(r_mean)) || crs(r_mean) == "") {
    crs(r_mean) <- "EPSG:4326"
  }
  
  return(r_mean)
}

#' Area-weighted global mean of a SpatRaster
area_weighted_mean <- function(r) {
  a  <- cellSize(r, unit = "km")
  v  <- values(r)[, 1]
  w  <- values(a)[, 1]
  ok <- !is.na(v) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(v[ok], w[ok])
}

#' Area-weighted zonal mean by 1-degree latitude band
zonal_mean_lat <- function(r, lat_res = 1.0) {
  areas  <- cellSize(r, unit = "km")
  lats   <- seq(-89.5, 89.5, by = lat_res)
  result <- data.frame(lat = lats, tau_mean = NA_real_)
  for (j in seq_along(lats)) {
    b_ext <- ext(-180, 180,
                 lats[j] - lat_res / 2,
                 lats[j] + lat_res / 2)
    r_b   <- crop(r,     b_ext)
    a_b   <- crop(areas, b_ext)
    vals  <- values(r_b)[, 1]
    wts   <- values(a_b)[, 1]
    ok    <- !is.na(vals) & !is.na(wts) & wts > 0
    if (any(ok))
      result$tau_mean[j] <- weighted.mean(vals[ok], wts[ok])
  }
  return(result)
}

# =============================================================================
# 2. Derive tau per model (with caching)
# =============================================================================

cat("\n=== Deriving tau per CMIP6 model ===\n\n")

tau_rasters <- list()

for (i in seq_len(nrow(MODELS))) {
  
  model_id   <- MODELS$model_id[i]
  csoil_var  <- MODELS$csoil_var[i]
  cache_file <- file.path(OUTPUT_DIR, sprintf("cmip6_tau_%s.tif", model_id))
  
  if (file.exists(cache_file)) {
    cat(sprintf("[CACHED] %s\n", model_id))
    tau_rasters[[model_id]] <- rast(cache_file)
    next
  }
  
  cat(sprintf("[COMPUTING] %s  (%s | %s | csoil=%s)\n",
              model_id, MODELS$lsm[i], MODELS$institution[i], csoil_var))
  
  tryCatch({
    
    csoil   <- load_cmip6_mean(model_id, csoil_var, YEAR_START, YEAR_END)
    rh_raw  <- load_cmip6_mean(model_id, "rh",      YEAR_START, YEAR_END)
    rh_yr   <- rh_raw * SEC_PER_YEAR
    
    target_r <- rast(TARGET_EXT, resolution = TARGET_RES, crs = TARGET_CRS)
    csoil_1  <- project(csoil,  target_r, method = "bilinear")
    rh_1     <- project(rh_yr,  target_r, method = "bilinear")
    
    rh_1[rh_1 <= 0] <- NA
    tau <- csoil_1 / rh_1
    tau[tau < 1]     <- NA
    tau[tau > 10000] <- NA
    
    names(tau) <- model_id
    
    writeRaster(tau, cache_file, overwrite = TRUE)
    tau_rasters[[model_id]] <- tau
    
    cat(sprintf("  tau range: %.1f - %.1f yr  (area-weighted mean: %.1f yr)\n",
                global(tau, "min", na.rm = TRUE)[[1]],
                global(tau, "max", na.rm = TRUE)[[1]],
                area_weighted_mean(tau)))
    
  }, error = function(e) {
    message(sprintf("  ERROR for %s: %s", model_id, e$message))
  })
}

# =============================================================================
# 3. Ensemble statistics
# =============================================================================

cat("\n=== Computing ensemble statistics ===\n\n")

tau_stack   <- rast(tau_rasters)
tau_ensmean <- mean(tau_stack,  na.rm = TRUE)
tau_enssd   <- stdev(tau_stack, na.rm = TRUE)

writeRaster(tau_ensmean,
            file.path(OUTPUT_DIR, "cmip6_tau_ensemble_mean.tif"),
            overwrite = TRUE)
writeRaster(tau_enssd,
            file.path(OUTPUT_DIR, "cmip6_tau_ensemble_sd.tif"),
            overwrite = TRUE)

cat(sprintf("Ensemble area-weighted global tau: %.1f yr\n",
            area_weighted_mean(tau_ensmean)))

# =============================================================================
# 4. Bias map: RF M7 minus CMIP6 ensemble mean
# =============================================================================

cat("\n=== Computing bias map ===\n\n")

rf_1 <- NULL
if (file.exists(RF_M7_PATH)) {
  
  rf_m7    <- rast(RF_M7_PATH)
  target_r <- rast(TARGET_EXT, resolution = TARGET_RES, crs = TARGET_CRS)
  rf_1     <- project(rf_m7, target_r, method = "bilinear")
  
  bias        <- rf_1 - tau_ensmean
  names(bias) <- "RF_M7_minus_CMIP6_ensemble"
  
  writeRaster(bias,
              file.path(OUTPUT_DIR, "cmip6_tau_bias.tif"),
              overwrite = TRUE)
  
  cat(sprintf("Area-weighted mean bias (RF - CMIP6): %.1f yr\n",
              area_weighted_mean(bias)))
  
} else {
  warning("RF M7 raster not found at: ", RF_M7_PATH,
          "\nSkipping bias map. Run step 14 first.")
}

# =============================================================================
# 5. Zonal means + summary + dispersion stats
# =============================================================================

cat("\n=== Computing zonal means ===\n\n")

zonal_list <- list()

for (model_id in names(tau_rasters)) {
  cat(sprintf("  Zonal means: %s\n", model_id))
  zm <- zonal_mean_lat(tau_rasters[[model_id]])
  zm$source <- model_id
  zm$type   <- "CMIP6"
  zonal_list[[model_id]] <- zm
}

zm_ens        <- zonal_mean_lat(tau_ensmean)
zm_ens$source <- "CMIP6_ensemble"
zm_ens$type   <- "CMIP6_ensemble"
zonal_list[["CMIP6_ensemble"]] <- zm_ens

if (!is.null(rf_1)) {
  cat("  Zonal means: RF M7\n")
  zm_rf        <- zonal_mean_lat(rf_1)
  zm_rf$source <- "RF_M7"
  zm_rf$type   <- "RF"
  zonal_list[["RF_M7"]] <- zm_rf
}

zonal_df <- bind_rows(zonal_list)
write.csv(zonal_df,
          file.path(OUTPUT_DIR, "cmip6_zonal_means.csv"),
          row.names = FALSE)

# Global mean summary (area-weighted)
summary_df <- bind_rows(lapply(names(tau_rasters), function(m) {
  data.frame(
    model              = m,
    lsm                = MODELS$lsm[MODELS$model_id == m],
    global_mean_tau_yr = round(area_weighted_mean(tau_rasters[[m]]), 1)
  )
}))
if (!is.null(rf_1)) {
  summary_df <- bind_rows(summary_df, data.frame(
    model = "RF_M7", lsm = "RandomForest",
    global_mean_tau_yr = round(area_weighted_mean(rf_1), 1)
  ))
}
write.csv(summary_df,
          file.path(OUTPUT_DIR, "cmip6_tau_summary.csv"),
          row.names = FALSE)
cat("\nGlobal mean tau summary (area-weighted):\n")
print(summary_df)

# --- NEW: dispersion statistics ----------------------------------------------
# Magnitude-free measures of spatial spread in tau — use these for any
# quantitative variance comparison between products in the manuscript.

cat("\n=== Computing dispersion statistics ===\n\n")

dispersion_one <- function(r, label) {
  v <- values(r)[, 1]
  v <- v[!is.na(v) & v > 0]
  p05 <- as.numeric(quantile(v, 0.05))
  p95 <- as.numeric(quantile(v, 0.95))
  data.frame(
    source        = label,
    mean          = round(mean(v), 1),
    sd            = round(sd(v), 1),
    cv            = round(sd(v) / mean(v), 3),
    p05           = round(p05, 1),
    p95           = round(p95, 1),
    p95_p05_ratio = round(p95 / p05, 2),
    log_sd        = round(sd(log(v)), 3),
    stringsAsFactors = FALSE
  )
}

dispersion_df <- bind_rows(lapply(names(tau_rasters), function(m) {
  dispersion_one(tau_rasters[[m]], m)
}))
if (!is.null(rf_1)) {
  dispersion_df <- bind_rows(dispersion_df, dispersion_one(rf_1, "RF_M7"))
}
write.csv(dispersion_df,
          file.path(OUTPUT_DIR, "cmip6_dispersion_stats.csv"),
          row.names = FALSE)
cat("Dispersion statistics:\n")
print(dispersion_df)

# =============================================================================
# 6. Plot A - Latitudinal profile
# =============================================================================

cat("\n=== Plot A: Latitudinal profile ===\n\n")

cmip6_spread <- zonal_df %>%
  filter(type == "CMIP6") %>%
  group_by(lat) %>%
  summarise(tau_min = min(tau_mean, na.rm = TRUE),
            tau_max = max(tau_mean, na.rm = TRUE), .groups = "drop")

model_colors <- setNames(
  RColorBrewer::brewer.pal(6, "Set2"),
  MODELS$model_id
)

p_profile <- ggplot() +
  geom_ribbon(data = cmip6_spread,
              aes(x = lat, ymin = tau_min, ymax = tau_max),
              fill = "steelblue", alpha = 0.18) +
  geom_line(data = filter(zonal_df, type == "CMIP6"),
            aes(x = lat, y = tau_mean, color = source),
            linewidth = 0.5, alpha = 0.75) +
  geom_line(data = filter(zonal_df, type == "CMIP6_ensemble"),
            aes(x = lat, y = tau_mean),
            color = "steelblue4", linewidth = 1.0, linetype = "dashed") +
  geom_line(data = filter(zonal_df, type == "RF"),
            aes(x = lat, y = tau_mean),
            color = "firebrick", linewidth = 1.3) +
  scale_color_manual(values = model_colors, name = "CMIP6 model") +
  scale_y_log10(
    breaks = c(1, 5, 10, 50, 100, 500, 1000),
    labels = c("1", "5", "10", "50", "100", "500", "1000")
  ) +
  coord_flip() +
  labs(
    x        = "Latitude",
    y        = expression(paste("Mean residence time ", tau, " (years, log scale)")),
    title    = "Latitudinal profile of soil carbon turnover time",
    subtitle = paste0("RF M7 (red) vs. CMIP6 ensemble ",
                      YEAR_START, "\u2013", YEAR_END,
                      " (blue); shading = model spread")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(PLOT_DIR, "fig_latitudinal_profile.png"),
       p_profile, width = 6, height = 8, dpi = 300)
cat("Saved: fig_latitudinal_profile.png\n")

# =============================================================================
# 7. Plot B - Latitude x model heatmap  (log-ratio tau / mean)
# =============================================================================
#
# Each column normalised to its own global latitudinal mean (tau / tau_bar),
# then log10-transformed: 0 = product's mean, +1 = 10x mean, -1 = 0.1x mean.
# This preserves magnitude information around each product's centre and is
# the version to cite in the manuscript when discussing relative structure.
# =============================================================================

cat("\n=== Plot B: Latitude x model heatmap (log-ratio) ===\n\n")

label_map <- c(
  "CMIP6_ensemble" = "Ensemble\nmean",
  "CESM2"          = "CESM2\n(CLM5)",
  "UKESM1-0-LL"    = "UKESM1\n(JULES)",
  "IPSL-CM6A-LR"   = "IPSL\n(ORCHIDEE)",
  "MPI-ESM1-2-LR"  = "MPI\n(JSBACH)",
  "CanESM5"        = "CanESM5\n(CLASS)",
  "MIROC-ES2L"     = "MIROC\n(VISIT-e)",
  "GAP"            = "",
  "RF_M7"          = "RF M7\n(this study)"
)
band_order <- names(label_map)

gap_df <- data.frame(
  lat      = seq(-89.5, 89.5, by = 1.0),
  source   = "GAP",
  tau_mean = NA_real_,
  type     = "gap"
)

source_means <- zonal_df %>%
  filter(source %in% setdiff(band_order, "GAP"), !is.na(tau_mean)) %>%
  group_by(source) %>%
  summarise(source_mean = mean(tau_mean, na.rm = TRUE), .groups = "drop")

heatmap_df <- zonal_df %>%
  filter(source %in% setdiff(band_order, "GAP")) %>%
  left_join(source_means, by = "source") %>%
  mutate(tau_rel = tau_mean / source_mean) %>%
  bind_rows(gap_df %>% mutate(source_mean = NA, tau_rel = NA)) %>%
  mutate(
    source      = factor(source, levels = band_order),
    log_tau_rel = log10(tau_rel)
  )

p_heatmap <- ggplot(heatmap_df,
                    aes(x = source, y = lat, fill = log_tau_rel)) +
  geom_tile(width = 0.9) +
  scale_fill_gradientn(
    colours  = viridis::turbo(256),
    limits   = c(-1.5, 1.5),
    na.value = "white",
    name     = expression(tau / bar(tau)),
    breaks   = c(-1, 0, 1),
    labels   = c("0.1x", "mean", "10x")
  ) +
  geom_hline(yintercept = c(-66.5, -23.5, 0, 23.5, 66.5),
             linetype = "dotted", color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = label_map) +
  scale_y_continuous(
    breaks = c(-90, -60, -30, 0, 30, 60, 90),
    labels = c("90\u00b0S", "60\u00b0S", "30\u00b0S",
               "0\u00b0",
               "30\u00b0N", "60\u00b0N", "90\u00b0N"),
    expand = c(0, 0)
  ) +
  labs(
    x        = NULL,
    y        = "Latitude",
    title    = "Soil carbon turnover time by latitude and model",
    subtitle = paste0("Each product normalised to its own latitudinal mean ",
                      "(centre = mean, red = longer, blue = shorter).\n",
                      "CMIP6 historical ", YEAR_START, "\u2013", YEAR_END,
                      " vs. RF M7 (this study).")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x       = element_text(size = 8, lineheight = 0.9),
    axis.text.y       = element_text(size = 9),
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(3, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, color = "grey40"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(PLOT_DIR, "fig_latitude_heatmap.png"),
       p_heatmap, width = 9, height = 8, dpi = 300)
cat("Saved: fig_latitude_heatmap.png\n")

# =============================================================================
# 7b. Plot B2 - Latitude x model heatmap  (min-max normalised)
# =============================================================================
#
# Complementary view of Plot B. Each column min-max normalised to [0, 1]:
# 0 = that product's shortest zonal tau, 1 = its longest. Shows WHERE within
# a product's own range the extremes sit and how widely they are spread
# across latitude. Note: by construction every column spans [0, 1], so this
# plot cannot be used to compare absolute spread between products -- use
# cmip6_dispersion_stats.csv for that.
# =============================================================================

cat("\n=== Plot B2: Latitude x model heatmap (min-max normalised) ===\n\n")

source_stats <- zonal_df %>%
  filter(source %in% setdiff(band_order, "GAP"), !is.na(tau_mean)) %>%
  group_by(source) %>%
  summarise(
    src_min = min(tau_mean, na.rm = TRUE),
    src_max = max(tau_mean, na.rm = TRUE),
    .groups = "drop"
  )

heatmap_df_rel <- zonal_df %>%
  filter(source %in% setdiff(band_order, "GAP")) %>%
  left_join(source_stats, by = "source") %>%
  mutate(tau_norm = (tau_mean - src_min) / (src_max - src_min)) %>%
  bind_rows(gap_df %>% mutate(src_min = NA, src_max = NA, tau_norm = NA)) %>%
  mutate(source = factor(source, levels = band_order))

p_heatmap_rel <- ggplot(heatmap_df_rel,
                        aes(x = source, y = lat, fill = tau_norm)) +
  geom_tile(width = 0.9) +
  scale_fill_gradientn(
    colours  = viridis::turbo(256),
    limits   = c(0, 1),
    na.value = "white",
    name     = expression(tau[norm]),
    breaks   = c(0, 0.5, 1),
    labels   = c("min", "mid", "max")
  ) +
  geom_hline(yintercept = c(-66.5, -23.5, 0, 23.5, 66.5),
             linetype = "dotted", color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = label_map) +
  scale_y_continuous(
    breaks = c(-90, -60, -30, 0, 30, 60, 90),
    labels = c("90\u00b0S", "60\u00b0S", "30\u00b0S",
               "0\u00b0",
               "30\u00b0N", "60\u00b0N", "90\u00b0N"),
    expand = c(0, 0)
  ) +
  labs(
    x        = NULL,
    y        = "Latitude",
    title    = "Relative soil carbon turnover time by latitude and model",
    subtitle = paste0("Each product min-max normalised (0 = shortest ",
                      "\u03c4, 1 = longest \u03c4 for that product). ",
                      "Shows distribution within each product's own range.")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x       = element_text(size = 8, lineheight = 0.9),
    axis.text.y       = element_text(size = 9),
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(3, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, color = "grey40"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(PLOT_DIR, "fig_latitude_heatmap_relative.png"),
       p_heatmap_rel, width = 9, height = 8, dpi = 300)
cat("Saved: fig_latitude_heatmap_relative.png\n")

# =============================================================================
# 8. Plot C - Per-model bias maps (RF M7 minus each CMIP6 model)
# =============================================================================

cat("\n=== Plot C: Per-model bias maps ===\n\n")

if (!is.null(rf_1) && length(tau_rasters) > 0) {
  
  world      <- ne_countries(scale = "medium", returnclass = "sf")
  crs_robin  <- "+proj=robin +lon_0=0 +datum=WGS84"
  BIAS_LIMIT <- 500
  
  # Pre-project world outline once
  world_robin <- st_transform(world, crs_robin)
  
  #' Build a bias-map panel.
  #' IMPORTANT: the raster is projected to Robinson BEFORE being converted
  #' to a data.frame, so the x/y columns are in Robinson meters (matching
  #' coord_sf(crs = crs_robin)). If you skip this, the raster ends up drawn
  #' at the origin of a ~34 million-meter-wide plot and the map looks empty.
  bias_to_panel <- function(bias_r, title_str, limit = BIAS_LIMIT) {
    
    # 1. Clamp extreme values in raster space
    bias_r[bias_r >  limit] <-  limit
    bias_r[bias_r < -limit] <- -limit
    
    # 2. Project raster to Robinson (resamples to a regular Robinson grid)
    bias_proj <- project(bias_r, crs_robin, method = "bilinear")
    
    # 3. Convert to data.frame (x, y now in Robinson meters)
    df <- as.data.frame(bias_proj, xy = TRUE, na.rm = TRUE)
    colnames(df) <- c("x", "y", "bias")
    
    ggplot() +
      geom_raster(data = df, aes(x = x, y = y, fill = bias)) +
      geom_sf(data = world_robin, fill = NA, color = "grey30",
              linewidth = 0.15, inherit.aes = FALSE) +
      scale_fill_distiller(
        palette   = "RdBu",
        limits    = c(-limit, limit),
        direction = -1,
        na.value  = "white",
        name      = expression(Delta * tau ~ (yr)),
        guide     = guide_colorbar(barheight = unit(3, "cm"),
                                   barwidth  = unit(0.35, "cm"))
      ) +
      coord_sf(crs = crs_robin, expand = FALSE, datum = NA) +
      labs(title = title_str) +
      theme_void(base_size = 9) +
      theme(
        plot.title       = element_text(size = 8, face = "bold",
                                        hjust = 0.5, margin = margin(b = 2)),
        legend.position  = "right",
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 7),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "grey92", color = NA)
      )
  }
  
  panels <- list()
  for (model_id in names(tau_rasters)) {
    lsm_label     <- MODELS$lsm[MODELS$model_id == model_id]
    title_str     <- sprintf("%s\n(%s)", model_id, lsm_label)
    bias_i        <- rf_1 - tau_rasters[[model_id]]
    names(bias_i) <- "bias"
    writeRaster(bias_i,
                file.path(OUTPUT_DIR,
                          sprintf("cmip6_tau_bias_%s.tif", model_id)),
                overwrite = TRUE)
    panels[[model_id]] <- bias_to_panel(bias_i, title_str)
    cat(sprintf("  Panel ready: %s\n", model_id))
  }
  panels[["ensemble"]] <- bias_to_panel(
    bias, "Ensemble mean\n(all 6 models)"
  )
  
  p_bias_panel <- wrap_plots(panels, ncol = 4) +
    plot_annotation(
      title    = expression(paste("RF M7 ", tau, " minus CMIP6 ", tau,
                                  " (years)")),
      subtitle = paste0("Positive (red) = RF predicts longer residence time ",
                        "than ESM; negative (blue) = RF shorter.\n",
                        "CMIP6 historical mean ",
                        YEAR_START, "\u2013", YEAR_END, "."),
      theme = theme(
        plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle   = element_text(size = 8, hjust = 0.5, color = "grey40"),
        plot.background = element_rect(fill = "white", color = NA)
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
  
  ggsave(file.path(PLOT_DIR, "fig_bias_maps_panel.png"),
         p_bias_panel, width = 16, height = 8, dpi = 300)
  cat("Saved: fig_bias_maps_panel.png\n")
  
} else {
  warning("Skipping Plot C: RF M7 raster or CMIP6 tau rasters not available.")
}

# =============================================================================
# 9. Plot D - Biome-level relative tau barplot
# =============================================================================
#
# tau / global_mean per product, compared across four broad latitudinal zones:
# Tropical / Temperate / Boreal / Polar. (Arid zone is NOT implemented; this
# is a pragmatic latitude-band proxy. Replace with a Koppen-Geiger raster
# if one is available in ./spatialized_layers/ to get a true biome split.)
# =============================================================================

cat("\n=== Plot D: Biome-level relative tau barplot ===\n\n")

if (!is.null(rf_1) && length(tau_rasters) > 0) {
  
  target_r <- rast(TARGET_EXT, resolution = TARGET_RES, crs = TARGET_CRS)
  lat_r    <- init(target_r, "y")
  
  # Latitude-band proxy for biome zones (no Arid class)
  zone_r <- classify(lat_r, rbind(
    c(-90,   -66.5, 5),  # Polar S
    c(-66.5, -50,   4),  # Boreal S
    c(-50,   -23.5, 3),  # Temperate S
    c(-23.5,  23.5, 1),  # Tropical
    c( 23.5,  50,   3),  # Temperate N
    c( 50,    66.5, 4),  # Boreal N
    c( 66.5,  90,   5)   # Polar N
  ), include.lowest = TRUE)
  
  zone_labels <- c("1" = "Tropical", "3" = "Temperate",
                   "4" = "Boreal",   "5" = "Polar")
  
  zonal_biome <- function(tau_r, zone_r, label) {
    areas  <- cellSize(tau_r, unit = "km")
    zones  <- as.integer(values(zone_r)[, 1])
    tau_v  <- values(tau_r)[, 1]
    area_v <- values(areas)[, 1]
    
    result <- data.frame()
    for (z in c(1, 3, 4, 5)) {
      idx <- which(zones == z & !is.na(tau_v) &
                     !is.na(area_v) & area_v > 0)
      if (length(idx) < 10) next
      wmean <- weighted.mean(tau_v[idx], area_v[idx], na.rm = TRUE)
      result <- rbind(result, data.frame(
        source = label,
        zone   = zone_labels[as.character(z)],
        tau    = wmean
      ))
    }
    return(result)
  }
  
  biome_list <- list()
  for (model_id in names(tau_rasters)) {
    biome_list[[model_id]] <- zonal_biome(
      tau_rasters[[model_id]], zone_r, model_id
    )
  }
  biome_list[["RF_M7"]] <- zonal_biome(rf_1, zone_r, "RF_M7")
  biome_df <- bind_rows(biome_list)
  
  # Normalise within each product: tau_rel = tau / area-weighted global mean
  global_means <- bind_rows(lapply(names(tau_rasters), function(m) {
    data.frame(source = m, global_mean = area_weighted_mean(tau_rasters[[m]]))
  }))
  global_means <- bind_rows(global_means, data.frame(
    source = "RF_M7", global_mean = area_weighted_mean(rf_1)
  ))
  
  biome_norm <- biome_df %>%
    left_join(global_means, by = "source") %>%
    mutate(
      tau_rel = tau / global_mean,
      label   = recode(source,
                       "RF_M7"         = "RF M7\n(this study)",
                       "CESM2"         = "CESM2\n(CLM5)",
                       "UKESM1-0-LL"   = "UKESM1\n(JULES)",
                       "IPSL-CM6A-LR"  = "IPSL\n(ORCHIDEE)",
                       "MPI-ESM1-2-LR" = "MPI\n(JSBACH)",
                       "CanESM5"       = "CanESM5\n(CLASS)",
                       "MIROC-ES2L"    = "MIROC\n(VISIT-e)"
      ),
      zone = factor(zone, levels = c("Tropical", "Temperate",
                                     "Boreal",   "Polar"))
    )
  
  zone_colors <- c(
    "Tropical"  = "#d62828",
    "Temperate" = "#f77f00",
    "Boreal"    = "#4361ee",
    "Polar"     = "#7209b7"
  )
  
  model_order <- c("RF M7\n(this study)",
                   "CESM2\n(CLM5)",    "UKESM1\n(JULES)",
                   "IPSL\n(ORCHIDEE)", "MPI\n(JSBACH)",
                   "CanESM5\n(CLASS)", "MIROC\n(VISIT-e)")
  
  biome_norm <- biome_norm %>%
    mutate(label = factor(label, levels = model_order))
  
  # Clean highlight: rect first (background), then a single geom_col
  p_biome <- ggplot(biome_norm, aes(x = label, y = tau_rel, fill = zone)) +
    annotate("rect",
             xmin = 0.5, xmax = 1.5,
             ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey40", linewidth = 0.6) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.75, alpha = 0.9) +
    scale_fill_manual(values = zone_colors, name = "Biome zone") +
    scale_y_continuous(
      breaks = c(0.25, 0.5, 1, 2, 4, 8),
      trans  = "log2",
      labels = function(x) ifelse(x == 1, "1 (mean)", as.character(x))
    ) +
    labs(
      x        = NULL,
      y        = expression(paste("Relative ", tau,
                                  "  (product mean = 1, log2 scale)")),
      title    = expression(paste("Relative spatial structure of soil ",
                                  "carbon ", tau, " across biome zones")),
      subtitle = paste0(
        "Each product normalised to its own global mean to remove depth/",
        "pool-definition differences.\n",
        "Values > 1 = longer than average; < 1 = shorter than average. ",
        "Grey band = RF M7 (this study)."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x        = element_text(size = 8, lineheight = 0.9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right",
      plot.background    = element_rect(fill = "white", color = NA),
      plot.title         = element_text(face = "bold", size = 11),
      plot.subtitle      = element_text(size = 8, color = "grey40")
    )
  
  ggsave(file.path(PLOT_DIR, "fig_biome_relative_tau.png"),
         p_biome, width = 10, height = 6, dpi = 300)
  cat("Saved: fig_biome_relative_tau.png\n")
  
  write.csv(biome_norm %>%
              select(source, zone, tau, global_mean, tau_rel),
            file.path(OUTPUT_DIR, "cmip6_biome_relative_tau.csv"),
            row.names = FALSE)
  
} else {
  warning("Skipping Plot D: RF M7 or CMIP6 rasters not available.")
}

# =============================================================================
# 10. Done
# =============================================================================

cat("\n\u2550\u2550\u2550 Step 20 complete \u2550\u2550\u2550\n")
cat(sprintf("Models processed : %d / 6\n", length(tau_rasters)))
cat(sprintf("Outputs          : %s\n", OUTPUT_DIR))
cat(sprintf("Plots            : %s\n", PLOT_DIR))
cat("\nOutputs summary:\n")
cat("  Rasters : cmip6_tau_[model].tif, ensemble_mean/sd, bias maps\n")
cat("  Tables  : cmip6_tau_summary.csv, cmip6_zonal_means.csv,\n")
cat("            cmip6_dispersion_stats.csv, cmip6_biome_relative_tau.csv\n")
cat("  Plots   : fig_latitudinal_profile.png, fig_latitude_heatmap.png,\n")
cat("            fig_latitude_heatmap_relative.png, fig_bias_maps_panel.png,\n")
cat("            fig_biome_relative_tau.png\n")