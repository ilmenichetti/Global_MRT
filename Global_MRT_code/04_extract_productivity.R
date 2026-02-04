# =============================================================================
# Step 4: Extract MODIS NPP productivity data
# Load annual NPP rasters, extract at points, calculate equilibrium diagnostics
# using only years BEFORE sampling date
# =============================================================================

library(terra)
library(dplyr)
library(ggplot2)

# Source config and utilities
source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# Configuration
NPP_DIR <- "./Global_MRT_code/spatialized_layers/productivity/MODIS_NPP"
NPP_YEARS <- 2001:2020
TREND_PCT_THRESHOLD <- 2  # Flag if |trend| > 2% of mean per year

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║     Step 4: Extract MODIS NPP Productivity Data               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

extract_year <- function(date_str) {
  # Extract year from site_obsdate (handles multiple formats)
  if (is.na(date_str) || date_str == "" || date_str == "NULL" || date_str == "0") {
    return(NA_integer_)
  }
  
  # 4-digit year
  if (grepl("^\\d{4}$", date_str)) {
    yr <- as.integer(date_str)
    if (yr < 1900 || yr > 2025) return(NA_integer_)
    return(yr)
  }
  
  # Full date (YYYY-MM-DD)
  if (grepl("^\\d{4}-", date_str)) {
    yr <- as.integer(substr(date_str, 1, 4))
    if (yr < 1900 || yr > 2025) return(NA_integer_)
    return(yr)
  }
  
  # Short date (YY-MM-DD like "12-09-14")
  if (grepl("^\\d{2}-\\d{2}-\\d{2}$", date_str)) {
    yy <- as.integer(substr(date_str, 1, 2))
    yr <- ifelse(yy <= 25, 2000L + yy, 1900L + yy)
    return(yr)
  }
  
  return(NA_integer_)
}

load_npp_stack <- function(npp_dir = NPP_DIR, years = NPP_YEARS) {
  # Load all annual NPP rasters into a stack
  log_step("Loading annual NPP rasters")
  
  npp_files <- c()
  loaded_years <- c()
  
  for (year in years) {
    # Try different naming patterns
    patterns <- c(
      sprintf("NPP_%d_global.tif", year),
      sprintf("NPP_%d_0.1deg.tif", year),
      sprintf("NPP_%d.tif", year)
    )
    
    found <- FALSE
    for (pattern in patterns) {
      filepath <- file.path(npp_dir, pattern)
      if (file.exists(filepath)) {
        npp_files <- c(npp_files, filepath)
        loaded_years <- c(loaded_years, year)
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      warning(sprintf("  NPP file not found for year %d", year))
    }
  }
  
  if (length(npp_files) == 0) {
    stop("No NPP files found in ", npp_dir)
  }
  
  cat(sprintf("  Found %d annual NPP files\n", length(npp_files)))
  
  npp_stack <- rast(npp_files)
  names(npp_stack) <- paste0("NPP_", loaded_years)
  
  cat(sprintf("  Stack dimensions: %d x %d x %d layers\n", 
              nrow(npp_stack), ncol(npp_stack), nlyr(npp_stack)))
  cat(sprintf("  Years covered: %d - %d\n", min(loaded_years), max(loaded_years)))
  
  return(npp_stack)
}

extract_npp_at_points <- function(npp_stack, points_data) {
  # Extract NPP values at point locations for all years
  log_step("Extracting NPP at point locations")
  
  pts <- terra::vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %d points across %d years...\n", 
              nrow(points_data), nlyr(npp_stack)))
  
  npp_values <- terra::extract(npp_stack, pts, ID = FALSE)
  npp_matrix <- as.matrix(npp_values)
  
  cat(sprintf("  ✓ Extracted %d values\n", sum(!is.na(npp_matrix))))
  cat(sprintf("  Coverage: %.1f%%\n", 
              100 * sum(!is.na(npp_matrix)) / length(npp_matrix)))
  
  return(npp_matrix)
}

calculate_npp_diagnostics_by_date <- function(npp_matrix, obs_years, npp_years = NPP_YEARS) {
  # Calculate temporal statistics using only years BEFORE sampling
  log_step("Calculating NPP diagnostics (date-aware)")
  
  n_sites <- nrow(npp_matrix)
  
  cat(sprintf("  Processing %d sites\n", n_sites))
  cat(sprintf("  NPP years available: %d - %d\n", min(npp_years), max(npp_years)))
  
  # Initialize output vectors
  npp_mean <- rep(NA_real_, n_sites)
  npp_sd <- rep(NA_real_, n_sites)
  npp_cv <- rep(NA_real_, n_sites)
  npp_trend <- rep(NA_real_, n_sites)
  npp_trend_pct <- rep(NA_real_, n_sites)
  npp_n_years <- rep(NA_integer_, n_sites)
  npp_years_used <- rep(NA_character_, n_sites)
  
  flag_no_date <- rep(FALSE, n_sites)
  flag_pre_npp <- rep(FALSE, n_sites)
  flag_insufficient_years <- rep(FALSE, n_sites)
  
  # Helper function to calculate stats for a subset of years
  calc_stats <- function(npp_vals, years_used) {
    valid_idx <- !is.na(npp_vals)
    n_valid <- sum(valid_idx)
    
    if (n_valid < 3) return(NULL)
    
    npp_valid <- npp_vals[valid_idx]
    years_valid <- years_used[valid_idx]
    
    mean_val <- mean(npp_valid)
    sd_val <- sd(npp_valid)
    cv_val <- sd_val / mean_val
    
    # Linear trend
    if (n_valid >= 3) {
      trend_fit <- lm(npp_valid ~ years_valid)
      trend_val <- coef(trend_fit)[2]
      trend_pct_val <- 100 * (trend_val / mean_val)
    } else {
      trend_val <- NA_real_
      trend_pct_val <- NA_real_
    }
    
    list(
      mean = mean_val,
      sd = sd_val,
      cv = cv_val,
      trend = trend_val,
      trend_pct = trend_pct_val,
      n_years = n_valid,
      years_used = paste(sort(years_valid), collapse = ",")
    )
  }
  
  # Process each site
  for (i in 1:n_sites) {
    if (i %% 50000 == 0) {
      cat(sprintf("  Progress: %d / %d (%.1f%%)\n", i, n_sites, 100 * i / n_sites))
    }
    
    obs_yr <- obs_years[i]
    site_npp <- npp_matrix[i, ]
    
    # Case 1: No observation date
    if (is.na(obs_yr)) {
      flag_no_date[i] <- TRUE
      stats <- calc_stats(site_npp, npp_years)
      if (!is.null(stats)) {
        npp_mean[i] <- stats$mean
        npp_sd[i] <- stats$sd
        npp_cv[i] <- stats$cv
        npp_trend[i] <- stats$trend
        npp_trend_pct[i] <- stats$trend_pct
        npp_n_years[i] <- stats$n_years
        npp_years_used[i] <- stats$years_used
      }
      next
    }
    
    # Case 2: Sampled before NPP era
    if (obs_yr < min(npp_years)) {
      flag_pre_npp[i] <- TRUE
      next
    }
    
    # Case 3: Use only years before observation
    years_before <- npp_years[npp_years < obs_yr]
    
    if (length(years_before) < 3) {
      flag_insufficient_years[i] <- TRUE
      next
    }
    
    npp_before <- site_npp[npp_years < obs_yr]
    stats <- calc_stats(npp_before, years_before)
    
    if (!is.null(stats)) {
      npp_mean[i] <- stats$mean
      npp_sd[i] <- stats$sd
      npp_cv[i] <- stats$cv
      npp_trend[i] <- stats$trend
      npp_trend_pct[i] <- stats$trend_pct
      npp_n_years[i] <- stats$n_years
      npp_years_used[i] <- stats$years_used
    } else {
      flag_insufficient_years[i] <- TRUE
    }
  }
  
  cat("  ✓ Diagnostics calculated\n")
  
  # Report summary
  cat(sprintf("\n  Sites with NPP data: %d (%.1f%%)\n", 
              sum(!is.na(npp_mean)), 
              100 * sum(!is.na(npp_mean)) / n_sites))
  cat(sprintf("  Sites flagged - no date: %d\n", sum(flag_no_date)))
  cat(sprintf("  Sites flagged - pre-NPP era: %d\n", sum(flag_pre_npp)))
  cat(sprintf("  Sites flagged - insufficient years: %d\n", sum(flag_insufficient_years)))
  
  return(data.frame(
    npp_mean = npp_mean,
    npp_sd = npp_sd,
    npp_cv = npp_cv,
    npp_trend = npp_trend,
    npp_trend_pct = npp_trend_pct,
    npp_n_years = npp_n_years,
    npp_years_used = npp_years_used,
    npp_flag_no_date = flag_no_date,
    npp_flag_pre_npp_era = flag_pre_npp,
    npp_flag_insufficient = flag_insufficient_years
  ))
}

# -----------------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------------

log_step("Loading climate data")
data <- readRDS("./Global_MRT_code/outputs/03_with_climate.rds")

cat(sprintf("  Loaded %d sites\n", nrow(data)))

# Extract observation years
log_step("Extracting observation years")
data$obs_year <- sapply(data$site_obsdate, extract_year)

cat(sprintf("  Sites with valid year: %d (%.1f%%)\n",
            sum(!is.na(data$obs_year)),
            100 * sum(!is.na(data$obs_year)) / nrow(data)))
cat(sprintf("  Year range: %d - %d\n",
            min(data$obs_year, na.rm = TRUE),
            max(data$obs_year, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# Extract NPP data
# -----------------------------------------------------------------------------

# Load NPP rasters
npp_stack <- load_npp_stack()

# Extract at points
npp_matrix <- extract_npp_at_points(npp_stack, data)

# Calculate diagnostics
npp_diagnostics <- calculate_npp_diagnostics_by_date(
  npp_matrix, 
  data$obs_year, 
  NPP_YEARS
)

# Add to data
data <- cbind(data, npp_diagnostics)

# -----------------------------------------------------------------------------
# Calculate equilibrium flags
# -----------------------------------------------------------------------------

log_step("Calculating equilibrium flags")

# High CV flag
data$npp_flag_high_cv <- !is.na(data$npp_cv) & data$npp_cv > 0.3

# Significant trend flag
data$npp_flag_trend <- !is.na(data$npp_trend_pct) & abs(data$npp_trend_pct) > TREND_PCT_THRESHOLD

# Trend direction (for sites with trend data)
data$npp_trend_direction <- NA_character_
data$npp_trend_direction[data$npp_trend_pct > TREND_PCT_THRESHOLD] <- "increasing"
data$npp_trend_direction[data$npp_trend_pct < -TREND_PCT_THRESHOLD] <- "decreasing"
data$npp_trend_direction[abs(data$npp_trend_pct) <= TREND_PCT_THRESHOLD & !is.na(data$npp_trend_pct)] <- "stable"

# Combined non-equilibrium flag
data$npp_non_equilibrium <- (
  data$npp_flag_no_date | 
    data$npp_flag_pre_npp_era | 
    data$npp_flag_insufficient | 
    data$npp_flag_high_cv | 
    data$npp_flag_trend
)

cat(sprintf("  Sites in equilibrium: %d (%.1f%%)\n",
            sum(!data$npp_non_equilibrium, na.rm = TRUE),
            100 * sum(!data$npp_non_equilibrium, na.rm = TRUE) / nrow(data)))

# -----------------------------------------------------------------------------
# Save output
# -----------------------------------------------------------------------------

log_step("Saving results")

output_path <- "./Global_MRT_code/outputs/04_with_productivity.rds"
saveRDS(data, output_path)

cat(sprintf("  ✓ Saved: %s\n", output_path))
cat(sprintf("  Output contains %d sites with %d variables\n", nrow(data), ncol(data)))

# -----------------------------------------------------------------------------
# Create diagnostic plots
# -----------------------------------------------------------------------------

log_step("Creating diagnostic plots")

plot_dir <- "./Global_MRT_code/plots/04_extract_productivity"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Get world map
world <- map_data("world")

# Subset to sites with NPP data
sites_with_npp <- data %>%
  filter(!is.na(npp_mean)) %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees, 
           npp_mean, npp_trend_direction, npp_non_equilibrium)

# Plot 1: NPP coverage map
p1 <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray70", linewidth = 0.1) +
  geom_point(data = sites_with_npp,
             aes(x = longitude_decimal_degrees, y = latitude_decimal_degrees,
                 color = npp_mean),
             size = 0.3, alpha = 0.6) +
  scale_color_viridis_c(name = "NPP\n(kg C/m²/yr)", option = "viridis") +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
  labs(title = "Sites with MODIS NPP data",
       subtitle = sprintf("n = %s sites", format(nrow(sites_with_npp), big.mark = ",")),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid = element_line(color = "gray95"))

ggsave(file.path(plot_dir, "01_npp_coverage_map.png"), p1,
       width = 12, height = 6, dpi = 150)
cat(sprintf("  ✓ Saved: %s\n", file.path(plot_dir, "01_npp_coverage_map.png")))

# Plot 2: Trend direction map
sites_with_trend <- sites_with_npp %>%
  filter(!is.na(npp_trend_direction))

p2 <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray70", linewidth = 0.1) +
  geom_point(data = sites_with_trend,
             aes(x = longitude_decimal_degrees, y = latitude_decimal_degrees,
                 color = npp_trend_direction),
             size = 0.3, alpha = 0.6) +
  scale_color_manual(name = "NPP Trend",
                     values = c("increasing" = "forestgreen",
                                "stable" = "gray50",
                                "decreasing" = "firebrick")) +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
  labs(title = "NPP trend direction at sample sites",
       subtitle = sprintf("n = %s sites with trend data", 
                          format(nrow(sites_with_trend), big.mark = ",")),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid = element_line(color = "gray95"))

ggsave(file.path(plot_dir, "02_npp_trend_map.png"), p2,
       width = 12, height = 6, dpi = 150)
cat(sprintf("  ✓ Saved: %s\n", file.path(plot_dir, "02_npp_trend_map.png")))

# Plot 3: NPP distribution
p3 <- ggplot(sites_with_npp, aes(x = npp_mean)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(sites_with_npp$npp_mean, na.rm = TRUE),
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of mean NPP across sites",
       subtitle = sprintf("Median = %.1f kg C/m²/yr (red dashed line)",
                          median(sites_with_npp$npp_mean, na.rm = TRUE)),
       x = "Mean NPP (kg C/m²/year)", y = "Number of sites") +
  theme_minimal()

ggsave(file.path(plot_dir, "03_npp_distribution.png"), p3,
       width = 8, height = 5, dpi = 150)
cat(sprintf("  ✓ Saved: %s\n", file.path(plot_dir, "03_npp_distribution.png")))

# -----------------------------------------------------------------------------
# Generate summary tables
# -----------------------------------------------------------------------------

log_step("Creating summary tables")

# Variable coverage table
npp_vars <- c("npp_mean", "npp_sd", "npp_cv", "npp_trend", "npp_trend_pct",
              "npp_trend_direction", "npp_n_years", "npp_years_used",
              "npp_flag_no_date", "npp_flag_pre_npp_era", "npp_flag_insufficient",
              "npp_flag_high_cv", "npp_flag_trend", "npp_non_equilibrium")

coverage_table <- data.frame(
  variable = npp_vars,
  n_total = nrow(data),
  n_valid = sapply(npp_vars, function(v) {
    if (v %in% names(data)) sum(!is.na(data[[v]])) else NA
  }),
  stringsAsFactors = FALSE
)

coverage_table$pct_valid <- round(100 * coverage_table$n_valid / coverage_table$n_total, 1)
coverage_table$n_missing <- coverage_table$n_total - coverage_table$n_valid

coverage_table$mean_value <- sapply(npp_vars, function(v) {
  if (v %in% names(data) && is.numeric(data[[v]])) {
    round(mean(data[[v]], na.rm = TRUE), 3)
  } else {
    NA
  }
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                    NPP Variable Coverage Summary                         ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

print(coverage_table, row.names = FALSE)

write.csv(coverage_table, file.path(plot_dir, "04_npp_coverage_table.csv"), 
          row.names = FALSE)
cat(sprintf("\n  ✓ Saved: %s\n", file.path(plot_dir, "04_npp_coverage_table.csv")))

# Equilibrium flag summary
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                    Equilibrium Flag Summary                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

flag_summary <- data.frame(
  flag = c("No observation date", "Sampled pre-2001", "Insufficient years (<3)",
           "High CV (>0.3)", "Significant trend", "Non-equilibrium (combined)"),
  n_true = c(
    sum(data$npp_flag_no_date, na.rm = TRUE),
    sum(data$npp_flag_pre_npp_era, na.rm = TRUE),
    sum(data$npp_flag_insufficient, na.rm = TRUE),
    sum(data$npp_flag_high_cv, na.rm = TRUE),
    sum(data$npp_flag_trend, na.rm = TRUE),
    sum(data$npp_non_equilibrium, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
flag_summary$pct <- round(100 * flag_summary$n_true / nrow(data), 1)

print(flag_summary, row.names = FALSE)

cat("\n  Trend direction breakdown:\n")
print(table(data$npp_trend_direction, useNA = "ifany"))

write.csv(flag_summary, file.path(plot_dir, "05_equilibrium_flags_summary.csv"),
          row.names = FALSE)
cat(sprintf("\n  ✓ Saved: %s\n", file.path(plot_dir, "05_equilibrium_flags_summary.csv")))

log_step("Step 4 complete")