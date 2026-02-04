# Load required packages
library(ggplot2)
library(maps)
library(dplyr)

# =============================================================================
# Step 2.5: Standardize all soil data to 0-20cm depth interval
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║     Step 2.5: Standardize to 0-20cm Depth Interval            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

most_common <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

depth_weight_horizons <- function(horizons, target_top = 0, target_bot = 20) {
  # Depth-weight soil properties to a target interval (0-20cm)
  # For each variable, calculate weighted average based on horizon overlap
  
  agg_vars <- c(
    "organic_carbon",
    "total_nitrogen_ncs",
    "ph_h2o",
    "clay_total",
    "silt_total",
    "sand_total",
    "bulk_density_oven_dry"
  )
  
  result <- list()
  
  for (var in agg_vars) {
    # Skip if variable doesn't exist
    if (!(var %in% names(horizons))) {
      result[[var]] <- NA_real_
      next
    }
    
    values <- horizons[[var]]
    if (all(is.na(values))) {
      result[[var]] <- NA_real_
      next
    }
    
    # Find horizons that overlap with target interval
    hzn_top <- horizons$hzn_top
    hzn_bot <- horizons$hzn_bot
    overlaps <- (hzn_bot > target_top) & (hzn_top < target_bot)
    
    if (sum(overlaps & !is.na(values)) == 0) {
      result[[var]] <- NA_real_
      next
    }
    
    # Calculate weighted average
    weighted_sum <- 0
    total_weight <- 0
    
    for (i in which(overlaps & !is.na(values))) {
      overlap_top <- max(hzn_top[i], target_top)
      overlap_bot <- min(hzn_bot[i], target_bot)
      overlap_thickness <- overlap_bot - overlap_top
      
      weighted_sum <- weighted_sum + (values[i] * overlap_thickness)
      total_weight <- total_weight + overlap_thickness
    }
    
    result[[var]] <- if (total_weight > 0) weighted_sum / total_weight else NA_real_
  }
  
  return(as.data.frame(result))
}

# -----------------------------------------------------------------------------
# Main processing function
# -----------------------------------------------------------------------------

standardize_to_0_20cm <- function() {
  log_step("Standardizing soil data to 0-20cm depth interval")
  
  if (!exists("data")) {
    stop("Data not loaded. Run 02_load_filter.R first.")
  }
  
  # Report input data
  cat(sprintf("  Input: %s records from %s sites\n",
              format(nrow(data), big.mark = ","),
              format(n_distinct(data$dsiteid), big.mark = ",")))
  
  sites <- unique(data$dsiteid[!is.na(data$dsiteid)])
  cat(sprintf("  Processing %s unique sites...\n", format(length(sites), big.mark = ",")))
  
  n_na <- sum(is.na(data$dsiteid))
  if (n_na > 0) {
    cat(sprintf("  Note: Excluding %d records with missing dsiteid\n", n_na))
  }
  
  # Process each site
  results <- list()
  
  for (i in seq_along(sites)) {
    site_id <- sites[i]
    site_data <- data %>% filter(dsiteid == site_id)
    
    # Extract site-level metadata
    result_row <- data.frame(
      dsiteid = site_id,
      longitude_decimal_degrees = site_data$longitude_decimal_degrees[1],
      latitude_decimal_degrees = site_data$latitude_decimal_degrees[1],
      site_obsdate = most_common(site_data$site_obsdate),
      source_db = most_common(site_data$source_db),
      confidence_degree = as.character(tryCatch({
        conf_vals <- suppressWarnings(as.numeric(site_data$confidence_degree))
        if (all(is.na(conf_vals))) {
          most_common(site_data$confidence_degree)
        } else {
          round(mean(conf_vals, na.rm = TRUE), 2)
        }
      }, error = function(e) {
        most_common(site_data$confidence_degree)
      }))
    )
    
    # Depth-weight soil variables to 0-20cm
    weighted <- depth_weight_horizons(site_data, target_top = 0, target_bot = 20)
    result_row <- cbind(result_row, weighted)
    results[[i]] <- result_row
    
    # Progress updates
    if (i %% 10000 == 0 || i == length(sites)) {
      cat(sprintf("  Progress: %d / %d sites\n", i, length(sites)))
    }
  }
  
  standardized <- bind_rows(results)
  cat(sprintf("  ✓ Standardized to %s sites\n", format(nrow(standardized), big.mark = ",")))
  
  # Report data availability
  cat("\n  Data availability in 0-20cm interval:\n")
  agg_vars <- c(
    "organic_carbon", "total_nitrogen_ncs", "ph_h2o",
    "clay_total", "silt_total", "sand_total", "bulk_density_oven_dry"
  )
  
  for (var in agg_vars) {
    n_valid <- sum(!is.na(standardized[[var]]))
    pct <- 100 * n_valid / nrow(standardized)
    if (n_valid > 0) {
      mean_val <- mean(standardized[[var]], na.rm = TRUE)
      cat(sprintf("    %s: %d sites (%.1f%%), mean = %.2f\n",
                  var, n_valid, pct, mean_val))
    } else {
      cat(sprintf("    %s: NO DATA\n", var))
    }
  }
  
  return(standardized)
}

# -----------------------------------------------------------------------------
# Run standardization
# -----------------------------------------------------------------------------

data_0_20cm <- load_or_run(
  "02_5_standardized_0_20cm.rds",
  standardize_to_0_20cm
)

log_step("Standardization complete")
cat(sprintf("\nOutput saved: ./outputs/02_5_standardized_0_20cm.rds\n"))
cat(sprintf("Ready for climate extraction and spatialization.\n\n"))

# -----------------------------------------------------------------------------
# Inspect output
# -----------------------------------------------------------------------------

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                    Output Inspection                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("[1] Dimensions:\n")
cat(sprintf("  Rows: %s\n", format(nrow(data_0_20cm), big.mark = ",")))
cat(sprintf("  Columns: %d\n\n", ncol(data_0_20cm)))

cat("[2] Column names:\n")
print(names(data_0_20cm))

cat("\n[3] Data types:\n")
str(data_0_20cm, max.level = 1)

cat("\n[4] First 5 rows:\n")
print(head(data_0_20cm, 5))

cat("\n[5] Metadata summary:\n")
cat(sprintf("  Unique source_db: %d\n", n_distinct(data_0_20cm$source_db)))
cat(sprintf("  Source databases: %s\n", paste(unique(data_0_20cm$source_db), collapse = ", ")))
cat(sprintf("  Year range: %s\n", paste(range(data_0_20cm$site_obsdate, na.rm = TRUE), collapse = " - ")))
cat(sprintf("  Mean confidence: %.2f\n\n", mean(as.numeric(data_0_20cm$confidence_degree), na.rm = TRUE)))

cat("[6] Soil variable completeness (0-20cm):\n")
soil_vars <- c("organic_carbon", "total_nitrogen_ncs", "ph_h2o", "clay_total", 
               "silt_total", "sand_total", "bulk_density_oven_dry")
for (var in soil_vars) {
  n_valid <- sum(!is.na(data_0_20cm[[var]]))
  pct <- 100 * n_valid / nrow(data_0_20cm)
  if (n_valid > 0) {
    mean_val <- mean(data_0_20cm[[var]], na.rm = TRUE)
    sd_val <- sd(data_0_20cm[[var]], na.rm = TRUE)
    cat(sprintf("  %-25s: %d sites (%.1f%%), mean = %.2f, sd = %.2f\n",
                var, n_valid, pct, mean_val, sd_val))
  } else {
    cat(sprintf("  %-25s: NO DATA\n", var))
  }
}

cat("\n")

# -----------------------------------------------------------------------------
# Generate diagnostic plots
# -----------------------------------------------------------------------------

dir.create("./Global_MRT_code/plots/02_5_standardization", recursive = TRUE, showWarnings = FALSE)

# Histograms
pdf("./Global_MRT_code/plots/02_5_standardization/histograms.pdf", width = 10, height = 8)
par(mfrow = c(3, 3))
for (var in soil_vars) {
  vals <- data_0_20cm[[var]][!is.na(data_0_20cm[[var]])]
  if (length(vals) > 0) {
    hist(vals, main = var, xlab = var, col = "steelblue", breaks = 50)
  }
}
dev.off()
cat("Saved: ./Global_MRT_code/plots/02_5_standardization/histograms.pdf\n")

# Global maps
world <- map_data("world")

pdf("./Global_MRT_code/plots/02_5_standardization/global_maps.pdf", width = 12, height = 8)
for (var in soil_vars) {
  df <- data_0_20cm %>%
    filter(!is.na(.data[[var]])) %>%
    rename(lon = longitude_decimal_degrees, lat = latitude_decimal_degrees)
  
  if (nrow(df) > 0) {
    p <- ggplot() +
      geom_polygon(data = world, aes(x = long, y = lat, group = group),
                   fill = "grey90", color = "grey70", linewidth = 0.1) +
      geom_point(data = df, aes(x = lon, y = lat, color = .data[[var]]),
                 size = 0.3, alpha = 0.5) +
      scale_color_viridis_c(option = "plasma") +
      coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
      labs(title = paste0(var, " (0-20cm)"),
           subtitle = sprintf("n = %s sites", format(nrow(df), big.mark = ",")),
           x = "Longitude", y = "Latitude") +
      theme_minimal() +
      theme(legend.position = "bottom")
    print(p)
  }
}
dev.off()
cat("Saved: ./Global_MRT_code/plots/02_5_standardization/global_maps.pdf\n")

# Site density map
pdf("./Global_MRT_code/plots/02_5_standardization/site_density.pdf", width = 12, height = 8)
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey70", linewidth = 0.1) +
  geom_point(data = data_0_20cm, 
             aes(x = longitude_decimal_degrees, y = latitude_decimal_degrees),
             size = 0.1, alpha = 0.2, color = "darkblue") +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
  labs(title = "All sites (0-20cm standardized)",
       subtitle = sprintf("n = %s sites", format(nrow(data_0_20cm), big.mark = ",")),
       x = "Longitude", y = "Latitude") +
  theme_minimal()
dev.off()
cat("Saved: ./Global_MRT_code/plots/02_5_standardization/site_density.pdf\n")