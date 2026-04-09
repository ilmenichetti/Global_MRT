################################################################################
# Step 17: Global Visualization Plots
# Creates publication-quality visualizations of the SOC dataset
#
# Updates:
#   - MRT maps with log-scale palette (better use of color range)
#   - log(MRT) histogram colored by mean latitude per bin
################################################################################

library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(scales)
library(dplyr)  
library(paletteer)

# Set up paths
output_dir <- "./Global_MRT_code/plots/step_17_visualizations"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load processed SOC data (depth-standardized to 0-20cm)
cat("Loading depth-standardized SOC data...\n")
soc_data <- readRDS("./Global_MRT_code/outputs/12b_model_ready.rds")

names(soc_data)

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

################################################################################
# PLOT 1: Global distribution of SOC sampling points
################################################################################
cat("Creating global SOC sampling points map...\n")
# Prepare data for plotting
plot_data <- data.frame(
  lon = soc_data$longitude_decimal_degrees,
  lat = soc_data$latitude_decimal_degrees,
  soc = soc_data$SOC_stock_g_m2
)
plot_data = na.omit(plot_data)
plot_data_sf <- st_as_sf(plot_data, coords = c("lon", "lat"), crs = 4326)
cat(sprintf("Plotting %s sampling points...\n", 
            format(nrow(plot_data), big.mark = ",")))
# Create the map
p1 <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70", size = 0.2) +
  geom_sf(data = plot_data_sf, color = "#2E86AB", size = 0.8, alpha = 0.6) +
  # No coord_sf - uses plain lon/lat like p2
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_text(size = 8)
  ) +
  labs(
    title = "Global Distribution of Soil Organic Carbon Sampling Points",
    subtitle = sprintf("n = %s observations from %s unique sites", 
                       format(nrow(plot_data), big.mark = ","),
                       format(length(unique(paste(plot_data$lon, plot_data$lat))), 
                              big.mark = ","))
  )

print(p1)

# Save the plot
ggsave(
  filename = file.path(output_dir, "global_SOC_sampling_points.png"),
  plot = p1,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("✓ Global sampling points map saved\n")

################################################################################
# PLOT 2: Point density heatmap
################################################################################

cat("Creating point density heatmap...\n")

# Extract coordinates from sf geometry
coords <- st_coordinates(plot_data_sf)
plot_data_hex <- data.frame(
  soc = plot_data_sf$soc,
  lon = coords[, 1],
  lat = coords[, 2]
)

p2 <- ggplot(plot_data_hex, aes(x = lon, y = lat)) +
  stat_binhex(aes(fill = stat(count)), bins = 80, color = NA, alpha = 0.8) +
  geom_sf(data = world, fill = NA, color = "grey40", size = 0.3, 
          inherit.aes = FALSE) +
  scale_fill_viridis_c(
    name = "Observation\nDensity",
    option = "turbo",
    trans = "log10"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 9)
  ) +
  labs(
    title = "Global Distribution of SOC Observations",
    subtitle = "Sampling density by geographic region (log-scale)",
    caption = sprintf("Total observations: %s", format(nrow(plot_data_hex), big.mark = ","))
  )

print(p2)


ggsave(
  filename = file.path(output_dir, "global_SOC_density_heatmap.png"),
  plot = p2,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("✓ Density heatmap saved\n")



################################################################################
# PLOT 4: Sampling coverage by latitude
################################################################################

cat("Creating latitudinal sampling coverage plot...\n")

p4 <- ggplot(plot_data, aes(x = lat)) +
  geom_density(fill = "#2E86AB", alpha = 0.6, color = "#2E86AB", linewidth = 0.8) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  ) +
  labs(
    title = "Latitudinal Distribution of SOC Observations",
    x = "Latitude (°)",
    y = "Density",
    subtitle = "Kernel density estimate"
  )  # Optional: rotates to match your original orientation

ggsave(
  filename = file.path(output_dir, "latitudinal_sampling_coverage.png"),
  plot = p4,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)
cat("✓ Latitudinal coverage plot saved\n")


################################################################################
# Summary statistics
################################################################################

cat("\n=== Global Sampling Summary ===\n")
cat(sprintf("Total observations: %s\n", 
            format(nrow(plot_data), big.mark = ",")))
cat(sprintf("Unique locations: %s\n", 
            format(length(unique(paste(plot_data$lon, plot_data$lat))), 
                   big.mark = ",")))
cat(sprintf("Latitude range: %.2f° to %.2f°\n", 
            min(plot_data$lat), max(plot_data$lat)))
cat(sprintf("Longitude range: %.2f° to %.2f°\n", 
            min(plot_data$lon), max(plot_data$lon)))
cat(sprintf("SOC stock range: %.2f to %.2f kg/m² (0-20cm)\n",
            min(plot_data$soc, na.rm = TRUE),
            max(plot_data$soc, na.rm = TRUE)))
cat(sprintf("Mean SOC stock: %.2f kg/m²\n",
            mean(plot_data$soc, na.rm = TRUE)))
cat(sprintf("Median SOC stock: %.2f kg/m²\n",
            median(plot_data$soc, na.rm = TRUE)))

################################################################################
# PLOT 5: Temporal distribution of sampling dates
################################################################################

cat("Creating temporal sampling distribution plot...\n")

if ("Year_of_Observation" %in% names(soc_data)) {
  temporal_data <- data.frame(
    year = soc_data$Year_of_Observation,
    database = soc_data$Database
  )
  temporal_data <- temporal_data[!is.na(temporal_data$year), ]
  
  # Overall temporal distribution
  p5a <- ggplot(temporal_data, aes(x = year)) +
    geom_histogram(binwidth = 5, fill = "#2E86AB", alpha = 0.7, color = "white") +
    scale_x_continuous(breaks = seq(1900, 2030, by = 20)) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(
      title = "Temporal Distribution of SOC Observations",
      subtitle = "5-year bins",
      x = "Year of Observation",
      y = "Number of Observations"
    )
  
  ggsave(
    filename = file.path(output_dir, "temporal_distribution_overall.png"),
    plot = p5a,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  # Temporal distribution by database
  p5b <- ggplot(temporal_data, aes(x = year, fill = database)) +
    geom_histogram(binwidth = 5, alpha = 0.7, color = "white", position = "stack") +
    scale_x_continuous(breaks = seq(1900, 2030, by = 20)) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    scale_fill_viridis_d(option = "turbo", name = "Database") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "Temporal Distribution by Data Source",
      subtitle = "5-year bins, stacked by database",
      x = "Year of Observation",
      y = "Number of Observations"
    )
  
  ggsave(
    filename = file.path(output_dir, "temporal_distribution_by_database.png"),
    plot = p5b,
    width = 12,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✓ Temporal distribution plots saved\n")
} else {
  cat("! Year_of_Observation not found - skipping temporal plots\n")
}

################################################################################
# PLOT 6: SOC vs Climate Space Coverage
################################################################################

cat("Creating SOC vs climate space coverage plot...\n")

# Load climate data (from step 03)
if (file.exists("./outputs/step_03_climate/soc_with_climate.rds")) {
  climate_data <- readRDS("./outputs/step_03_climate/soc_with_climate.rds")
  
  # Merge with SOC stock data
  climate_soc <- merge(
    climate_data[, c("Sample_ID", "MAT", "MAP")],
    soc_data[, c("Sample_ID", "SOC_stock_0_20cm_kgm2")],
    by = "Sample_ID"
  )
  
  climate_soc <- climate_soc[complete.cases(climate_soc), ]
  
  # Limit to reasonable ranges for visualization
  climate_soc <- climate_soc[
    climate_soc$MAT > -30 & climate_soc$MAT < 35 &
      climate_soc$MAP > 0 & climate_soc$MAP < 4000,
  ]
  
  # SOC vs MAT and MAP
  p6a <- ggplot(climate_soc, aes(x = MAT, y = SOC_stock_0_20cm_kgm2)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(
      option = "inferno",
      trans = "log10",
      name = "Count",
      labels = label_number()
    ) +
    scale_y_log10(labels = label_number()) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "SOC Stock vs Mean Annual Temperature",
      x = "Mean Annual Temperature (°C)",
      y = "SOC Stock (kg/m², log scale)",
      subtitle = sprintf("n = %s", format(nrow(climate_soc), big.mark = ","))
    )
  
  ggsave(
    filename = file.path(output_dir, "SOC_vs_MAT.png"),
    plot = p6a,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  p6b <- ggplot(climate_soc, aes(x = MAP, y = SOC_stock_0_20cm_kgm2)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(
      option = "inferno",
      trans = "log10",
      name = "Count",
      labels = label_number()
    ) +
    scale_x_log10(labels = label_number()) +
    scale_y_log10(labels = label_number()) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "SOC Stock vs Mean Annual Precipitation",
      x = "Mean Annual Precipitation (mm, log scale)",
      y = "SOC Stock (kg/m², log scale)",
      subtitle = sprintf("n = %s", format(nrow(climate_soc), big.mark = ","))
    )
  
  ggsave(
    filename = file.path(output_dir, "SOC_vs_MAP.png"),
    plot = p6b,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  # Climate space coverage (MAT vs MAP)
  p6c <- ggplot(climate_soc, aes(x = MAP, y = MAT)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(
      option = "mako",
      trans = "log10",
      name = "Observation\nCount",
      labels = label_number()
    ) +
    scale_x_log10(labels = label_number()) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "Climate Space Coverage",
      x = "Mean Annual Precipitation (mm, log scale)",
      y = "Mean Annual Temperature (°C)",
      subtitle = "Sampling distribution across climate gradients"
    )
  
  ggsave(
    filename = file.path(output_dir, "climate_space_coverage.png"),
    plot = p6c,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  cat("✓ Climate space coverage plots saved\n")
} else {
  cat("! Climate data not found - skipping climate space plots\n")
}

################################################################################
# PLOT 7: Data source contributions
################################################################################

cat("Creating data source contribution plots...\n")

if ("Database" %in% names(soc_data)) {
  # Count by database
  database_counts <- as.data.frame(table(soc_data$Database))
  names(database_counts) <- c("Database", "Count")
  database_counts <- database_counts[order(database_counts$Count, decreasing = TRUE), ]
  database_counts$Percentage <- 100 * database_counts$Count / sum(database_counts$Count)
  
  # Bar plot of contributions
  p7a <- ggplot(database_counts, aes(x = reorder(Database, Count), y = Count)) +
    geom_col(fill = "#2E86AB", alpha = 0.7) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)", 
                                  format(Count, big.mark = ","),
                                  Percentage)),
              hjust = -0.1, size = 3) +
    scale_y_continuous(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0, 0.15))
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.y = element_text(size = 9)
    ) +
    labs(
      title = "Data Source Contributions",
      x = NULL,
      y = "Number of Observations",
      subtitle = sprintf("Total: %s observations from %s databases",
                         format(sum(database_counts$Count), big.mark = ","),
                         nrow(database_counts))
    )
  
  ggsave(
    filename = file.path(output_dir, "database_contributions_bar.png"),
    plot = p7a,
    width = 10,
    height = max(6, nrow(database_counts) * 0.4),
    dpi = 300,
    bg = "white"
  )
  
  # Pie chart for top contributors
  database_counts$Label <- ifelse(
    database_counts$Percentage > 2,
    sprintf("%s\n%.1f%%", database_counts$Database, database_counts$Percentage),
    ""
  )
  
  p7b <- ggplot(database_counts, aes(x = "", y = Count, fill = Database)) +
    geom_col(color = "white", size = 0.5) +
    geom_text(aes(label = Label), 
              position = position_stack(vjust = 0.5),
              size = 3) +
    coord_polar(theta = "y") +
    scale_fill_viridis_d(option = "turbo", name = "Database") +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right"
    ) +
    labs(
      title = "Proportional Data Source Contributions",
      subtitle = "Only sources >2% labeled"
    )
  
  ggsave(
    filename = file.path(output_dir, "database_contributions_pie.png"),
    plot = p7b,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  # Geographic distribution by database
  if (all(c("Longitude", "Latitude") %in% names(soc_data))) {
    # Select top 6 databases for clearer visualization
    top_databases <- database_counts$Database[1:min(6, nrow(database_counts))]
    
    plot_data_db <- soc_data[soc_data$Database %in% top_databases, ]
    
    p7c <- ggplot() +
      geom_sf(data = world, fill = "grey95", color = "grey70", size = 0.2) +
      geom_point(
        data = plot_data_db,
        aes(x = Longitude, y = Latitude, color = Database),
        size = 0.3,
        alpha = 0.4
      ) +
      scale_color_viridis_d(option = "turbo", name = "Database") +
      coord_sf(crs = "+proj=robin") +
      theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "grey80", size = 0.2),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_text(size = 8),
        legend.position = "right"
      ) +
      labs(
        title = "Geographic Distribution by Data Source",
        subtitle = sprintf("Top %d databases shown", length(top_databases))
      )
    
    ggsave(
      filename = file.path(output_dir, "database_geographic_distribution.png"),
      plot = p7c,
      width = 14,
      height = 8,
      dpi = 300,
      bg = "white"
    )
  }
  
  # Print summary table
  cat("\n=== Database Contributions ===\n")
  print(database_counts, row.names = FALSE)
  
  cat("\n✓ Data source contribution plots saved\n")
} else {
  cat("! Database field not found - skipping data source plots\n")
}


################################################################################
# PLOT 8: MRT prediction maps with log-scale palette
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  MRT PREDICTION MAPS (log-scale palette)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

PRED_DIR <- "./Global_MRT_code/outputs/MRT_predictions"
VIS_RES  <- 0.5  # aggregation for faster plotting

# Load borders
sf_use_s2(FALSE)
world_borders <- ne_countries(scale = 50, returnclass = "sf")
world_borders_geom <- st_geometry(world_borders)

# Helper: aggregate raster for visualization
prep_raster <- function(r, agg_res = VIS_RES) {
  agg_factor <- round(agg_res / res(r)[1])
  if (agg_factor > 1) {
    r_agg <- aggregate(r, fact = agg_factor, fun = "mean", na.rm = TRUE)
  } else {
    r_agg <- r
  }
  return(r_agg)
}

# Find prediction TIFs
tif_files <- list.files(PRED_DIR, pattern = "^MRT_M[0-9].*\\.tif$", full.names = TRUE)
cat("Found", length(tif_files), "prediction maps\n")

if (length(tif_files) > 0) {
  
  # Determine global range from the M7 (full model) raster
  m7_file <- grep("M7_full", tif_files, value = TRUE)
  if (length(m7_file) > 0) {
    m7_vals <- values(rast(m7_file[1]), na.rm = TRUE)
    q01 <- quantile(m7_vals, 0.01)
    q99 <- quantile(m7_vals, 0.99)
    cat(sprintf("  M7 value range: %.1f - %.1f years\n", min(m7_vals), max(m7_vals)))
    cat(sprintf("  Using p1-p99 for color scale: %.1f - %.1f years\n", q01, q99))
    rm(m7_vals)
  } else {
    q01 <- 5; q99 <- 80
  }
  
  # Log-scale breaks: evenly spaced in log, labeled in years
  log_min <- log(max(q01, 1))
  log_max <- log(q99)
  n_cols  <- 100
  log_breaks <- seq(log_min, log_max, length.out = n_cols + 1)
  mrt_breaks <- exp(log_breaks)
  
  # Color palette
  cols_mrt <- hcl.colors(n_cols, "Spectral", rev = TRUE)
  
  # Legend tick marks at nice round numbers
  legend_vals <- c(5, 10, 15, 20, 30, 50, 75, 100, 150, 200)
  legend_vals <- legend_vals[legend_vals >= exp(log_min) & legend_vals <= exp(log_max)]
  legend_at   <- log(legend_vals)
  
  # --- Individual maps (log-scale) ---
  cat("\nCreating individual MRT maps (log-scale)...\n")
  
  for (tif_file in tif_files) {
    
    model_name <- gsub("\\.tif$", "", basename(tif_file))
    cat("  Processing", model_name, "...\n")
    
    r <- rast(tif_file)
    r_agg <- prep_raster(r)
    
    # Clamp and log-transform
    r_clamped <- clamp(r_agg, lower = exp(log_min), upper = exp(log_max))
    r_log <- log(r_clamped)
    
    out_png <- file.path(output_dir, paste0(model_name, "_map_log.png"))
    
    png(out_png, width = 2000, height = 1000, res = 150)
    
    plot(r_log,
         main = paste0(model_name, " — Mean Residence Time"),
         col = cols_mrt,
         range = c(log_min, log_max),
         mar = c(2, 2, 2, 5),
         axes = TRUE,
         legend = FALSE)
    
    plot(world_borders_geom, add = TRUE, col = NA, border = "grey30", lwd = 0.3)
    
    # Custom legend with years labels
    legend_y <- seq(-50, 60, length.out = n_cols)
    for (i in seq_len(n_cols)) {
      rect(175, legend_y[i], 180, legend_y[min(i + 1, n_cols)],
           col = cols_mrt[i], border = NA, xpd = TRUE)
    }
    # Tick labels
    for (j in seq_along(legend_vals)) {
      y_pos <- -50 + (legend_at[j] - log_min) / (log_max - log_min) * 110
      text(184, y_pos, labels = legend_vals[j], cex = 0.65, xpd = TRUE, adj = 0)
      segments(180, y_pos, 181.5, y_pos, xpd = TRUE, lwd = 0.5)
    }
    text(182, 68, "MRT (yr)", cex = 0.7, xpd = TRUE, font = 2)
    
    dev.off()
    cat("    ✓", basename(out_png), "\n")
  }
  
  # --- Multi-panel comparison (log-scale) ---
  cat("\n  Creating multi-panel comparison (log-scale)...\n")
  
  png(file.path(output_dir, "MRT_all_models_comparison_log.png"),
      width = 2800, height = 1400, res = 150)
  
  n_models <- length(tif_files)
  n_cols_panel <- 4
  n_rows_panel <- ceiling(n_models / n_cols_panel)
  
  par(mfrow = c(n_rows_panel, n_cols_panel),
      mar = c(0.2, 0.2, 1.5, 0.2), oma = c(0, 0, 2, 0))
  
  for (tif_file in tif_files) {
    model_label <- gsub("^MRT_", "", gsub("\\.tif$", "", basename(tif_file)))
    
    r <- rast(tif_file)
    r_agg <- prep_raster(r)
    r_clamped <- clamp(r_agg, lower = exp(log_min), upper = exp(log_max))
    r_log <- log(r_clamped)
    
    plot(r_log,
         main = model_label,
         col = cols_mrt,
         range = c(log_min, log_max),
         axes = FALSE,
         legend = FALSE,
         mar = c(0.1, 0.1, 1, 0.1))
    
    plot(world_borders_geom, add = TRUE, col = NA, border = "grey40", lwd = 0.3)
  }
  
  mtext("MRT Predictions — Log-scale Comparison", outer = TRUE,
        cex = 1.2, font = 2, line = 0.3)
  
  dev.off()
  cat("  ✓ MRT_all_models_comparison_log.png\n")
  
  cat("\n✓ All MRT maps saved\n")
  
} else {
  cat("! No MRT prediction files found in ", PRED_DIR, "\n")
}


################################################################################
# PLOT 9: log(MRT) histogram colored by mean latitude per bin
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  log(MRT) HISTOGRAM COLORED BY LATITUDE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Use M7 full model raster
m7_file <- list.files(PRED_DIR, pattern = "MRT_M7_full\\.tif$", full.names = TRUE)

if (length(m7_file) > 0) {
  
  cat("Loading M7 raster and extracting values + coordinates...\n")
  r <- rast(m7_file[1])
  
  # Extract all values and their coordinates
  vals_all <- values(r, na.rm = FALSE)[, 1]
  coords_all <- xyFromCell(r, 1:ncell(r))
  
  # Keep non-NA
  keep <- !is.na(vals_all) & vals_all > 0
  df_hist <- data.frame(
    mrt     = vals_all[keep],
    log_mrt = log(vals_all[keep]),
    lat     = coords_all[keep, 2]
  )
  cat(sprintf("  Non-NA pixels: %s\n", format(nrow(df_hist), big.mark = ",")))
  
  # Compute bin statistics
  n_bins <- 120
  df_hist$bin <- cut(df_hist$log_mrt, breaks = n_bins)
  
  bin_stats <- df_hist %>%
    group_by(bin) %>%
    summarise(
      mean_lat = mean(abs(lat)),
      count    = n(),
      .groups  = "drop"
    ) %>%
    filter(!is.na(bin))
  
  # Extract bin midpoints from factor labels
  bin_stats$mid <- sapply(as.character(bin_stats$bin), function(b) {
    nums <- as.numeric(regmatches(b, gregexpr("-?[0-9.]+", b))[[1]])
    mean(nums)
  })
  
  # MRT in years for the secondary x-axis labels
  mrt_ticks <- c(5, 10, 15, 20, 30, 50, 75, 100, 150)
  log_ticks <- log(mrt_ticks)
  
  # Build the plot
  p_lat_hist <- ggplot(bin_stats, aes(x = mid, y = count, fill = mean_lat)) +
    geom_col(width = diff(range(bin_stats$mid)) / n_bins * 0.95, color = NA) +
    scale_fill_gradientn(
      colours = rev(paletteer::paletteer_c("ggthemes::Temperature Diverging", n = 256)),
      name    = "Mean distance\nfrom equator (°)",
      limits  = c(0, 65),
      oob     = scales::squish
    ) +
    scale_x_continuous(
      name   = "MRT (years)",
      breaks = log_ticks,
      labels = mrt_ticks,
      sec.axis = sec_axis(~ ., name = "log(MRT)", labels = function(x) round(x, 1))
    ) +
    scale_y_continuous(
      labels = label_comma(),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      plot.title          = element_text(size = 14, face = "bold"),
      plot.subtitle       = element_text(size = 10, colour = "grey40"),
      legend.position     = "right",
      axis.title          = element_text(size = 11),
      axis.text           = element_text(size = 9)
    ) +
    labs(
      title    = "Global Distribution of Predicted SOC Mean Residence Time",
      subtitle = sprintf(
        "Bin colour = mean latitude of pixels in each bin · n = %s pixels",
        format(nrow(df_hist), big.mark = ",")),
      y = "Number of pixels"
    )
  
  print(p_lat_hist)
  
  ggsave(
    filename = file.path(output_dir, "MRT_histogram_by_latitude.png"),
    plot = p_lat_hist,
    width = 12,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✓ MRT_histogram_by_latitude.png saved\n")
  
  # --- Also save a summary table ---
  cat("\n  MRT distribution summary:\n")
  cat(sprintf("    Median: %.1f years\n", median(df_hist$mrt)))
  cat(sprintf("    Mean:   %.1f years\n", mean(df_hist$mrt)))
  cat(sprintf("    p5-p95: %.1f - %.1f years\n",
              quantile(df_hist$mrt, 0.05),
              quantile(df_hist$mrt, 0.95)))
  cat(sprintf("    p1-p99: %.1f - %.1f years\n",
              quantile(df_hist$mrt, 0.01),
              quantile(df_hist$mrt, 0.99)))
  
} else {
  cat("! M7 prediction raster not found — skipping latitude histogram\n")
}


################################################################################
cat("\n✓ Step 17 complete - all visualizations saved to:\n")
cat(sprintf("  %s\n", output_dir))
################################################################################
