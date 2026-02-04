################################################################################
# Step 17: Global Visualization Plots
# Creates publication-quality visualizations of the SOC dataset
################################################################################

library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(scales)
library(dplyr)  

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

cat("\n✓ Step 17 complete - all visualizations saved to:\n")
cat(sprintf("  %s\n", output_dir))