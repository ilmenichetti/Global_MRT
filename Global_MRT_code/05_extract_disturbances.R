# =============================================================================
# Step 5: Extract disturbance data (Hansen forest loss + MODIS burned area)
# Date-aware: flags disturbances that occurred BEFORE sampling
# =============================================================================

library(terra)
library(dplyr)

# Source config and utilities
source("./Global_MRT_code/00_config.R")
source("./Global_MRT_code/01_utils.R")

# --- Configuration ---
HANSEN_DIR <- "./Global_MRT_code//spatialized_layers/disturbances/Hansen_forest_loss"
BURNED_DIR <- "./Global_MRT_code//spatialized_layers/disturbances/MODIS_burned"
BURNED_YEARS <- 2001:2020

# Hansen loss year encoding: 0 = no loss, 1 = 2001, 2 = 2002, ..., 23 = 2023
HANSEN_YEAR_OFFSET <- 2000  # Add this to lossyear value to get actual year

# =============================================================================
# Helper Functions
# =============================================================================

#' Load Hansen forest loss raster
#' @param hansen_dir Directory containing Hansen_lossyear_global.tif
#' @return SpatRaster with loss year values
load_hansen_raster <- function(hansen_dir = HANSEN_DIR) {
  log_step("Loading Hansen forest loss raster")
  
  # Try different naming patterns
  patterns <- c(
    "Hansen_lossyear_global.tif",
    "Hansen_lossyear.tif",
    "lossyear_global.tif",
    "lossyear.tif"
  )
  
  for (pattern in patterns) {
    filepath <- file.path(hansen_dir, pattern)
    if (file.exists(filepath)) {
      cat(sprintf("  Found: %s\n", filepath))
      r <- rast(filepath)
      cat(sprintf("  Dimensions: %d x %d\n", nrow(r), ncol(r)))
      
      # Check value range
      vals <- values(r, na.rm = TRUE)
      cat(sprintf("  Value range: %d to %d\n", min(vals), max(vals)))
      cat(sprintf("  Loss year range: %d to %d (encoded as 1-23)\n", 
                  min(vals[vals > 0]), max(vals[vals > 0])))
      
      return(r)
    }
  }
  
  # List what's actually in the directory
  files <- list.files(hansen_dir, pattern = "\\.tif$", full.names = FALSE)
  stop("Hansen raster not found. Files in directory: ", paste(files, collapse = ", "))
}


#' Load MODIS burned area rasters into a stack
#' @param burned_dir Directory containing Burned_YYYY_global.tif files
#' @param years Vector of years to load
#' @return SpatRaster stack with one layer per year
load_burned_stack <- function(burned_dir = BURNED_DIR, years = BURNED_YEARS) {
  log_step("Loading MODIS burned area rasters")
  
  burned_files <- c()
  loaded_years <- c()
  
  for (year in years) {
    # Try different naming patterns
    patterns <- c(
      sprintf("Burned_%d_global.tif", year),
      sprintf("Burned_%d.tif", year),
      sprintf("burned_%d_global.tif", year),
      sprintf("burned_%d.tif", year)
    )
    
    found <- FALSE
    for (pattern in patterns) {
      filepath <- file.path(burned_dir, pattern)
      if (file.exists(filepath)) {
        burned_files <- c(burned_files, filepath)
        loaded_years <- c(loaded_years, year)
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      warning(sprintf("  Burned area file not found for year %d", year))
    }
  }
  
  if (length(burned_files) == 0) {
    # List what's actually in the directory
    files <- list.files(burned_dir, pattern = "\\.tif$", full.names = FALSE)
    stop("No burned area files found. Files in directory: ", paste(files, collapse = ", "))
  }
  
  cat(sprintf("  Found %d annual burned area files\n", length(burned_files)))
  
  # Load as stack
  burned_stack <- rast(burned_files)
  names(burned_stack) <- paste0("Burned_", loaded_years)
  
  cat(sprintf("  Stack dimensions: %d x %d x %d layers\n", 
              nrow(burned_stack), ncol(burned_stack), nlyr(burned_stack)))
  cat(sprintf("  Years covered: %d - %d\n", min(loaded_years), max(loaded_years)))
  
  return(burned_stack)
}


#' Extract Hansen forest loss at point locations
#' @param hansen_rast SpatRaster with loss year
#' @param points_data Data frame with coordinates
#' @return Vector of loss years (0 = no loss)
extract_hansen_at_points <- function(hansen_rast, points_data) {
  log_step("Extracting Hansen forest loss at points")
  
  pts <- terra::vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %d points...\n", nrow(points_data)))
  
  loss_vals <- terra::extract(hansen_rast, pts, ID = FALSE)[, 1]
  
  # Handle NoData
  loss_vals[loss_vals < 0] <- NA
  loss_vals[loss_vals > 30] <- NA  # Invalid values
  
  n_with_loss <- sum(loss_vals > 0, na.rm = TRUE)
  cat(sprintf("  ✓ Sites with forest loss: %d (%.1f%%)\n", 
              n_with_loss, 100 * n_with_loss / length(loss_vals)))
  
  return(loss_vals)
}


#' Extract MODIS burned area at point locations for all years
#' @param burned_stack SpatRaster stack of annual burned area
#' @param points_data Data frame with coordinates
#' @return Matrix of burned values (sites x years), 1 = burned
extract_burned_at_points <- function(burned_stack, points_data) {
  log_step("Extracting MODIS burned area at points")
  
  pts <- terra::vect(
    points_data,
    geom = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
    crs = "EPSG:4326"
  )
  
  cat(sprintf("  Extracting for %d points across %d years...\n", 
              nrow(points_data), nlyr(burned_stack)))
  
  burned_vals <- terra::extract(burned_stack, pts, ID = FALSE)
  burned_matrix <- as.matrix(burned_vals)
  
  # Handle NoData (negative values)
  burned_matrix[burned_matrix < 0] <- NA
  
  n_ever_burned <- sum(rowSums(burned_matrix, na.rm = TRUE) > 0)
  cat(sprintf("  ✓ Sites burned at least once: %d (%.1f%%)\n",
              n_ever_burned, 100 * n_ever_burned / nrow(burned_matrix)))
  
  return(burned_matrix)
}


#' Calculate disturbance diagnostics using observation year
#' @param hansen_loss Vector of Hansen loss years (encoded 1-23)
#' @param burned_matrix Matrix of burned values (sites x years)
#' @param obs_years Vector of observation years
#' @param burned_years Vector of years for burned columns
#' @return Data frame with disturbance diagnostics
calculate_disturbance_diagnostics <- function(hansen_loss, burned_matrix, 
                                              obs_years, burned_years = BURNED_YEARS) {
  log_step("Calculating disturbance diagnostics (date-aware)")
  
  n_sites <- length(hansen_loss)
  cat(sprintf("  Processing %d sites\n", n_sites))
  
  # --- Hansen forest loss ---
  # Convert encoded loss year (1-23) to actual year (2001-2023)
  hansen_loss_year <- rep(NA_integer_, n_sites)
  has_loss <- hansen_loss > 0 & !is.na(hansen_loss)
  hansen_loss_year[has_loss] <- hansen_loss[has_loss] + HANSEN_YEAR_OFFSET
  
  # Flag: forest loss occurred BEFORE observation
  hansen_loss_before_obs <- rep(FALSE, n_sites)
  for (i in 1:n_sites) {
    if (!is.na(hansen_loss_year[i]) && !is.na(obs_years[i])) {
      hansen_loss_before_obs[i] <- hansen_loss_year[i] < obs_years[i]
    }
  }
  
  # Flag: forest loss occurred (any time)
  hansen_any_loss <- hansen_loss > 0 & !is.na(hansen_loss)
  
  cat(sprintf("\n  === Hansen Forest Loss ===\n"))
  cat(sprintf("  Sites with any forest loss: %d (%.1f%%)\n",
              sum(hansen_any_loss), 100 * mean(hansen_any_loss)))
  cat(sprintf("  Sites with loss BEFORE sampling: %d (%.1f%%)\n",
              sum(hansen_loss_before_obs), 100 * mean(hansen_loss_before_obs)))
  
  # --- MODIS Burned Area ---
  # Total burn count (all years)
  burn_count_total <- rowSums(burned_matrix, na.rm = TRUE)
  
  # Burn count BEFORE observation year
  burn_count_before_obs <- rep(0L, n_sites)
  burned_before_obs <- rep(FALSE, n_sites)
  burn_years_before <- rep(NA_character_, n_sites)
  
  for (i in 1:n_sites) {
    obs_yr <- obs_years[i]
    
    if (is.na(obs_yr)) {
      # No obs date: use total burn count
      burn_count_before_obs[i] <- burn_count_total[i]
      burned_before_obs[i] <- burn_count_total[i] > 0
      next
    }
    
    # Get years before observation
    years_before_idx <- which(burned_years < obs_yr)
    
    if (length(years_before_idx) == 0) {
      burn_count_before_obs[i] <- 0L
      burned_before_obs[i] <- FALSE
      next
    }
    
    # Count burns before obs
    burns_before <- burned_matrix[i, years_before_idx]
    burn_count_before_obs[i] <- sum(burns_before, na.rm = TRUE)
    burned_before_obs[i] <- burn_count_before_obs[i] > 0
    
    # Record which years burned
    if (burned_before_obs[i]) {
      burn_yrs <- burned_years[years_before_idx][burns_before == 1 & !is.na(burns_before)]
      if (length(burn_yrs) > 0) {
        burn_years_before[i] <- paste(burn_yrs, collapse = ",")
      }
    }
  }
  
  cat(sprintf("\n  === MODIS Burned Area ===\n"))
  cat(sprintf("  Sites burned at least once (any time): %d (%.1f%%)\n",
              sum(burn_count_total > 0), 100 * mean(burn_count_total > 0)))
  cat(sprintf("  Sites burned BEFORE sampling: %d (%.1f%%)\n",
              sum(burned_before_obs), 100 * mean(burned_before_obs)))
  cat(sprintf("  Mean burn count (before sampling): %.2f\n",
              mean(burn_count_before_obs)))
  
  # --- Combined disturbance flag ---
  # Non-equilibrium if EITHER forest loss OR fire before sampling
  disturbed_before_obs <- hansen_loss_before_obs | burned_before_obs
  
  cat(sprintf("\n  === Combined Disturbance ===\n"))
  cat(sprintf("  Sites disturbed BEFORE sampling: %d (%.1f%%)\n",
              sum(disturbed_before_obs), 100 * mean(disturbed_before_obs)))
  cat(sprintf("    - Forest loss only: %d\n",
              sum(hansen_loss_before_obs & !burned_before_obs)))
  cat(sprintf("    - Fire only: %d\n",
              sum(!hansen_loss_before_obs & burned_before_obs)))
  cat(sprintf("    - Both: %d\n",
              sum(hansen_loss_before_obs & burned_before_obs)))
  
  # --- Create output data frame ---
  diagnostics <- data.frame(
    # Hansen forest loss
    hansen_loss_year = hansen_loss_year,
    hansen_any_loss = hansen_any_loss,
    hansen_loss_before_obs = hansen_loss_before_obs,
    
    # MODIS burned area
    burn_count_total = burn_count_total,
    burn_count_before_obs = burn_count_before_obs,
    burned_before_obs = burned_before_obs,
    burn_years_before = burn_years_before,
    
    # Combined
    disturbed_before_obs = disturbed_before_obs
  )
  
  return(diagnostics)
}


# =============================================================================
# Main Execution
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║         Step 5: Disturbance Extraction (Forest + Fire)        ║\n")
cat("║         Date-aware: flags disturbances BEFORE sampling        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# --- Step 1: Load input data ---
log_step("Loading input data")

input_file <- file.path(PATHS$output_dir, "04_with_productivity.rds")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, 
       "\nRun 04_extract_productivity.R first")
}

data <- readRDS(input_file)
cat(sprintf("  Loaded %d observations\n", nrow(data)))

# Check for obs_year column
if (!"obs_year" %in% names(data)) {
  stop("obs_year column not found. Run 04_extract_productivity.R first.")
}


# --- Step 2: Get unique site locations ---
log_step("Identifying unique locations")

sites <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees, obs_year)

cat(sprintf("  %d unique site locations\n", nrow(sites)))


# --- Step 3: Load Hansen forest loss ---
hansen_rast <- load_hansen_raster(HANSEN_DIR)

# Check if raster is valid
if (inherits(try(nrow(hansen_rast), silent = TRUE), "try-error")) {
  cat("  Reloading Hansen raster...\n")
  hansen_rast <- load_hansen_raster(HANSEN_DIR)
}


# --- Step 4: Load MODIS burned area stack ---
burned_stack <- load_burned_stack(BURNED_DIR, BURNED_YEARS)

# Check if stack is valid
if (inherits(try(nlyr(burned_stack), silent = TRUE), "try-error")) {
  cat("  Reloading burned area stack...\n")
  burned_stack <- load_burned_stack(BURNED_DIR, BURNED_YEARS)
}


# --- Step 5: Extract Hansen at points ---
hansen_loss <- load_or_run(
  "05_hansen_extracted.rds",
  function() {
    extract_hansen_at_points(hansen_rast, sites)
  }
)


# --- Step 6: Extract burned area at points ---
burned_matrix <- load_or_run(
  "05_burned_matrix.rds",
  function() {
    extract_burned_at_points(burned_stack, sites)
  }
)


# --- Step 7: Calculate diagnostics ---
disturbance_diagnostics <- load_or_run(
  "05_disturbance_diagnostics.rds",
  function() {
    calculate_disturbance_diagnostics(hansen_loss, burned_matrix, 
                                      sites$obs_year, BURNED_YEARS)
  }
)

# Add site IDs for joining
disturbance_diagnostics$dsiteid <- sites$dsiteid


# --- Step 8: Merge with main dataset ---
log_step("Merging disturbance data with main dataset")

data <- data %>%
  left_join(disturbance_diagnostics, by = "dsiteid")

cat(sprintf("  ✓ Added %d disturbance variables\n", ncol(disturbance_diagnostics) - 1))


# --- Step 9: Save output ---
log_step("Saving output")

output_file <- file.path(PATHS$output_dir, "05_with_disturbance.rds")
saveRDS(data, output_file)
cat(sprintf("  ✓ Saved: %s\n", output_file))


# =============================================================================
# Diagnostics
# =============================================================================

log_step("Generating diagnostics")

diag_dir <- "./plots/05_extract_disturbance"
if (!dir.exists(diag_dir)) {
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
}

# --- 1. Data coverage table ---
log_step("Creating data coverage table")

dist_vars <- c("hansen_loss_year", "hansen_any_loss", "hansen_loss_before_obs",
               "burn_count_total", "burn_count_before_obs", "burned_before_obs",
               "disturbed_before_obs")

coverage_df <- data.frame(
  variable = dist_vars,
  n_total = nrow(data),
  n_valid = sapply(dist_vars, function(v) sum(!is.na(data[[v]]))),
  pct_valid = sapply(dist_vars, function(v) round(100 * mean(!is.na(data[[v]])), 2))
)

# For logical columns, add TRUE counts
coverage_df$n_true <- sapply(dist_vars, function(v) {
  if (is.logical(data[[v]])) sum(data[[v]], na.rm = TRUE) else NA
})
coverage_df$pct_true <- sapply(dist_vars, function(v) {
  if (is.logical(data[[v]])) round(100 * mean(data[[v]], na.rm = TRUE), 2) else NA
})

write.csv(coverage_df, file.path(diag_dir, "disturbance_data_coverage.csv"), row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "disturbance_data_coverage.csv")))


# --- 2. Map: Forest loss before sampling ---
log_step("Creating forest loss map")

sites_for_plot <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees,
           hansen_loss_before_obs, hansen_loss_year)

png(file.path(diag_dir, "map_hansen_loss.png"), 
    width = 2400, height = 1200, res = 150)

par(mar = c(2, 2, 3, 1))

plot(c(-180, 180), c(-90, 90), type = "n", 
     xlab = "", ylab = "", asp = 1,
     main = "Hansen Forest Loss (before sampling)")

if (requireNamespace("maps", quietly = TRUE)) {
  maps::map("world", add = TRUE, col = "gray90", fill = TRUE, border = "gray70")
}

# Sites without loss (gray)
no_loss <- sites_for_plot[!sites_for_plot$hansen_loss_before_obs | 
                            is.na(sites_for_plot$hansen_loss_before_obs), ]
points(no_loss$longitude_decimal_degrees, no_loss$latitude_decimal_degrees,
       pch = 16, cex = 0.3, col = rgb(0.5, 0.5, 0.5, 0.3))

# Sites with loss (red)
with_loss <- sites_for_plot[sites_for_plot$hansen_loss_before_obs == TRUE & 
                              !is.na(sites_for_plot$hansen_loss_before_obs), ]
points(with_loss$longitude_decimal_degrees, with_loss$latitude_decimal_degrees,
       pch = 16, cex = 0.5, col = rgb(0.8, 0, 0, 0.7))

legend("bottomleft",
       legend = c(sprintf("Forest loss before sampling (n=%d)", nrow(with_loss)),
                  sprintf("No loss (n=%d)", nrow(no_loss))),
       pch = 16, col = c("firebrick", "gray50"),
       bg = "white", cex = 0.9)

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "map_hansen_loss.png")))


# --- 3. Map: Fire before sampling ---
log_step("Creating fire map")

sites_for_plot <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees,
           burned_before_obs, burn_count_before_obs)

png(file.path(diag_dir, "map_burned_area.png"), 
    width = 2400, height = 1200, res = 150)

par(mar = c(2, 2, 3, 1))

plot(c(-180, 180), c(-90, 90), type = "n", 
     xlab = "", ylab = "", asp = 1,
     main = "MODIS Burned Area (before sampling)")

if (requireNamespace("maps", quietly = TRUE)) {
  maps::map("world", add = TRUE, col = "gray90", fill = TRUE, border = "gray70")
}

# Sites without fire (gray)
no_fire <- sites_for_plot[!sites_for_plot$burned_before_obs | 
                            is.na(sites_for_plot$burned_before_obs), ]
points(no_fire$longitude_decimal_degrees, no_fire$latitude_decimal_degrees,
       pch = 16, cex = 0.3, col = rgb(0.5, 0.5, 0.5, 0.3))

# Sites with fire - color by burn count
with_fire <- sites_for_plot[sites_for_plot$burned_before_obs == TRUE & 
                              !is.na(sites_for_plot$burned_before_obs), ]

if (nrow(with_fire) > 0) {
  fire_colors <- colorRampPalette(c("orange", "red", "darkred"))(10)
  burn_idx <- pmin(with_fire$burn_count_before_obs, 10)
  points(with_fire$longitude_decimal_degrees, with_fire$latitude_decimal_degrees,
         pch = 16, cex = 0.5, col = fire_colors[burn_idx])
}

legend("bottomleft",
       legend = c(sprintf("Burned before sampling (n=%d)", nrow(with_fire)),
                  sprintf("Not burned (n=%d)", nrow(no_fire))),
       pch = 16, col = c("red", "gray50"),
       bg = "white", cex = 0.9)

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "map_burned_area.png")))


# --- 4. Map: Combined disturbance ---
log_step("Creating combined disturbance map")

sites_for_plot <- data %>%
  distinct(dsiteid, longitude_decimal_degrees, latitude_decimal_degrees,
           hansen_loss_before_obs, burned_before_obs, disturbed_before_obs)

png(file.path(diag_dir, "map_combined_disturbance.png"), 
    width = 2400, height = 1200, res = 150)

par(mar = c(2, 2, 3, 1))

plot(c(-180, 180), c(-90, 90), type = "n", 
     xlab = "", ylab = "", asp = 1,
     main = "Combined Disturbance (before sampling)")

if (requireNamespace("maps", quietly = TRUE)) {
  maps::map("world", add = TRUE, col = "gray90", fill = TRUE, border = "gray70")
}

# Undisturbed (gray)
undisturbed <- sites_for_plot[!sites_for_plot$disturbed_before_obs | 
                                is.na(sites_for_plot$disturbed_before_obs), ]
points(undisturbed$longitude_decimal_degrees, undisturbed$latitude_decimal_degrees,
       pch = 16, cex = 0.3, col = rgb(0.5, 0.5, 0.5, 0.3))

# Forest loss only (red)
forest_only <- sites_for_plot[sites_for_plot$hansen_loss_before_obs & 
                                !sites_for_plot$burned_before_obs &
                                !is.na(sites_for_plot$hansen_loss_before_obs), ]
if (nrow(forest_only) > 0) {
  points(forest_only$longitude_decimal_degrees, forest_only$latitude_decimal_degrees,
         pch = 16, cex = 0.5, col = "firebrick")
}

# Fire only (orange)
fire_only <- sites_for_plot[!sites_for_plot$hansen_loss_before_obs & 
                              sites_for_plot$burned_before_obs &
                              !is.na(sites_for_plot$burned_before_obs), ]
if (nrow(fire_only) > 0) {
  points(fire_only$longitude_decimal_degrees, fire_only$latitude_decimal_degrees,
         pch = 16, cex = 0.5, col = "orange")
}

# Both (purple)
both <- sites_for_plot[sites_for_plot$hansen_loss_before_obs & 
                         sites_for_plot$burned_before_obs &
                         !is.na(sites_for_plot$hansen_loss_before_obs) &
                         !is.na(sites_for_plot$burned_before_obs), ]
if (nrow(both) > 0) {
  points(both$longitude_decimal_degrees, both$latitude_decimal_degrees,
         pch = 16, cex = 0.6, col = "purple")
}

legend("bottomleft",
       legend = c(sprintf("Undisturbed (n=%d)", nrow(undisturbed)),
                  sprintf("Forest loss only (n=%d)", nrow(forest_only)),
                  sprintf("Fire only (n=%d)", nrow(fire_only)),
                  sprintf("Both (n=%d)", nrow(both))),
       pch = 16, col = c("gray50", "firebrick", "orange", "purple"),
       bg = "white", cex = 0.9)

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "map_combined_disturbance.png")))


# --- 5. Histogram: Forest loss years ---
png(file.path(diag_dir, "hist_hansen_loss_years.png"), 
    width = 1200, height = 800, res = 150)

par(mar = c(5, 4, 3, 2))

loss_years <- data$hansen_loss_year[!is.na(data$hansen_loss_year)]
if (length(loss_years) > 0) {
  hist(loss_years, 
       breaks = seq(2000.5, 2024.5, by = 1),
       main = "Distribution of Forest Loss Years",
       xlab = "Year of forest loss", ylab = "Number of sites",
       col = "firebrick", border = "white")
} else {
  plot(1, type = "n", main = "No forest loss data available")
}

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "hist_hansen_loss_years.png")))


# --- 6. Histogram: Burn counts ---
png(file.path(diag_dir, "hist_burn_counts.png"), 
    width = 1200, height = 800, res = 150)

par(mar = c(5, 4, 3, 2))

burn_counts <- data$burn_count_before_obs[data$burn_count_before_obs > 0]
if (length(burn_counts) > 0) {
  hist(burn_counts, 
       breaks = seq(0.5, max(burn_counts) + 0.5, by = 1),
       main = "Distribution of Burn Counts (before sampling)",
       xlab = "Number of burn events", ylab = "Number of sites",
       col = "orange", border = "white")
} else {
  plot(1, type = "n", main = "No fire data available")
}

dev.off()
cat(sprintf("  ✓ Saved: %s\n", file.path(diag_dir, "hist_burn_counts.png")))

# --- Summary ---
log_step("Disturbance extraction complete")