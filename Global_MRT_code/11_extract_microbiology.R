# =============================================================================
# Step 11: Extract microbial predictor variables
# 
# Variables: fungal_proportion (Yu 2022), AM/EcM root colonization (Barceló 2023),
#            AM/EcM richness & endemism (SPUN 2025)
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(readr)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
RASTER_DIR   <- file.path(PIPELINE_DIR, "spatialized_layers/microbial")
MICROBIAL_DATA_DIR <- file.path(PIPELINE_DIR, "microbial_data")

TARGET_RES <- 0.1  # Target resolution in degrees for global predictions

dir.create(RASTER_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  MICROBIAL VARIABLE EXTRACTION (Step 11)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Pipeline directory:", PIPELINE_DIR, "\n")
cat("Microbial data:    ", MICROBIAL_DATA_DIR, "\n")
cat("Output rasters:    ", RASTER_DIR, "\n")
cat("Target resolution: ", TARGET_RES, "degrees\n\n")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

find_input_file <- function() {
  # Find most recent pipeline output file
  candidates <- c(
    file.path(OUTPUT_DIR, "10_with_bd_filled.rds"),
    file.path(OUTPUT_DIR, "09_with_landcover.rds"),
    file.path(OUTPUT_DIR, "08_with_soilclass.rds"),
    file.path(OUTPUT_DIR, "07_with_koppen.rds")
  )
  
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  
  # Fallback: find latest numbered file
  available <- list.files(OUTPUT_DIR, pattern = "^[0-9]+.*\\.rds$", full.names = TRUE)
  if (length(available) > 0) return(sort(available, decreasing = TRUE)[1])
  
  stop("No input file found in ", OUTPUT_DIR)
}

aggregate_to_target <- function(r, target_res = TARGET_RES) {
  # Aggregate or resample raster to target resolution
  current_res <- res(r)[1]
  
  template <- rast(
    xmin = -180, xmax = 180,
    ymin = -90, ymax = 90,
    res = target_res,
    crs = "EPSG:4326"
  )
  
  if (current_res >= target_res) {
    return(resample(r, template, method = "bilinear"))
  } else {
    agg_factor <- max(1, round(target_res / current_res))
    r_agg <- aggregate(r, fact = agg_factor, fun = "mean", na.rm = TRUE)
    return(resample(r_agg, template, method = "bilinear"))
  }
}

extract_with_coverage <- function(r, pts, var_name) {
  # Extract values at points and report coverage
  vals <- terra::extract(r, pts, ID = FALSE)[[1]]
  coverage <- sum(!is.na(vals)) / length(vals) * 100
  cat(sprintf("    %-25s: %.1f%% coverage\n", var_name, coverage))
  return(vals)
}

# -----------------------------------------------------------------------------
# Section 1: Fungal proportion (Yu et al. 2022)
# -----------------------------------------------------------------------------

extract_fungal_proportion <- function(obs_points = NULL) {
  cat("─── 1. Fungal Proportion (Yu et al. 2022) ───\n\n")
  
  fb_file <- file.path(MICROBIAL_DATA_DIR, "fungal_biomass/pro_fun_mean_boot.tif")
  
  if (!file.exists(fb_file)) {
    cat("  ✗ File not found:", fb_file, "\n\n")
    return(list(points = NULL, raster = NULL))
  }
  
  cat("  Loading:", basename(fb_file), "\n")
  fb_rast <- rast(fb_file)
  names(fb_rast) <- "fungal_proportion"
  cat("  Original resolution:", round(res(fb_rast)[1], 5), "degrees\n")
  
  # Extract at observation points
  pts_result <- NULL
  if (!is.null(obs_points)) {
    cat("  Extracting at observation points...\n")
    pts_result <- data.frame(
      fungal_proportion = extract_with_coverage(fb_rast, obs_points, "fungal_proportion")
    )
  }
  
  # Aggregate to target resolution
  cat("  Aggregating to", TARGET_RES, "degree for predictions...\n")
  fb_agg <- aggregate_to_target(fb_rast)
  
  # Save global raster
  out_file <- file.path(RASTER_DIR, "fungal_proportion_0.1deg.tif")
  writeRaster(fb_agg, out_file, overwrite = TRUE)
  cat("  ✓ Saved:", out_file, "\n\n")
  
  return(list(
    points = pts_result,
    raster = fb_agg,
    source = "Yu et al. 2022 ESSD",
    doi = "10.5194/essd-14-4339-2022"
  ))
}

# -----------------------------------------------------------------------------
# Section 2: AM/EcM root colonization (Barceló et al. 2023)
# -----------------------------------------------------------------------------

extract_barcelo_mycorrhiza <- function(obs_points = NULL) {
  cat("─── 2. AM/EcM Root Colonization (Barceló et al. 2023) ───\n\n")
  
  am_dir <- file.path(MICROBIAL_DATA_DIR, "AM_data")
  
  if (!dir.exists(am_dir)) {
    cat("  ✗ Directory not found:", am_dir, "\n\n")
    return(list(points = NULL, raster = NULL))
  }
  
  files_to_extract <- list(
    AM_roots_colonized = "AM_roots_colonized.tif",
    EcM_roots_colonized = "EcM_roots_colonized.tif"
  )
  
  rast_list <- list()
  pts_list <- list()
  
  for (var_name in names(files_to_extract)) {
    fpath <- file.path(am_dir, files_to_extract[[var_name]])
    
    if (!file.exists(fpath)) {
      cat("  ✗", var_name, "not found\n")
      next
    }
    
    cat("  Loading:", files_to_extract[[var_name]], "\n")
    r <- rast(fpath)
    names(r) <- var_name
    
    # Extract at points
    if (!is.null(obs_points)) {
      pts_list[[var_name]] <- extract_with_coverage(r, obs_points, var_name)
    }
    
    # Aggregate for predictions
    cat("    Aggregating to", TARGET_RES, "degree...\n")
    r_agg <- aggregate_to_target(r)
    rast_list[[var_name]] <- r_agg
  }
  
  if (length(rast_list) == 0) {
    return(list(points = NULL, raster = NULL))
  }
  
  # Combine rasters
  barc_stack <- rast(rast_list)
  
  # Calculate EcM:AM ratio
  if (all(c("AM_roots_colonized", "EcM_roots_colonized") %in% names(barc_stack))) {
    cat("  Calculating EcM:AM root ratio...\n")
    ecm_am_ratio <- barc_stack[["EcM_roots_colonized"]] / 
      (barc_stack[["AM_roots_colonized"]] + 0.01)
    names(ecm_am_ratio) <- "EcM_AM_root_ratio"
    barc_stack <- c(barc_stack, ecm_am_ratio)
    
    if (!is.null(obs_points) && length(pts_list) == 2) {
      pts_list[["EcM_AM_root_ratio"]] <- pts_list[["EcM_roots_colonized"]] / 
        (pts_list[["AM_roots_colonized"]] + 0.01)
    }
  }
  
  # Save global raster stack
  out_file <- file.path(RASTER_DIR, "mycorrhiza_barcelo_0.1deg.tif")
  writeRaster(barc_stack, out_file, overwrite = TRUE)
  cat("  ✓ Saved:", out_file, "\n")
  cat("    Layers:", paste(names(barc_stack), collapse = ", "), "\n\n")
  
  pts_result <- if (length(pts_list) > 0) as.data.frame(pts_list) else NULL
  
  return(list(
    points = pts_result,
    raster = barc_stack,
    source = "Barceló et al. 2023 Scientific Data",
    doi = "10.1038/s41597-022-01913-2"
  ))
}

# -----------------------------------------------------------------------------
# Section 3: Mycorrhizal richness & endemism (SPUN 2025)
# -----------------------------------------------------------------------------

extract_spun_mycorrhiza <- function(obs_points = NULL) {
  cat("─── 3. Mycorrhizal Richness (SPUN 2025) ───\n\n")
  
  spun_dir <- file.path(MICROBIAL_DATA_DIR, "SPUN")
  
  if (!dir.exists(spun_dir)) {
    cat("  ✗ Directory not found:", spun_dir, "\n\n")
    return(list(points = NULL, raster = NULL))
  }
  
  # RWR = Rarity-Weighted Richness (endemism metric)
  files_to_extract <- list(
    AM_richness  = "AM_fungi/AM_Fungi_Richness_Predicted.tif",
    EcM_richness = "EcM_fungi/EcM_Fungi_Richness_Predicted.tif",
    AM_endemism  = "AM_fungi/AM_Fungi_RWR_Empirical_Predicted.tif",
    EcM_endemism = "EcM_fungi/EcM_Fungi_RWR_Empirical_Predicted.tif"
  )
  
  rast_list <- list()
  pts_list <- list()
  
  for (var_name in names(files_to_extract)) {
    fpath <- file.path(spun_dir, files_to_extract[[var_name]])
    
    if (!file.exists(fpath)) {
      cat("  ✗", var_name, "not found at", files_to_extract[[var_name]], "\n")
      next
    }
    
    cat("  Loading:", files_to_extract[[var_name]], "\n")
    r <- rast(fpath)
    names(r) <- var_name
    cat("    Original resolution:", round(res(r)[1], 5), "degrees\n")
    
    # Extract at points
    if (!is.null(obs_points)) {
      pts_list[[var_name]] <- extract_with_coverage(r, obs_points, var_name)
    }
    
    # Aggregate for predictions
    cat("    Aggregating to", TARGET_RES, "degree...\n")
    r_agg <- aggregate_to_target(r)
    rast_list[[var_name]] <- r_agg
  }
  
  if (length(rast_list) == 0) {
    return(list(points = NULL, raster = NULL))
  }
  
  # Combine rasters
  spun_stack <- rast(rast_list)
  
  # Calculate EcM:AM richness ratio
  if (all(c("AM_richness", "EcM_richness") %in% names(spun_stack))) {
    cat("  Calculating EcM:AM richness ratio...\n")
    ecm_am_ratio <- spun_stack[["EcM_richness"]] / 
      (spun_stack[["AM_richness"]] + 0.1)
    names(ecm_am_ratio) <- "EcM_AM_richness_ratio"
    spun_stack <- c(spun_stack, ecm_am_ratio)
    
    if (!is.null(obs_points) && all(c("AM_richness", "EcM_richness") %in% names(pts_list))) {
      pts_list[["EcM_AM_richness_ratio"]] <- pts_list[["EcM_richness"]] / 
        (pts_list[["AM_richness"]] + 0.1)
    }
  }
  
  # Save global raster stack
  out_file <- file.path(RASTER_DIR, "mycorrhiza_spun_0.1deg.tif")
  writeRaster(spun_stack, out_file, overwrite = TRUE)
  cat("  ✓ Saved:", out_file, "\n")
  cat("    Layers:", paste(names(spun_stack), collapse = ", "), "\n\n")
  
  pts_result <- if (length(pts_list) > 0) as.data.frame(pts_list) else NULL
  
  return(list(
    points = pts_result,
    raster = spun_stack,
    source = "Van Nuland et al. 2025 Nature (SPUN Underground Atlas)",
    doi = "10.1038/s41586-025-09277-4"
  ))
}

# =============================================================================
# Main execution
# =============================================================================

cat("Finding input file...\n")
obs_file <- find_input_file()
cat("Input file:", basename(obs_file), "\n\n")

# Load observation data
cat("Loading observation data...\n")
obs_df <- readRDS(obs_file)
cat("  Loaded", format(nrow(obs_df), big.mark = ","), "observations\n")
cat("  Variables:", ncol(obs_df), "\n\n")

# Find coordinate columns
coord_patterns <- list(
  lon = c("longitude_decimal_degrees", "longitude", "lon", "long", "x", "coords.x1"),
  lat = c("latitude_decimal_degrees", "latitude", "lat", "y", "coords.x2")
)

lon_col <- NULL
lat_col <- NULL

for (cand in coord_patterns$lon) {
  if (cand %in% names(obs_df)) {
    lon_col <- cand
    break
  }
}

for (cand in coord_patterns$lat) {
  if (cand %in% names(obs_df)) {
    lat_col <- cand
    break
  }
}

if (is.null(lon_col) || is.null(lat_col)) {
  cat("Available columns:\n")
  cat(paste(names(obs_df)[1:min(30, ncol(obs_df))], collapse = ", "), "\n")
  if (ncol(obs_df) > 30) cat("... and", ncol(obs_df) - 30, "more\n")
  stop("Cannot find coordinate columns")
}

cat("  Coordinates:", lon_col, "/", lat_col, "\n\n")

# Create spatial points
obs_points <- vect(obs_df, geom = c(lon_col, lat_col), crs = "EPSG:4326")

# -----------------------------------------------------------------------------
# Run extractions
# -----------------------------------------------------------------------------

results <- list()
results$fungal  <- extract_fungal_proportion(obs_points)
results$barcelo <- extract_barcelo_mycorrhiza(obs_points)
results$spun    <- extract_spun_mycorrhiza(obs_points)

# -----------------------------------------------------------------------------
# Combine point extractions
# -----------------------------------------------------------------------------

cat("═══════════════════════════════════════════════════════════════\n")
cat("  COMBINING RESULTS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

pt_dfs <- list()
for (nm in names(results)) {
  if (!is.null(results[[nm]]$points)) {
    pt_dfs[[nm]] <- results[[nm]]$points
  }
}

if (length(pt_dfs) > 0) {
  microbial_pts <- bind_cols(pt_dfs)
  cat("Extracted", ncol(microbial_pts), "microbial variables:\n")
  cat(" ", paste(names(microbial_pts), collapse = "\n  "), "\n\n")
  
  # Add to original data
  obs_with_microbial <- bind_cols(obs_df, microbial_pts)
  
  # Save
  out_file <- file.path(OUTPUT_DIR, "11_with_microbial.rds")
  saveRDS(obs_with_microbial, out_file)
  cat("✓ Saved point data:", out_file, "\n")
  cat("  ", format(nrow(obs_with_microbial), big.mark = ","), "observations,",
      ncol(obs_with_microbial), "variables\n\n")
  
  # Summary statistics
  cat("Variable summary:\n")
  for (v in names(microbial_pts)) {
    vals <- microbial_pts[[v]]
    n_valid <- sum(!is.na(vals))
    if (n_valid > 0) {
      cat(sprintf("  %-25s: n=%7d (%.1f%%), mean=%.3f, range=[%.3f, %.3f]\n",
                  v, n_valid, n_valid/length(vals)*100,
                  mean(vals, na.rm = TRUE),
                  min(vals, na.rm = TRUE),
                  max(vals, na.rm = TRUE)))
    } else {
      cat(sprintf("  %-25s: no valid data\n", v))
    }
  }
} else {
  cat("⚠ No point extractions succeeded!\n")
}

# -----------------------------------------------------------------------------
# Combine global rasters for predictions
# -----------------------------------------------------------------------------

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  GLOBAL PREDICTION RASTERS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

raster_files <- list.files(RASTER_DIR, pattern = "\\.tif$", full.names = TRUE)
raster_files <- raster_files[!grepl("combined", raster_files)]

cat("Individual rasters saved to:", RASTER_DIR, "\n")
for (f in raster_files) {
  r <- rast(f)
  cat(sprintf("  %s: %d layer(s)\n", basename(f), nlyr(r)))
}

# Create combined stack
if (length(raster_files) > 0) {
  cat("\nCreating combined raster stack...\n")
  
  all_rasts <- lapply(raster_files, rast)
  
  template <- rast(
    xmin = -180, xmax = 180,
    ymin = -90, ymax = 90,
    res = TARGET_RES,
    crs = "EPSG:4326"
  )
  
  aligned <- lapply(all_rasts, function(r) {
    resample(r, template, method = "bilinear")
  })
  
  combined <- rast(aligned)
  combined_file <- file.path(RASTER_DIR, "microbial_combined_0.1deg.tif")
  writeRaster(combined, combined_file, overwrite = TRUE)
  
  cat("✓ Saved combined stack:", combined_file, "\n")
  cat("  Layers (", nlyr(combined), "):\n", sep = "")
  for (nm in names(combined)) {
    cat("    -", nm, "\n")
  }
}

# -----------------------------------------------------------------------------
# Data sources and citations
# -----------------------------------------------------------------------------

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  DATA SOURCES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

citations <- data.frame(
  dataset = c("Fungal proportion", "AM/EcM root colonization", "AM/EcM richness & endemism"),
  source = c(
    "Yu et al. 2022 Earth System Science Data",
    "Barceló et al. 2023 Scientific Data",
    "Van Nuland et al. 2025 Nature (SPUN Underground Atlas)"
  ),
  variables = c(
    "fungal_proportion",
    "AM_roots_colonized, EcM_roots_colonized, EcM_AM_root_ratio",
    "AM_richness, EcM_richness, AM_endemism, EcM_endemism, EcM_AM_richness_ratio"
  ),
  doi = c(
    "10.5194/essd-14-4339-2022",
    "10.1038/s41597-022-01913-2",
    "10.1038/s41586-025-09277-4"
  ),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(citations)) {
  cat(citations$dataset[i], "\n")
  cat("  Source:", citations$source[i], "\n")
  cat("  DOI:", citations$doi[i], "\n")
  cat("  Variables:", citations$variables[i], "\n\n")
}

write_csv(citations, file.path(OUTPUT_DIR, "11_microbial_citations.csv"))
cat("✓ Saved citations:", file.path(OUTPUT_DIR, "11_microbial_citations.csv"), "\n")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  STEP 11 COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")