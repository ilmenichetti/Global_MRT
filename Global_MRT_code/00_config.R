# =============================================================================
# Configuration: paths, packages, and global settings
# =============================================================================

# --- Paths (all relative to Global_MRT_code/) ---
PATHS <- list(
  # Input data (one level up)
  data_dir           = "../Datasets/Global_SOC",
  
  # Downloaded/cached data (in Global_MRT_code/)
  temp_dir           = "./Global_MRT_code/temp",
  era5_dir           = "./Global_MRT_code/era5_data",
  soilgrids_dir      = "./Global_MRT_code/soilgrids",
  koppen_dir         = "./Global_MRT_code/koppen_data",
  hwsd_dir           = "./Global_MRT_code/hwsd_data",
  
  # Output data (in Global_MRT_code/)
  output_dir         = "./Global_MRT_code/outputs",
  raster_output_dir  = "./Global_MRT_code/spatialized_layers",  # organized by data source
  
  # Specific file references (within downloaded directories)
  koppen_file        = "./koppen_data/koppen_geiger_0p1.tif",
  hwsd_raster        = "./hwsd_data/hwsd.bil",
  hwsd_lookup        = "./hwsd_data/hwsd_lookup.csv"
)

# --- API credentials (set via environment variables for security) ---
CDS_KEY  <- Sys.getenv("CDS_API_KEY", "fba6dc92-aa6e-49cb-98db-d8060a941fd1")
CDS_USER <- "cds"

# --- Processing parameters ---
PARAMS <- list(
  # Mineral soil threshold
  oc_organic_threshold = 120,  # g/kg - soils above this are organic
  
  # Depth standardization
  target_depth = 20,           # cm - standardize OC to this depth
  min_depth    = 10,           # cm - minimum sample depth to include
  max_depth    = 30,           # cm - maximum sample depth to include
  
  # Climate data period
  climate_years = 1999:2020,
  
  # Extraction batch sizes
  batch_size_era5 = 5000,
  batch_size_soilgrids = 5000,
  
  # Random forest settings
  rf_ntree = 500,
  rf_seed  = 42,
  
  # SPATIAL RASTERIZATION PARAMETERS
  spatial_resolution = 0.5,    # degrees (lon/lat)
  spatial_crs = "EPSG:4326",   # WGS84
  global_extent = terra::ext(-180, 180, -90, 90),  # global bounding box
  raster_datatype = "FLT4S"    # float32, good for soil data
)

# --- Required packages ---
required_packages <- c(
  "dplyr", "tidyr", "purrr",    # data manipulation
  "terra", "sf",                 # spatial
  "ecmwfr",                      # ERA5 download
  "elevatr",                     # elevation data
  "geodata",                     # land cover
  "randomForest", "caret",       # modeling
  "ggplot2"                      # visualization
)

# --- Load packages ---
load_packages <- function() {
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' required. Install with: install.packages('%s')", pkg, pkg))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  message("All packages loaded.")
}

# --- Create directories ---
setup_directories <- function() {
  # Point workflow directories
  dirs <- c(
    PATHS$output_dir, 
    PATHS$temp_dir, 
    PATHS$era5_dir, 
    PATHS$soilgrids_dir,
    PATHS$koppen_dir,
    PATHS$hwsd_dir
  )
  
  # Spatial raster directories (organized by data source)
  raster_subdirs <- file.path(
    PATHS$raster_output_dir,
    c("climate", "soilgrids", "topography", "landcover", "soilclass", "productivity")
  )
  dirs <- c(dirs, raster_subdirs)
  
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("Created: ", d)
    } else {
      message("Exists: ", d)
    }
  }
  
  cat("\nDirectory structure ready.\n")
}



message("Configuration loaded. Paths and parameters set.")

