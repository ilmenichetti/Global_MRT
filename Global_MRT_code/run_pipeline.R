# =============================================================================
# Main Pipeline: Soil Carbon Data Processing
# =============================================================================

# --- Setup ---
rm(list = ls())

# TODO: filter for duplicates based on a similarity score
# TODO: Define a process based temperature and moisture model to substitute the RM model
# TODO: Write an alternative data analysis path, point-based


# --- Data preparation ---
source("./Global_MRT_code/00_config.R") 
source("./Global_MRT_code/01_utils.R")
source("./Global_MRT_code/02_load_filter.R")
source("./Global_MRT_code/02_5_standardize_to_0_20cm.R")
source("./Global_MRT_code/03_extract_climate.R")
source("./Global_MRT_code/04_extract_productivity.R")
source("./Global_MRT_code/05_extract_disturbances.R")
source("./Global_MRT_code/06_extract_soilgrids.R") 
source("./Global_MRT_code/07_extract_covariates.R")
source("./Global_MRT_code/08_extract_soil_class.R")
source("./Global_MRT_code/09_extract_landcover.R")
source("./Global_MRT_code/10_gap_fill_bd.R")
source("./Global_MRT_code/11_extract_microbiology.R")
source("./Global_MRT_code/12_MRT_calculation.R")
source("./Global_MRT_code/12_1_add_missing_variables.R") #

# --- Data analysis, raster-based ---
# This pathway works with inferred geographical rasters
source("./Global_MRT_code/13_model_fitting.R")
source("./Global_MRT_code/14_model_extrapolation.R")
source("./Global_MRT_code/15_map_visualization.R")

# --- Data analysis, point-based pathway ---
# This pathway works with as much measured variables as possible, on point clouds

git remote add origin https://github.com/ilmenichetti/Soil_biodiversity_task_force/tree/main/Global_MRT.git