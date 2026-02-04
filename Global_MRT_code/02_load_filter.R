# =============================================================================
# Step 2: Load and filter soil data to mineral soils
# =============================================================================

#' Load raw soil chemistry data and filter to mineral soils
#' The limit for organic soils comes from the configuration file, object PARAMS$oc_organic_threshold
#' @return Filtered dataframe with mineral soils only
load_and_filter_mineral_soils <- function() {
  log_step("Loading soil chemistry data")
  
  # Load primary dataset
  data_path <- file.path(PATHS$data_dir, "sol_chem.pnts_horizons.rds")
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s", data_path))
  }
  
  sol_chem <- readRDS(data_path)
  cat(sprintf("  Loaded: %s rows, %s columns\n", 
              format(nrow(sol_chem), big.mark = ","),
              ncol(sol_chem)))
  
  # Filter to mineral soils (OC < threshold)
  mineral_soils <- sol_chem %>%
    filter(organic_carbon < PARAMS$oc_organic_threshold | is.na(organic_carbon))
  
  n_removed <- nrow(sol_chem) - nrow(mineral_soils)
  cat(sprintf("  Removed %s organic soil records (OC >= %d g/kg)\n",
              format(n_removed, big.mark = ","), 
              PARAMS$oc_organic_threshold))
  cat(sprintf("  Retained %s mineral soil records\n",
              format(nrow(mineral_soils), big.mark = ",")))
  
  # Get unique sites for reference
  n_sites <- n_distinct(mineral_soils$dsiteid)
  cat(sprintf("  Unique sites: %s\n", format(n_sites, big.mark = ","))
  )
  
  return(mineral_soils)
}



# --- Main execution ---

data <- load_or_run(
  "02_mineral_soils_filtered.rds",
  load_and_filter_mineral_soils
)
