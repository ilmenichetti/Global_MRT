################################################################################
# 12_calculate_MRT.R
#
# Calculate Mean Residence Time (MRT) of soil organic carbon
# MRT = SOC_stock / C_input
#
# References used for NPP allocation rates:
#
#   Mokany K, Raison RJ, Prokushkin AS (2006) Critical analysis of root:shoot
#     ratios in terrestrial biomes. Global Change Biology 12: 84-96.
#     DOI: 10.1111/j.1365-2486.2005.001043.x
#
#   Malhi Y et al. (2011) The allocation of ecosystem net primary productivity
#     in tropical forests. Phil Trans R Soc B 366: 3225-3245.
#     DOI: 10.1098/rstb.2011.0062
#
#   Litton CM, Giardina CP (2008) Below-ground carbon flux and partitioning:
#     global patterns and response to temperature. Functional Ecology 22: 941-954.
#     DOI: 10.1111/j.1365-2435.2008.01479.x

# Input:  ./Global_MRT_code/outputs/11_with_microbial.rds
# Output: ./Global_MRT_code/outputs/12_with_MRT.rds
#
# Author: Lorenzo
# Date: 2026-01-07
################################################################################

library(dplyr)
library(readr)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")

INPUT_FILE  <- file.path(OUTPUT_DIR, "11_with_microbial.rds")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "12_with_MRT.rds")

cat("═══════════════════════════════════════════════════════════════\n")
cat("  MEAN RESIDENCE TIME CALCULATION (Step 12)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
soil_data <- readRDS(INPUT_FILE)
cat("  Loaded", format(nrow(soil_data), big.mark = ","), "observations\n")
cat("  Variables:", ncol(soil_data), "\n\n")

# =============================================================================
# DEFINE ALLOCATION FRACTIONS BY LAND COVER
# =============================================================================
#
# Belowground NPP fraction (BNPP/NPP) represents the proportion of total NPP
# that enters the soil as carbon input (roots, root exudates, mycorrhizal transfer).
#
# Values derived from:
#   - Mokany et al. (2006): Root:shoot ratios across biomes
#   - Malhi et al. (2011): Tropical forest NPP allocation (~30-40% belowground)
#   - Litton & Giardina (2008): Temperature effects on belowground C flux
#
# Conversion from R:S to BNPP fraction:
#   If R:S = 0.25, then root fraction = 0.25/(1+0.25) = 0.20
#   Adding root turnover, exudates, mycorrhizal transfer increases this to ~0.30-0.40
#
# Land cover classes from ESA WorldCover 2021:
#   10 = Tree cover
#   20 = Shrubland
#   30 = Grassland
#   40 = Cropland
#   50 = Built-up
#   60 = Bare/sparse vegetation
#   70 = Snow and ice
#   80 = Permanent water bodies
#   90 = Herbaceous wetland
#   95 = Mangroves
#   100 = Moss and lichen
# =============================================================================

# Allocation lookup table
# BNPP_fraction = fraction of NPP entering soil as C input
allocation_fractions <- data.frame(
  landcover_code = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100),
  landcover_name = c("Tree cover", "Shrubland", "Grassland", "Cropland",
                     "Built-up", "Bare/sparse", "Snow/ice", "Water",
                     "Herbaceous wetland", "Mangroves", "Moss/lichen"),
  bnpp_fraction = c(
    0.35,   # Tree cover: R:S ~0.25, plus turnover/exudates (Mokany 2006, Malhi 2011)
    0.45,   # Shrubland: higher R:S than trees (Mokany 2006)
    0.55,   # Grassland: R:S ~0.9, high root turnover (Mokany 2006)
    0.30,   # Cropland: annual systems, lower allocation (Bolinder et al. 2007)
    0.20,   # Built-up: minimal vegetation, low estimate
    0.30,   # Bare/sparse: conservative estimate for sparse vegetation
    NA,     # Snow/ice: no NPP
    NA,     # Water: aquatic systems excluded
    0.50,   # Herbaceous wetland: high belowground allocation (Saunders et al. 2006)
    0.40,   # Mangroves: moderate allocation (Komiyama et al. 2008)
    0.60    # Moss/lichen: high belowground in tundra systems (Campioli et al. 2009)
  ),
  stringsAsFactors = FALSE
)

cat("Allocation fractions by land cover:\n")
print(allocation_fractions[!is.na(allocation_fractions$bnpp_fraction), ], row.names = FALSE)
cat("\n")

# =============================================================================
# CALCULATE SOC STOCK
# =============================================================================
#
# SOC stock (g C/m²) for 0-20 cm depth:
#   SOC_stock = OC × BD × depth × 10
#
# Where:
#   OC = organic carbon concentration (g/kg)
#   BD = bulk density (g/cm³)
#   depth = soil depth (cm) = 20
#   10 = unit conversion factor (g/kg × g/cm³ × cm × 10 = g/m²)
#
# Derivation:
#   OC (g C / kg soil) × BD (g soil / cm³) × depth (cm) × (1 kg/1000 g) × (10000 cm²/m²)
#   = OC × BD × depth × 10 g C/m²
# =============================================================================

cat("Calculating SOC stocks...\n")

# Use gap-filled bulk density (bd_filled) from step 10
soil_data <- soil_data %>%
  mutate(
    # SOC stock in g C/m² for 0-20 cm
    SOC_stock_g_m2 = organic_carbon * bd_filled * 20 * 10
  )

cat("  SOC stock (g C/m²) summary:\n")
print(summary(soil_data$SOC_stock_g_m2))
cat("\n")

# =============================================================================
# CALCULATE C INPUT FROM NPP
# =============================================================================
#
# NPP from MODIS MOD17A3 is in raw units (scale factor 0.1)
# True NPP (g C/m²/yr) = MODIS_NPP / 10
#
# C input to soil = NPP × BNPP_fraction
#
# This represents total belowground C flux including:
#   - Fine root production and turnover
#   - Coarse root production
#   - Root exudates
#   - Mycorrhizal C transfer
#   - Aboveground litter is NOT included (would add ~10-20% more)
#
# Note: This is a simplification. In reality, soil C input also includes
# aboveground litterfall, but belowground inputs dominate in many systems
# and are more directly linked to soil C stabilization.
# =============================================================================

cat("Calculating C inputs...\n")

# Check which land cover column exists
lc_col <- NULL
if ("lc_majority" %in% names(soil_data)) {
  lc_col <- "lc_majority"
} else if ("landcover" %in% names(soil_data)) {
  lc_col <- "landcover"
} else if ("lc_code" %in% names(soil_data)) {
  lc_col <- "lc_code"
} else {
  # Find any column with landcover in name
  
  lc_candidates <- grep("landcover|lc_", names(soil_data), value = TRUE, ignore.case = TRUE)
  if (length(lc_candidates) > 0) {
    lc_col <- lc_candidates[1]
  }
}

if (is.null(lc_col)) {
  cat("  WARNING: No land cover column found. Using global average allocation (0.40)\n")
  soil_data$bnpp_fraction <- 0.40
} else {
  cat("  Using land cover column:", lc_col, "\n")
  
  # Merge allocation fractions
  soil_data <- soil_data %>%
    left_join(
      allocation_fractions %>% select(landcover_code, bnpp_fraction),
      by = setNames("landcover_code", lc_col)
    )
  
  # Fill missing with global average
  n_missing <- sum(is.na(soil_data$bnpp_fraction))
  if (n_missing > 0) {
    cat("  Filling", n_missing, "missing allocation values with global average (0.40)\n")
    soil_data$bnpp_fraction[is.na(soil_data$bnpp_fraction)] <- 0.40
  }
}

# Calculate NPP in correct units and C input
soil_data <- soil_data %>%
  mutate(
    # NPP in g C/m²/yr (MODIS scale factor correction)
    NPP_g_C_m2_yr = npp_mean / 10,
    
    # Belowground C input to soil (g C/m²/yr)
    C_input_g_m2_yr = NPP_g_C_m2_yr * bnpp_fraction
  )

cat("\n  NPP (g C/m²/yr) summary:\n")
print(summary(soil_data$NPP_g_C_m2_yr))

cat("\n  C input (g C/m²/yr) summary:\n")
print(summary(soil_data$C_input_g_m2_yr))
cat("\n")

# =============================================================================
# CALCULATE MEAN RESIDENCE TIME
# =============================================================================
#
# MRT (years) = SOC_stock / C_input
#
# Interpretation:
#   MRT represents the average time a carbon atom spends in the soil pool
#   before being respired or leached. It is the inverse of the decomposition
#   rate constant (k) in first-order decay models: MRT = 1/k
#
# Assumptions:
#   - Steady state
#   - Homogeneous pool
#   - The calculated MRT is an "apparent" or "bulk" MRT 
#
# =============================================================================

cat("Calculating Mean Residence Time...\n")

soil_data <- soil_data %>%
  mutate(
    # MRT in years
    MRT_years = SOC_stock_g_m2 / C_input_g_m2_yr
  )

# =============================================================================
# QUALITY CONTROL
# =============================================================================
#
# Flag potentially unreliable MRT values:
#   - Very low MRT (<1 year): likely measurement error or non-steady state
#   - Very high MRT (>1000 years): possible in permafrost, but may indicate
#     data issues (very low NPP, very high SOC)
#   - Negative or infinite values: calculation errors
# =============================================================================

cat("\nQuality control...\n")

soil_data <- soil_data %>%
  mutate(
    # QC flags
    MRT_QC = case_when(
      is.na(MRT_years) ~ "missing_data",
      is.infinite(MRT_years) ~ "infinite",
      MRT_years < 0 ~ "negative",
      MRT_years < 1 ~ "very_low",
      MRT_years > 1000 ~ "very_high",
      TRUE ~ "valid"
    )
  )

cat("  MRT quality flags:\n")
print(table(soil_data$MRT_QC, useNA = "ifany"))

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  MRT SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Overall summary (valid values only)
valid_mrt <- soil_data$MRT_years[soil_data$MRT_QC == "valid"]

cat("MRT (years) - valid observations only:\n")
cat("  n =", length(valid_mrt), "\n")
cat("  Mean:", round(mean(valid_mrt, na.rm = TRUE), 1), "\n")
cat("  Median:", round(median(valid_mrt, na.rm = TRUE), 1), "\n")
cat("  SD:", round(sd(valid_mrt, na.rm = TRUE), 1), "\n")
cat("  Range:", round(min(valid_mrt, na.rm = TRUE), 1), "-", 
    round(max(valid_mrt, na.rm = TRUE), 1), "\n")
cat("  IQR:", round(quantile(valid_mrt, 0.25, na.rm = TRUE), 1), "-",
    round(quantile(valid_mrt, 0.75, na.rm = TRUE), 1), "\n")

# Summary by Köppen main group (if available)
if ("koppen_main_group" %in% names(soil_data)) {
  cat("\nMRT by Köppen main climate group (median, IQR):\n")
  
  mrt_by_koppen <- soil_data %>%
    filter(MRT_QC == "valid") %>%
    group_by(koppen_main_group) %>%
    summarise(
      n = n(),
      median_MRT = round(median(MRT_years, na.rm = TRUE), 1),
      q25 = round(quantile(MRT_years, 0.25, na.rm = TRUE), 1),
      q75 = round(quantile(MRT_years, 0.75, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  
  print(as.data.frame(mrt_by_koppen), row.names = FALSE)
}

# =============================================================================
# SAVE OUTPUT
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  SAVING OUTPUT\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

saveRDS(soil_data, OUTPUT_FILE)
cat("✓ Saved:", OUTPUT_FILE, "\n")
cat("  Observations:", format(nrow(soil_data), big.mark = ","), "\n")
cat("  Variables:", ncol(soil_data), "\n")

# List new variables added
new_vars <- c("SOC_stock_g_m2", "NPP_g_C_m2_yr", "bnpp_fraction", 
              "C_input_g_m2_yr", "MRT_years", "MRT_QC")
cat("\n  New variables added:\n")
for (v in new_vars) {
  if (v %in% names(soil_data)) {
    cat("    -", v, "\n")
  }
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  STEP 12 COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")