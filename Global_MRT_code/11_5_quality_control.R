################################################################################
# 11_5_quality_control.R
#
# Quality control filters applied AFTER all covariate extractions (steps 03-11)
# and BEFORE MRT calculation (step 12)
#
# Filters:
#   1. OC == 0 → set to NA (missing coded as zero; 83k profiles)
#   2. Remove exact coordinate + OC duplicates (cross-database overlap; ~102k)
#   3. Implausible BD (< 0.1 or > 2.5) → set to NA (re-gap-filled in step 12)
#   4. Remove profiles with no OC (cannot compute MTT)
#   5. Flag BD = 1.600 plateau (possible default value in source DB)
#   6. Flag low-precision coordinates (integer-only)
#
# Input:  ./Global_MRT_code/outputs/11_with_microbial.rds
# Output: ./Global_MRT_code/outputs/11_5_quality_controlled.rds
#         ./Global_MRT_code/plots/11_5_qc/  (diagnostics)
#
# NOTE: Step 12_MRT_calculation.R must be updated to read 11_5 instead of 11
################################################################################

library(dplyr)
library(ggplot2)

# =============================================================================
# CONFIGURATION
# =============================================================================

PIPELINE_DIR <- "./Global_MRT_code"
OUTPUT_DIR   <- file.path(PIPELINE_DIR, "outputs")
PLOT_DIR     <- file.path(PIPELINE_DIR, "plots", "11_5_qc")
INPUT_FILE   <- file.path(OUTPUT_DIR, "11_with_microbial.rds")
OUTPUT_FILE  <- file.path(OUTPUT_DIR, "11_5_quality_controlled.rds")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  DATA QUALITY CONTROL (Step 11.5)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

d <- readRDS(INPUT_FILE)
n_start <- nrow(d)
cat(sprintf("Input: %s profiles, %d variables\n\n",
            format(n_start, big.mark = ","), ncol(d)))

# Removal log: track each filter
removal_log <- data.frame(
  step = character(),
  description = character(),
  n_removed = integer(),
  n_remaining = integer(),
  stringsAsFactors = FALSE
)

log_filter <- function(log, step, desc, n_before, n_after) {
  rbind(log, data.frame(
    step = step,
    description = desc,
    n_removed = n_before - n_after,
    n_remaining = n_after,
    stringsAsFactors = FALSE
  ))
}

removal_log <- log_filter(removal_log, "0_start", "Input data", n_start, n_start)


# =============================================================================
# PRE-FILTER DIAGNOSTICS
# =============================================================================

cat("─── Pre-filter diagnostics ───\n\n")

cat(sprintf("  Total profiles: %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("  OC valid: %s  |  OC == 0: %s  |  OC NA: %s\n",
            format(sum(!is.na(d$organic_carbon) & d$organic_carbon > 0, na.rm = TRUE), big.mark = ","),
            format(sum(d$organic_carbon == 0, na.rm = TRUE), big.mark = ","),
            format(sum(is.na(d$organic_carbon)), big.mark = ",")))

if ("bulk_density_oven_dry" %in% names(d)) {
  bd <- d$bulk_density_oven_dry
  cat(sprintf("  BD valid: %s  |  BD NA: %s  |  BD==1.6: %s\n",
              format(sum(!is.na(bd)), big.mark = ","),
              format(sum(is.na(bd)), big.mark = ","),
              format(sum(bd == 1.6, na.rm = TRUE), big.mark = ",")))
}
cat("\n")


# =============================================================================
# FILTER 1: OC == 0 → NA (missing coded as zero)
# =============================================================================
#
# Rationale: True OC = 0 in the 0-20cm mineral soil layer is physically
# implausible. In WoSIS, zero values typically indicate "not measured" in
# legacy datasets that used 0 as a missing-data code.
# We set these to NA rather than removing the profile, because the profile
# may still carry useful covariate information if BD is available from
# gap-filling. However, profiles with OC = NA will be excluded from MTT
# calculation in step 12.
# =============================================================================

cat("─── Filter 1: OC == 0 → NA ───\n")

n_oc_zero <- sum(d$organic_carbon == 0, na.rm = TRUE)
cat(sprintf("  Profiles with OC == 0: %s (%.1f%% of total)\n",
            format(n_oc_zero, big.mark = ","),
            100 * n_oc_zero / nrow(d)))

# Diagnostic: which databases contribute OC == 0?
if ("source_db" %in% names(d)) {
  oc_zero_by_db <- d %>%
    filter(!is.na(organic_carbon)) %>%
    group_by(source_db) %>%
    summarise(
      n_total = n(),
      n_zero = sum(organic_carbon == 0),
      pct_zero = round(100 * n_zero / n_total, 1),
      .groups = "drop"
    ) %>%
    filter(n_zero > 0) %>%
    arrange(desc(n_zero))
  
  cat("\n  OC == 0 by source database (top 10):\n")
  print(head(as.data.frame(oc_zero_by_db), 10), row.names = FALSE)
  write.csv(oc_zero_by_db, file.path(PLOT_DIR, "oc_zero_by_source.csv"), row.names = FALSE)
}

d$organic_carbon[!is.na(d$organic_carbon) & d$organic_carbon == 0] <- NA

cat(sprintf("\n  → Set %s OC values to NA\n\n", format(n_oc_zero, big.mark = ",")))


# =============================================================================
# FILTER 2: Remove exact duplicates (coord + OC + BD)
# =============================================================================
#
# Rationale: WoSIS integrates many national/regional databases that overlap.
# The same profile can appear in multiple source databases with identical
# coordinates and measurements. Retaining duplicates inflates sample size
# and biases spatial CV (same data in train and test).
# We use coord + OC + BD as the composite key — profiles at the same
# location with different OC values are legitimate (different sampling dates
# or depths in the original data, now standardised to 0-20cm).
# =============================================================================

cat("─── Filter 2: Remove exact duplicates ───\n")

n_before <- nrow(d)

# Composite key: coordinates + organic carbon + bulk density
# Use NA-safe string representation
d$dup_key <- paste(
  d$longitude_decimal_degrees,
  d$latitude_decimal_degrees,
  ifelse(is.na(d$organic_carbon), "NA_oc", d$organic_carbon),
  ifelse(is.na(d$bulk_density_oven_dry), "NA_bd", d$bulk_density_oven_dry),
  sep = "|"
)

is_dup <- duplicated(d$dup_key)
n_dupes <- sum(is_dup)

cat(sprintf("  Exact duplicates (coord + OC + BD): %s (%.1f%%)\n",
            format(n_dupes, big.mark = ","),
            100 * n_dupes / nrow(d)))

# Which databases lose the most?
if ("source_db" %in% names(d)) {
  dup_by_source <- d %>%
    mutate(is_dup = is_dup) %>%
    group_by(source_db) %>%
    summarise(
      n_total = n(),
      n_dup = sum(is_dup),
      pct_dup = round(100 * n_dup / n_total, 1),
      .groups = "drop"
    ) %>%
    filter(n_dup > 0) %>%
    arrange(desc(n_dup))
  
  cat("\n  Duplicates by source (top 10):\n")
  print(head(as.data.frame(dup_by_source), 10), row.names = FALSE)
  write.csv(dup_by_source, file.path(PLOT_DIR, "duplicates_by_source.csv"), row.names = FALSE)
}

d <- d[!is_dup, ]
d$dup_key <- NULL

cat(sprintf("\n  → Removed %s duplicates\n", format(n_dupes, big.mark = ",")))
cat(sprintf("  → Remaining: %s profiles\n\n", format(nrow(d), big.mark = ",")))

removal_log <- log_filter(removal_log, "2_duplicates", "Exact coord+OC+BD duplicates",
                          n_before, nrow(d))


# =============================================================================
# FILTER 3: Implausible bulk density → NA
# =============================================================================
#
# Rationale: BD < 0.1 is below any known mineral soil. BD > 2.5 approaches
# mineral particle density. These are measurement or unit errors.
# Setting to NA lets the gap-filling in step 10 handle them.
# BD = 1.600 appears suspiciously often (plateau at median and p75),
# suggesting a default fill value in one or more source databases.
# We flag but do not remove these.
# =============================================================================

cat("─── Filter 3: Implausible bulk density → NA ───\n")

if ("bulk_density_oven_dry" %in% names(d)) {
  bd <- d$bulk_density_oven_dry
  
  n_bd_low  <- sum(bd < 0.1, na.rm = TRUE)
  n_bd_high <- sum(bd > 2.5, na.rm = TRUE)
  n_bd_plateau <- sum(bd == 1.6, na.rm = TRUE)
  n_bd_valid <- sum(!is.na(bd))
  
  cat(sprintf("  BD < 0.1: %d profiles\n", n_bd_low))
  cat(sprintf("  BD > 2.5: %d profiles\n", n_bd_high))
  cat(sprintf("  BD == 1.600 exactly: %d (%.1f%% of measured BD)\n",
              n_bd_plateau, 100 * n_bd_plateau / n_bd_valid))
  
  # Investigate the plateau: is it from one source?
  if ("source_db" %in% names(d)) {
    plateau_by_db <- d %>%
      filter(!is.na(bulk_density_oven_dry)) %>%
      group_by(source_db) %>%
      summarise(
        n_bd = n(),
        n_plateau = sum(bulk_density_oven_dry == 1.6),
        pct_plateau = round(100 * n_plateau / n_bd, 1),
        .groups = "drop"
      ) %>%
      filter(n_plateau > 0) %>%
      arrange(desc(n_plateau))
    
    cat("\n  BD == 1.600 by source (top 5):\n")
    print(head(as.data.frame(plateau_by_db), 5), row.names = FALSE)
    write.csv(plateau_by_db, file.path(PLOT_DIR, "bd_plateau_by_source.csv"), row.names = FALSE)
  }
  
  # Set implausible to NA
  implausible <- !is.na(bd) & (bd < 0.1 | bd > 2.5)
  d$bulk_density_oven_dry[implausible] <- NA
  
  # Flag plateau (retain but mark)
  d$bd_plateau_flag <- !is.na(d$bulk_density_oven_dry) & d$bulk_density_oven_dry == 1.6
  
  cat(sprintf("\n  → Set %d implausible BD to NA\n", sum(implausible)))
  cat(sprintf("  → Flagged %d BD plateau profiles\n\n", sum(d$bd_plateau_flag)))
}


# =============================================================================
# FILTER 4: Remove profiles with no OC
# =============================================================================
#
# Rationale: Without OC, SOC stock and MTT cannot be computed. These profiles
# contribute nothing to the modelling target.
# =============================================================================

cat("─── Filter 4: Remove profiles with no OC ───\n")

n_before <- nrow(d)
n_oc_na <- sum(is.na(d$organic_carbon))

cat(sprintf("  Profiles with OC = NA: %s (%.1f%%)\n",
            format(n_oc_na, big.mark = ","),
            100 * n_oc_na / nrow(d)))

d <- d[!is.na(d$organic_carbon), ]

cat(sprintf("  → Removed %s profiles\n", format(n_before - nrow(d), big.mark = ",")))
cat(sprintf("  → Remaining: %s profiles\n\n", format(nrow(d), big.mark = ",")))

removal_log <- log_filter(removal_log, "4_no_oc", "No OC value (including former OC==0)",
                          n_before, nrow(d))


# =============================================================================
# FILTER 5: Flag low-precision coordinates
# =============================================================================

cat("─── Filter 5: Flag low-precision coordinates ───\n")

d$coord_integer_flag <- (
  d$longitude_decimal_degrees == round(d$longitude_decimal_degrees) &
  d$latitude_decimal_degrees == round(d$latitude_decimal_degrees)
)

cat(sprintf("  Integer-only coordinates: %d (%.2f%%)\n",
            sum(d$coord_integer_flag),
            100 * mean(d$coord_integer_flag)))
cat("  → Flagged (not removed)\n\n")


# =============================================================================
# FILTER 6: Extreme OC values
# =============================================================================
#
# Rationale: Very low OC (< 0.5 g/kg) in 0-20cm mineral soil is unusual and
# may reflect analytical issues. Very high OC near the 120 g/kg threshold
# may be organic soils that barely passed the filter.
# We flag but do not remove.
# =============================================================================

cat("─── Filter 6: Flag extreme OC values ───\n")

d$oc_low_flag <- d$organic_carbon < 0.5
d$oc_high_flag <- d$organic_carbon > 100

cat(sprintf("  OC < 0.5 g/kg: %d\n", sum(d$oc_low_flag, na.rm = TRUE)))
cat(sprintf("  OC > 100 g/kg: %d\n", sum(d$oc_high_flag, na.rm = TRUE)))
cat("  → Flagged (not removed)\n\n")


# =============================================================================
# SUMMARY
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  QC SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

removal_log <- log_filter(removal_log, "final", "Final dataset", nrow(d), nrow(d))

cat(sprintf("  Start:     %s profiles\n", format(n_start, big.mark = ",")))
cat(sprintf("  OC → NA:   %s (OC==0, not removed but excluded from MTT)\n",
            format(n_oc_zero, big.mark = ",")))
cat(sprintf("  Removed:   %s duplicates\n", format(n_dupes, big.mark = ",")))
cat(sprintf("  Removed:   %s profiles with no OC\n",
            format(n_start - n_dupes - nrow(d), big.mark = ",")))
cat(sprintf("  Final:     %s profiles (%.1f%% of input)\n",
            format(nrow(d), big.mark = ","),
            100 * nrow(d) / n_start))

cat(sprintf("\n  Flags retained in data:\n"))
cat(sprintf("    bd_plateau_flag:    %d profiles\n", sum(d$bd_plateau_flag, na.rm = TRUE)))
cat(sprintf("    coord_integer_flag: %d profiles\n", sum(d$coord_integer_flag)))
cat(sprintf("    oc_low_flag:        %d profiles\n", sum(d$oc_low_flag, na.rm = TRUE)))
cat(sprintf("    oc_high_flag:       %d profiles\n", sum(d$oc_high_flag, na.rm = TRUE)))

# Save removal log
write.csv(removal_log, file.path(PLOT_DIR, "removal_log.csv"), row.names = FALSE)
cat(sprintf("\n  Removal log: %s\n", file.path(PLOT_DIR, "removal_log.csv")))


# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================

cat("\n  Creating diagnostic plots...\n")

# --- OC distribution before/after ---
d_before <- readRDS(INPUT_FILE)

png(file.path(PLOT_DIR, "oc_distribution_comparison.png"),
    width = 1400, height = 600, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

hist(d_before$organic_carbon[d_before$organic_carbon > 0], breaks = 100,
     col = "grey70", border = NA,
     main = "OC distribution (before QC)", xlab = "OC (g/kg)")
abline(v = median(d_before$organic_carbon[d_before$organic_carbon > 0], na.rm = TRUE),
       col = "red", lwd = 2, lty = 2)

hist(d$organic_carbon, breaks = 100,
     col = "steelblue", border = NA,
     main = "OC distribution (after QC)", xlab = "OC (g/kg)")
abline(v = median(d$organic_carbon, na.rm = TRUE),
       col = "red", lwd = 2, lty = 2)

dev.off()
cat("  ✓ oc_distribution_comparison.png\n")

# --- BD distribution with plateau highlighted ---
if ("bulk_density_oven_dry" %in% names(d)) {
  png(file.path(PLOT_DIR, "bd_distribution.png"),
      width = 800, height = 600, res = 150)
  
  bd_vals <- d$bulk_density_oven_dry[!is.na(d$bulk_density_oven_dry)]
  hist(bd_vals, breaks = 100, col = "steelblue", border = NA,
       main = "Measured BD after QC", xlab = "BD (g/cm³)")
  abline(v = 1.6, col = "red", lwd = 2, lty = 2)
  text(1.65, par("usr")[4] * 0.9, "1.600\nplateau", col = "red", cex = 0.8, adj = 0)
  
  dev.off()
  cat("  ✓ bd_distribution.png\n")
}

# --- Source database breakdown after QC ---
if ("source_db" %in% names(d)) {
  src_after <- d %>%
    count(source_db, name = "n") %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    arrange(desc(n))
  
  write.csv(src_after, file.path(PLOT_DIR, "source_db_after_qc.csv"), row.names = FALSE)
  cat("  ✓ source_db_after_qc.csv\n")
}

# --- Spatial coverage before/after ---
png(file.path(PLOT_DIR, "spatial_coverage_comparison.png"),
    width = 2000, height = 800, res = 150)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))

plot(d_before$longitude_decimal_degrees, d_before$latitude_decimal_degrees,
     pch = ".", col = rgb(0, 0, 0, 0.1), asp = 1,
     main = sprintf("Before QC (n = %s)", format(nrow(d_before), big.mark = ",")),
     xlab = "", ylab = "", xlim = c(-180, 180), ylim = c(-60, 85))

plot(d$longitude_decimal_degrees, d$latitude_decimal_degrees,
     pch = ".", col = rgb(0, 0, 0, 0.1), asp = 1,
     main = sprintf("After QC (n = %s)", format(nrow(d), big.mark = ",")),
     xlab = "", ylab = "", xlim = c(-180, 180), ylim = c(-60, 85))

dev.off()
cat("  ✓ spatial_coverage_comparison.png\n")

rm(d_before)


# =============================================================================
# SAVE OUTPUT
# =============================================================================

cat("\n")
saveRDS(d, OUTPUT_FILE)
cat(sprintf("✓ Saved: %s\n", OUTPUT_FILE))
cat(sprintf("  %s profiles, %d variables\n\n", format(nrow(d), big.mark = ","), ncol(d)))

cat("═══════════════════════════════════════════════════════════════\n")
cat("  STEP 11.5 COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Next: update step 12 INPUT_FILE to read '11_5_quality_controlled.rds'\n")
cat("  Then re-run steps 12 → 18\n")
