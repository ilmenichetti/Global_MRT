# =============================================================================
# Verification script: Test config, utilities, and directory structure
# Run this FIRST to verify everything is set up correctly
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║          Pipeline Setup Verification                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# --- Test 1: Load configuration ---
cat("[1] Loading configuration...\n")
tryCatch({
  source("./Global_MRT_code/00_config.R")
  cat("  ✓ Config loaded successfully\n")
}, error = function(e) {
  cat("  ✗ Error loading config:\n")
  cat(sprintf("    %s\n", e$message))
  stop("Config loading failed")
})

# --- Test 2: Load utilities ---
cat("\n[2] Loading utilities...\n")
tryCatch({
  source("./Global_MRT_code/01_utils.R")
  cat("  ✓ Utilities loaded successfully\n")
}, error = function(e) {
  cat("  ✗ Error loading utilities:\n")
  cat(sprintf("    %s\n", e$message))
  stop("Utilities loading failed")
})

# --- Test 3: Load packages ---
cat("\n[3] Loading required packages...\n")
tryCatch({
  load_packages()
  cat("  ✓ All packages loaded\n")
}, error = function(e) {
  cat("  ✗ Error loading packages:\n")
  cat(sprintf("    %s\n", e$message))
  cat("  Install missing packages with: install.packages('package_name')\n")
})

# --- Test 4: Verify input data exists ---
cat("\n[4] Checking input data directory...\n")
if (dir.exists(PATHS$data_dir)) {
  files_in_data <- list.files(PATHS$data_dir)
  cat(sprintf("  ✓ Found: %s\n", PATHS$data_dir))
  cat(sprintf("    Contents (%d items):\n", length(files_in_data)))
  for (f in head(files_in_data, 5)) {
    cat(sprintf("      - %s\n", f))
  }
  if (length(files_in_data) > 5) {
    cat(sprintf("      ... and %d more\n", length(files_in_data) - 5))
  }
} else {
  cat(sprintf("  ✗ NOT FOUND: %s\n", PATHS$data_dir))
  cat("  Make sure you're running from Global_MRT_code/ directory\n")
}

# --- Test 5: Create directory structure ---
cat("\n[5] Creating output directories...\n")
setup_directories()

# --- Test 6: Verify directory structure ---
cat("\n[6] Verifying directory structure...\n")
key_dirs <- c(
  PATHS$output_dir,
  PATHS$raster_output_dir,
  file.path(PATHS$raster_output_dir, "climate"),
  file.path(PATHS$raster_output_dir, "soilgrids"),
  PATHS$temp_dir,
  PATHS$era5_dir
)

all_exist <- TRUE
for (d in key_dirs) {
  exists <- dir.exists(d)
  status <- if (exists) "✓" else "✗"
  cat(sprintf("  %s %s\n", status, d))
  if (!exists) all_exist <- FALSE
}

# --- Test 7: Test utilities ---
cat("\n[7] Testing utility functions...\n")
tryCatch({
  grid <- create_global_grid(0.5)
  cat(sprintf("  ✓ Global grid created: %d x %d cells\n", 
              nrow(grid), ncol(grid)))
  cat(sprintf("    Resolution: 0.5° (expected ~720 x 360)\n")
  )
}, error = function(e) {
  cat(sprintf("  ✗ Error creating grid: %s\n", e$message))
})

# --- Summary ---
cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                  Verification Summary                         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

if (all_exist) {
  cat("✓ All checks passed! Ready to proceed.\n\n")
  cat("Next step: Run climate extraction module (03_extract_climate.R)\n")
} else {
  cat("✗ Some issues detected. Please fix above before proceeding.\n")
}

cat("\nKey paths:\n")
cat(sprintf("  Input data:     %s\n", PATHS$data_dir))
cat(sprintf("  Point output:   %s\n", PATHS$output_dir))
cat(sprintf("  Raster output:  %s\n", PATHS$raster_output_dir))
cat(sprintf("  Working dir:    %s\n", getwd()))
cat("\n")
cat("Verification complete.\n")
