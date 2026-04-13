# --- 1. Model comparison table (exact values) ---
models <- readRDS("./Global_MRT_code/outputs/13_rf_models.rds")
comparison <- read.csv("./Global_MRT_code/outputs/13_model_comparison.csv")

# or rebuild:
cat("=== MODEL COMPARISON ===\n")
for (m in names(models)) {
  met <- models[[m]]$metrics
  cat(sprintf("%-30s  n_pred=%d  n_obs=%d  R2_cv=%.4f  RMSE_cv=%.4f  MAE_cv=%.4f  R2_oob=%.4f\n",
              models[[m]]$config$name,
              met$n_predictors, met$n_obs,
              met$R2_cv, met$RMSE_cv, met$MAE_cv, met$R2_oob))
}

# --- 2. Variable importance top 20 (from M7) ---
cat("\n=== VARIABLE IMPORTANCE (M7, top 20) ===\n")
imp <- sort(models$M7_full$importance, decreasing = TRUE)
imp_df <- data.frame(variable = names(imp), importance = unname(imp))
imp_df$group <- ifelse(imp_df$variable %in% models$M7_full$config$vars, "", "")
# Just print top 20
print(head(imp_df, 20))

# --- 3. MTT distribution summary ---
cat("\n=== MTT DISTRIBUTION (M7 raster) ===\n")
library(terra)
r <- rast("./Global_MRT_code/outputs/MRT_predictions/MRT_M7_full.tif")
v <- values(r, na.rm = TRUE)
cat("n pixels:", length(v), "\n")
cat("Quantiles:\n")
print(round(quantile(v, probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)), 2))
cat("Mean:", round(mean(v), 2), "\n")
cat("SD:", round(sd(v), 2), "\n")

# --- 4. Sample sizes at each pipeline step ---
cat("\n=== SAMPLE SIZES ===\n")
for (f in c("02_mineral_soils_filtered.rds",
            "02_5_standardized_0_20cm.rds",
            "12_with_MRT.rds",
            "12b_model_ready.rds")) {
  fp <- file.path("./Global_MRT_code/outputs", f)
  if (file.exists(fp)) {
    d <- readRDS(fp)
    cat(sprintf("%-35s  rows=%s  cols=%d\n", f, format(nrow(d), big.mark=","), ncol(d)))
    if ("MRT_QC" %in% names(d)) {
      cat("  MRT_QC:\n")
      print(table(d$MRT_QC, useNA = "ifany"))
    }
  }
}

# --- 5. BD gap-fill source distribution ---
cat("\n=== BD SOURCE ===\n")
d <- readRDS("./Global_MRT_code/outputs/10_with_bd_filled.rds")
print(table(d$bd_source, useNA = "ifany"))

# --- 6. Group-level importance sums ---
cat("\n=== GROUP IMPORTANCE (M7) ===\n")
var_groups <- readRDS("./Global_MRT_code/outputs/13_var_groups.rds")
imp_all <- models$M7_full$importance
for (grp in names(var_groups)) {
  vars_in <- intersect(names(imp_all), var_groups[[grp]])
  cat(sprintf("%-12s  n_vars=%d  total_imp=%.4f  mean_imp=%.4f\n",
              grp, length(vars_in), sum(imp_all[vars_in]), mean(imp_all[vars_in])))
}



# =============================================================================
# ALE MAGNITUDES — exact values from saved results
# =============================================================================

cat("=== ALE MAGNITUDES ===\n\n")

ale_results <- readRDS("./Global_MRT_code/outputs/18_ale_results.rds")

# Check structure
cat("Structure of ale_results:\n")
str(ale_results, max.level = 1)
cat("\n")

# If it's a list of per-variable results:
if (is.list(ale_results) && !is.data.frame(ale_results)) {
  for (var in names(ale_results)) {
    a <- ale_results[[var]]
    if (is.data.frame(a)) {
      # Find ALE effect column (might be named differently)
      ale_col <- grep("ale|effect", names(a), ignore.case = TRUE, value = TRUE)
      if (length(ale_col) > 0) {
        vals <- a[[ale_col[1]]]
        cat(sprintf("%-28s  range=%.3f  min=%.3f  max=%.3f\n",
                    var, diff(range(vals, na.rm = TRUE)),
                    min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
      }
    }
  }
}

# If it's a single data frame (long format):
if (is.data.frame(ale_results)) {
  cat("Columns:", paste(names(ale_results), collapse = ", "), "\n\n")
  
  # Try to find the right columns
  var_col <- grep("variable|var", names(ale_results), ignore.case = TRUE, value = TRUE)[1]
  ale_col <- grep("ale|effect", names(ale_results), ignore.case = TRUE, value = TRUE)[1]
  band_col <- grep("band|lat", names(ale_results), ignore.case = TRUE, value = TRUE)[1]
  bin_col <- grep("mid|bin|x", names(ale_results), ignore.case = TRUE, value = TRUE)[1]
  
  cat("Detected columns:", var_col, ale_col, band_col, bin_col, "\n\n")
  
  # Global ALE ranges per variable
  if (!is.na(band_col)) {
    global <- ale_results[ale_results[[band_col]] == "Global", ]
  } else {
    global <- ale_results
  }
  
  ale_summary <- aggregate(
    global[[ale_col]],
    by = list(variable = global[[var_col]]),
    FUN = function(x) c(range = diff(range(x, na.rm = TRUE)),
                        min = min(x, na.rm = TRUE),
                        max = max(x, na.rm = TRUE))
  )
  ale_summary <- do.call(data.frame, ale_summary)
  names(ale_summary) <- c("variable", "range", "min", "max")
  ale_summary <- ale_summary[order(-ale_summary$range), ]
  
  cat("ALE effect range (years) — all variables, ranked:\n")
  print(ale_summary, row.names = FALSE, digits = 3)
  
  # Top 10: show the actual x-axis values at key ALE points
  cat("\n\nTop 10 — detailed response:\n")
  top10 <- head(ale_summary$variable, 10)
  
  for (var in top10) {
    sub <- global[global[[var_col]] == var, ]
    sub <- sub[order(sub[[bin_col]]), ]
    cat(sprintf("\n--- %s ---\n", var))
    cat(sprintf("  x range: %.3f to %.3f\n", min(sub[[bin_col]]), max(sub[[bin_col]])))
    cat(sprintf("  ALE range: %.3f years (min=%.3f at x=%.3f, max=%.3f at x=%.3f)\n",
                diff(range(sub[[ale_col]])),
                min(sub[[ale_col]]), sub[[bin_col]][which.min(sub[[ale_col]])],
                max(sub[[ale_col]]), sub[[bin_col]][which.max(sub[[ale_col]])]))
  }
}


# =============================================================================
# Quick post-pipeline QC sanity check
# Run after step 18 (or whenever) to confirm filtering propagated
# =============================================================================

cat("═══ POST-PIPELINE QC CHECK ═══\n\n")

# The model-ready file is what step 13 actually trained on
d <- readRDS("./Global_MRT_code/outputs/12b_model_ready.rds")

cat(sprintf("Profiles in model-ready data: %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("OC == 0:  %d  (should be 0)\n", sum(d$organic_carbon == 0, na.rm = TRUE)))
cat(sprintf("OC == NA: %d  (should be 0)\n", sum(is.na(d$organic_carbon))))
cat(sprintf("Valid MTT: %s  (should be >> 100k)\n", 
            format(sum(d$MRT_QC == "valid", na.rm = TRUE), big.mark = ",")))

cat(sprintf("\n✓ Filtering %s\n",
            ifelse(nrow(d) < 300000 & sum(d$organic_carbon == 0, na.rm = TRUE) == 0,
                   "CONFIRMED — clean data",
                   "FAILED — still seeing old data")))
