# --- 1. Model comparison ---
models <- readRDS("./Global_MRT_code/outputs/13_rf_models.rds")
cat("=== MODEL COMPARISON ===\n")
for (m in names(models)) {
  met <- models[[m]]$metrics
  cat(sprintf("%-30s  n_pred=%d  n_obs=%d  R2_cv=%.4f  RMSE_cv=%.4f  MAE_cv=%.4f  R2_oob=%.4f\n",
              models[[m]]$config$name,
              met$n_predictors, met$n_obs,
              met$R2_cv, met$RMSE_cv, met$MAE_cv, met$R2_oob))
}

# --- 2. Variable importance top 20 ---
cat("\n=== VARIABLE IMPORTANCE (M7, top 20) ===\n")
imp <- sort(models$M7_full$importance, decreasing = TRUE)
print(data.frame(variable = names(imp), importance = unname(imp))[1:20, ])

# --- 3. Group-level importance ---
cat("\n=== GROUP IMPORTANCE ===\n")
var_groups <- readRDS("./Global_MRT_code/outputs/13_var_groups.rds")
imp_all <- models$M7_full$importance
for (grp in names(var_groups)) {
  vars_in <- intersect(names(imp_all), var_groups[[grp]])
  cat(sprintf("%-12s  n_vars=%d  total=%.4f  mean=%.4f\n",
              grp, length(vars_in), sum(imp_all[vars_in]), mean(imp_all[vars_in])))
}

# --- 4. MTT raster distribution ---
cat("\n=== MTT DISTRIBUTION (M7 raster) ===\n")
library(terra)
r <- rast("./Global_MRT_code/outputs/MRT_predictions/MRT_M7_full.tif")
v <- values(r, na.rm = TRUE)
cat("n pixels:", length(v), "\n")
print(round(quantile(v, probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1)), 2))
cat("Mean:", round(mean(v), 2), "  SD:", round(sd(v), 2), "\n")

# --- 5. Sample sizes ---
cat("\n=== SAMPLE SIZES ===\n")
for (f in c("11_5_quality_controlled.rds", "12_with_MRT.rds", "12b_model_ready.rds")) {
  fp <- file.path("./Global_MRT_code/outputs", f)
  if (file.exists(fp)) {
    d <- readRDS(fp)
    cat(sprintf("%-35s  rows=%s\n", f, format(nrow(d), big.mark = ",")))
    if ("MRT_QC" %in% names(d)) print(table(d$MRT_QC, useNA = "ifany"))
  }
}

# --- 6. BD source ---
cat("\n=== BD SOURCE ===\n")
d <- readRDS("./Global_MRT_code/outputs/12b_model_ready.rds")
if ("bd_source" %in% names(d)) print(table(d$bd_source, useNA = "ifany"))

# --- 7. ALE top 10 ---
cat("\n=== ALE RANKING ===\n")
ale <- readRDS("./Global_MRT_code/outputs/18_ale_results.rds")
global_ale <- ale[ale$lat_band == "Global", ]
ale_range <- aggregate(ale_effect ~ variable, data = global_ale,
                       FUN = function(x) diff(range(x, na.rm = TRUE)))
ale_range <- ale_range[order(-ale_range$ale_effect), ]
names(ale_range)[2] <- "range_years"
print(ale_range)

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
