# =============================================================================
# Shared peatland exclusion filter  (MAIN default = drop peat, GPM 2.0 source)
# =============================================================================
# Sourced by the fitting scripts (13, 13c, 13e, 13g) so the peat-excluded scope
# is applied identically everywhere and cannot drift between them.
#
# Env-controlled (defaults = the production MAIN, i.e. peat excluded):
#   MRT_13C_PEAT    = drop (DEFAULT, MAIN) | keep
#   MRT_PEAT_SOURCE = gpm  (DEFAULT, MAIN) | histosol   (modal-Histosol cross-check)
#
# The gpm source requires outputs/08b_peat_flag.rds (built by 08b_extract_peat.R),
# joined to the observation table by dsiteid. See the decision record
# manuscript/decisions/2026-07-08_peat_exclusion_plan.md.
# =============================================================================

apply_peat_filter <- function(soil_data, output_dir,
                              id_col = "dsiteid", class_col = "soilclass_wrb_name",
                              peat_mode = NULL, peat_source = NULL) {
  # peat_mode/peat_source default to the env vars (production MAIN); callers with an
  # explicit intent (e.g. 13j's per-convention sweep) can override them directly.
  if (is.null(peat_mode))   peat_mode   <- tolower(Sys.getenv("MRT_13C_PEAT",   "drop"))
  if (is.null(peat_source)) peat_source <- tolower(Sys.getenv("MRT_PEAT_SOURCE", "gpm"))
  stopifnot(peat_mode %in% c("keep", "drop"),
            peat_source %in% c("gpm", "histosol"))
  if (peat_mode != "drop") return(soil_data)

  n0 <- nrow(soil_data)
  if (peat_source == "gpm") {
    flag_file <- file.path(output_dir, "08b_peat_flag.rds")
    if (!file.exists(flag_file))
      stop("GPM peat flag not found: ", flag_file, "\nRun 08b_extract_peat.R first.")
    pf      <- readRDS(flag_file)
    is_peat <- soil_data[[id_col]] %in% pf[[id_col]][pf$peat_flag_gpm]
    src_lbl <- "GPM 2.0 (both tiers)"
  } else {
    is_peat <- !is.na(soil_data[[class_col]]) & soil_data[[class_col]] == "Histosols"
    src_lbl <- "modal-Histosol"
  }
  out <- soil_data[!is_peat, ]
  cat(sprintf("Peatland filter [%s]: dropped %d rows (%.2f%%); %d remain\n",
              src_lbl, n0 - nrow(out), 100 * (n0 - nrow(out)) / n0, nrow(out)))
  out
}
