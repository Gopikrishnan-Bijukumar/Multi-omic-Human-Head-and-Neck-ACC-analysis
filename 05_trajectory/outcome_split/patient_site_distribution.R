#!/usr/bin/env Rscript
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# patient_site_distribution.R
#
# Builds the per-patient clinical/site table for the outcome-split cohorts.
# READ-ONLY with respect to the pipeline; writes one CSV.
#
# WHY THIS EXISTS
# ---------------
# Anatomical site is almost perfectly separated between the two outcome arms:
# of the 19 patients with a recorded site, no site appears in both cohorts.
# We do NOT model site and make no site claim -- ACC is regarded as molecularly
# site-independent, and 19 patients across 10 sites supports no claim in either
# direction. But the distribution is REPORTED rather than left implicit, so a
# reader meets it as a disclosed characteristic of the cohort rather than as
# something they discovered for themselves.
#
# NOTE ON LABELS: `BOT` and `BaseofTongue` are recorded as distinct strings but
# are almost certainly the same anatomical site (base of tongue). They are NOT
# merged here -- this table reports what the h5ad actually contains. The
# harmonised view is provided as a separate column so the ambiguity is visible
# rather than silently resolved. If they are the same site, base of tongue is
# the only site represented in both arms (2 Good vs 2 Poor).
#
# Usage:  Rscript patient_site_distribution.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr) })

BASE_DIR <- ACC_DATA_ROOT
OUT_CSV  <- file.path(BASE_DIR, "outputs", "outcome_split",
                      "patient_site_distribution.csv")

# Same string, two spellings. Recorded, not silently merged (see header).
SITE_ALIASES <- c("BOT" = "BaseofTongue")

rows <- lapply(c("good", "poor"), function(co) {
  f <- file.path(BASE_DIR, "outputs", "outcome_split", co, "monocle3_analysis",
                 "tables", sprintf("cds_per_cell_metadata_%s.csv", co))
  if (!file.exists(f)) stop("missing input: ", f)
  d <- read.csv(f)
  need <- c("batch", "outcome_group", "years_of_survival", "site", "annot_5")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("missing columns in ", basename(f), ": ",
                         paste(miss, collapse = ", "))
  d %>%
    group_by(patient = batch, outcome_group, years_of_survival, site) %>%
    summarise(n_malignant_cells = dplyr::n(),
              n_annot5_classes  = dplyr::n_distinct(annot_5),
              .groups = "drop")
})
tab <- bind_rows(rows)

# One row per patient. If a patient ever spanned two sites this would surface as
# a duplicate rather than being averaged away.
dup <- tab$patient[duplicated(tab$patient)]
if (length(dup)) stop("patient(s) with >1 site row: ", paste(unique(dup), collapse = ", "))

tab <- tab %>%
  mutate(site_recorded   = ifelse(is.na(site) | site == "", NA_character_, site),
         site_harmonised = ifelse(!is.na(site_recorded) & site_recorded %in% names(SITE_ALIASES),
                                  unname(SITE_ALIASES[site_recorded]), site_recorded),
         site_known      = !is.na(site_recorded)) %>%
  select(patient, outcome_group, years_of_survival, site_recorded, site_harmonised,
         site_known, n_malignant_cells, n_annot5_classes) %>%
  arrange(outcome_group, site_harmonised, patient)

write.csv(tab, OUT_CSV, row.names = FALSE)
message("Saved: ", OUT_CSV, " (", nrow(tab), " patients)")

# --- Console summary ---------------------------------------------------------
message("\n=== patients per arm ===")
print(as.data.frame(table(tab$outcome_group)), row.names = FALSE)

message("\n=== site (as recorded) x outcome ===")
print(table(tab$site_recorded, tab$outcome_group, useNA = "ifany"))

message("\n=== site (harmonised: BOT == BaseofTongue) x outcome ===")
h <- table(tab$site_harmonised, tab$outcome_group, useNA = "ifany")
print(h)

shared <- rownames(h)[rowSums(h > 0) > 1 & !is.na(rownames(h))]
if (length(shared)) {
  message("\n  Sites represented in BOTH arms after harmonisation: ",
          paste(shared, collapse = ", "))
  message("  -> the only place site and outcome are separable in this cohort.")
} else {
  message("\n  NO site is represented in both arms, even after harmonisation:")
  message("  site and outcome are completely confounded here. Reported, not modelled.")
}
message("\n  Patients with no recorded site: ", sum(!tab$site_known),
        " (Good ", sum(!tab$site_known & tab$outcome_group == "Good"),
        ", Poor ", sum(!tab$site_known & tab$outcome_group == "Poor"), ")")
