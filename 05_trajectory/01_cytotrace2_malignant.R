#!/usr/bin/env Rscript
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# cytotrace2_malignant_analysis.R
#
# CytoTRACE2 differentiation/stemness analysis of malignant cells from a
# single-cell RNA-seq dataset of adenoid cystic carcinoma (ACC) tumors.
#
# We use CytoTRACE2 (digitalcytometry/cytotrace2), NOT the original CytoTRACE
# package: original CytoTRACE has significant install friction under R 4.4,
# while CytoTRACE2 is the actively-maintained successor, natively handles
# heterogeneous cell-type mixtures in a single run, and produces both a
# continuous stemness score (CytoTRACE2_Score) and a discrete potency category
# (CytoTRACE2_Potency).
#
# Cell-type identity is taken from the `annot_5` column of colData(sce) -- the
# FINEST annotation level in this dataset (26 levels), and the authoritative
# filter used everywhere else in this script. NOT `prect_1`, and no longer
# `annot_6`.
#
# WHY annot_5 AND NOT annot_6 (changed 2026-08-26): annot_6 is a coarsening of
# annot_5 that merges the two POLES of the trajectory this pipeline infers.
# Its dominant malignant label, "Epithelial cells - Tumor" (57,938 cells, 80%
# of the malignant compartment), is the union of three annot_5 classes:
#     Myoepithelial cells - Tumor        29,177   <- one pole
#     Epithelial cells - Tumor           26,840   <- the other pole
#     Epithelial cells - Basal - Tumor    1,921
# At annot_6 the myoepithelial <-> ductal axis therefore has no annotation to
# stand on, and the root population cannot be named -- it was only ever visible
# through the marker-derived score_myoepithelial proxy. annot_5 is strictly
# hierarchical over annot_6 (no annot_5 class splits across two annot_6
# classes), so this is a pure relabelling: the SAME 71,984 cells are selected.
# `annot_6` is still carried through to every output alongside annot_5, for
# continuity with the reviewed results and for crosswalking to the superseded
# run (see SUPERSEDED_annot6_RUN.md).
#
# The malignant compartment is cross-tabulated against `prect_2` and written to
# qc/celltype_vs_cancer_profile_crosstab.csv. Read that as a CONSISTENCY log,
# NOT as evidence: prect_2 is a 1:1 deterministic recode of the annot_5
# malignant flag (every annot_5 class is 100% "Cancer profile" or 100% "Normal
# profile", zero mixing across all 26 classes), so citing it as independent
# confirmation of malignancy would be circular. See ANALYSIS_NOTES.md 6.7.
#
# RECODE: annot_5 labels the 10,197 actively-dividing tumour cells "Actively
# Dividing cells", with no "- Tumor" tag, even though annot_6 calls the very
# same cells "ADC - Tumor" and annot_7_c2l puts them in "Tumor cells". They are
# tumour. RECODE_ANNOT5 (section 2) renames them to "Actively Dividing cells -
# Tumor" in-script, immediately after the obs read. The source h5ad is never
# modified.
#
# Outcome grouping is controlled by the single constant OUTCOME_COL (see
# section 2 below), currently `clinical_outcome` ("Good"/"Poor") -> 13 Good /
# 11 Poor patients. Verified against inprogress_3.h5ad: `clinical_outcome` is
# identical to the `outcome` column for all 24 patients, and differs from the
# previously-used `model_outcome` (17 Good / 7 Poor) on exactly 4 patients --
# P15, P17, P20, P24 -- which are Good under `model_outcome`
# but Poor under `clinical_outcome`. Whichever column OUTCOME_COL names is
# copied once into an internal `outcome_group` column that the rest of this
# script uses, so switching groupings is a one-line edit.
#
# ----------------------------------------------------------------------------
# Downstream contract: a separate Monocle3 trajectory script consumes two
# cache artifacts written by this script verbatim:
#   1. outputs/h5ad_files/cache/malignant_subset_counts_meta.rds
#   2. outputs/cytotrace_analysis/tables/cytotrace2_scores_malignant.csv
# Do not rename/relocate these without updating that script too.
# =============================================================================

## ---- 0. Smoke-test toggle -------------------------------------------------
SMOKE_TEST <- FALSE  # set FALSE for the full production run
SMOKE_N    <- 5000   # approx. number of malignant cells to keep when smoke-testing

## ---- 1. Package install (idempotent) --------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
cran_pkgs <- c("HiClimR", "RSpectra")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
if (!requireNamespace("CytoTRACE2", quietly = TRUE)) {
  remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r", dependencies = TRUE, upgrade = "never")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(Matrix)
  library(CytoTRACE2); library(DESeq2); library(ggplot2); library(patchwork)
  library(cowplot); library(dplyr); library(tidyr); library(tibble)
  library(svglite); library(ragg); library(ggrastr); library(ggrepel)
  library(viridis); library(RColorBrewer)
})
# CytoTRACE2 (Depends/Imports) attaches `plyr`, which masks several dplyr
# verbs (count, arrange, desc, mutate, rename, summarise, summarize, slice)
# on the search path. Both packages are required by CytoTRACE2's namespace
# so neither can be cleanly detached; instead force the dplyr versions to
# win by binding them directly in the global environment.
for (.fn in c("count", "arrange", "desc", "mutate", "rename",
              "summarise", "summarize", "slice")) {
  assign(.fn, get(.fn, envir = asNamespace("dplyr")), envir = .GlobalEnv)
}
rm(.fn)
set.seed(42)

## ---- 2. Paths ---------------------------------------------------------------
BASE_DIR <- ACC_DATA_ROOT
H5AD_PATH <- file.path(BASE_DIR, "outputs/h5ad_files/inprogress_3.h5ad")

# Single switch controlling the Good/Poor grouping used for every outcome
# comparison below (see script header for how the candidate columns differ).
# The named column is copied once into `outcome_group`; nothing downstream
# refers to the source column name directly.
OUTCOME_COL <- "clinical_outcome"   # was "model_outcome"

# Smoke-test runs are written to a completely separate output root so they can
# never overwrite production results or -- more dangerously -- the shared cache
# RDS that monocle3_malignant_trajectory.R consumes. Before this split, setting
# SMOKE_TEST <- TRUE would clobber both with a ~5k-cell subsample.
OUT_DIR   <- file.path(BASE_DIR, if (SMOKE_TEST) "outputs/smoke_test/cytotrace_analysis"
                                 else            "outputs/cytotrace_analysis")
CACHE_DIR <- file.path(BASE_DIR, if (SMOKE_TEST) "outputs/smoke_test/cache"
                                 else            "outputs/h5ad_files/cache")

for (d in c(OUT_DIR,
            file.path(OUT_DIR, "qc"),
            file.path(OUT_DIR, "umap"),
            file.path(OUT_DIR, "patient_level"),
            file.path(OUT_DIR, "markers"),
            file.path(OUT_DIR, "by_celltype"),
            file.path(OUT_DIR, "tables"),
            CACHE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ---- 3. Plotting helper -----------------------------------------------------
save_plot <- function(p, name, subdir = "", width = 8, height = 6, dpi = 600) {
  out_dir <- file.path(OUT_DIR, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), plot = p, width = width, height = height,
         dpi = dpi, units = "in", bg = "white", device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot = p, width = width, height = height,
         units = "in", bg = "white", device = svglite::svglite)
  message("  Saved: ", subdir, "/", name)
}

## ---- 3b. Shared scVI-UMAP panel style ---------------------------------------
# These panels sit side by side in the same figure:
#   umap/umap_cytotrace_score_malignant       (this script)
#   umap/umap_cytotrace_potency_malignant     (this script)
#   trajectory_original_umap/trajectory_by_*  (monocle3_malignant_trajectory.R)
# They draw the SAME cells at the SAME scVI UMAP coordinates, so any rendering
# difference between them is a visual artefact that makes the comparison
# equivocal. Everything controlling how a cell is drawn is therefore pinned here
# rather than left to a ggplot default.
#
# IMPORTANT: an identical copy of this block lives in
# monocle3_malignant_trajectory.R (section 7a). Each script is deliberately
# standalone -- the pipeline runs them as separate Rscript invocations -- so
# edit BOTH or the panels silently drift apart again.
#
# The `pt_stroke` pin is the one that actually bit us: geom_point()'s default
# stroke = 0.5 adds ~0.95 mm to the drawn diameter, so `size = 0.4` here and
# `size = 0.4, stroke = 0` in the trajectory script produced dots of 2.08 mm vs
# 1.14 mm -- a ~3.3x difference in area from two call sites that both said 0.4.
UMAP_PANEL <- list(
  pt_size    = 0.4,
  pt_stroke  = 0.5,
  pt_alpha   = 0.6,
  raster_dpi = 600,
  width      = 7,      # canvas, inches -- passed to save_plot()
  height     = 6.5,
  panel_w    = 5.6,    # drawn data area, inches -- see force_panel_size()
  panel_h    = 4.8,
  legend_h   = 0.45    # min height of the bottom legend strip, inches
)

# Scatter on the original scVI UMAP coordinates. The colour SCALE is left to the
# caller so each panel keeps its own (viridis-c for continuous, turbo-d for
# potency); only the geometry and the frame are shared.
umap_panel <- function(df, color_col, title, legend_title = color_col,
                       x_col = "X_umap_1", y_col = "X_umap_2") {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                      color = .data[[color_col]])) +
    geom_point(size = UMAP_PANEL$pt_size, stroke = UMAP_PANEL$pt_stroke,
               alpha = UMAP_PANEL$pt_alpha, na.rm = TRUE) +
    labs(title = title, x = "scVI UMAP 1", y = "scVI UMAP 2",
         color = legend_title) +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "bottom")
  # A bottom colourbar defaults to a short bar, which crowds the tick labels
  # into each other ("0.10.20.3..."). Discrete legends are left alone -- keys
  # widened to match would push the potency legend past the canvas edge.
  if (is.numeric(df[[color_col]])) {
    p <- p + guides(color = guide_colourbar(theme = theme(
      legend.key.width  = grid::unit(2.2, "in"),
      legend.key.height = grid::unit(0.16, "in"))))
  }
  ggrastr::rasterise(p, dpi = UMAP_PANEL$raster_dpi)
}

# Pin the panel (the data area, excluding axes/title/legend) to an exact size.
# Without this, a continuous colourbar and a 4-level discrete legend consume
# different amounts of the canvas, so the same canvas size still yields
# different plot areas -- and identical dots then cover different fractions of
# the panel. Equivalent to egg::set_panel_size(), inlined because neither egg
# nor ggh4x is installed here. ggsave() draws the returned gtable directly.
force_panel_size <- function(p, w_in = UMAP_PANEL$panel_w,
                             h_in = UMAP_PANEL$panel_h) {
  # Both ggplotGrob() and convertHeight() need an open graphics device for font
  # metrics. Under Rscript there is none, so R silently opens the DEFAULT device
  # -- pdf() -- which dumps an Rplots.pdf into the working directory (the script
  # dir, since run_production_pipeline.sh cd's there). Open a throwaway device
  # up front instead, before the first call that would trigger that.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  g <- ggplot2::ggplotGrob(p)
  is_panel <- grepl("^panel", g$layout$name)
  g$widths[unique(g$layout$l[is_panel])]  <- grid::unit(w_in, "in")
  g$heights[unique(g$layout$t[is_panel])] <- grid::unit(h_in, "in")

  # Pin the bottom legend strip as well. A horizontal colourbar (~0.37 in) and a
  # row of discrete keys (~0.18 in) are different heights, which would push the
  # panel up or down between figures even with the panel itself pinned -- the
  # panels would then be the same size but not aligned. Only ever grows the
  # strip, so a tall legend (one key per sample) is never squeezed.
  gb <- which(g$layout$name == "guide-box-bottom")
  if (length(gb) == 1) {
    r <- g$layout$t[gb]
    nat <- grid::convertHeight(g$heights[r], "in", valueOnly = TRUE)
    g$heights[r] <- grid::unit(max(nat, UMAP_PANEL$legend_h), "in")
  }
  g
}

# Build + save an scVI-UMAP panel at the pinned panel size. The canvas is sized
# to whatever the pinned panel plus that panel's own decorations actually need,
# rather than fixed: a legend with one key per sample is far taller than a
# colourbar, and a fixed canvas would clip it. The DATA AREA is identical across
# every panel either way, which is the property the figure comparison rests on.
save_umap_panel <- function(p, name, subdir) {
  g <- force_panel_size(p)
  # Throwaway device again, same reason as in force_panel_size(). The small pad
  # absorbs the metric difference between this device and the ragg/svglite
  # devices that actually render the file -- it is outer whitespace and does not
  # touch the panel.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  w <- sum(grid::convertWidth(g$widths, "in", valueOnly = TRUE)) + 0.08
  h <- sum(grid::convertHeight(g$heights, "in", valueOnly = TRUE)) + 0.08
  save_plot(g, name, subdir, width = max(w, UMAP_PANEL$width),
            height = max(h, UMAP_PANEL$height))
}

# ---- Annotation level -------------------------------------------------------
# annot_5 is the finest annotation available (26 levels) and the authoritative
# cell-type column for this pipeline. annot_6 (19 levels) is a coarsening of it
# that merges the trajectory's two poles; it is retained as a passenger column
# only. See the header for the full rationale.
ANNOT_COL        <- "annot_5"
ANNOT_COL_LEGACY <- "annot_6"

# The actively-dividing tumour cells carry no "- Tumor" tag at annot_5, though
# annot_6 labels the identical cells "ADC - Tumor". Applied immediately after
# the obs read, before any filtering, and asserted to hit exactly this many
# cells so a silent upstream relabelling can never pass unnoticed.
RECODE_ANNOT5 <- c("Actively Dividing cells" = "Actively Dividing cells - Tumor")
RECODE_ANNOT5_EXPECTED_N <- 10197L

# The six malignant annot_5 classes. This selects exactly the same 71,984 cells
# the previous four-class annot_6 filter did -- annot_5 is strictly hierarchical
# over annot_6 -- so every EXPECTED_* size gate downstream stays valid.
MALIGNANT_TYPES <- c(
  "Myoepithelial cells - Tumor",                    # 29,177
  "Epithelial cells - Tumor",                       # 26,840
  "Actively Dividing cells - Tumor",                # 10,197  (post-recode)
  "Actively Dividing Myoepithelial cells - Tumor",  #  3,573
  "Epithelial cells - Basal - Tumor",               #  1,921
  "Basal cells - Tumor"                             #    276  (Good-only)
)
EXPECTED_MALIGNANT_N <- 71984L

# Coarser compartment grouping used for the "by cell type" figure. All 26
# annot_5 levels must appear or the hard stop in section 10 fires.
#
# Label strings copied from the validated annot_5 vocabulary in
# ../processed_celltype_composition_barplots.R:94-125. Do NOT retype them --
# annot_5 uses "Macrophage - M2" / "Macrophage - AP/TAM" (SINGULAR) where
# annot_6 used "Macrophages - ...", and "Muscle cells - Tongue" where annot_6
# had "Muscle cells".
COMPARTMENT_MAP <- c(
  # malignant (6)
  "Myoepithelial cells - Tumor"                   = "Malignant",
  "Epithelial cells - Tumor"                      = "Malignant",
  "Actively Dividing cells - Tumor"               = "Malignant",
  "Actively Dividing Myoepithelial cells - Tumor" = "Malignant",
  "Epithelial cells - Basal - Tumor"              = "Malignant",
  "Basal cells - Tumor"                           = "Malignant",
  # non-malignant epithelium (6) -- annot_6 collapsed all of these to one
  # "Epithelial cells" label (4,897 cells)
  "Epithelial cells - Basal"                      = "Normal Epithelial",
  "Epithelial cells - EA"                         = "Normal Epithelial",
  "Epithelial cells - Ductal"                     = "Normal Epithelial",
  "Epithelial cells - Glandular"                  = "Normal Epithelial",
  "Epithelial cells - Secretory"                  = "Normal Epithelial",
  "Epithelial cells - Ciliated"                   = "Normal Epithelial",
  # fibroblast (3)
  "Fibroblast cells"                              = "Fibroblast",
  "Fibroblast cells - CAF"                        = "Fibroblast",
  "Fibroblast cells - Inflammatory"               = "Fibroblast",
  # macrophage (2) -- SINGULAR at annot_5
  "Macrophage - M2"                               = "Macrophage",
  "Macrophage - AP/TAM"                           = "Macrophage",
  # T (2)
  "CD4+ T cells"                                  = "T cells",
  "CD8+ T cells"                                  = "T cells",
  # B / plasma (3)
  "B cells"                                       = "B/Plasma cells",
  "Memory B cells"                                = "B/Plasma cells",
  "Plasma cells"                                  = "B/Plasma cells",
  # singletons (4)
  "Endothelial cells"                             = "Endothelial",
  "Mural cells"                                   = "Mural",
  "Muscle cells - Tongue"                         = "Muscle",
  "Mast cells"                                    = "Mast"
)
stopifnot(length(COMPARTMENT_MAP) == 26L)

# Fixed display order for the six malignant classes, arranged along the
# myoepithelial -> ductal axis so trajectory and violin panels read left to
# right as the axis itself. Used wherever annot_5 is a discrete colour key.
MALIGNANT_LEVELS <- c(
  "Myoepithelial cells - Tumor",
  "Actively Dividing Myoepithelial cells - Tumor",
  "Basal cells - Tumor",
  "Epithelial cells - Basal - Tumor",
  "Actively Dividing cells - Tumor",
  "Epithelial cells - Tumor"
)
MALIGNANT_PALETTE <- c(
  "Myoepithelial cells - Tumor"                   = "#1B6CA8",
  "Actively Dividing Myoepithelial cells - Tumor" = "#5FA8D3",
  "Basal cells - Tumor"                           = "#7B5AA6",
  "Epithelial cells - Basal - Tumor"              = "#C77DA5",
  "Actively Dividing cells - Tumor"               = "#E8A33D",
  "Epithelial cells - Tumor"                      = "#C1272D"
)

OUTCOME_COLORS <- c("Good" = "forestgreen", "Poor" = "firebrick3")

cliffs_delta <- function(x, y) {
  d <- outer(x, y, FUN = function(a, b) sign(a - b))
  mean(d)
}

message("=== [1/14] Loading h5ad (reader='R', use_hdf5=TRUE) ===")

## ---- 4. Load data -------------------------------------------------------
sce <- zellkonverter::readH5AD(
  file = H5AD_PATH,
  use_hdf5 = TRUE, reader = "R"
)
message("Loaded sce: ", nrow(sce), " genes x ", ncol(sce), " cells")

obs_df <- as.data.frame(colData(sce))
obs_df$cell_id <- colnames(sce)

## ---- 4a. Annotation level: annot_5, with the ADC tumour-tag recode ---------
if (!ANNOT_COL %in% colnames(obs_df))
  stop("ANNOT_COL '", ANNOT_COL, "' not found in colData(sce). Annotation-like ",
       "columns present: ",
       paste(grep("^annot_|^prect_", colnames(obs_df), value = TRUE), collapse = ", "))

obs_df[[ANNOT_COL]] <- as.character(obs_df[[ANNOT_COL]])
if (ANNOT_COL_LEGACY %in% colnames(obs_df))
  obs_df[[ANNOT_COL_LEGACY]] <- as.character(obs_df[[ANNOT_COL_LEGACY]])

# Apply the "- Tumor" tag to the actively-dividing tumour cells. Asserted on an
# exact count: if an upstream re-annotation ever changes this population, the
# run stops here rather than silently selecting a different cell set.
.n_recode <- sum(obs_df[[ANNOT_COL]] %in% names(RECODE_ANNOT5))
if (.n_recode != RECODE_ANNOT5_EXPECTED_N)
  stop("RECODE_ANNOT5 matched ", .n_recode, " cells but expected ",
       RECODE_ANNOT5_EXPECTED_N, ". The upstream ", ANNOT_COL, " annotation has ",
       "changed -- re-verify the malignant class list before proceeding.")
if (any(RECODE_ANNOT5 %in% unique(obs_df[[ANNOT_COL]])))
  stop("RECODE_ANNOT5 target name(s) already exist in ", ANNOT_COL,
       "; the recode would merge two distinct classes.")
.hit <- obs_df[[ANNOT_COL]] %in% names(RECODE_ANNOT5)
obs_df[[ANNOT_COL]][.hit] <- RECODE_ANNOT5[obs_df[[ANNOT_COL]][.hit]]
message("Recoded ", .n_recode, " cells: '", names(RECODE_ANNOT5)[1], "' -> '",
        RECODE_ANNOT5[[1]], "' (annot_6 calls the same cells 'ADC - Tumor')")

# Every annot_5 level must be mapped, checked here rather than 40 minutes later
# at the compartment-assignment stop in section 10.
.unmapped <- setdiff(unique(obs_df[[ANNOT_COL]]),
                     c(MALIGNANT_TYPES, names(COMPARTMENT_MAP)))
if (length(.unmapped) > 0)
  stop("COMPARTMENT_MAP is missing an entry for the following ", ANNOT_COL,
       " value(s): ", paste(.unmapped, collapse = ", "))
message("Annotation level: ", ANNOT_COL, " (",
        length(unique(obs_df[[ANNOT_COL]])), " levels), all mapped")

# Resolve OUTCOME_COL -> internal `outcome_group` once, here, so every
# downstream summary/plot/test reads from one canonical column.
if (!OUTCOME_COL %in% colnames(obs_df)) {
  stop("OUTCOME_COL '", OUTCOME_COL, "' not found in colData(sce). ",
       "Outcome-like columns present: ",
       paste(grep("outcome", colnames(obs_df), value = TRUE, ignore.case = TRUE),
             collapse = ", "))
}
obs_df$outcome_group <- as.character(obs_df[[OUTCOME_COL]])

.n_missing_outcome <- sum(is.na(obs_df$outcome_group))
if (.n_missing_outcome > 0)
  warning(.n_missing_outcome, " cell(s) have NA in outcome column '", OUTCOME_COL, "'.")

# Log the patient-level split up front: this is the headline design parameter
# of the whole analysis and should be visible at the top of every run log.
.patient_outcome <- unique(obs_df[, c("batch", "outcome_group")])
message("Outcome grouping from '", OUTCOME_COL, "': ",
        paste(names(table(.patient_outcome$outcome_group)),
              table(.patient_outcome$outcome_group), sep = " = ", collapse = ", "),
        " (", nrow(.patient_outcome), " patients total)")
if (nrow(.patient_outcome) != length(unique(obs_df$batch)))
  stop("At least one patient (batch) maps to more than one '", OUTCOME_COL,
       "' value; outcome must be constant within a patient.")

# KNOWN_OUTCOME_MISMATCH -- surfaced, not corrected.
#
# Cell barcodes carry a 'good_'/'poor_' prefix from an earlier processing pass.
# For exactly one of the 24 patients that prefix disagrees with the outcome
# column this pipeline actually keys off:
#
#     P07   barcodes 'poor_...'   clinical_outcome = 'Good'
#              years_of_survival = 3   (every other Good = 5, every Poor = 2)
#              3,776 cells total / 2,770 malignant = 8.4% of the Good cohort
#
# OUTCOME_COL wins, so P07 is analysed as Good -- unchanged behaviour, this
# is not a regression. It is logged so a reviewer meets it here rather than
# discovering it in a barcode. Do NOT "fix" it by trusting the prefix without
# checking the source clinical records.
.prefix <- sub("_.*$", "", obs_df$cell_id)
.mismatch <- !is.na(obs_df$outcome_group) &
             tolower(.prefix) != tolower(obs_df$outcome_group)
if (any(.mismatch)) {
  .mm <- unique(obs_df$batch[.mismatch])
  for (.p in .mm) {
    .i <- obs_df$batch == .p
    message("KNOWN_OUTCOME_MISMATCH: patient ", .p,
            "  barcode_prefix='", unique(.prefix[.i])[1],
            "'  ", OUTCOME_COL, "='", unique(obs_df$outcome_group[.i])[1], "'",
            if ("years_of_survival" %in% colnames(obs_df))
              paste0("  years_of_survival=", unique(obs_df$years_of_survival[.i])[1])
            else "",
            "  n_cells=", sum(.i),
            "  -> analysed as ", unique(obs_df$outcome_group[.i])[1])
  }
  message("KNOWN_OUTCOME_MISMATCH: ", length(.mm), " of ",
          length(unique(obs_df$batch)), " patients affected. ",
          "This is expected and non-fatal; see ANALYSIS_NOTES.md.")
} else {
  message("Barcode prefix agrees with ", OUTCOME_COL, " for all patients.")
}

# NOTE: on this machine, realizing assay(sce, "counts") via zellkonverter's
# HDF5Array/DelayedArray-backed sparse reader throws "H5Dread() returned an
# error" (both on arbitrary column subsets and on a full realize) -- this
# reproduces even outside this script and is a read-path issue specific to
# the H5SparseMatrix -> DelayedArray -> h5mread chain in this package stack,
# not a problem with the underlying HDF5 file (plain rhdf5::h5read() of the
# same datasets succeeds instantly). colData/reducedDim access via
# zellkonverter is unaffected and used normally above/below.
#
# Workaround: read the raw CSR triplet (/layers/counts/{data,indices,indptr})
# directly with rhdf5 and reassemble it in R as a genes x cells dgCMatrix.
# AnnData stores counts as an n_obs x n_genes CSR matrix; its indices are not
# guaranteed sorted within a row on this file, so we build via the (i, j, x)
# triplet constructor (which sorts/validates) rather than assuming a directly
# transposable CSR->CSC layout. Values are cross-checked below against
# obs$total_counts / obs$n_genes (Pearson r > 0.98; the genes here are the
# pre-subset top-10,000 HVGs, so raw sums are a fraction of the full-
# transcriptome obs totals, as expected -- see script header).
message("  Reconstructing counts matrix from raw HDF5 CSR triplet (H5Dread workaround)...")
suppressPackageStartupMessages(library(rhdf5))
h5_data    <- rhdf5::h5read(H5AD_PATH, "/layers/counts/data")
h5_indices <- rhdf5::h5read(H5AD_PATH, "/layers/counts/indices")
h5_indptr  <- rhdf5::h5read(H5AD_PATH, "/layers/counts/indptr")
nnz_per_cell <- diff(h5_indptr)
counts_full <- Matrix::sparseMatrix(
  i = as.integer(h5_indices) + 1L,
  j = rep.int(seq_along(nnz_per_cell), times = nnz_per_cell),
  x = as.double(h5_data),
  dims = c(nrow(sce), ncol(sce)),
  dimnames = list(rownames(sce), colnames(sce))
)
rm(h5_data, h5_indices, h5_indptr, nnz_per_cell); gc()

.cs_check <- Matrix::colSums(counts_full)
.corr_check <- suppressWarnings(cor(.cs_check, colData(sce)$total_counts))
message("  Sanity check: cor(colSums(counts_full), obs$total_counts) = ", round(.corr_check, 4))
if (is.na(.corr_check) || .corr_check < 0.9) {
  stop("Reconstructed counts matrix failed sanity check against obs$total_counts (cor = ",
       .corr_check, "). Aborting rather than proceeding on possibly-misaligned data.")
}
rm(.cs_check, .corr_check)
message("  counts_full: ", nrow(counts_full), " genes x ", ncol(counts_full), " cells (", class(counts_full)[1], ")")

## ---- 5. QC / malignant filter -----------------------------------------------
message("=== [2/14] QC + malignant filtering ===")

keep_qc        <- is.na(colData(sce)$marked_cells) | colData(sce)$marked_cells != "marked_removal"
# NB: read the annotation from obs_df, NOT colData(sce) -- obs_df is where
# RECODE_ANNOT5 was applied (section 4a). obs_df is built from colData(sce) and
# is in the same cell order by construction, so the index is interchangeable.
stopifnot(identical(obs_df$cell_id, colnames(sce)))
keep_malignant <- obs_df[[ANNOT_COL]] %in% MALIGNANT_TYPES
malignant_idx  <- which(keep_qc & keep_malignant)

message("QC-passed cells: ", sum(keep_qc), " / ", ncol(sce))
message("Malignant (QC-passed) cells: ", length(malignant_idx),
        " (raw ", ANNOT_COL, " malignant count, pre-QC: ", sum(keep_malignant), ")")

# The malignant compartment must be the same 71,984 cells the annot_6 filter
# selected -- annot_5 is strictly hierarchical over annot_6, so any deviation
# means the upstream annotation changed, not that the migration is working.
if (sum(keep_malignant) != EXPECTED_MALIGNANT_N)
  stop("Expected ", EXPECTED_MALIGNANT_N, " malignant cells from ", ANNOT_COL,
       " but got ", sum(keep_malignant), ". MALIGNANT_TYPES and the upstream ",
       "annotation are out of sync -- do not proceed.")

# Per-class breakdown, so the split annot_6 used to hide is in the run log.
print(table(obs_df[[ANNOT_COL]][keep_malignant]))

# Cross-tab annot_5 x prect_2 as a CONSISTENCY log (prect_2 is NOT used for
# filtering, and is NOT independent evidence -- it is a 1:1 deterministic recode
# of the malignant flag: every annot_5 class is 100% "Cancer profile" or 100%
# "Normal profile". Agreement here confirms bookkeeping, not biology.)
if ("prect_2" %in% colnames(obs_df)) {
  crosstab <- table(annot_5 = obs_df[[ANNOT_COL]], prect_2 = obs_df$prect_2)
  write.csv(as.data.frame.matrix(crosstab),
            file.path(OUT_DIR, "qc", "celltype_vs_cancer_profile_crosstab.csv"))
  n_cancer_profile <- sum(obs_df$prect_2 == "Cancer profile", na.rm = TRUE)
  message("Consistency log (NOT independent evidence): ", ANNOT_COL,
          "-based malignant (pre-QC) = ", sum(keep_malignant),
          " vs prect_2 == 'Cancer profile' = ", n_cancer_profile)
} else {
  message("WARNING: prect_2 column not found; skipping cross-tab sanity log.")
}

# annot_5 -> annot_6 crosswalk, plus the assertion that the coarsening really is
# strictly hierarchical. If any annot_5 class ever splits across two annot_6
# classes, every "annot_6 result" would stop being re-derivable from annot_5 and
# the crosswalk to the superseded run would be invalid -- so this is a stop(),
# not a warning.
if (ANNOT_COL_LEGACY %in% colnames(obs_df)) {
  .xt <- table(obs_df[[ANNOT_COL]], obs_df[[ANNOT_COL_LEGACY]])
  .split <- rownames(.xt)[rowSums(.xt > 0) > 1]
  if (length(.split) > 0)
    stop("annot_5 -> annot_6 is NOT strictly hierarchical; these annot_5 ",
         "class(es) map to more than one annot_6 class: ",
         paste(.split, collapse = ", "),
         ". Results can no longer be crosswalked to the superseded annot_6 run.")
  .cw <- as.data.frame(.xt, stringsAsFactors = FALSE)
  names(.cw) <- c("annot_5", "annot_6", "n_cells")
  .cw <- .cw[.cw$n_cells > 0, ]
  .cw$is_malignant <- .cw$annot_5 %in% MALIGNANT_TYPES
  .cw <- .cw[order(-.cw$n_cells), ]
  write.csv(.cw, file.path(OUT_DIR, "qc", "annot5_to_annot6_crosswalk.csv"),
            row.names = FALSE)
  message("Saved qc/annot5_to_annot6_crosswalk.csv (", nrow(.cw),
          " non-empty pairs); hierarchy assertion passed.")
}

## ---- 6. Smoke-test subsampling ----------------------------------------------
if (SMOKE_TEST) {
  message("=== SMOKE_TEST = TRUE : subsampling to ~", SMOKE_N, " malignant cells ===")
  mal_meta <- obs_df[malignant_idx, ]
  frac <- min(1, SMOKE_N / nrow(mal_meta))
  set.seed(42)
  mal_sub <- mal_meta %>%
    group_by(batch) %>%
    slice_sample(prop = frac) %>%
    ungroup()
  malignant_idx <- match(mal_sub$cell_id, obs_df$cell_id)
  message("Smoke-test malignant subset: ", length(malignant_idx), " cells")

  # proportionally subsample the all-cell-type QC-passed set as well, using the
  # SAME downsampling fraction applied to the malignant compartment above
  qc_meta <- obs_df[keep_qc, ]
  set.seed(42)
  qc_sub <- qc_meta %>%
    group_by(batch) %>%
    slice_sample(prop = frac) %>%
    ungroup()
  keep_qc_idx <- match(qc_sub$cell_id, obs_df$cell_id)
  keep_qc <- rep(FALSE, ncol(sce))
  keep_qc[keep_qc_idx] <- TRUE
  message("Smoke-test all-QC-passed subset: ", sum(keep_qc), " cells")
}

## ---- 7. Sample-level QC figures (all QC-passed cells) -----------------------
message("=== [3/14] Sample-level QC figures ===")

qc_df <- obs_df[keep_qc, ]

sample_order <- qc_df %>%
  count(batch, outcome_group) %>%
  arrange(outcome_group, desc(n)) %>%
  pull(batch)
qc_df$batch <- factor(qc_df$batch, levels = sample_order)

# 01: cells per sample by outcome
p01 <- qc_df %>%
  count(batch, outcome_group) %>%
  ggplot(aes(x = batch, y = n, fill = outcome_group)) +
  geom_col() +
  scale_fill_manual(values = OUTCOME_COLORS) +
  coord_flip() +
  labs(title = "Cells per sample (QC-passed)", x = "Sample (batch)", y = "N cells",
       fill = "Outcome") +
  theme_cowplot(font_size = 11)
save_plot(p01, "01_cells_per_sample_by_outcome", "qc", width = 7, height = 8)

# 02: QC metrics by outcome (3-panel)
qc_long <- qc_df %>%
  select(outcome_group, n_genes, total_counts, pct_counts_mt) %>%
  pivot_longer(cols = c(n_genes, total_counts, pct_counts_mt), names_to = "metric", values_to = "value")

metric_labels <- c(n_genes = "N genes", total_counts = "Total counts", pct_counts_mt = "% mito counts")
panels02 <- lapply(names(metric_labels), function(m) {
  ggplot(filter(qc_long, metric == m), aes(x = outcome_group, y = value, fill = outcome_group)) +
    geom_violin(alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = OUTCOME_COLORS) +
    labs(title = metric_labels[[m]], x = NULL, y = NULL, fill = "Outcome") +
    theme_cowplot(font_size = 11)
})
p02 <- wrap_plots(panels02, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p02 <- p02 + plot_annotation(
  title = "QC metrics by outcome (cell-level distributions)",
  subtitle = "Note: significance statistics for these metrics are computed at the patient level elsewhere in this pipeline, not shown here."
)
save_plot(p02, "02_qc_metrics_by_outcome", "qc", width = 13, height = 5.5)

# 03: cell-type proportions by sample (stacked, fill)
# 26 annot_5 levels (was 19 at annot_6), so the legend is taller and the
# canvas widens to keep the bars readable.
n_types <- length(unique(qc_df[[ANNOT_COL]]))
type_pal <- viridis::viridis(n_types, option = "turbo")
p03 <- qc_df %>%
  count(batch, .data[[ANNOT_COL]]) %>%
  ggplot(aes(x = batch, y = n, fill = .data[[ANNOT_COL]])) +
  geom_col(position = "fill") +
  scale_fill_manual(values = type_pal) +
  coord_flip() +
  labs(title = "Cell-type proportions by sample", x = "Sample (batch, ordered by outcome)",
       y = "Proportion", fill = ANNOT_COL) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_cowplot(font_size = 11)
save_plot(p03, "03_celltype_proportions_by_sample", "qc", width = 12, height = 8.5)

# 04: malignant fraction by outcome
mal_frac_df <- obs_df[keep_qc, ] %>%
  mutate(is_malignant = .data[[ANNOT_COL]] %in% MALIGNANT_TYPES) %>%
  group_by(batch, outcome_group) %>%
  summarise(n_total = n(), n_malignant = sum(is_malignant), .groups = "drop") %>%
  mutate(malignant_fraction = n_malignant / n_total)

wt04 <- wilcox.test(malignant_fraction ~ outcome_group, data = mal_frac_df, exact = FALSE)
p04 <- ggplot(mal_frac_df, aes(x = outcome_group, y = malignant_fraction, color = outcome_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.1, size = 2) +
  scale_color_manual(values = OUTCOME_COLORS) +
  annotate("text", x = 1.5, y = max(mal_frac_df$malignant_fraction) * 1.05,
           label = paste0("Wilcoxon p = ", signif(wt04$p.value, 3)), size = 4) +
  labs(title = "Per-patient malignant fraction by outcome", x = "Outcome", y = "Malignant fraction",
       color = "Outcome") +
  theme_cowplot(font_size = 12)
save_plot(p04, "04_malignant_fraction_by_outcome", "qc", width = 6, height = 6)

# 05: UMAP by site, faceted by outcome
umap_coords <- reducedDim(sce, "X_umap")
umap_df <- obs_df
umap_df$UMAP_1 <- umap_coords[, 1]
umap_df$UMAP_2 <- umap_coords[, 2]
umap_df <- umap_df[keep_qc, ]

p05 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = site)) +
  ggrastr::rasterise(geom_point(size = 0.3, alpha = 0.4), dpi = 600) +
  facet_wrap(~outcome_group) +
  labs(title = "UMAP by site, faceted by outcome", color = "Site") +
  theme_cowplot(font_size = 12)
save_plot(p05, "05_umap_by_site_outcome", "qc", width = 12, height = 6)

## ---- 8. CytoTRACE2 run A: malignant only -------------------------------------
message("=== [4/14] CytoTRACE2 run A: malignant cells only ===")

# NOTE on ncores: CytoTRACE2's internal parallelism (parallelize_models /
# parallelize_smoothing) uses parallel::mclapply-style forking. Under R's
# copy-on-write + refcounting GC, forked workers touching a multi-GB shared
# matrix during GC can each end up privately copying large chunks of it,
# so peak RSS can grow much faster than a naive "N workers x per-worker-size"
# estimate. On this 60GB-RAM/8GB-swap machine (shared with a concurrent
# Monocle3 job), ncores=20 on the ~95k-cell run OOM-killed the process
# (`sendMaster(...) : ignoring SIGPIPE signal`, a forked-worker-died symptom).
# Capped at 6 cores here as a deliberate deviation from the ncores=20 shown
# in the task's example call, to keep peak memory well within budget; this
# trades some wall-clock time for reliability.
N_CORES <- min(6, parallel::detectCores())

counts_malignant <- as.matrix(counts_full[, malignant_idx])
colnames(counts_malignant) <- obs_df$cell_id[malignant_idx]

scores_malignant <- CytoTRACE2::cytotrace2(
  input = counts_malignant, species = "human", is_seurat = FALSE, slot_type = "counts",
  batch_size = 10000, smooth_batch_size = 1000,
  parallelize_models = TRUE, parallelize_smoothing = TRUE, ncores = N_CORES, seed = 14
)

mal_meta <- obs_df[malignant_idx, ]
rownames(mal_meta) <- mal_meta$cell_id
mal_meta <- cbind(mal_meta, scores_malignant[rownames(mal_meta), , drop = FALSE])

message("CytoTRACE2 run A complete: ", nrow(mal_meta), " malignant cells scored")
# Free the dense malignant-only matrix now; it is cheaply reconstructed later
# (Marker DE section, cache-writing section) from the still-resident sparse
# counts_full rather than being held in memory for the entire run-B step.
rm(counts_malignant); gc()

## ---- 9. Cell-cycle independence diagnostic (malignant CytoTRACE2 vs cell cycle) ----
message("=== [5/14] Cell-cycle independence diagnostic (CytoTRACE2 vs S/G2M score) ===")

# CytoTRACE2 is deliberately scored on raw counts (see script header), which
# means this DE pipeline never sees the scVI batch correction or cell-cycle
# regression applied upstream in Scanpy. Before trusting any downstream DE
# result, check whether the stemness score is simply tracking proliferation.
cc_cor_S   <- suppressWarnings(cor(mal_meta$CytoTRACE2_Score, mal_meta$S_score,
                                    use = "pairwise.complete.obs", method = "pearson"))
cc_cor_G2M <- suppressWarnings(cor(mal_meta$CytoTRACE2_Score, mal_meta$G2M_score,
                                    use = "pairwise.complete.obs", method = "pearson"))

cc_cor_df <- data.frame(
  metric    = c("CytoTRACE2_Score_vs_S_score", "CytoTRACE2_Score_vs_G2M_score"),
  pearson_r = c(cc_cor_S, cc_cor_G2M),
  n_cells   = nrow(mal_meta)
)
write.csv(cc_cor_df, file.path(OUT_DIR, "qc", "cytotrace_cellcycle_correlation.csv"), row.names = FALSE)

message("  Pearson r(CytoTRACE2_Score, S_score)   = ", round(cc_cor_S, 4))
message("  Pearson r(CytoTRACE2_Score, G2M_score) = ", round(cc_cor_G2M, 4))
if (!is.na(cc_cor_S) && abs(cc_cor_S) > 0.5) {
  warning("CytoTRACE2_Score correlates strongly with S_score (r = ", round(cc_cor_S, 3),
          "); stemness signal may partly reflect proliferation rather than true differentiation state.")
}
if (!is.na(cc_cor_G2M) && abs(cc_cor_G2M) > 0.5) {
  warning("CytoTRACE2_Score correlates strongly with G2M_score (r = ", round(cc_cor_G2M, 3),
          "); stemness signal may partly reflect proliferation rather than true differentiation state.")
}

p_cc_S <- ggplot(mal_meta, aes(x = S_score, y = CytoTRACE2_Score)) +
  ggrastr::rasterise(geom_point(size = 0.4, alpha = 0.5, color = "steelblue"), dpi = 600) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  annotate("text", x = min(mal_meta$S_score, na.rm = TRUE),
           y = max(mal_meta$CytoTRACE2_Score, na.rm = TRUE),
           label = paste0("r = ", round(cc_cor_S, 3)), hjust = 0, vjust = 1, size = 4.5) +
  labs(title = "CytoTRACE2 score vs S-phase score (malignant cells)",
       x = "S_score", y = "CytoTRACE2_Score") +
  theme_cowplot(font_size = 12)
save_plot(p_cc_S, "cytotrace_vs_S_score_malignant", "qc", width = 7, height = 6)

p_cc_G2M <- ggplot(mal_meta, aes(x = G2M_score, y = CytoTRACE2_Score)) +
  ggrastr::rasterise(geom_point(size = 0.4, alpha = 0.5, color = "darkorange"), dpi = 600) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  annotate("text", x = min(mal_meta$G2M_score, na.rm = TRUE),
           y = max(mal_meta$CytoTRACE2_Score, na.rm = TRUE),
           label = paste0("r = ", round(cc_cor_G2M, 3)), hjust = 0, vjust = 1, size = 4.5) +
  labs(title = "CytoTRACE2 score vs G2M-phase score (malignant cells)",
       x = "G2M_score", y = "CytoTRACE2_Score") +
  theme_cowplot(font_size = 12)
save_plot(p_cc_G2M, "cytotrace_vs_G2M_score_malignant", "qc", width = 7, height = 6)

## ---- 10. CytoTRACE2 run B: all QC-passed cells (all 26 annot_5 types) --------
message("=== [6/14] CytoTRACE2 run B: all QC-passed cells (all cell types) ===")

counts_all <- as.matrix(counts_full[, which(keep_qc)])
colnames(counts_all) <- obs_df$cell_id[which(keep_qc)]

scores_all <- CytoTRACE2::cytotrace2(
  input = counts_all, species = "human", is_seurat = FALSE, slot_type = "counts",
  batch_size = 10000, smooth_batch_size = 1000,
  parallelize_models = TRUE, parallelize_smoothing = TRUE, ncores = N_CORES, seed = 14
)

all_meta <- obs_df[which(keep_qc), ]
rownames(all_meta) <- all_meta$cell_id
all_meta <- cbind(all_meta, scores_all[rownames(all_meta), , drop = FALSE])
all_meta$compartment <- ifelse(all_meta[[ANNOT_COL]] %in% MALIGNANT_TYPES, "Malignant",
                                COMPARTMENT_MAP[as.character(all_meta[[ANNOT_COL]])])

# Fail loudly rather than silently dropping unmapped cell types (NA compartment)
# from the by-cell-type figures below.
unmapped_annot <- unique(all_meta[[ANNOT_COL]][is.na(all_meta$compartment)])
if (length(unmapped_annot) > 0) {
  stop("COMPARTMENT_MAP is missing an entry for the following ", ANNOT_COL, " value(s): ",
       paste(unmapped_annot, collapse = ", "),
       ". Add them to COMPARTMENT_MAP (or MALIGNANT_TYPES) before proceeding -- ",
       "otherwise these cells would silently be dropped (NA compartment) from the ",
       "by-cell-type figures.")
}

message("CytoTRACE2 run B complete: ", nrow(all_meta), " cells scored (all types)")
rm(counts_all); gc()

## ---- 11. UMAP figures: CytoTRACE score/potency on malignant cells -----------
message("=== [7/14] UMAP figures (malignant CytoTRACE2 score/potency) ===")

mal_umap <- umap_coords[malignant_idx, , drop = FALSE]
mal_meta$UMAP_1 <- mal_umap[, 1]
mal_meta$UMAP_2 <- mal_umap[, 2]

# Geometry (dot size/stroke/alpha, theme, canvas, panel size) comes from the
# shared UMAP_PANEL block in section 3b so these two panels are directly
# comparable with trajectory_original_umap/* from the Monocle3 script. Only the
# colour scale is per-panel.
p_umap_score <- umap_panel(mal_meta, "CytoTRACE2_Score",
                           "CytoTRACE2 score - malignant cells",
                           "CytoTRACE2 Score",
                           x_col = "UMAP_1", y_col = "UMAP_2") +
  scale_color_viridis_c(option = "viridis")
save_umap_panel(p_umap_score, "umap_cytotrace_score_malignant", "umap")

p_umap_potency <- umap_panel(mal_meta, "CytoTRACE2_Potency",
                             "CytoTRACE2 potency - malignant cells",
                             "CytoTRACE2 Potency",
                             x_col = "UMAP_1", y_col = "UMAP_2") +
  scale_color_viridis_d(option = "turbo")
save_umap_panel(p_umap_potency, "umap_cytotrace_potency_malignant", "umap")

## ---- 12. Patient-level summary table -----------------------------------------
message("=== [8/14] Patient-level summary table ===")

patient_counts <- mal_meta %>% count(batch, name = "n_cells_malignant")
low_n_patients <- patient_counts %>% filter(n_cells_malignant < 10) %>% pull(batch)
if (length(low_n_patients) > 0) {
  warning("Dropping patients with <10 malignant cells: ", paste(low_n_patients, collapse = ", "))
}
keep_patients <- patient_counts %>% filter(n_cells_malignant >= 10) %>% pull(batch)

mal_meta_kept <- mal_meta %>% filter(batch %in% keep_patients)

# GLOBAL quantile thresholds pooled across all malignant cells
q_top10 <- quantile(mal_meta_kept$CytoTRACE2_Score, probs = 0.90, na.rm = TRUE)
q_top25 <- quantile(mal_meta_kept$CytoTRACE2_Score, probs = 0.75, na.rm = TRUE)
q_top33 <- quantile(mal_meta_kept$CytoTRACE2_Score, probs = 0.6667, na.rm = TRUE)

patient_level <- mal_meta_kept %>%
  group_by(batch, outcome_group) %>%
  summarise(
    n_cells_malignant = n(),
    mean_cytotrace_score = mean(CytoTRACE2_Score, na.rm = TRUE),
    median_cytotrace_score = median(CytoTRACE2_Score, na.rm = TRUE),
    frac_high_cytotrace_top10 = mean(CytoTRACE2_Score >= q_top10, na.rm = TRUE),
    frac_high_cytotrace_top25 = mean(CytoTRACE2_Score >= q_top25, na.rm = TRUE),
    frac_high_cytotrace_top33 = mean(CytoTRACE2_Score >= q_top33, na.rm = TRUE),
    frac_G1  = mean(phase == "G1", na.rm = TRUE),
    frac_S   = mean(phase == "S", na.rm = TRUE),
    frac_G2M = mean(phase == "G2M", na.rm = TRUE),
    .groups = "drop"
  )

# FILENAME IS STAGE-QUALIFIED ON PURPOSE. This file and
# monocle3_malignant_trajectory.R's patient_level_summary_table.csv used to
# share a name, differing only by output tree (cytotrace_analysis/ vs
# monocle3_analysis/). Any time the two trees were copied or bundled together
# for review, one silently shadowed the other and a reviewer opened a
# CytoTRACE2-only table believing it carried the pseudotime/tertile/cluster
# columns stage 2 advertises. The names can no longer collide.
write.csv(patient_level, file.path(OUT_DIR, "patient_level", "patient_level_summary_cytotrace2.csv"), row.names = FALSE)
message("Patient-level table: ", nrow(patient_level), " patients (",
        sum(patient_level$outcome_group == "Good"), " Good, ",
        sum(patient_level$outcome_group == "Poor"), " Poor)")
message("  Saved: patient_level/patient_level_summary_cytotrace2.csv")

## ---- 13. Statistics: Wilcoxon + Cliff's delta, and LOO sensitivity ----------
message("=== [9/14] Statistics: Wilcoxon + Cliff's delta ===")

metrics <- c("mean_cytotrace_score", "median_cytotrace_score",
             "frac_high_cytotrace_top10", "frac_high_cytotrace_top25", "frac_high_cytotrace_top33")

run_wilcox_battery <- function(df, metrics) {
  good_vals <- df %>% filter(outcome_group == "Good")
  poor_vals <- df %>% filter(outcome_group == "Poor")
  purrr_map <- lapply(metrics, function(m) {
    x <- good_vals[[m]]; y <- poor_vals[[m]]
    wt <- tryCatch(wilcox.test(df[[m]] ~ df$outcome_group, exact = FALSE), error = function(e) NULL)
    data.frame(
      metric = m,
      W = if (!is.null(wt)) unname(wt$statistic) else NA,
      p_value = if (!is.null(wt)) wt$p.value else NA,
      cliffs_delta = cliffs_delta(x, y),
      n_good = length(x),
      n_poor = length(y)
    )
  })
  do.call(rbind, purrr_map)
}

wilcoxon_results <- run_wilcox_battery(patient_level, metrics)
# Stage-qualified for the same reason as the summary table above. These five
# CytoTRACE2 tests are also carried, unchanged, into stage 2's combined
# wilcoxon_results_patient_level.csv under test_family = "cytotrace2", so the
# monocle3 file is a superset and no reviewer has to open two files.
write.csv(wilcoxon_results, file.path(OUT_DIR, "patient_level", "wilcoxon_results_patient_level_cytotrace2.csv"), row.names = FALSE)
message("  Saved: patient_level/wilcoxon_results_patient_level_cytotrace2.csv")

# Leave-one-out: drop each Poor patient one at a time
poor_patients <- patient_level %>% filter(outcome_group == "Poor") %>% pull(batch)
loo_list <- lapply(poor_patients, function(pt) {
  df_loo <- patient_level %>% filter(batch != pt)
  res <- run_wilcox_battery(df_loo, metrics)
  res$dropped_patient <- pt
  res
})
loo_results <- do.call(rbind, loo_list)
# Stage-qualified: stage 2 writes its own leave-one-Poor-patient-out table, on
# different metrics and with a different column schema. Sharing a basename across
# the two output trees is the exact defect the Round 2 review caught elsewhere.
write.csv(loo_results, file.path(OUT_DIR, "patient_level", "loo_sensitivity_poor_dropout_cytotrace2.csv"), row.names = FALSE)
message("  Saved: patient_level/loo_sensitivity_poor_dropout_cytotrace2.csv")

## ---- 14. Patient-level figures -----------------------------------------------
message("=== [10/14] Patient-level figures ===")

p_med <- patient_level %>%
  ggplot(aes(x = outcome_group, y = median_cytotrace_score, color = outcome_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.1, size = 2.5) +
  scale_color_manual(values = OUTCOME_COLORS) +
  annotate("text", x = 1.5, y = max(patient_level$median_cytotrace_score) * 1.03,
           label = paste0("Wilcoxon p = ",
                           signif(wilcoxon_results$p_value[wilcoxon_results$metric == "median_cytotrace_score"], 3)),
           size = 4) +
  labs(title = "Patient-level median CytoTRACE2 score by outcome", x = "Outcome",
       y = "Median CytoTRACE2 score", color = "Outcome") +
  theme_cowplot(font_size = 12)
save_plot(p_med, "patient_median_cytotrace_by_outcome", "patient_level", width = 6, height = 6)

thresh_metrics <- c("frac_high_cytotrace_top10" = "Top 10%", "frac_high_cytotrace_top25" = "Top 25%",
                     "frac_high_cytotrace_top33" = "Top 33%")
panels_thresh <- lapply(names(thresh_metrics), function(m) {
  pval <- wilcoxon_results$p_value[wilcoxon_results$metric == m]
  ggplot(patient_level, aes(x = outcome_group, y = .data[[m]], color = outcome_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.1, size = 2) +
    scale_color_manual(values = OUTCOME_COLORS) +
    labs(title = paste0(thresh_metrics[[m]], " (p = ", signif(pval, 3), ")"), x = NULL,
         y = "Fraction high-CytoTRACE", color = "Outcome") +
    theme_cowplot(font_size = 11)
})
p_thresh <- wrap_plots(panels_thresh, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_thresh <- p_thresh + plot_annotation(title = "Patient-level fraction of high-CytoTRACE malignant cells")
save_plot(p_thresh, "patient_frac_high_cytotrace_thresholds", "patient_level", width = 13, height = 5.5)

## ---- 15. Marker DE: CytoTRACE-high vs -low malignant cells (per-patient pseudobulk DESeq2) ----
message("=== [11/14] Marker DE: CytoTRACE-high vs -low (per-patient pseudobulk DESeq2) ===")

# --- Per-patient (LOCAL) quantile thresholding ------------------------------
# Each cell's High/Low/Mid label is relative to its OWN patient's malignant
# score distribution, not the pooled/global distribution (unlike Section 12's
# separate frac_high_cytotrace_* summary table, which is intentionally left as
# a global/pooled comparison against a shared cohort-wide reference).
mal_meta <- mal_meta %>%
  group_by(batch) %>%
  mutate(
    .q33_local = quantile(CytoTRACE2_Score, probs = 1/3, na.rm = TRUE),
    .q67_local = quantile(CytoTRACE2_Score, probs = 2/3, na.rm = TRUE),
    cytotrace_group = case_when(
      CytoTRACE2_Score > .q67_local ~ "High",
      CytoTRACE2_Score < .q33_local ~ "Low",
      TRUE ~ "Mid"
    )
  ) %>%
  ungroup() %>%
  select(-.q33_local, -.q67_local) %>%
  as.data.frame()
rownames(mal_meta) <- mal_meta$cell_id

message("  Per-patient local cytotrace_group counts:")
print(table(mal_meta$batch, mal_meta$cytotrace_group))

# --- Per-patient pseudobulk DESeq2 (replaces per-cell Seurat FindMarkers) ---
# Cells from the same patient are correlated, not independent replicates; a
# per-cell Wilcoxon test across the pooled cohort (the old approach) is
# pseudoreplication. Instead, SUM raw counts per (patient, group) into
# pseudobulk samples and run DESeq2 with a patient-blocked design, so each
# patient contributes one paired vote per group.

# Reconstruct the malignant-cell sparse count matrix (also reused verbatim by
# the cache-artifacts section below).
counts_malignant_sparse <- counts_full[, malignant_idx]
colnames(counts_malignant_sparse) <- obs_df$cell_id[malignant_idx]

hl_meta <- mal_meta %>% filter(cytotrace_group %in% c("High", "Low"))

patient_hl_counts <- hl_meta %>%
  count(batch, cytotrace_group) %>%
  tidyr::pivot_wider(names_from = cytotrace_group, values_from = n, values_fill = 0)
if (!"High" %in% colnames(patient_hl_counts)) patient_hl_counts$High <- 0
if (!"Low"  %in% colnames(patient_hl_counts)) patient_hl_counts$Low  <- 0

MIN_CELLS_PER_GROUP <- 10
patient_hl_counts <- patient_hl_counts %>%
  mutate(pass = High >= MIN_CELLS_PER_GROUP & Low >= MIN_CELLS_PER_GROUP)

dropped_patients_pb <- patient_hl_counts %>% filter(!pass) %>% pull(batch)
kept_patients_pb    <- patient_hl_counts %>% filter(pass)  %>% pull(batch)

if (length(dropped_patients_pb) > 0) {
  message("  Dropping ", length(dropped_patients_pb),
          " patient(s) with <", MIN_CELLS_PER_GROUP,
          " malignant cells in local-High or local-Low: ",
          paste(dropped_patients_pb, collapse = ", "))
}
message("  Patients retained for pseudobulk DESeq2: ", length(kept_patients_pb),
        " / ", nrow(patient_hl_counts))

if (length(kept_patients_pb) < 2) {
  warning("Fewer than 2 patients have >= ", MIN_CELLS_PER_GROUP,
          " malignant cells in BOTH local High and Low groups; a patient-blocked ",
          "DESeq2 design (~ patient + group) cannot be fit (no cross-patient ",
          "replication). Skipping pseudobulk DE step -- no CSV/volcano plot written.")
} else {
  hl_meta_kept <- hl_meta %>% filter(batch %in% kept_patients_pb)
  hl_meta_kept$pb_sample <- paste(hl_meta_kept$batch, hl_meta_kept$cytotrace_group, sep = "_")
  pb_samples <- sort(unique(hl_meta_kept$pb_sample))

  # Sum raw counts per (patient, group) via a sparse cells x samples indicator
  # matrix multiply (same aggregation idiom as the sibling PyDESeq2 notebook's
  # pseudobulk_by_sample(): sparse indicator matrix @ counts).
  cell_ids_hl <- hl_meta_kept$cell_id
  counts_hl <- counts_malignant_sparse[, cell_ids_hl, drop = FALSE]  # genes x cells

  ind_mat <- Matrix::sparseMatrix(
    i = seq_along(cell_ids_hl),
    j = match(hl_meta_kept$pb_sample, pb_samples),
    x = 1,
    dims = c(length(cell_ids_hl), length(pb_samples)),
    dimnames = list(cell_ids_hl, pb_samples)
  )
  pb_counts <- as.matrix(counts_hl %*% ind_mat)   # genes x pseudobulk-samples
  colnames(pb_counts) <- pb_samples
  storage.mode(pb_counts) <- "integer"

  pb_coldata <- data.frame(
    row.names = pb_samples,
    patient   = factor(sub("_(High|Low)$", "", pb_samples)),
    group     = factor(sub("^.*_(High|Low)$", "\\1", pb_samples), levels = c("Low", "High"))
  )
  # "Low" is the explicit reference level; the contrast is also passed
  # explicitly to results() below so this choice is documented in both places.

  message("  Pseudobulk design: ", nrow(pb_coldata), " samples (",
          sum(pb_coldata$group == "High"), " High, ", sum(pb_coldata$group == "Low"), " Low) from ",
          length(unique(pb_coldata$patient)), " patients")
  # Full-rank note: by construction every kept patient contributes exactly one
  # High and one Low pseudobulk sample (a balanced paired design, analogous to
  # a paired t-test). With k patients, design ~patient + group has k + 1 free
  # parameters estimated from n = 2k samples, which is full rank whenever
  # k >= 1; the guard above already requires k >= 2 (k = 1 would confound
  # patient and group entirely).

  # Gene-level low-count filter (no equivalent existed for the old FindMarkers
  # path; mirrors the documented-but-commented-out convention in the sibling
  # PyDESeq2 notebook: total count >= 10 AND detected/nonzero in >= 2 samples).
  gene_keep <- (Matrix::rowSums(pb_counts) >= 10) & (Matrix::rowSums(pb_counts > 0) >= 2)
  message("  Genes passing low-count filter: ", sum(gene_keep), " / ", nrow(pb_counts))
  pb_counts <- pb_counts[gene_keep, , drop = FALSE]

  de_res_pb <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = pb_counts, colData = pb_coldata,
                                           design = ~ patient + group)
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, contrast = c("group", "High", "Low"))
    as.data.frame(res)
  }, error = function(e) {
    message("Pseudobulk DESeq2 DE failed: ", conditionMessage(e))
    data.frame()
  })

  if (nrow(de_res_pb) > 0) {
    de_res_pb$gene <- rownames(de_res_pb)
    de_res_pb <- de_res_pb %>% arrange(padj)
    write.csv(de_res_pb, file.path(OUT_DIR, "markers", "de_cytotrace_high_vs_low_pseudobulk_deseq2.csv"),
              row.names = FALSE)

    de_res_pb_plot <- de_res_pb %>%
      filter(!is.na(padj)) %>%
      mutate(neglog10padj = -log10(pmax(padj, .Machine$double.xmin)))
    top_genes_pb <- de_res_pb_plot %>% arrange(desc(abs(log2FoldChange))) %>% head(15)

    p_volcano_pb <- ggplot(de_res_pb_plot, aes(x = log2FoldChange, y = neglog10padj)) +
      ggrastr::rasterise(geom_point(alpha = 0.5, size = 1), dpi = 600) +
      ggrepel::geom_text_repel(data = top_genes_pb, aes(label = gene), size = 3, max.overlaps = 30) +
      labs(title = "CytoTRACE-high vs -low malignant cells: pseudobulk DESeq2 DE genes",
           subtitle = paste0(nrow(pb_coldata), " pseudobulk samples from ",
                              length(unique(pb_coldata$patient)), " patients"),
           x = "log2FoldChange (High vs Low)", y = "-log10(BH-adjusted p-value)") +
      theme_cowplot(font_size = 12)
    save_plot(p_volcano_pb, "de_cytotrace_high_vs_low_pseudobulk_deseq2_volcano", "markers", width = 8, height = 7)
  } else {
    message("WARNING: pseudobulk DESeq2 DE step produced no results; skipping volcano plot and CSV.")
  }
}

## ---- 16. By-cell-type figures (run B) ----------------------------------------
message("=== [12/14] By-cell-type figures (run B) ===")

n_types_b <- length(unique(all_meta[[ANNOT_COL]]))
type_pal_b <- viridis::viridis(n_types_b, option = "turbo")

celltype_order <- all_meta %>% group_by(.data[[ANNOT_COL]]) %>%
  summarise(med = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
  arrange(med) %>% pull(1)
all_meta[[ANNOT_COL]] <- factor(all_meta[[ANNOT_COL]], levels = celltype_order)

p_by_type <- ggplot(all_meta, aes(x = .data[[ANNOT_COL]], y = CytoTRACE2_Score,
                                   fill = .data[[ANNOT_COL]])) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8) +
  scale_fill_manual(values = type_pal_b) +
  coord_flip() +
  labs(title = paste0("CytoTRACE2 score by cell type (", ANNOT_COL, ")"),
       x = NULL, y = "CytoTRACE2 score") +
  guides(fill = "none") +
  theme_cowplot(font_size = 11)
save_plot(p_by_type, "cytotrace_score_by_annot5", "by_celltype", width = 10, height = 9)

compartment_order <- all_meta %>% group_by(compartment) %>% summarise(med = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
  arrange(med) %>% pull(compartment)
all_meta$compartment <- factor(all_meta$compartment, levels = compartment_order)
n_comp <- length(unique(all_meta$compartment))
comp_pal <- viridis::viridis(n_comp, option = "turbo")

p_by_compartment <- ggplot(all_meta, aes(x = compartment, y = CytoTRACE2_Score, fill = compartment)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8) +
  scale_fill_manual(values = comp_pal) +
  coord_flip() +
  labs(title = "CytoTRACE2 score by compartment", x = NULL, y = "CytoTRACE2 score") +
  guides(fill = "none") +
  theme_cowplot(font_size = 11)
save_plot(p_by_compartment, "cytotrace_score_by_compartment", "by_celltype", width = 7, height = 6)

all_meta_out <- all_meta %>% select(cell_id, batch, outcome_group,
                                     any_of(c(ANNOT_COL, ANNOT_COL_LEGACY)),
                                     compartment, CytoTRACE2_Score, CytoTRACE2_Potency)
write.csv(all_meta_out, file.path(OUT_DIR, "by_celltype", "cytotrace2_scores_all_celltypes.csv"), row.names = FALSE)

## ---- 17. Cache artifacts for downstream Monocle3 script ----------------------
message("=== [13/14] Writing cache artifacts for downstream Monocle3 script ===")

# Carry the resolved `outcome_group` AND every candidate source column that
# exists, so the cache stays usable if OUTCOME_COL is later switched without
# paying for a full CytoTRACE2 re-run.
cache_meta <- mal_meta %>%
  select(cell_id, batch, outcome_group,
         any_of(c("clinical_outcome", "model_outcome", "outcome")),
         any_of(c(ANNOT_COL, ANNOT_COL_LEGACY)), phase, S_score, G2M_score,
         cytotrace_group, n_genes, total_counts, pct_counts_mt, UMAP_1, UMAP_2) %>%
  rename(X_umap_1 = UMAP_1, X_umap_2 = UMAP_2) %>%
  as.data.frame()
rownames(cache_meta) <- cache_meta$cell_id
cache_meta$cell_id <- NULL

# Reuse the sparse malignant-cell matrix built for the Marker DE step above
# (already genes x malignant-cells, correctly named) rather than re-densifying.
counts_dgc <- as(counts_malignant_sparse, "CsparseMatrix")

saveRDS(list(counts = counts_dgc, meta = cache_meta),
        file = file.path(CACHE_DIR, "malignant_subset_counts_meta.rds"))

# Per-cell CytoTRACE2 run A scores. Row identifier = cell barcode, written as an
# explicit `cell_id` column (row.names = FALSE) for robust downstream parsing.
cyto_out <- mal_meta %>% select(cell_id, CytoTRACE2_Score, CytoTRACE2_Potency, cytotrace_group,
                                 any_of(setdiff(colnames(scores_malignant), c("CytoTRACE2_Score", "CytoTRACE2_Potency"))))
write.csv(cyto_out, file.path(OUT_DIR, "tables", "cytotrace2_scores_malignant.csv"), row.names = FALSE)

# Provenance: record which outcome column produced these results, so a set of
# output files can never be misattributed to the wrong grouping after the fact.
run_config <- data.frame(
  key = c("outcome_col", "smoke_test", "n_patients", "n_patients_good",
          "n_patients_poor", "n_malignant_cells_scored", "run_timestamp"),
  value = c(OUTCOME_COL, as.character(SMOKE_TEST),
            as.character(nrow(.patient_outcome)),
            as.character(sum(.patient_outcome$outcome_group == "Good", na.rm = TRUE)),
            as.character(sum(.patient_outcome$outcome_group == "Poor", na.rm = TRUE)),
            as.character(nrow(mal_meta)),
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  stringsAsFactors = FALSE
)
write.csv(run_config, file.path(OUT_DIR, "tables", "run_config.csv"), row.names = FALSE)

message("=== [14/14] DONE ===")
message("Outcome column: ", OUTCOME_COL, " (",
        sum(.patient_outcome$outcome_group == "Good", na.rm = TRUE), " Good / ",
        sum(.patient_outcome$outcome_group == "Poor", na.rm = TRUE), " Poor patients)")
message("SMOKE_TEST = ", SMOKE_TEST)
message("Malignant cells scored (run A): ", nrow(mal_meta))
message("All QC-passed cells scored (run B): ", nrow(all_meta))
message("Output root: ", OUT_DIR)
