#!/usr/bin/env Rscript
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# cytotrace2_cohort.R
#
# CytoTRACE2 stemness/potency scoring of the MALIGNANT (tumour) compartment of
# an adenoid cystic carcinoma scRNA-seq cohort, run SEPARATELY INSIDE ONE
# OUTCOME COHORT (Good or Poor) rather than across the pooled 24-patient set.
#
#   Usage:  Rscript cytotrace2_cohort.R good
#           Rscript cytotrace2_cohort.R poor
#
# WHY THIS EXISTS
#   The pooled analysis (../cytotrace2_malignant_analysis.R ->
#   ../monocle3_malignant_trajectory.R) scores and embeds all 71,984 malignant
#   cells together. This script instead re-fits CytoTRACE2 *within* one cohort,
#   so the scores answer "where does this cell sit within its own cohort" and
#   the downstream Monocle3 trajectory is built on that cohort's cells alone.
#   The pooled scripts and all pooled outputs are untouched by this run.
#
# INPUT: the shared malignant-cell cache written by the pooled stage 1,
#   outputs/h5ad_files/cache/malignant_subset_counts_meta.rds
#   It already holds exactly the QC-passed malignant cells, so this script
#   never reads the 8 GB h5ad. It is opened READ-ONLY and never rewritten.
#
# OUTPUT: outputs/outcome_split/<cohort>/ ...
#   Every basename carries a _good / _poor suffix. Directory separation alone
#   does not satisfy the no-duplicate-basename rule (see ANALYSIS_NOTES.md S7):
#   that rule is about basenames, because bundling two trees for review makes
#   one silently shadow the other.
#
# DIFFERENCES FROM THE POOLED STAGE 1 (source line ranges refer to it):
#   dropped  262-338  h5ad load            -> replaced by a cache subset
#   dropped  340-388  QC / malignant filter-> the cache is already filtered
#   dropped  390-482  all-cell-type QC figs-> need cell types the cache lacks
#   dropped  571-601  CytoTRACE2 run B     -> removed from scope
#   dropped  673-753  Wilcoxon / Cliff's / LOO and its one-box figures
#                       -> degenerate inside a single cohort: wilcox.test()
#                          errors on a single-level grouping factor
#   dropped  755-905  pseudobulk marker DE -> out of scope for this subproject
#   dropped  906-943  by-cell-type figures (run B)
#   carried  42-76, 108-226, 484-517, 519-569, 603-626, 628-671, 945-986
# =============================================================================

## ---- 0. Cohort argument ------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) < 1 || !.args[1] %in% c("good", "poor")) {
  stop("Usage: Rscript cytotrace2_cohort.R <good|poor>\n",
       "  The cohort MUST be given explicitly -- there is no default, so a\n",
       "  mistyped or missing argument can never silently write one cohort's\n",
       "  results into the other's tree. It also makes the running cohort\n",
       "  visible in `ps -eo args`.")
}
COHORT       <- .args[1]                                   # "good" / "poor"
COHORT_LABEL <- c(good = "Good", poor = "Poor")[[COHORT]]  # value in the meta
SUFFIX       <- paste0("_", COHORT)

# Expected cohort size, verified against the pooled per-cell export. A run that
# lands far from this is reading the wrong cache or the wrong column.
EXPECTED_CELLS    <- c(good = 32991L, poor = 38993L)[[COHORT]]
EXPECTED_PATIENTS <- c(good = 13L,    poor = 11L)[[COHORT]]

message("=== CytoTRACE2, ", COHORT_LABEL, "-outcome cohort (tumour cells only) ===")

## ---- 1. Package install (idempotent) ----------------------------------------
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
  remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r",
                          dependencies = TRUE, upgrade = "never")
}
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix)
  library(CytoTRACE2); library(ggplot2); library(patchwork)
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

# Read-only. Written by the pooled stage 1; this script must never write here.
POOLED_CACHE <- file.path(BASE_DIR, "outputs", "h5ad_files", "cache",
                          "malignant_subset_counts_meta.rds")

OUT_DIR   <- file.path(BASE_DIR, "outputs", "outcome_split", COHORT,
                       "cytotrace_analysis")
CACHE_DIR <- file.path(BASE_DIR, "outputs", "outcome_split", COHORT, "cache")

for (d in c(OUT_DIR,
            file.path(OUT_DIR, "qc"),
            file.path(OUT_DIR, "umap"),
            file.path(OUT_DIR, "patient_level"),
            file.path(OUT_DIR, "tables"),
            CACHE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ---- 3. Plotting / table helpers --------------------------------------------
# The cohort suffix is applied HERE, once, rather than at ~30 call sites, so
# every call site below can be copied verbatim from the pooled script and still
# produce a collision-free basename.
save_plot <- function(p, name, subdir = "", width = 8, height = 6, dpi = 600) {
  name    <- paste0(name, SUFFIX)
  out_dir <- file.path(OUT_DIR, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), plot = p, width = width, height = height,
         dpi = dpi, units = "in", bg = "white", device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot = p, width = width, height = height,
         units = "in", bg = "white", device = svglite::svglite)
  message("  Saved: ", subdir, "/", name)
}

write_tab <- function(df, name, subdir = "") {
  out_dir <- file.path(OUT_DIR, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fname <- paste0(name, SUFFIX, ".csv")
  write.csv(df, file.path(out_dir, fname), row.names = FALSE)
  message("  Saved: ", subdir, "/", fname)
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
# IMPORTANT: identical copies of this block now live in FOUR files:
#   ../cytotrace2_malignant_analysis.R      (section 3b)
#   ../monocle3_malignant_trajectory.R      (section 7a)
#   ./cytotrace2_cohort.R                   (this block)
#   ./monocle3_cohort.R                     (section 7a)
# Each script is deliberately standalone -- the pipeline runs them as separate
# Rscript invocations -- so edit ALL FOUR or the panels silently drift apart
# again. OUTCOME_SPLIT_NOTES.md carries the awk-range + diff recipe.
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

OUTCOME_COLORS <- c("Good" = "forestgreen", "Poor" = "firebrick3")
COHORT_COLOR   <- OUTCOME_COLORS[[COHORT_LABEL]]

# ---- Annotation level -------------------------------------------------------
# annot_5 is the annotation of record. It is INHERITED from the pooled cache --
# this script never reads the h5ad, so RECODE_ANNOT5 has already been applied
# upstream in ../cytotrace2_malignant_analysis.R section 4a. annot_6 rides along
# as a passenger column for crosswalking to the superseded run.
ANNOT_COL        <- "annot_5"
ANNOT_COL_LEGACY <- "annot_6"
MALIGNANT_TYPES <- c(
  "Myoepithelial cells - Tumor",
  "Epithelial cells - Tumor",
  "Actively Dividing cells - Tumor",
  "Actively Dividing Myoepithelial cells - Tumor",
  "Epithelial cells - Basal - Tumor",
  "Basal cells - Tumor"
)
# Display order along the myoepithelial -> ductal axis. Duplicated verbatim in
# monocle3_cohort.R and both pooled scripts -- edit all four.
MALIGNANT_LEVELS <- c(
  "Myoepithelial cells - Tumor",
  "Actively Dividing Myoepithelial cells - Tumor",
  "Basal cells - Tumor",
  "Epithelial cells - Basal - Tumor",
  "Actively Dividing cells - Tumor",
  "Epithelial cells - Tumor"
)
# Poor has no myoepithelial reserve in these patients: with <= 60 myoepithelial
# tumour cells, align_cds(alignment_k = 20) cannot find honest mutual nearest
# neighbours for them in that compartment and will manufacture some. Recorded,
# not corrected -- see OUTCOME_SPLIT_NOTES.md section 8.
MYO_SPARSE_MAX_CELLS <- 60L

# =============================================================================
# 4. Load the cohort subset from the shared malignant cache
# =============================================================================
message("=== [1/7] Loading cohort subset from the malignant cache ===")

if (!file.exists(POOLED_CACHE))
  stop("Malignant cache not found: ", POOLED_CACHE,
       "\nRun the pooled ../cytotrace2_malignant_analysis.R first -- this script ",
       "subsets its cache and never reads the h5ad itself.")

.cached <- readRDS(POOLED_CACHE)
counts_full <- .cached$counts
meta_full   <- .cached$meta
rm(.cached); gc()

# Guard the PARENT cache at the pooled constant, where it is actually correct:
# this cache IS the full malignant compartment, so a smoke-sized one means a
# companion script clobbered the shared artifact. That is a bug to fix at the
# source, not something to work around by re-reading the 8 GB h5ad here.
EXPECTED_MIN_MALIGNANT_CELLS_FULL <- 50000
if (ncol(counts_full) < EXPECTED_MIN_MALIGNANT_CELLS_FULL)
  stop("The shared malignant cache holds only ", ncol(counts_full), " cells, well ",
       "below the expected full compartment (~71,984). This looks like a ",
       "smoke-sized subset written by a companion script sharing this path. ",
       "Refusing to derive a cohort from it -- regenerate the cache instead.")

stopifnot(identical(colnames(counts_full), rownames(meta_full)))
message("  Parent cache: ", ncol(counts_full), " malignant cells x ",
        nrow(counts_full), " genes")

if (!"outcome_group" %in% colnames(meta_full))
  stop("The cache meta carries no `outcome_group` column. Columns present: ",
       paste(colnames(meta_full), collapse = ", "))

keep <- which(meta_full$outcome_group == COHORT_LABEL)
if (length(keep) == 0)
  stop("No cells with outcome_group == '", COHORT_LABEL, "'. Values present: ",
       paste(sort(unique(meta_full$outcome_group)), collapse = ", "))

counts_cohort <- counts_full[, keep, drop = FALSE]
mal_meta      <- meta_full[keep, , drop = FALSE]
rm(counts_full, meta_full); gc()

# The cache deliberately stores the cell barcode as the ROWNAME and drops the
# `cell_id` column (pooled stage 1, lines 951-960). Restore it immediately:
# every downstream block here indexes by mal_meta$cell_id, and without this the
# patient table, the score join and the cache write all fail or silently
# misalign.
mal_meta$cell_id <- rownames(mal_meta)

# Factor columns (batch, annot_6, ...) are levelled on the FULL 24-patient
# dataset. After subsetting to one cohort the other cohort's levels survive with
# zero cells. That is not cosmetic: the cache written at the end of this script
# feeds scran::modelGeneVar(block = batch) in stage 2, which splits by block and
# would be handed empty blocks. Drop unused levels once, here, so the cohort
# cache is clean at the source rather than patched downstream.
mal_meta <- droplevels(mal_meta)

stopifnot(identical(colnames(counts_cohort), rownames(mal_meta)))
stopifnot(identical(colnames(counts_cohort), mal_meta$cell_id))

## ---- 4a. Annotation gate: annot_5 must have come through the pooled cache ---
if (!ANNOT_COL %in% colnames(mal_meta))
  stop("The pooled cache carries no '", ANNOT_COL, "' column -- it predates the ",
       "annot_5 migration. Re-run ../cytotrace2_malignant_analysis.R (full stage 1, ",
       "NOT --from=2) before running this script. Columns present: ",
       paste(colnames(mal_meta), collapse = ", "))

mal_meta[[ANNOT_COL]] <- as.character(mal_meta[[ANNOT_COL]])
.stray <- setdiff(unique(mal_meta[[ANNOT_COL]]), MALIGNANT_TYPES)
if (length(.stray) > 0)
  stop("Non-malignant ", ANNOT_COL, " value(s) in the cohort cache: ",
       paste(.stray, collapse = ", "),
       ". Either the recode did not run upstream or the cache is stale.")

# Fixed level order, so every figure orders the classes along the axis rather
# than alphabetically. intersect() keeps only the classes this cohort actually
# has -- Poor has no `Basal cells - Tumor` at all (276 cells, all Good).
mal_meta[[ANNOT_COL]] <- factor(mal_meta[[ANNOT_COL]],
                                levels = intersect(MALIGNANT_LEVELS,
                                                   unique(mal_meta[[ANNOT_COL]])))
.a5_tab <- table(mal_meta[[ANNOT_COL]])
message("Annotation level ", ANNOT_COL, ": ", nlevels(mal_meta[[ANNOT_COL]]),
        " malignant class(es) present in the ", COHORT_LABEL, " cohort")
print(.a5_tab)

# Which patients are myoepithelial-sparse? This is the population align_cds has
# nothing honest to match, so it belongs in the run log and the run config.
.myo_by_pt <- table(as.character(mal_meta$batch)[
                      mal_meta[[ANNOT_COL]] == "Myoepithelial cells - Tumor"])
.all_pt <- unique(as.character(mal_meta$batch))
.myo_n  <- setNames(rep(0L, length(.all_pt)), .all_pt)
.myo_n[names(.myo_by_pt)] <- as.integer(.myo_by_pt)
MYO_SPARSE_PATIENTS <- sort(names(.myo_n)[.myo_n <= MYO_SPARSE_MAX_CELLS])
if (length(MYO_SPARSE_PATIENTS) > 0) {
  message("MYO_SPARSE: ", length(MYO_SPARSE_PATIENTS), " patient(s) with <= ",
          MYO_SPARSE_MAX_CELLS, " myoepithelial tumour cells: ",
          paste0(MYO_SPARSE_PATIENTS, " (n=", .myo_n[MYO_SPARSE_PATIENTS], ")",
                 collapse = ", "),
          " -- MNN alignment has no honest neighbours for these in that compartment.")
} else {
  message("MYO_SPARSE: none (every patient has > ", MYO_SPARSE_MAX_CELLS,
          " myoepithelial tumour cells).")
}

# Exactly one outcome level, and the resolved column still agrees with its
# source -- a cohort that quietly contains both groups would invalidate the
# whole point of the split.
if (length(unique(mal_meta$outcome_group)) != 1)
  stop("Cohort subset contains ", length(unique(mal_meta$outcome_group)),
       " outcome levels; expected exactly 1.")
if ("clinical_outcome" %in% colnames(mal_meta) &&
    !all(mal_meta$outcome_group == mal_meta$clinical_outcome))
  stop("`outcome_group` and `clinical_outcome` disagree in the cohort subset; ",
       "the cache was written with a different OUTCOME_COL than assumed here.")

n_cells    <- ncol(counts_cohort)
n_patients <- length(unique(mal_meta$batch))
message("  ", COHORT_LABEL, " cohort: ", n_cells, " tumour cells, ",
        n_patients, " patients")
.per_pt <- sort(table(mal_meta$batch))
message("  Cells per patient: ",
        paste(names(.per_pt), as.integer(.per_pt), sep = "=", collapse = ", "))

.chk <- function(actual, expected, what) {
  if (actual == expected) return(invisible(NULL))
  rel <- abs(actual - expected) / expected
  msg <- paste0(what, " = ", actual, ", expected ", expected,
                " (", round(100 * rel, 1), "% off)")
  if (rel > 0.05)
    stop("Cohort composition is far from expectation: ", msg,
         ". Refusing to run -- check the cache and OUTCOME_COL before ",
         "spending an hour of compute on the wrong cells.")
  warning("Cohort composition differs from expectation: ", msg)
}
.chk(n_cells,    EXPECTED_CELLS,    "cohort cell count")
.chk(n_patients, EXPECTED_PATIENTS, "cohort patient count")

# The cache carries a `cytotrace_group` computed from the POOLED run's scores.
# It is stale the moment we re-fit, and silently wrong rather than missing, so
# drop it now and recompute from this cohort's own scores further down.
mal_meta$cytotrace_group <- NULL

# =============================================================================
# 5. CytoTRACE2 on this cohort's tumour cells
# =============================================================================
message("=== [2/7] CytoTRACE2: ", COHORT_LABEL, " cohort, tumour cells only ===")

# NOTE on ncores: CytoTRACE2's internal parallelism (parallelize_models /
# parallelize_smoothing) uses parallel::mclapply-style forking. Under R's
# copy-on-write + refcounting GC, forked workers touching a multi-GB shared
# matrix during GC can each end up privately copying large chunks of it,
# so peak RSS can grow much faster than a naive "N workers x per-worker-size"
# estimate. On this 60GB-RAM/8GB-swap machine, ncores=20 on the ~95k-cell
# pooled run OOM-killed the process. Capped at 6 here, identical to the pooled
# run: every CytoTRACE2 parameter below is deliberately unchanged, so the only
# difference between this run and the pooled one is the input cell set.
N_CORES <- min(6, parallel::detectCores())

counts_dense <- as.matrix(counts_cohort)

scores_cohort <- CytoTRACE2::cytotrace2(
  input = counts_dense, species = "human", is_seurat = FALSE, slot_type = "counts",
  batch_size = 10000, smooth_batch_size = 1000,
  parallelize_models = TRUE, parallelize_smoothing = TRUE, ncores = N_CORES, seed = 14
)

mal_meta <- cbind(mal_meta, scores_cohort[rownames(mal_meta), , drop = FALSE])
message("CytoTRACE2 complete: ", nrow(mal_meta), " ", COHORT_LABEL,
        " tumour cells scored")
message("  Score range: ", round(min(mal_meta$CytoTRACE2_Score, na.rm = TRUE), 4),
        " to ", round(max(mal_meta$CytoTRACE2_Score, na.rm = TRUE), 4),
        ", median ", round(median(mal_meta$CytoTRACE2_Score, na.rm = TRUE), 4))
message("  Potency levels: ",
        paste(names(table(mal_meta$CytoTRACE2_Potency)),
              table(mal_meta$CytoTRACE2_Potency), sep = "=", collapse = ", "))

rm(counts_dense); gc()

# =============================================================================
# 6. Cell-cycle independence diagnostic
# =============================================================================
message("=== [3/7] Cell-cycle independence diagnostic (CytoTRACE2 vs S/G2M) ===")

# CytoTRACE2 is scored on raw counts, so it never sees the scVI batch
# correction or cell-cycle regression applied upstream in Scanpy. Before
# trusting the rooting that stage 2 builds on this score, check whether the
# stemness signal is simply tracking proliferation.
cc_cor_S   <- suppressWarnings(cor(mal_meta$CytoTRACE2_Score, mal_meta$S_score,
                                    use = "pairwise.complete.obs", method = "pearson"))
cc_cor_G2M <- suppressWarnings(cor(mal_meta$CytoTRACE2_Score, mal_meta$G2M_score,
                                    use = "pairwise.complete.obs", method = "pearson"))

write_tab(data.frame(
  cohort    = COHORT_LABEL,
  metric    = c("CytoTRACE2_Score_vs_S_score", "CytoTRACE2_Score_vs_G2M_score"),
  pearson_r = c(cc_cor_S, cc_cor_G2M),
  n_cells   = nrow(mal_meta)
), "cytotrace_cellcycle_correlation", "qc")

message("  Pearson r(CytoTRACE2_Score, S_score)   = ", round(cc_cor_S, 4))
message("  Pearson r(CytoTRACE2_Score, G2M_score) = ", round(cc_cor_G2M, 4))
if (!is.na(cc_cor_S) && abs(cc_cor_S) > 0.5)
  warning("CytoTRACE2_Score correlates strongly with S_score (r = ", round(cc_cor_S, 3),
          "); stemness signal may partly reflect proliferation rather than true ",
          "differentiation state.")
if (!is.na(cc_cor_G2M) && abs(cc_cor_G2M) > 0.5)
  warning("CytoTRACE2_Score correlates strongly with G2M_score (r = ", round(cc_cor_G2M, 3),
          "); stemness signal may partly reflect proliferation rather than true ",
          "differentiation state.")

p_cc_S <- ggplot(mal_meta, aes(x = S_score, y = CytoTRACE2_Score)) +
  ggrastr::rasterise(geom_point(size = 0.4, alpha = 0.5, color = "steelblue"), dpi = 600) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  annotate("text", x = min(mal_meta$S_score, na.rm = TRUE),
           y = max(mal_meta$CytoTRACE2_Score, na.rm = TRUE),
           label = paste0("r = ", round(cc_cor_S, 3)), hjust = 0, vjust = 1, size = 4.5) +
  labs(title = paste0("CytoTRACE2 score vs S-phase score (", COHORT_LABEL,
                      " cohort, tumour cells)"),
       x = "S_score", y = "CytoTRACE2_Score") +
  theme_cowplot(font_size = 12)
save_plot(p_cc_S, "cytotrace_vs_S_score", "qc", width = 7, height = 6)

p_cc_G2M <- ggplot(mal_meta, aes(x = G2M_score, y = CytoTRACE2_Score)) +
  ggrastr::rasterise(geom_point(size = 0.4, alpha = 0.5, color = "darkorange"), dpi = 600) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  annotate("text", x = min(mal_meta$G2M_score, na.rm = TRUE),
           y = max(mal_meta$CytoTRACE2_Score, na.rm = TRUE),
           label = paste0("r = ", round(cc_cor_G2M, 3)), hjust = 0, vjust = 1, size = 4.5) +
  labs(title = paste0("CytoTRACE2 score vs G2M-phase score (", COHORT_LABEL,
                      " cohort, tumour cells)"),
       x = "G2M_score", y = "CytoTRACE2_Score") +
  theme_cowplot(font_size = 12)
save_plot(p_cc_G2M, "cytotrace_vs_G2M_score", "qc", width = 7, height = 6)

# =============================================================================
# 7. UMAP figures on the original scVI coordinates
# =============================================================================
message("=== [4/7] UMAP figures (CytoTRACE2 score / potency) ===")

# The cache stores the scVI UMAP as X_umap_1 / X_umap_2 (pooled stage 1 renames
# UMAP_1/UMAP_2 on the way in, line 956). These are the FULL-COHORT coordinates
# computed once upstream, simply restricted to this cohort's cells -- they are
# not recomputed here, so the two cohorts' panels are in the same coordinate
# space and can be laid side by side.
stopifnot(all(c("X_umap_1", "X_umap_2") %in% colnames(mal_meta)))

p_umap_score <- umap_panel(mal_meta, "CytoTRACE2_Score",
                           paste0("CytoTRACE2 score - ", COHORT_LABEL,
                                  " cohort tumour cells"),
                           "CytoTRACE2 Score") +
  scale_color_viridis_c(option = "viridis")
save_umap_panel(p_umap_score, "umap_cytotrace_score", "umap")

p_umap_potency <- umap_panel(mal_meta, "CytoTRACE2_Potency",
                             paste0("CytoTRACE2 potency - ", COHORT_LABEL,
                                    " cohort tumour cells"),
                             "CytoTRACE2 Potency") +
  scale_color_viridis_d(option = "turbo")
save_umap_panel(p_umap_potency, "umap_cytotrace_potency", "umap")

# =============================================================================
# 8. Patient-level summary + within-cohort patient figures
# =============================================================================
message("=== [5/7] Patient-level summary table ===")

patient_counts <- mal_meta %>% count(batch, name = "n_cells_malignant")
low_n_patients <- patient_counts %>% filter(n_cells_malignant < 10) %>% pull(batch)
if (length(low_n_patients) > 0)
  warning("Dropping patients with <10 tumour cells: ",
          paste(low_n_patients, collapse = ", "))
keep_patients <- patient_counts %>% filter(n_cells_malignant >= 10) %>% pull(batch)

mal_meta_kept <- mal_meta %>% filter(batch %in% keep_patients)

# WITHIN-COHORT quantile thresholds: pooled across this cohort's cells only.
# These are NOT comparable to the pooled run's thresholds or to the other
# cohort's -- by construction every cohort puts 10% of its own cells above its
# own top-10% cut. See OUTCOME_SPLIT_NOTES.md.
q_top10 <- quantile(mal_meta_kept$CytoTRACE2_Score, probs = 0.90,   na.rm = TRUE)
q_top25 <- quantile(mal_meta_kept$CytoTRACE2_Score, probs = 0.75,   na.rm = TRUE)
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
    frac_G1  = mean(phase == "G1",  na.rm = TRUE),
    frac_S   = mean(phase == "S",   na.rm = TRUE),
    frac_G2M = mean(phase == "G2M", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_cytotrace_score))

write_tab(patient_level, "patient_level_summary_cytotrace2", "patient_level")
message("  ", nrow(patient_level), " patients in the ", COHORT_LABEL, " cohort")

# Ranked patient medians. Coloured by the cohort's own colour so the Good and
# Poor twins of this figure are distinguishable at a glance.
p_rank <- ggplot(patient_level,
                 aes(x = reorder(batch, median_cytotrace_score),
                     y = median_cytotrace_score)) +
  geom_segment(aes(xend = batch, y = min(patient_level$median_cytotrace_score),
                   yend = median_cytotrace_score),
               color = "grey70", linewidth = 0.5) +
  geom_point(size = 3.2, color = COHORT_COLOR) +
  coord_flip() +
  labs(title = paste0("Median CytoTRACE2 score per patient - ", COHORT_LABEL,
                      " cohort"),
       subtitle = "Within-cohort ranking only; scores are re-fit inside this cohort and do not carry a cross-cohort level",
       x = "Patient (batch)", y = "Median CytoTRACE2 score") +
  theme_bw(base_size = 12)
save_plot(p_rank, "patient_median_cytotrace_ranked", "patient_level",
          width = 7.5, height = max(4, 0.35 * nrow(patient_level)))

# Cell-level distribution per patient, ordered by that patient's median.
.ord <- patient_level$batch[order(patient_level$median_cytotrace_score)]
p_viol <- ggplot(mal_meta_kept %>% mutate(batch = factor(batch, levels = .ord)),
                 aes(x = batch, y = CytoTRACE2_Score)) +
  geom_violin(fill = COHORT_COLOR, alpha = 0.45, color = NA, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  coord_flip() +
  labs(title = paste0("Per-cell CytoTRACE2 score by patient - ", COHORT_LABEL,
                      " cohort"),
       x = "Patient (batch)", y = "CytoTRACE2 score") +
  theme_bw(base_size = 12)
save_plot(p_viol, "patient_cytotrace_violin", "patient_level",
          width = 7.5, height = max(4.5, 0.4 * nrow(patient_level)))

# =============================================================================
# 9. Recompute cytotrace_group from THIS cohort's scores
# =============================================================================
message("=== [6/7] Recomputing per-patient cytotrace_group ===")

# Same rule as the pooled stage 1 (lines 762-781): each cell's High/Low/Mid
# label is relative to its OWN patient's score distribution. Recomputed rather
# than inherited -- the cache's copy was derived from the pooled fit and would
# be silently wrong here, which is worse than absent.
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
stopifnot(identical(colnames(counts_cohort), rownames(mal_meta)))

message("  Per-patient local cytotrace_group counts:")
print(table(mal_meta$batch, mal_meta$cytotrace_group))

# =============================================================================
# 10. Cache artifacts for the cohort Monocle3 stage
# =============================================================================
message("=== [7/7] Writing cohort cache artifacts ===")

# Same schema as the pooled cache (cell barcode as rowname, no cell_id column),
# so monocle3_cohort.R's load path is identical to the pooled script's.
cache_meta <- mal_meta %>%
  select(cell_id, batch, outcome_group,
         any_of(c("clinical_outcome", "model_outcome", "outcome")),
         any_of(c(ANNOT_COL, ANNOT_COL_LEGACY)), phase, S_score, G2M_score,
         cytotrace_group, n_genes, total_counts, pct_counts_mt,
         X_umap_1, X_umap_2) %>%
  as.data.frame()
rownames(cache_meta) <- cache_meta$cell_id
cache_meta$cell_id <- NULL

counts_dgc <- as(counts_cohort, "CsparseMatrix")
stopifnot(identical(colnames(counts_dgc), rownames(cache_meta)))

saveRDS(list(counts = counts_dgc, meta = cache_meta),
        file = file.path(CACHE_DIR,
                         paste0("malignant_subset_counts_meta", SUFFIX, ".rds")))
message("  Saved: ", CACHE_DIR, "/malignant_subset_counts_meta", SUFFIX, ".rds")

# Per-cell scores. Row identifier written as an explicit `cell_id` column so
# stage 2's header-sniffing loader resolves it the same way it does the pooled
# file.
cyto_out <- mal_meta %>%
  select(cell_id, CytoTRACE2_Score, CytoTRACE2_Potency, cytotrace_group,
         any_of(setdiff(colnames(scores_cohort),
                        c("CytoTRACE2_Score", "CytoTRACE2_Potency"))))
write_tab(cyto_out, "cytotrace2_scores", "tables")

# --- How much did re-fitting inside the cohort actually change the scores? ---
# One join and one correlation against the pooled run's per-cell scores. This
# is the single number that says whether "cohort-refit CytoTRACE2" is a
# materially different quantity from the pooled score, and it costs nothing.
# It is a DIAGNOSTIC, not a cross-cohort comparison.
POOLED_SCORES <- file.path(BASE_DIR, "outputs", "cytotrace_analysis", "tables",
                           "cytotrace2_scores_malignant.csv")
rho_refit_vs_pooled <- NA_real_
if (file.exists(POOLED_SCORES)) {
  .pool <- read.csv(POOLED_SCORES, stringsAsFactors = FALSE)
  .m <- match(mal_meta$cell_id, .pool$cell_id)
  if (sum(!is.na(.m)) >= 100) {
    rho_refit_vs_pooled <- suppressWarnings(stats::cor(
      mal_meta$CytoTRACE2_Score, .pool$CytoTRACE2_Score[.m],
      use = "complete.obs", method = "spearman"))
    message("  Spearman rho(cohort-refit score, pooled score) = ",
            round(rho_refit_vs_pooled, 4), " over ", sum(!is.na(.m)), " cells")
    if (!is.na(rho_refit_vs_pooled) && rho_refit_vs_pooled < 0.8)
      message("  NOTE: rho < 0.8 -- CytoTRACE2 is more cohort-dependent than ",
              "assumed. Every within-cohort potency statement downstream is ",
              "about THIS cohort's internal ranking only.")
  }
  rm(.pool); gc()
} else {
  message("  NOTE: pooled score file not found -- rho(refit, pooled) not computed.")
}

run_config <- data.frame(
  key = c("cohort", "outcome_label", "n_tumour_cells", "n_patients",
          "n_genes", "cytotrace2_ncores", "cytotrace2_seed",
          "rho_refit_vs_pooled_cytotrace", "source_cache", "run_timestamp",
          "annot_col", "n_annot5_classes_present", "annot5_class_counts",
          "myo_sparse_patients"),
  value = c(COHORT, COHORT_LABEL, as.character(n_cells), as.character(n_patients),
            as.character(nrow(counts_cohort)), as.character(N_CORES), "14",
            as.character(round(rho_refit_vs_pooled, 4)), POOLED_CACHE,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            ANNOT_COL, as.character(length(.a5_tab)),
            paste0(names(.a5_tab), "=", as.integer(.a5_tab), collapse = ";"),
            if (length(MYO_SPARSE_PATIENTS))
              paste0(MYO_SPARSE_PATIENTS, "(n=", .myo_n[MYO_SPARSE_PATIENTS], ")",
                     collapse = ";") else "none"),
  stringsAsFactors = FALSE
)
write_tab(run_config, "run_config_cytotrace2", "tables")

message("=== DONE: CytoTRACE2 complete for the ", COHORT_LABEL, " cohort ===")
message("  ", n_cells, " tumour cells, ", n_patients, " patients")
message("  Next: Rscript monocle3_cohort.R ", COHORT)
