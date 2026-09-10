#!/usr/bin/env Rscript
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

#
# monocle3_cohort.R
#
# Monocle3 trajectory / pseudotime analysis of the MALIGNANT (tumour)
# compartment of an adenoid cystic carcinoma scRNA-seq cohort, run SEPARATELY
# INSIDE ONE OUTCOME COHORT (Good or Poor) rather than across the pooled set.
#
#   Usage:  Rscript monocle3_cohort.R good
#           Rscript monocle3_cohort.R poor
#
# Stage 2 of the two-stage outcome-split pipeline:
#   cytotrace2_cohort.R <cohort>   ->   monocle3_cohort.R <cohort>
# It reads that cohort's cache and cohort-refit CytoTRACE2 scores, builds the
# cohort's OWN embedding and trajectory, roots it, and orders the cells.
#
# SCOPE: malignant cells only (annot_5 %in% MALIGNANT_TYPES; formerly annot_6
#   "ADC - Tumor","ADMEC - Tumor","Basal cells - Tumor")), after excluding
#   cells flagged marked_cells == "marked_removal" -- both filters are already
#   applied in the cache this reads, so they are not re-applied here.
#
# ROOTING: unchanged from the pooled script. The root is the principal-graph
#   node closest to the Monocle3 cluster with the HIGHEST mean CytoTRACE2
#   score, restricted to the DOMINANT principal-graph partition (order_cells()
#   only assigns finite pseudotime inside the root's own partition, so rooting
#   in a small residual partition would leave almost every cell undefined).
#   Highest CytoTRACE2 = most stem-like, so pseudotime increases toward more
#   differentiated states with no downstream reversal. The rule NEVER reads
#   outcome_group or batch. Inside a single cohort that outcome-blindness is
#   automatic -- outcome_group is constant -- but the rule is kept identical
#   to the pooled one so the two runs differ only in their input cell set.
#   Two alternative orientations (myoepithelial-program-based, and a
#   cluster-independent CytoTRACE2 top-5% centroid) are computed alongside as
#   a sensitivity check: order_cells() is cheap once learn_graph() has run, and
#   without them the root is chosen on CytoTRACE2 and validated by nothing.
#
# WHAT IS AND IS NOT COMPARABLE ACROSS THE TWO COHORT RUNS
#   Each cohort gets its OWN HVG panel, PCA, MNN alignment, UMAP, Leiden
#   clusters, principal graph and pseudotime. Everything in a re-derived
#   coordinate system is therefore MEANINGLESS across the two runs: UMAP axes,
#   cluster IDs, partition IDs, principal-graph node IDs, pseudotime VALUES,
#   Moran's I ranks and gene-module numbers. "Good cluster 4" and "Poor
#   cluster 4" are unrelated objects that happen to share an integer.
#   See OUTCOME_SPLIT_NOTES.md before comparing anything between the trees.
#
# DIFFERENCES FROM THE POOLED monocle3_malignant_trajectory.R (its line
# numbers): the h5ad counts fallback (493-518) is DELETED and replaced by a
# stop(); align_cds gains a retry ladder; the patient-level Wilcoxon/BH
# battery (1823-1919), focal-cluster selection (1921-1949), the per-patient
# outcome figures (1951-2004), the Poor-patient LOO (2005-2027) and the
# structural LOO (2073-2229) are all DROPPED -- every one of them is either
# degenerate on a single-level outcome factor or out of scope. Outcome-coloured
# panels are re-pointed at annot_6 or at patient.
#
# ---------------------------------------------------------------------------

suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
set.seed(42)

# The source h5ad's sparse layers (e.g. layers/counts) are written with the
# HDF5 "lzf" compression filter (a non-standard, dynamically-loaded filter
# used by h5py). The system HDF5 library used by R's rhdf5/HDF5Array does not
# ship this filter, so H5Dread() fails on any read of such a dataset unless
# HDF5's dynamic filter plugin search path is pointed at a filter plugin that
# provides it. The Bioconductor package `rhdf5filters` already bundles a
# compatible libH5Zlzf.so; pointing HDF5_PLUGIN_PATH at it (before any HDF5
# read happens) resolves this with no system/root changes required. This must
# be set before `rhdf5`/`HDF5Array`/`zellkonverter` perform any read.
if (requireNamespace("rhdf5filters", quietly = TRUE)) {
  Sys.setenv(HDF5_PLUGIN_PATH = rhdf5filters::hdf5_plugin_path())
}

# =============================================================================
# 0. Constants / toggles
# =============================================================================

BASE_DIR <- ACC_DATA_ROOT

## ---- Cohort argument --------------------------------------------------------
# The cohort MUST be given explicitly -- there is no default, so a mistyped or
# missing argument can never silently write one cohort's trajectory into the
# other's tree. It also makes the running cohort visible in `ps -eo args`,
# which is how the standing liveness check tells the two stages apart.
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) < 1 || !.args[1] %in% c("good", "poor")) {
  stop("Usage: Rscript monocle3_cohort.R <good|poor>")
}
COHORT       <- .args[1]
COHORT_LABEL <- c(good = "Good", poor = "Poor")[[COHORT]]
SUFFIX       <- paste0("_", COHORT)

# Cohort sizes verified against the pooled per-cell export.
EXPECTED_MIN_CELLS_COHORT <- c(good = 30000L, poor = 35000L)[[COHORT]]

message("=== Monocle3 trajectory, ", COHORT_LABEL,
        "-outcome cohort (tumour cells only) ===")

# SMOKE_TEST is retained only so the blocks copied verbatim from the pooled
# script still reference a defined symbol. This subproject has no smoke tree.
SMOKE_TEST <- FALSE
SMOKE_N    <- 5000

# Must match the OUTCOME_COL used by the pooled cytotrace2_malignant_analysis.R,
# which wrote the parent cache that cytotrace2_cohort.R subset.
OUTCOME_COL <- "clinical_outcome"

OUT_DIR <- file.path(BASE_DIR, "outputs", "outcome_split", COHORT,
                     "monocle3_analysis")

# Written by cytotrace2_cohort.R for THIS cohort -- not the shared pooled cache.
CACHE_PATH <- file.path(BASE_DIR, "outputs", "outcome_split", COHORT, "cache",
                        paste0("malignant_subset_counts_meta", SUFFIX, ".rds"))
H5AD_PATH  <- file.path(BASE_DIR, "outputs", "h5ad_files", "inprogress_3.h5ad")
# Cohort-refit scores, also from cytotrace2_cohort.R. The rooting below is
# therefore driven by this cohort's own CytoTRACE2 fit, not the pooled one.
CYTO_CSV_PATH <- file.path(BASE_DIR, "outputs", "outcome_split", COHORT,
                           "cytotrace_analysis", "tables",
                           paste0("cytotrace2_scores", SUFFIX, ".csv"))

# ---- Annotation level -------------------------------------------------------
# annot_5 is the annotation of record; it arrives via the cohort cache, already
# recoded upstream. annot_6 is carried as a passenger column only.
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
# cytotrace2_cohort.R and both pooled scripts -- edit all four.
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
# The myoepithelial-lineage classes, for the root-composition summary.
MYO_LINEAGE_TYPES <- c("Myoepithelial cells - Tumor",
                       "Actively Dividing Myoepithelial cells - Tumor",
                       "Basal cells - Tumor",
                       "Epithelial cells - Basal - Tumor")
# Minimum cells before a class gets its own program-score facet. Kept from the
# annot_6 version, but the decision is now RECORDED rather than silent -- see
# program_scores_facet_inclusion_<c>.csv.
MIN_CELLS_PER_ANNOT5_FACET <- 200L
.MALIGNANT_TYPES_LEGACY <- c("Epithelial cells - Tumor", "ADC - Tumor",
                      "ADMEC - Tumor", "Basal cells - Tumor")

# Matches the convention already established in cytotrace2_malignant_analysis.R
# (line 124), so outcome colors are consistent across the whole pipeline.
OUTCOME_COLORS <- c("Good" = "forestgreen", "Poor" = "firebrick3")
COHORT_COLOR   <- OUTCOME_COLORS[[COHORT_LABEL]]

RECOMPUTE_MALIGNANT_HVGS <- TRUE
N_HVGS <- 2000

# Highly variable genes are selected with the per-patient mean-variance trend
# fitted SEPARATELY per patient (scran::modelGeneVar(block = batch)) and the
# resulting per-block statistics combined, rather than from the pooled malignant
# cells. In tumor scRNA-seq, patient-specific programs (CNV-driven expression,
# clone-private modules) are frequently the single largest source of variance;
# selecting HVGs from the pooled data therefore risks building the trajectory on
# features that encode patient identity rather than shared biology. Blocking
# removes the between-patient component from the variance decomposition, so the
# selected genes are those variable *within* patients. Set to FALSE only to
# reproduce the pre-revision (pooled) behaviour.
BLOCK_HVGS_ON_PATIENT <- TRUE

# Prespecified marker programs. These are fixed in advance and used for (a) the
# independent myoepithelial-based root orientation in Step 5, (b) per-cell
# program scores exported with the metadata, and (c) the smoothed
# program-vs-pseudotime figure that carries the "continuum" claim. Defined once
# here and re-read by cluster_characterization.R so the two scripts can never
# drift apart.
#
# The gene space of this dataset is a 10,000-gene HVG subset chosen upstream in
# Scanpy, not the full transcriptome; every gene below was checked to be present
# in it. Scoring silently tolerates absent genes (they are dropped, and the
# count actually used is logged), so a future change of gene space degrades the
# scores rather than aborting the run.
PROGRAM_GENES <- list(
  myoepithelial = c("TP63", "ACTA2", "TAGLN", "CNN1", "MYH11"),
  ductal_notch  = c("KIT", "KRT7", "KRT19", "GABRP", "ELF5",
                    "NOTCH3", "HEY1", "HES1")
)

# Fraction of top-scoring cells used by the "cytotrace_top5pct_centroid" root
# rule (Step 5). Cluster-independent by construction.
ROOT_TOP_FRAC <- 0.05

# Extra per-cell columns pulled directly from the h5ad `obs` table (obs only --
# the counts matrix is never touched, so this costs seconds rather than the
# minutes a full h5ad read would). These are QC / doublet / malignancy-profile
# fields that the companion script's cache RDS does not carry; pulling them here
# avoids a 40-minute re-run of cytotrace2_malignant_analysis.R purely to widen
# the cache. Any column absent from the h5ad is filled with NA and logged.
EXTRA_OBS_COLS <- c("prect_2", "doublet_score", "dd_score", "scrublet_score",
                    "solo_score", "doublet_consensus_2of3_str",
                    "triple_agree_doublet", "site", "years_of_survival",
                    "n_genes_by_counts", "total_counts_mt")

# NOTE: SMOKE_TEST / SMOKE_N are declared at the top of this section, above the
# path block that branches on them.

for (d in c("trajectory", "trajectory_original_umap", "pseudotime", "clusters",
            "gene_modules", "patient_level", "tables", "sensitivity")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# 1. Package install (idempotent) + library load
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c("batchelor")
for (pkg in bioc_pkgs)
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)

cran_pkgs <- c("assertthat", "RhpcBLASctl", "rsample", "slam", "proxy",
               "sf", "spdep", "ggdist", "pbmcapply")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Hard check for the known libudunits2-dev risk (units -> sf -> spdep -> monocle3).
if (!requireNamespace("units", quietly = TRUE) ||
    !requireNamespace("sf", quietly = TRUE)) {
  stop(
    "FATAL: required package 'units' and/or 'sf' failed to install.\n",
    "This is the known-risk failure mode: 'units' needs the system library\n",
    "libudunits2-dev (udunits2.h / -ludunits2), which is NOT installed on\n",
    "this machine, and this session has no sudo access.\n",
    "ACTION REQUIRED: someone with root must run:\n",
    "  sudo apt-get install -y libudunits2-dev\n",
    "before monocle3's dependency chain (spdep -> sf -> units) can be ",
    "installed.\n",
    "Do NOT attempt to work around this with a stub/fake package."
  )
}

if (!requireNamespace("leidenbase", quietly = TRUE))
  remotes::install_github("cole-trapnell-lab/leidenbase")
if (!requireNamespace("speedglm", quietly = TRUE))
  remotes::install_github("cole-trapnell-lab/speedglm")
if (!requireNamespace("BPCells", quietly = TRUE))
  remotes::install_github("bnprks/BPCells", subdir = "r")
if (!requireNamespace("monocle3", quietly = TRUE))
  remotes::install_github("cole-trapnell-lab/monocle3")

# NOTE on a known transient failure mode: when this script performs a
# first-time install of the whole dependency chain in a single R session, a
# late GitHub install (e.g. monocle3 or one of its deps) can pull a newer
# `igraph` build than the version already loaded into memory earlier in the
# same session (e.g. via batchelor/scran). R then tries to unload/reload
# igraph to reconcile versions and fails with "Package 'igraph' version X
# cannot be unloaded" because it is already imported by other loaded
# namespaces. This is purely a same-session install-then-load artifact: once
# all packages are actually installed on disk, a fresh Rscript invocation of
# this file loads cleanly (the install lines above are then no-ops via the
# requireNamespace guards). If the block below fails with that specific
# error, simply re-running this script resolves it.
load_core_libraries <- function() {
  suppressPackageStartupMessages({
    library(monocle3)
    library(SingleCellExperiment)
    library(Matrix)
    library(scran)
    library(scuttle)
    library(zellkonverter)
    library(ggplot2)
    library(patchwork)
    library(cowplot)
    library(dplyr)
    library(tidyr)
    library(svglite)
    library(ragg)
    library(ggrastr)
    library(igraph)
  })
}

lib_load_result <- tryCatch({
  load_core_libraries()
  TRUE
}, error = function(e) {
  msg <- conditionMessage(e)
  if (grepl("cannot be unloaded", msg, fixed = TRUE) ||
      grepl("igraph", msg, fixed = TRUE)) {
    message("Caught a same-session package-version-reconciliation error ",
            "while loading libraries (see NOTE above): ", msg)
    message("Retrying library load once (this typically self-resolves once ",
            "installs from this same session have settled)...")
    Sys.sleep(1)
    tryCatch({ load_core_libraries(); TRUE }, error = function(e2) {
      stop("Library loading failed on retry as well. This usually means a ",
           "package version was changed mid-session by the install steps ",
           "above. Please re-run this script as a fresh `Rscript` process ",
           "(all dependencies are now installed on disk, so the install ",
           "steps will be skipped and library loading will succeed). ",
           "Original error: ", msg, " / Retry error: ", conditionMessage(e2))
    })
  } else {
    stop(e)
  }
})
set.seed(42)

message("=== All packages loaded OK ===")

# =============================================================================
# Helper functions
# =============================================================================

# Inject the cohort suffix into an output basename. Applied here, once, so the
# ~12 write.csv/saveRDS call sites copied verbatim from the pooled script still
# produce collision-free basenames -- directory separation alone does not
# satisfy the no-duplicate-basename rule (ANALYSIS_NOTES.md S7).
sfx <- function(fname) sub("\\.([A-Za-z0-9]+)$", paste0(SUFFIX, ".\\1"), fname)

save_plot <- function(p, name, subdir = "", width = 8, height = 6, dpi = 600) {
  name    <- paste0(name, SUFFIX)
  out_dir <- file.path(OUT_DIR, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), plot = p, width = width,
         height = height, dpi = dpi, units = "in", bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot = p, width = width,
         height = height, units = "in", bg = "white", device = svglite::svglite)
  message("  Saved: ", subdir, "/", name)
}
score_program <- function(logmat, genes, label = "") {
  # Mean of per-gene z-scores (z-scored ACROSS CELLS) over the genes of a
  # prespecified program. Z-scoring first prevents a single high-expression
  # gene (e.g. ACTA2) from dominating the average purely through its scale.
  # Genes absent from the gene space, or with zero variance across cells, are
  # dropped and reported rather than silently producing NaN.
  present <- intersect(genes, rownames(logmat))
  missing <- setdiff(genes, present)
  if (length(missing) > 0)
    message("    [", label, "] not in gene space, dropped: ",
            paste(missing, collapse = ", "))
  if (length(present) == 0) {
    message("    [", label, "] WARNING: no genes available -- score is all NA")
    return(rep(NA_real_, ncol(logmat)))
  }
  sub <- as.matrix(logmat[present, , drop = FALSE])
  mu  <- rowMeans(sub)
  sdv <- apply(sub, 1, stats::sd)
  keep <- is.finite(sdv) & sdv > 0
  if (!any(keep)) {
    message("    [", label, "] WARNING: all genes have zero variance -- score is all NA")
    return(rep(NA_real_, ncol(logmat)))
  }
  if (any(!keep))
    message("    [", label, "] zero-variance, dropped: ",
            paste(present[!keep], collapse = ", "))
  z <- (sub[keep, , drop = FALSE] - mu[keep]) / sdv[keep]
  message("    [", label, "] scored on ", sum(keep), "/", length(genes), " genes: ",
          paste(present[keep], collapse = ", "))
  colMeans(z)
}

read_h5ad_obs_subset <- function(h5ad_path, cell_ids, cols) {
  # Reads selected columns from an AnnData `obs` table WITHOUT touching X or
  # any layer, and returns them aligned to `cell_ids`. Used to recover QC /
  # doublet / malignancy-profile fields that the shared cache RDS does not
  # carry. Reading /obs only costs seconds against this 8 GB file.
  #
  # AnnData encodes a column either as a plain HDF5 dataset (numeric/boolean)
  # or, for categoricals, as a GROUP holding `categories` + 0-based `codes`
  # (with -1 meaning NA). Both forms are handled. Every failure mode --
  # missing file, missing column, length mismatch, unreadable index -- degrades
  # to an all-NA column with a logged warning rather than aborting a 40-minute
  # run over an optional annotation.
  out <- data.frame(row.names = cell_ids)
  fill_na <- function(o) { for (cl in cols) o[[cl]] <- NA; o }

  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    message("  WARNING: rhdf5 unavailable -- extra obs columns filled with NA")
    return(fill_na(out))
  }
  if (!file.exists(h5ad_path)) {
    message("  WARNING: h5ad not found (", h5ad_path,
            ") -- extra obs columns filled with NA")
    return(fill_na(out))
  }

  obs_index <- tryCatch(as.character(rhdf5::h5read(h5ad_path, "/obs/_index")),
                        error = function(e) NULL)
  if (is.null(obs_index)) {
    message("  WARNING: could not read /obs/_index -- extra obs columns filled with NA")
    return(fill_na(out))
  }
  m <- match(cell_ids, obs_index)
  n_unmatched <- sum(is.na(m))
  if (n_unmatched > 0)
    message("  WARNING: ", n_unmatched, " / ", length(cell_ids),
            " cell IDs not found in /obs/_index -- those entries will be NA")

  for (cl in cols) {
    v <- tryCatch({
      x <- rhdf5::h5read(h5ad_path, paste0("/obs/", cl))
      if (is.list(x) && all(c("categories", "codes") %in% names(x))) {
        cats  <- as.character(x$categories)
        codes <- as.integer(x$codes)
        ifelse(codes < 0L, NA_character_, cats[codes + 1L])
      } else {
        as.vector(x)
      }
    }, error = function(e) NULL)

    if (is.null(v)) {
      message("  WARNING: obs column '", cl, "' not readable -- filled with NA")
      out[[cl]] <- NA
    } else if (length(v) != length(obs_index)) {
      message("  WARNING: obs column '", cl, "' has length ", length(v),
              " but /obs/_index has ", length(obs_index), " -- filled with NA")
      out[[cl]] <- NA
    } else {
      out[[cl]] <- v[m]
    }
  }
  rhdf5::h5closeAll()
  out
}

# =============================================================================
# 2. Load the cohort's tumour cells (cohort cache only -- no h5ad fallback)
# =============================================================================

message("=== Step 2: Loading data (", COHORT_LABEL, " cohort) ===")

# THE h5ad FALLBACK PRESENT IN THE POOLED SCRIPT IS DELETED HERE, DELIBERATELY.
# In the pooled script, a missing/undersized cache falls back to reading the
# 8 GB h5ad and re-deriving the malignant subset. That branch is wrong twice
# over in this subproject: it would rebuild the POOLED compartment, not this
# cohort, and it costs hours while exiting 0. The only correct response to a
# missing cohort cache is to stop and re-run stage 1 for this cohort.
# H5AD_PATH is still used further down for the obs-only column read (seconds).
if (!file.exists(CACHE_PATH))
  stop("Cohort cache not found: ", CACHE_PATH,
       "\n  Run stage 1 first:  Rscript cytotrace2_cohort.R ", COHORT)

message("  Using cohort cache: ", CACHE_PATH)
cached <- readRDS(CACHE_PATH)
counts_mat <- cached$counts
meta_df    <- cached$meta
rm(cached); gc()

stopifnot(is(counts_mat, "dgCMatrix") || is(counts_mat, "Matrix"))
counts_mat <- as(counts_mat, "CsparseMatrix")
stopifnot(identical(colnames(counts_mat), rownames(meta_df)))

# Cohort-sized guard, replacing the pooled EXPECTED_MIN_MALIGNANT_CELLS_FULL.
# A cohort cache far below this means stage 1 wrote a partial or wrong subset.
if (ncol(counts_mat) < EXPECTED_MIN_CELLS_COHORT)
  stop("Cohort cache holds only ", ncol(counts_mat), " cells, below the ",
       "expected minimum of ", EXPECTED_MIN_CELLS_COHORT, " for the ",
       COHORT_LABEL, " cohort. Re-run stage 1 for this cohort rather than ",
       "proceeding on a partial cell set.")

message("  counts_mat dim: ", paste(dim(counts_mat), collapse = " x "),
        "  meta_df rows: ", nrow(meta_df))
stopifnot(identical(colnames(counts_mat), rownames(meta_df)))

required_meta_cols <- c("batch", OUTCOME_COL)
missing_cols <- setdiff(required_meta_cols, colnames(meta_df))
if (length(missing_cols) > 0)
  stop("meta_df is missing required column(s): ", paste(missing_cols, collapse = ", "),
       ". Outcome-like columns present: ",
       paste(grep("outcome", colnames(meta_df), value = TRUE, ignore.case = TRUE),
             collapse = ", "),
       ". If loading from the cache at ", CACHE_PATH,
       ", re-run cytotrace2_malignant_analysis.R with OUTCOME_COL = '", OUTCOME_COL,
       "' to regenerate it.")

# Resolve OUTCOME_COL -> internal `outcome_group` once, on whichever load path
# ran above, so every downstream summary/plot/test reads one canonical column.
meta_df$outcome_group <- as.character(meta_df[[OUTCOME_COL]])

.patient_outcome <- unique(meta_df[, c("batch", "outcome_group")])
if (nrow(.patient_outcome) != length(unique(meta_df$batch)))
  stop("At least one patient (batch) maps to more than one '", OUTCOME_COL,
       "' value; outcome must be constant within a patient.")
message("  Outcome grouping from '", OUTCOME_COL, "': ",
        paste(names(table(.patient_outcome$outcome_group)),
              table(.patient_outcome$outcome_group), sep = " = ", collapse = ", "),
        " (", nrow(.patient_outcome), " patients total)")

# Exactly one outcome level, by construction of the cohort cache. `outcome_group`
# is kept in colData all the same, so every export below keeps the pooled
# schema and the two runs' tables line up column-for-column.
if (length(unique(meta_df$outcome_group)) != 1)
  stop("The ", COHORT_LABEL, " cohort cache contains ",
       length(unique(meta_df$outcome_group)), " outcome levels; expected 1. ",
       "Stage 1 wrote the wrong cell set.")
if (unique(meta_df$outcome_group) != COHORT_LABEL)
  stop("Cohort cache holds outcome_group == '", unique(meta_df$outcome_group),
       "' but this run was invoked as '", COHORT, "'. Refusing to mislabel ",
       "an entire output tree.")
message("  Cohort confirmed: ", COHORT_LABEL, " -- ", ncol(counts_mat),
        " tumour cells, ", nrow(.patient_outcome), " patients")

# Factor columns arrive levelled on the FULL dataset. After subsetting to one
# cohort the unused levels linger and would print empty keys in every discrete
# legend (and empty facets), so drop them once, here.
meta_df <- droplevels(meta_df)

# =============================================================================
# 2b. Smoke-test subsample (stratified by batch)
# =============================================================================

if (SMOKE_TEST) {
  message("=== SMOKE_TEST = TRUE: subsampling to ~", SMOKE_N,
          " cells, stratified by batch ===")
  meta_df$.orig_idx <- seq_len(nrow(meta_df))
  frac <- min(1, SMOKE_N / nrow(meta_df))
  set.seed(42)
  sampled_idx <- unlist(lapply(split(meta_df$.orig_idx, meta_df$batch), function(idx) {
    n_take <- max(1, round(length(idx) * frac))
    if (n_take >= length(idx)) return(idx)
    sample(idx, n_take)
  }))
  sampled_idx <- sort(sampled_idx)
  meta_df$.orig_idx <- NULL
  meta_df <- meta_df[sampled_idx, , drop = FALSE]
  counts_mat <- counts_mat[, sampled_idx, drop = FALSE]
  message("  Post-subsample: ", ncol(counts_mat), " cells")
  stopifnot(identical(colnames(counts_mat), rownames(meta_df)))
}

# =============================================================================
# 3. Malignant-specific HVGs (patient-blocked) + prespecified program scores
# =============================================================================

message("=== Step 3: HVG selection (malignant-specific = ", RECOMPUTE_MALIGNANT_HVGS,
        ", patient-blocked = ", BLOCK_HVGS_ON_PATIENT, ") ===")

# The SCE is now built unconditionally: its logcounts are needed for the
# prespecified program scores below even when HVG recomputation is disabled.
sce_mal <- SingleCellExperiment(assays = list(counts = counts_mat),
                                 colData = meta_df)
sce_mal <- scuttle::logNormCounts(sce_mal)

if (RECOMPUTE_MALIGNANT_HVGS) {
  n_use <- min(N_HVGS, nrow(counts_mat))

  # Patient-blocked selection is the PRIMARY panel (see BLOCK_HVGS_ON_PATIENT).
  # modelGeneVar fits the mean-variance trend within each patient and combines
  # the per-block statistics, so between-patient variance no longer drives gene
  # ranking. Every patient in this cohort contributes >=500 malignant cells, so
  # per-block trend fitting is well conditioned.
  if (BLOCK_HVGS_ON_PATIENT) {
    .block <- colData(sce_mal)$batch
    message("  modelGeneVar(block = batch) over ", length(unique(.block)), " patients")
    dec <- scran::modelGeneVar(sce_mal, block = .block)
  } else {
    dec <- scran::modelGeneVar(sce_mal)
  }
  malignant_hvgs <- scran::getTopHVGs(dec, n = n_use)
  message("  Selected ", length(malignant_hvgs), " malignant-specific HVGs")

  # Companion unblocked panel, computed for COMPARISON ONLY -- it is never used
  # downstream. Its overlap with the blocked panel quantifies how much of the
  # original feature selection was driven by between-patient variance, which is
  # the question the revision asks, and does so without needing the previous
  # run's outputs.
  hvg_compare_path <- file.path(OUT_DIR, "tables", sfx("hvg_blocked_vs_unblocked.csv"))
  if (BLOCK_HVGS_ON_PATIENT) {
    dec_unblocked <- scran::modelGeneVar(sce_mal)
    hvgs_unblocked <- scran::getTopHVGs(dec_unblocked, n = n_use)
    n_overlap <- length(intersect(malignant_hvgs, hvgs_unblocked))
    message("  HVG panel comparison: ", n_overlap, " / ", n_use, " genes shared with ",
            "the unblocked panel (", round(100 * n_overlap / n_use, 1), "%); ",
            n_use - n_overlap, " genes are blocked-only")
    hvg_compare <- data.frame(
      gene = union(malignant_hvgs, hvgs_unblocked),
      stringsAsFactors = FALSE
    )
    hvg_compare$in_blocked_panel   <- hvg_compare$gene %in% malignant_hvgs
    hvg_compare$in_unblocked_panel <- hvg_compare$gene %in% hvgs_unblocked
    hvg_compare$rank_blocked   <- match(hvg_compare$gene, malignant_hvgs)
    hvg_compare$rank_unblocked <- match(hvg_compare$gene, hvgs_unblocked)
    hvg_compare <- hvg_compare[order(hvg_compare$rank_blocked,
                                     hvg_compare$rank_unblocked), ]
    write.csv(hvg_compare, hvg_compare_path, row.names = FALSE)
    message("  Saved: tables/hvg_blocked_vs_unblocked.csv (n_overlap = ", n_overlap,
            " / ", n_use, ")")
    rm(dec_unblocked, hvgs_unblocked, hvg_compare); gc()
  }
  rm(dec); gc()
} else {
  malignant_hvgs <- rownames(counts_mat)
  message("  Using all ", length(malignant_hvgs), " genes (HVG recompute disabled)")
}

# --- Prespecified program scores ---------------------------------------------
# Computed on log-normalised counts BEFORE any dimensionality reduction, so they
# are independent of the PCA/alignment/UMAP steps whose output they will later
# be used to orient and interpret. That independence is what makes the
# myoepithelial-based root in Step 5 a genuine sensitivity check on the
# CytoTRACE2-based root rather than a restatement of it.
message("  Computing prespecified program scores")
.logmat <- SummarizedExperiment::assay(sce_mal, "logcounts")
meta_df$score_myoepithelial <- score_program(.logmat, PROGRAM_GENES$myoepithelial,
                                              "myoepithelial")
meta_df$score_ductal_notch  <- score_program(.logmat, PROGRAM_GENES$ductal_notch,
                                              "ductal_notch")
message("  Program score correlation (myo vs ductal/NOTCH): ",
        round(stats::cor(meta_df$score_myoepithelial, meta_df$score_ductal_notch,
                         use = "complete.obs", method = "spearman"), 4))

rm(.logmat, sce_mal); gc()

# =============================================================================
# 4. Build cell_data_set and run standard Monocle3 pipeline
# =============================================================================

message("=== Step 4: Monocle3 pipeline (preprocess -> UMAP -> cluster -> graph) ===")

gene_metadata <- data.frame(gene_short_name = rownames(counts_mat),
                             row.names = rownames(counts_mat))
stopifnot(identical(colnames(counts_mat), rownames(meta_df)))

cds <- monocle3::new_cell_data_set(counts_mat, cell_metadata = meta_df,
                                    gene_metadata = gene_metadata)

n_dim <- min(50, ncol(cds) - 1, length(malignant_hvgs) - 1)
message("  preprocess_cds: num_dim = ", n_dim)
cds <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = n_dim,
                                 norm_method = "log", use_genes = malignant_hvgs)

# Batch/sample (patient) integration. Without this, malignant cells cluster
# strongly by patient (a well-known feature of tumor scRNA-seq, driven by
# patient-specific CNV/transcriptional programs), which fragments the UMAP
# graph into many disconnected per-patient partitions once
# learn_graph(use_partition = TRUE) is applied downstream. Since order_cells()
# only assigns finite pseudotime to cells reachable from the root within its
# own partition, unaligned data leaves most patients with zero or near-zero
# cells having any pseudotime at all -- observed directly in smoke testing
# (17 disconnected partitions, ~93/4989 cells with finite pseudotime,
# essentially confined to one patient). align_cds(alignment_group = "batch")
# uses batchelor's mutual-nearest-neighbors correction (batchelor is already
# an installed monocle3 dependency, so no new packages are required) to
# remove this per-patient artifact from the low-dimensional embedding before
# UMAP/clustering/graph learning, so the resulting trajectory reflects
# biological progression rather than patient identity.
# --- align_cds retry ladder (the one genuinely uncertain step) ---------------
# This cohort has 11-13 batches where the pooled run had 24. batchelor's
# reducedMNN has been observed to fail with a NaN batch magnitude at 23
# batches; that failure is not obviously batch-count-driven, so 11-13 is not
# provably safe either. Retry at progressively smaller alignment_k before
# giving up.
#
# The PCA fallback is here so the run finishes with an HONEST LABEL, not so the
# result gets used. Unaligned malignant cells cluster by patient and shatter
# learn_graph(use_partition = TRUE): smoke testing produced 17 partitions with
# ~93/4,989 cells having finite pseudotime. That output is useless AND exits 0,
# which is exactly why alignment_status is recorded and checked.
ALIGNMENT_K_LADDER <- c(20, 10, 5)
alignment_status <- NA_character_
for (.k in ALIGNMENT_K_LADDER) {
  message("  align_cds (batchelor MNN, alignment_k = ", .k, ") over ",
          length(unique(colData(cds)$batch)), " batches")
  .aligned <- tryCatch({
    monocle3::align_cds(cds, alignment_group = "batch", alignment_k = .k)
  }, error = function(e) {
    message("    FAILED at alignment_k = ", .k, ": ", conditionMessage(e))
    NULL
  })
  if (!is.null(.aligned)) {
    cds <- .aligned
    alignment_status <- paste0("OK_alignment_k_", .k)
    message("    align_cds succeeded at alignment_k = ", .k)
    break
  }
}
rm(.aligned)

if (is.na(alignment_status)) {
  alignment_status <- "FAILED_USED_PCA"
  message("  *** align_cds FAILED at every alignment_k in the ladder (",
          paste(ALIGNMENT_K_LADDER, collapse = ", "), ") ***")
  message("  *** Falling back to the raw PCA embedding. alignment_status = ",
          "FAILED_USED_PCA. DO NOT INTERPRET ANY TRAJECTORY RESULT FROM THIS ",
          "RUN -- see run_config", SUFFIX, ".csv. ***")
  message("  reduce_dimension (UMAP on PCA -- UNALIGNED)")
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP",
                                     preprocess_method = "PCA")
} else {
  message("  reduce_dimension (UMAP on the aligned embedding)")
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP",
                                     preprocess_method = "Aligned")
}

message("  cluster_cells")
cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")

message("  learn_graph")
cds <- monocle3::learn_graph(cds, use_partition = TRUE, close_loop = TRUE)

message("=== Step 4 COMPLETE: preprocess_cds / reduce_dimension / cluster_cells / ",
        "learn_graph all ran successfully on ", ncol(cds), " cells ===")

# Persist an intermediate checkpoint so a rerun (once the CytoTRACE2 file
# exists) doesn't have to redo Steps 1-4.
checkpoint_path <- file.path(OUT_DIR, "tables", sfx("cds_checkpoint_pre_rooting.rds"))
saveRDS(cds, checkpoint_path)
message("  Checkpoint saved: ", checkpoint_path)

# =============================================================================
# 5. Outcome-blind root selection (depends on companion CytoTRACE2 script)
# =============================================================================

message("=== Step 5: Outcome-blind trajectory rooting ===")

if (!file.exists(CYTO_CSV_PATH)) {
  stop("Required file not found: ", CYTO_CSV_PATH,
       ". Run cytotrace2_malignant_analysis.R first -- its output is required ",
       "here for outcome-blind trajectory rooting.")
}

cyto_raw <- read.csv(CYTO_CSV_PATH, nrows = 5)
message("  CytoTRACE2 CSV header: ", paste(colnames(cyto_raw), collapse = ", "))

# Determine whether the first column is a cell-id column (non-numeric) or
# already a proper CSV with rownames encoded as column 1.
first_col_is_id <- !is.numeric(cyto_raw[[1]])
if (first_col_is_id) {
  cyto_scores <- read.csv(CYTO_CSV_PATH, row.names = 1)
} else {
  cyto_scores <- read.csv(CYTO_CSV_PATH)
  id_col <- intersect(c("cell_id", "X", "barcode", "cell"), colnames(cyto_scores))
  if (length(id_col) == 0)
    stop("Could not identify a cell-id column in ", CYTO_CSV_PATH,
         ". Columns present: ", paste(colnames(cyto_scores), collapse = ", "))
  rownames(cyto_scores) <- cyto_scores[[id_col[1]]]
}

if (!"CytoTRACE2_Score" %in% colnames(cyto_scores))
  stop("Expected column 'CytoTRACE2_Score' not found in ", CYTO_CSV_PATH,
       ". Columns present: ", paste(colnames(cyto_scores), collapse = ", "))

stopifnot(all(colnames(cds) %in% rownames(cyto_scores)))
colData(cds)$CytoTRACE2_Score <- cyto_scores[colnames(cds), "CytoTRACE2_Score"]

get_root_principal_nodes <- function(cds, target_cluster) {
  cell_ids <- which(monocle3::clusters(cds, reduction_method = "UMAP") == target_cluster)
  closest_vertex <- as.matrix(
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  )[cell_ids, , drop = FALSE]
  igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex))))
  ]
}

.principal_node_positions <- function(cds) {
  # dp_mst is stored as dims x nodes; transpose to nodes x dims and guarantee
  # node names, which older monocle3 versions do not always set as colnames.
  np <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  if (is.null(rownames(np)))
    rownames(np) <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name
  np
}

graph_landmarks <- function(cds, root_nodes, n_endpoints = 3) {
  # The principal graph has hundreds of nodes, and monocle3's plot_cells()
  # defaults (label_roots / label_leaves / label_branch_points all TRUE) print a
  # number next to every one of them. On 72k cells that buries the biology under
  # numbering. This returns just the landmarks worth naming on a manuscript
  # figure: the root, and the terminal nodes that actually carry cells.
  #
  # "Relevant endpoint" = a degree-1 node (a true leaf of the principal tree)
  # that lies in the ROOT'S OWN connected component -- leaves in other
  # partitions serve cells with undefined pseudotime and would imply a
  # trajectory endpoint that no cell can reach from the root.
  #
  # Cell support is counted over each leaf's order-2 graph neighbourhood, not
  # over the node itself. Principal nodes are dense relative to the cell cloud,
  # so a genuine terminal state spreads its cells over the last few nodes of the
  # branch and the exact tip often owns almost none.
  #
  # SELECTION RULE, in two steps. First drop leaves with negligible cell support
  # -- a tip serving a handful of cells is graph noise, not a state. Then, among
  # the survivors, take the ones with the HIGHEST median pseudotime. Ranking by
  # cell count instead would label whichever lobes happen to be biggest,
  # including mid-trajectory ones, which is not what an oriented continuum's
  # "endpoints" means: the figure is claiming where the trajectory ENDS, so the
  # marks have to be its late terminal states. Every qualifying leaf is written
  # to the exported table with a `labelled` flag, so the ones left unmarked on
  # the figure are still inspectable.
  pg <- monocle3::principal_graph(cds)[["UMAP"]]
  np <- .principal_node_positions(cds)
  vnames <- igraph::V(pg)$name
  root_nodes <- intersect(root_nodes, vnames)
  if (length(root_nodes) == 0) return(NULL)

  comp <- igraph::components(pg)
  root_comp <- comp$membership[root_nodes[1]]
  in_root_comp <- names(comp$membership)[comp$membership == root_comp]

  cv <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)
  cell_node <- vnames[cv[, 1]]
  pt <- colData(cds)$pseudotime

  deg <- igraph::degree(pg)
  leaves <- setdiff(intersect(names(deg)[deg == 1], in_root_comp), root_nodes)

  summarise_node <- function(nd) {
    nbhd <- vnames[as.integer(igraph::ego(pg, order = 2, nodes = nd)[[1]])]
    sel <- cell_node %in% nbhd
    data.frame(node_id = nd,
               lm_x = unname(np[nd, 1]), lm_y = unname(np[nd, 2]),
               n_cells_at_node = sum(cell_node == nd),
               n_cells_within2 = sum(sel),
               median_pseudotime = if (any(sel & is.finite(pt)))
                 median(pt[sel & is.finite(pt)]) else NA_real_,
               stringsAsFactors = FALSE)
  }

  root_df <- do.call(rbind, lapply(root_nodes, summarise_node))
  root_df$role     <- "Root"
  root_df$label    <- if (nrow(root_df) == 1) "Root" else paste0("Root", seq_len(nrow(root_df)))
  root_df$labelled <- TRUE

  end_df <- NULL
  if (length(leaves) > 0) {
    n_connected <- sum(cell_node %in% in_root_comp)
    min_support <- max(50, 0.005 * n_connected)
    end_df <- do.call(rbind, lapply(leaves, summarise_node))
    end_df$role  <- "Endpoint"
    end_df <- end_df[end_df$n_cells_within2 >= min_support &
                       is.finite(end_df$median_pseudotime), , drop = FALSE]
    end_df <- end_df[order(-end_df$median_pseudotime), , drop = FALSE]
    message("  Landmarks: ", length(leaves), " leaves in the root component, ",
            nrow(end_df), " with >= ", round(min_support), " cells nearby, ",
            min(n_endpoints, nrow(end_df)), " labelled on the figures.")
    if (nrow(end_df) > 0) {
      end_df$labelled <- seq_len(nrow(end_df)) <= n_endpoints
      end_df$label <- NA_character_
      end_df$label[end_df$labelled] <- paste0("E", seq_len(sum(end_df$labelled)))
    } else {
      # Filtered down to nothing: drop it rather than rbind()ing a frame that is
      # missing the label/labelled columns root_df has.
      end_df <- NULL
    }
  }

  out <- rbind(root_df, end_df)
  rownames(out) <- NULL
  out[, c("node_id", "role", "label", "labelled", "lm_x", "lm_y",
          "n_cells_at_node", "n_cells_within2", "median_pseudotime")]
}

isTRUE_vec <- function(x) !is.na(x) & x

add_landmarks <- function(p, lm) {
  # No-op rather than error if the graph is degenerate -- a missing annotation is
  # recoverable, a failed figure job is not.
  if (is.null(lm) || nrow(lm) == 0) return(p)
  lm <- lm[isTRUE_vec(lm$labelled), , drop = FALSE]   # unlabelled leaves stay in the CSV only
  if (nrow(lm) == 0) return(p)
  root <- lm[lm$role == "Root", , drop = FALSE]
  ends <- lm[lm$role == "Endpoint", , drop = FALSE]
  # Constant aesthetics only (no aes() mapping for shape/fill): mapping them
  # would introduce a second discrete scale and collide with the scale
  # plot_cells() already set for color_cells_by.
  if (nrow(ends) > 0)
    p <- p + geom_point(data = ends, aes(x = lm_x, y = lm_y), inherit.aes = FALSE,
                        shape = 21, size = 2.6, stroke = 0.7,
                        fill = "black", colour = "white")
  p <- p + geom_point(data = root, aes(x = lm_x, y = lm_y), inherit.aes = FALSE,
                      shape = 24, size = 3.0, stroke = 0.9,
                      fill = "white", colour = "black")
  p + ggrepel::geom_text_repel(data = lm, aes(x = lm_x, y = lm_y, label = label),
                               inherit.aes = FALSE, size = 3.1, fontface = "bold",
                               colour = "black", bg.colour = "white", bg.r = 0.12,
                               min.segment.length = 0, segment.size = 0.3,
                               box.padding = 0.45, max.overlaps = Inf, seed = 42)
}

get_root_node_near_cells <- function(cds, cell_idx, eligible_idx) {
  # CLUSTER-INDEPENDENT root rule: take the UMAP centroid of a set of cells and
  # return the nearest principal-graph node. Candidate nodes are restricted to
  # those actually serving cells of the dominant partition, so the centroid of a
  # diffuse cell set cannot land on a node in a different partition (which would
  # leave almost every cell with undefined pseudotime).
  um <- SingleCellExperiment::reducedDims(cds)[["UMAP"]]
  centroid <- colMeans(um[cell_idx, , drop = FALSE])
  cv <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)
  vnames <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name
  allowed <- unique(vnames[cv[eligible_idx, 1]])
  np <- .principal_node_positions(cds)
  np <- np[rownames(np) %in% allowed, , drop = FALSE]
  if (nrow(np) == 0)
    stop("No principal-graph nodes found in the dominant partition -- cannot root.")
  d2 <- rowSums(sweep(np, 2, centroid, "-")^2)
  rownames(np)[which.min(d2)]
}

# Stable short slugs for the six malignant classes, used to build fixed column
# names. Duplicated verbatim in the pooled and cohort monocle3 scripts.
.annot5_slug <- function(x) {
  m <- c("Myoepithelial cells - Tumor"                   = "myo",
         "Actively Dividing Myoepithelial cells - Tumor" = "admec",
         "Basal cells - Tumor"                           = "basalT",
         "Epithelial cells - Basal - Tumor"              = "epibasal",
         "Actively Dividing cells - Tumor"               = "adc",
         "Epithelial cells - Tumor"                      = "epi")
  unname(m[as.character(x)])
}

pick_root <- function(cds, method, eligible_idx, cluster_assign) {
  # Returns a single-row description of the chosen root under one orientation
  # rule. Every rule is restricted to `eligible_idx` (the dominant partition),
  # and none of them reads outcome_group -- all three remain outcome-blind.
  cd <- colData(cds)
  if (method == "cytotrace_cluster_mean") {
    m <- tapply(cd$CytoTRACE2_Score[eligible_idx], cluster_assign[eligible_idx],
                mean, na.rm = TRUE)
    root_cluster <- names(which.max(m))
    nodes <- get_root_principal_nodes(cds, root_cluster)
  } else if (method == "myo_cluster_mean") {
    m <- tapply(cd$score_myoepithelial[eligible_idx], cluster_assign[eligible_idx],
                mean, na.rm = TRUE)
    root_cluster <- names(which.max(m))
    nodes <- get_root_principal_nodes(cds, root_cluster)
  } else if (method == "cytotrace_top5pct_centroid") {
    sc <- cd$CytoTRACE2_Score
    sc[-eligible_idx] <- NA_real_
    thr <- stats::quantile(sc, probs = 1 - ROOT_TOP_FRAC, na.rm = TRUE)
    top_idx <- which(!is.na(sc) & sc >= thr)
    root_cluster <- NA_character_
    nodes <- get_root_node_near_cells(cds, top_idx, eligible_idx)
  } else {
    stop("Unknown root method: ", method)
  }

  # The cells the rule actually used to place the root: the root cluster's cells
  # for the cluster-based rules, the top-scoring cells for the centroid rule.
  # Without this the centroid rule would report n = 0 and NA composition, which
  # would look like a failure rather than a rule that has no root cluster.
  root_idx <- if (is.na(root_cluster)) top_idx else which(cluster_assign == root_cluster)

  # --- annot_5 identity of the root ------------------------------------------
  # This is what the annot_6 -> annot_5 migration exists to produce. At annot_6
  # the root was ~90% "Epithelial cells - Tumor" in every cohort, which was
  # uninformative because that one label covered 80% of the compartment. At
  # annot_5 the myoepithelial and epithelial poles are distinct, so the root can
  # actually be named.
  #
  # The six fractions are emitted with FIXED names whether or not the class is
  # present, so the pooled and both cohort CSVs share one schema and stack.
  .a5   <- factor(as.character(cd[[ANNOT_COL]])[root_idx], levels = MALIGNANT_LEVELS)
  .tab  <- table(.a5)
  .frac <- if (length(root_idx)) as.numeric(.tab) / length(root_idx)
           else rep(NA_real_, length(MALIGNANT_LEVELS))
  names(.frac) <- paste0("root_frac_annot5_", .annot5_slug(MALIGNANT_LEVELS))

  # Patient concentration of the root. `Basal cells - Tumor` is 276 cells of
  # which 274 come from ONE patient, so a root that is 8% Basal-Tumor is a
  # statement about that patient, not about the Good cohort. These two columns
  # sit on the same row as the composition so the caveat cannot be read
  # separately from the claim.
  .pat <- table(as.character(cd$batch)[root_idx])

  data.frame(
    method = method,
    root_cluster = root_cluster,
    root_node = paste(nodes, collapse = ";"),
    root_n_cells = length(root_idx),
    root_mean_cytotrace = if (length(root_idx)) mean(cd$CytoTRACE2_Score[root_idx], na.rm = TRUE) else NA_real_,
    root_mean_myo_score = if (length(root_idx)) mean(cd$score_myoepithelial[root_idx], na.rm = TRUE) else NA_real_,
    root_frac_good = if (length(root_idx)) mean(cd$outcome_group[root_idx] == "Good") else NA_real_,
    root_frac_poor = if (length(root_idx)) mean(cd$outcome_group[root_idx] == "Poor") else NA_real_,
    root_top_annot5 = if (length(root_idx) && sum(.tab) > 0) names(which.max(.tab)) else NA_character_,
    root_top_annot5_frac = if (length(root_idx)) max(.frac, na.rm = TRUE) else NA_real_,
    root_frac_myo_lineage = if (length(root_idx))
      sum(.tab[intersect(names(.tab), MYO_LINEAGE_TYPES)]) / length(root_idx) else NA_real_,
    as.list(.frac),
    root_n_patients = length(.pat),
    root_top_patient = if (length(.pat)) names(which.max(.pat)) else NA_character_,
    root_top_patient_frac = if (length(.pat)) max(.pat) / length(root_idx) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# Restrict root-cluster candidates to the DOMINANT (largest) principal-graph
# partition. This uses only partition size -- a purely structural/topological
# property of the graph, never outcome_group or batch -- so it remains fully
# outcome-blind. It exists to avoid a degenerate rooting outcome: even after
# align_cds() integration collapses most per-patient fragmentation, a handful
# of small residual partitions can still exist (e.g. rare outlier cells), and
# the globally highest-mean-CytoTRACE2 cluster can by chance fall inside one
# of those tiny partitions rather than the dominant one. Since order_cells()
# only assigns finite pseudotime within the root's own partition, rooting in
# a tiny partition would make pseudotime undefined for almost all cells
# regardless of how well batch integration worked upstream. Restricting
# candidates to the dominant partition ensures the root -- and therefore
# finite pseudotime -- lands in the partition that actually contains the vast
# majority of cells (and, after alignment, spans essentially all patients).
partition_assign <- monocle3::partitions(cds, reduction_method = "UMAP")
dominant_partition <- names(sort(table(partition_assign), decreasing = TRUE))[1]
n_partitions <- length(unique(partition_assign))
message("  Principal-graph partitions: ", n_partitions,
        " total; dominant partition ", dominant_partition, " contains ",
        sum(partition_assign == dominant_partition), " / ", length(partition_assign),
        " cells (", round(100 * mean(partition_assign == dominant_partition), 1), "%)")
eligible_cells <- which(partition_assign == dominant_partition)

cluster_assign_all <- monocle3::clusters(cds, reduction_method = "UMAP")
cluster_mean_cyto <- tapply(colData(cds)$CytoTRACE2_Score[eligible_cells],
                             cluster_assign_all[eligible_cells],
                             mean, na.rm = TRUE)
root_cluster <- names(which.max(cluster_mean_cyto))
message("  Root cluster (highest mean CytoTRACE2 among clusters in the ",
        "dominant partition = most stem-like/potent): ",
        root_cluster, " (mean score = ", round(max(cluster_mean_cyto, na.rm = TRUE), 4), ")")

# --- Orientation is a CHOICE, not a proven direction --------------------------
# Rooting on CytoTRACE2 means CytoTRACE2 and pseudotime cannot disagree by
# construction, so their agreement is not independent evidence. Two additional
# orientations are therefore computed:
#
#   myo_cluster_mean            -- roots on the prespecified myoepithelial
#                                  program (Step 3), which is derived from
#                                  marker genes alone and never sees CytoTRACE2.
#                                  This is the INDEPENDENT check: if it recovers
#                                  the same broad myoepithelial -> ductal/NOTCH
#                                  ordering, the continuum is not an artefact of
#                                  how the root was chosen.
#   cytotrace_top5pct_centroid  -- same CytoTRACE2 signal but a
#                                  cluster-INDEPENDENT rule, addressing the
#                                  observed instability of the cluster-mean rule
#                                  (the root cluster changed in 5/5 leave-one-out
#                                  scenarios, because Leiden cluster boundaries
#                                  move between runs).
#
# The CytoTRACE2 cluster-mean rule remains PRIMARY; the other two contribute
# extra pseudotime columns and a sensitivity table, and change nothing else.
# What is compared is the broad ORDERING of the two programs along pseudotime --
# not numerical pseudotime values, which are not expected to match.
ROOT_METHODS <- c("cytotrace_cluster_mean", "myo_cluster_mean",
                  "cytotrace_top5pct_centroid")
PRIMARY_ROOT_METHOD <- "cytotrace_cluster_mean"
ROOT_PSEUDOTIME_COL <- c(cytotrace_cluster_mean     = "pseudotime",
                         myo_cluster_mean           = "pseudotime_root_myo",
                         cytotrace_top5pct_centroid = "pseudotime_root_top5")

# Diagnostic only (does not affect the analysis): make the root cluster's
# patient/outcome composition visible directly in the log, since the
# rooting procedure above is outcome-blind by construction but the resulting
# cluster can still happen to be concentrated in a handful of patients or
# skewed toward one outcome group -- worth surfacing for interpretation.
root_cluster_idx <- which(cluster_assign_all == root_cluster)
message("  Root cluster composition: ", length(root_cluster_idx), " cells (",
        round(100 * length(root_cluster_idx) / length(cluster_assign_all), 2),
        "% of malignant cells)")
root_cluster_outcome <- table(colData(cds)$outcome_group[root_cluster_idx])
message("  Root cluster outcome breakdown: ",
        paste(names(root_cluster_outcome), root_cluster_outcome, sep = "=", collapse = ", "))
root_cluster_patients <- sort(table(colData(cds)$batch[root_cluster_idx]), decreasing = TRUE)
top_n <- seq_len(min(5, length(root_cluster_patients)))
message("  Root cluster top patients: ",
        paste(names(root_cluster_patients)[top_n], root_cluster_patients[top_n],
              sep = "=", collapse = ", "))

root_nodes <- get_root_principal_nodes(cds, root_cluster)
message("  Root principal node(s): ", paste(root_nodes, collapse = ", "))

cds <- monocle3::order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_nodes)
colData(cds)$pseudotime <- monocle3::pseudotime(cds)

# --- Sensitivity orientations -------------------------------------------------
# order_cells() is cheap once learn_graph() has run, so each alternative
# orientation costs seconds. Each runs on a COPY; only the primary rooting above
# is retained on `cds` itself. A failure in any alternative is recorded and does
# not abort the run -- these are sensitivity checks, not the main result.
message("  --- Alternative root orientations (sensitivity) ---")
root_rows <- list()
root_rows[[PRIMARY_ROOT_METHOD]] <- cbind(
  pick_root(cds, PRIMARY_ROOT_METHOD, eligible_cells, cluster_assign_all),
  data.frame(is_primary = TRUE, error_message = NA_character_)
)

for (mth in setdiff(ROOT_METHODS, PRIMARY_ROOT_METHOD)) {
  pt_col <- ROOT_PSEUDOTIME_COL[[mth]]
  # The handler returns a value rather than assigning into `cds` from inside
  # the error function -- superassignment across a tryCatch handler is fragile
  # for S4 replacement functions like colData<-.
  alt <- tryCatch({
    spec <- pick_root(cds, mth, eligible_cells, cluster_assign_all)
    alt_nodes <- strsplit(spec$root_node, ";", fixed = TRUE)[[1]]
    message("    ", mth, ": root cluster = ",
            ifelse(is.na(spec$root_cluster), "<cluster-independent>", spec$root_cluster),
            ", node(s) = ", spec$root_node)
    cds_alt <- monocle3::order_cells(cds, reduction_method = "UMAP",
                                      root_pr_nodes = alt_nodes)
    pt <- monocle3::pseudotime(cds_alt)
    rm(cds_alt); gc()
    list(spec = cbind(spec, data.frame(is_primary = FALSE,
                                        error_message = NA_character_)),
         pt = pt)
  }, error = function(e) {
    message("    ", mth, ": FAILED -- ", conditionMessage(e))
    # Build the failure row from the PRIMARY row's schema rather than by hand.
    # pick_root() now emits six fixed root_frac_annot5_* columns plus the
    # patient-concentration columns; a hand-maintained NA row silently drifts
    # out of sync and then do.call(rbind, ...) fails at the very end of a
    # 25-minute run. Copy the schema, blank every value, keep method/error.
    .blank <- root_rows[[PRIMARY_ROOT_METHOD]][1, , drop = FALSE]
    .blank[] <- lapply(.blank, function(col) col[NA_integer_])
    .blank$method        <- mth
    .blank$is_primary    <- FALSE
    .blank$error_message <- conditionMessage(e)
    list(spec = .blank, pt = NULL)
  })
  colData(cds)[[pt_col]] <- if (is.null(alt$pt)) NA_real_ else alt$pt
  root_rows[[mth]] <- alt$spec
}

# The comparison that matters: does each orientation order the two prespecified
# programs the same way? A myoepithelial-like -> ductal/NOTCH continuum implies
# rho(pseudotime, myo) < 0 and rho(pseudotime, ductal) > 0. Agreement in SIGN
# across orientations is the claim; agreement in magnitude is not expected.
root_sens <- do.call(rbind, root_rows)
.cdf <- as.data.frame(colData(cds))
primary_pt <- .cdf[[ROOT_PSEUDOTIME_COL[[PRIMARY_ROOT_METHOD]]]]
.sp <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 10) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
}
root_sens$n_finite_pseudotime <- NA_integer_
root_sens$rho_vs_primary_pseudotime <- NA_real_
root_sens$rho_pseudotime_vs_myo <- NA_real_
root_sens$rho_pseudotime_vs_ductal_notch <- NA_real_
root_sens$rho_pseudotime_vs_cytotrace <- NA_real_
for (i in seq_len(nrow(root_sens))) {
  pt <- .cdf[[ROOT_PSEUDOTIME_COL[[root_sens$method[i]]]]]
  if (is.null(pt)) next
  root_sens$n_finite_pseudotime[i]            <- sum(is.finite(pt))
  root_sens$rho_vs_primary_pseudotime[i]      <- .sp(pt, primary_pt)
  root_sens$rho_pseudotime_vs_myo[i]          <- .sp(pt, .cdf$score_myoepithelial)
  root_sens$rho_pseudotime_vs_ductal_notch[i] <- .sp(pt, .cdf$score_ductal_notch)
  root_sens$rho_pseudotime_vs_cytotrace[i]    <- .sp(pt, .cdf$CytoTRACE2_Score)
}

.sign_ok <- function(v) {
  v <- v[is.finite(v)]
  length(v) > 1 && (all(v > 0) || all(v < 0))
}
root_sens$orientation_concordant <- .sign_ok(root_sens$rho_pseudotime_vs_myo) &&
  .sign_ok(root_sens$rho_pseudotime_vs_ductal_notch)

.rs_primary <- root_sens[which(root_sens$is_primary)[1], ]
SCOPE_LABEL <- COHORT
ROOT_COMP_PATH <- file.path(OUT_DIR, "tables", sfx("root_annot5_composition.csv"))

# =============================================================================
# ROOT COMPOSITION BY annot_5 -- the output this migration exists to produce
# =============================================================================
# Under annot_6 the root was ~90% "Epithelial cells - Tumor" in every run, which
# said nothing: that single label covered 80% of the malignant compartment. At
# annot_5 the myoepithelial and epithelial poles are separate classes, so "what
# IS the root population" becomes an answerable question -- and it is answered
# here for EVERY root rule, not just the primary one, because a result that
# holds only under the CytoTRACE2-derived root would be circular.
#
# Two denominators, deliberately:
#   frac_of_root           -- what the root is made of.
#   frac_of_class_in_scope -- what fraction of that whole class sits in the
#                             root. This is what separates "the root is
#                             myoepithelial" from "myoepithelial cells are
#                             simply everywhere in this cohort".
.scope_a5 <- factor(as.character(colData(cds)[[ANNOT_COL]]), levels = MALIGNANT_LEVELS)
.scope_n  <- table(.scope_a5)
.slugs    <- .annot5_slug(MALIGNANT_LEVELS)

root_comp <- do.call(rbind, lapply(seq_len(nrow(root_sens)), function(i) {
  r <- root_sens[i, ]
  fr <- vapply(.slugs, function(sl) {
    v <- r[[paste0("root_frac_annot5_", sl)]]
    if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  n_cells <- if (is.na(r$root_n_cells)) rep(NA_integer_, length(fr))
             else as.integer(round(fr * r$root_n_cells))
  data.frame(
    scope                  = SCOPE_LABEL,
    method                 = r$method,
    is_primary             = r$is_primary,
    root_cluster           = r$root_cluster,
    root_node              = r$root_node,
    root_n_cells           = r$root_n_cells,
    annot_5                = MALIGNANT_LEVELS,
    n_cells                = n_cells,
    frac_of_root           = as.numeric(fr),
    n_in_scope             = as.integer(.scope_n[MALIGNANT_LEVELS]),
    frac_of_class_in_scope = ifelse(as.integer(.scope_n[MALIGNANT_LEVELS]) > 0,
                                    n_cells / as.integer(.scope_n[MALIGNANT_LEVELS]),
                                    NA_real_),
    root_frac_myo_lineage  = r$root_frac_myo_lineage,
    root_n_patients        = r$root_n_patients,
    root_top_patient       = r$root_top_patient,
    root_top_patient_frac  = r$root_top_patient_frac,
    error_message          = r$error_message,
    stringsAsFactors = FALSE, row.names = NULL)
}))
write.csv(root_comp, ROOT_COMP_PATH, row.names = FALSE)
message("  Saved: tables/", basename(ROOT_COMP_PATH), " (",
        nrow(root_comp), " rows, ", nrow(root_sens), " root rules)")

# The single figure that answers the question by eye.
.rc_plot <- root_comp[!is.na(root_comp$frac_of_root) & root_comp$frac_of_root > 0, ]
if (nrow(.rc_plot) > 0) {
  # Drop levels with no cells in this scope. Keeping them (drop = FALSE) renders
  # a legend key with nothing filled in -- in the Poor cohort `Basal cells -
  # Tumor` has zero cells, and a blank swatch reads as a rendering bug rather
  # than as "absent here". The palette is keyed by NAME, so a class still gets
  # the same colour in both cohorts even when the legends differ in length.
  .rc_plot$annot_5 <- factor(.rc_plot$annot_5,
                             levels = intersect(MALIGNANT_LEVELS,
                                                unique(.rc_plot$annot_5)))
  .rc_plot$method_lab <- paste0(.rc_plot$method,
                                ifelse(.rc_plot$is_primary, "\n(primary)", ""))
  p_rc <- ggplot(.rc_plot, aes(x = method_lab, y = frac_of_root, fill = annot_5)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = MALIGNANT_PALETTE, drop = FALSE, name = ANNOT_COL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = paste0("Root population composition by ", ANNOT_COL,
                        " - ", SCOPE_LABEL),
         subtitle = paste0("One bar per root rule. The myoepithelial-score rule\n",
                           "never sees CytoTRACE2, so agreement across bars\n",
                           "is the non-circular result."),
         x = NULL, y = "Fraction of root cells") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 8))
  save_plot(p_rc, "root_annot5_composition", "clusters", width = 11, height = 5.5)
}

write.csv(root_sens, file.path(OUT_DIR, "tables", sfx("root_sensitivity_summary.csv")),
          row.names = FALSE)
message("  Saved: tables/root_sensitivity_summary.csv")
message("  Orientation concordance across ", nrow(root_sens), " root rules: ",
        ifelse(isTRUE(root_sens$orientation_concordant[1]),
               "CONCORDANT (all rules agree on the sign of both program correlations)",
               "DISCORDANT -- inspect root_sensitivity_summary.csv before interpreting"))
rm(.cdf); gc()

message("=== Step 5 COMPLETE: trajectory rooted (outcome-blind), ",
        length(ROOT_METHODS), " orientations evaluated ===")

# =============================================================================
# 6. Downstream analyses: pseudotime + graph_test + gene modules
# =============================================================================

message("=== Step 6: pseudotime + graph_test + gene modules ===")

pseudotime_vec <- monocle3::pseudotime(cds)
colData(cds)$pseudotime <- pseudotime_vec

n_cores_graph_test <- if (SMOKE_TEST) min(8, parallel::detectCores()) else 20
message("  graph_test (cores = ", n_cores_graph_test, ")")
gene_test_res <- monocle3::graph_test(cds, neighbor_graph = "principal_graph",
                                       cores = n_cores_graph_test)
gene_test_res$q_value <- p.adjust(gene_test_res$p_value, method = "BH")
gene_test_res$gene <- rownames(gene_test_res)

sig_genes <- rownames(subset(gene_test_res, q_value < 0.05))
message("  Significant genes (q < 0.05): ", length(sig_genes))

gene_module_df <- NULL
if (length(sig_genes) >= 10) {
  gene_module_df <- tryCatch(
    monocle3::find_gene_modules(cds[sig_genes, ], resolution = 1e-3),
    error = function(e) {
      message("  WARNING: find_gene_modules failed: ", conditionMessage(e))
      NULL
    }
  )
} else {
  message("  WARNING: fewer than 10 significant genes (n = ", length(sig_genes),
          ") -- skipping find_gene_modules (expected on a small smoke subsample).")
}

# =============================================================================
# 7. Figures
# =============================================================================

message("=== Step 7: Figures ===")

n_cells_plot <- ncol(cds)
rasterize_pt_layer <- function(p, dpi = 600) {
  # plot_cells() returns a ggplot; try to rasterize its GeomPoint layer(s)
  # post-hoc via ggrastr. Fall back to returning p unchanged if that fails.
  tryCatch({
    is_point_layer <- vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))
    if (any(is_point_layer)) {
      for (i in which(is_point_layer)) {
        p$layers[[i]] <- ggrastr::rasterise(p$layers[[i]], dpi = dpi)
      }
    }
    p
  }, error = function(e) p)
}

# Landmarks for the trajectory panels: root + the terminal nodes that carry
# cells. Computed once, exported, and reused by every panel so the marks are
# identical across figures and traceable to node IDs rather than to a legend.
graph_lm <- tryCatch(graph_landmarks(cds, root_nodes, n_endpoints = 3),
                     error = function(e) {
                       message("  WARNING: graph_landmarks failed (",
                               conditionMessage(e),
                               ") -- trajectory panels drawn without landmarks.")
                       NULL
                     })
if (!is.null(graph_lm) && nrow(graph_lm) > 0) {
  write.csv(graph_lm, file.path(OUT_DIR, "tables", sfx("principal_graph_landmarks.csv")),
            row.names = FALSE)
  .lab <- graph_lm[isTRUE_vec(graph_lm$labelled), ]
  message("  Saved: tables/principal_graph_landmarks.csv (", nrow(graph_lm),
          " nodes; marked on figures: ",
          paste(paste0(.lab$label, "=", .lab$node_id), collapse = ", "), ")")
} else {
  message("  NOTE: no principal-graph landmarks resolved -- panels drawn unlabelled.")
}

traj_plot <- function(color_by, title) {
  # Every label_* argument is pinned FALSE. monocle3's defaults number every
  # root, leaf and branch point of the principal graph, which on this graph
  # prints hundreds of numbers over the cells. The root and the relevant
  # endpoints are re-added explicitly by add_landmarks() below.
  p <- monocle3::plot_cells(cds, color_cells_by = color_by,
                             label_cell_groups = FALSE,
                             label_groups_by_cluster = FALSE,
                             label_principal_points = FALSE,
                             label_branch_points = FALSE,
                             label_leaves = FALSE,
                             label_roots = FALSE,
                             graph_label_size = 3,
                             cell_size = 0.4) +
    ggtitle(title) +
    theme(legend.position = "bottom")
  if (color_by == "outcome_group") {
    p <- p + scale_color_manual(values = OUTCOME_COLORS)
  }
  # A fixed, named palette in myoepithelial -> ductal order. Named (not
  # positional) so a class keeps its colour even though Poor has one fewer
  # class than Good -- that is what makes the two cohorts' panels comparable.
  if (color_by == ANNOT_COL) {
    p <- p + scale_color_manual(values = MALIGNANT_PALETTE, drop = FALSE)
  }
  # Rasterize the cells first, then add landmarks, so the marks and their text
  # stay vector in the SVG.
  add_landmarks(rasterize_pt_layer(p), graph_lm)
}

# `outcome_group` is constant inside a cohort, so the pooled script's
# trajectory_by_outcome panel would be a single flat colour. annot_6 (the
# tumour-subtype annotation) is the categorical split that still carries
# information here.
if (ANNOT_COL %in% colnames(colData(cds))) {
  colData(cds)[[ANNOT_COL]] <- factor(as.character(colData(cds)[[ANNOT_COL]]),
                                      levels = MALIGNANT_LEVELS)
  save_plot(traj_plot(ANNOT_COL, paste0("Trajectory: ", ANNOT_COL,
                                        " tumour subtype (", COHORT_LABEL,
                                        " cohort)")),
            "trajectory_by_annot5", "trajectory", width = 8.5, height = 6.5)
}
save_plot(traj_plot("batch", "Trajectory: sample (batch)"),
          "trajectory_by_sample", "trajectory", width = 8, height = 6.5)
save_plot(traj_plot("CytoTRACE2_Score", "Trajectory: CytoTRACE2 score"),
          "trajectory_by_cytotrace", "trajectory", width = 7, height = 6.5)
if ("phase" %in% colnames(colData(cds))) {
  save_plot(traj_plot("phase", "Trajectory: cell-cycle phase"),
            "trajectory_by_phase", "trajectory", width = 7, height = 6.5)
}
save_plot(traj_plot("pseudotime", "Trajectory: pseudotime"),
          "trajectory_by_pseudotime", "trajectory", width = 7, height = 6.5)
save_plot(traj_plot("cluster", "Trajectory: Monocle cluster"),
          "trajectory_by_cluster", "trajectory", width = 7, height = 6.5)
save_plot(traj_plot("partition", "Trajectory: principal-graph partition"),
          "trajectory_by_partition", "trajectory", width = 7, height = 6.5)
save_plot(traj_plot("score_myoepithelial", "Trajectory: myoepithelial program score"),
          "trajectory_by_myo_score", "trajectory", width = 7, height = 6.5)
save_plot(traj_plot("score_ductal_notch", "Trajectory: ductal/NOTCH program score"),
          "trajectory_by_ductal_score", "trajectory", width = 7, height = 6.5)

# =============================================================================
# 7b. Companion plots on ORIGINAL (scVI) UMAP coordinates
# =============================================================================
# Monocle3's own UMAP above is recomputed from malignant-specific HVGs/PCA/
# batch-aligned embedding and looks nothing like the original scVI-derived
# UMAP computed once on the full dataset upstream of this script. That
# original layout is preserved per-cell as X_umap_1/X_umap_2 (already present
# in colData(cds), passed through from the cache at cds construction). These
# companion plots use the SAME per-cell color values as the plots above, just
# at the original (x, y) positions, for direct visual comparison. The
# principal graph is intentionally NOT drawn here -- it was learned in
# Monocle3's own embedding and is not meaningful in this coordinate space.
message("=== Step 7b: Original scVI UMAP companion plots ===")

# ---- Shared scVI-UMAP panel style -------------------------------------------
# These panels sit side by side in the same figure:
#   trajectory_original_umap/trajectory_by_*  (this script)
#   umap/umap_cytotrace_score_malignant       (cytotrace2_malignant_analysis.R)
#   umap/umap_cytotrace_potency_malignant     (cytotrace2_malignant_analysis.R)
# They draw the SAME cells at the SAME scVI UMAP coordinates, so any rendering
# difference between them is a visual artefact that makes the comparison
# equivocal. Everything controlling how a cell is drawn is therefore pinned here
# rather than left to a ggplot default.
#
# IMPORTANT: an identical copy of this block lives in
# cytotrace2_malignant_analysis.R (section 3b). Each script is deliberately
# standalone -- the pipeline runs them as separate Rscript invocations -- so
# edit BOTH or the panels silently drift apart again.
#
# The `pt_stroke` pin is the one that actually bit us: geom_point()'s default
# stroke = 0.5 adds ~0.95 mm to the drawn diameter, so `size = 0.4, stroke = 0`
# here and `size = 0.4` (default stroke) in the CytoTRACE2 script produced dots
# of 1.14 mm vs 2.08 mm -- a ~3.3x difference in area from two call sites that
# both said 0.4.
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
# potency, OUTCOME_COLORS for outcome); only the geometry and the frame are
# shared.
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
  # widened to match would push the per-sample batch legend past the canvas.
  if (is.numeric(df[[color_col]])) {
    p <- p + guides(color = guide_colourbar(theme = theme(
      legend.key.width  = grid::unit(2.2, "in"),
      legend.key.height = grid::unit(0.16, "in"))))
  }
  ggrastr::rasterise(p, dpi = UMAP_PANEL$raster_dpi)
}

# Pin the panel (the data area, excluding axes/title/legend) to an exact size.
# Without this, a continuous colourbar and a multi-level discrete legend consume
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
# rather than fixed: trajectory_by_sample_original_umap has one legend key per
# sample and needs ~1.5 in of legend, which a fixed canvas would clip. The DATA
# AREA is identical across every panel either way, which is the property the
# figure comparison rests on.
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

old_coord_df <- as.data.frame(colData(cds)) %>%
  mutate(cell_id = colnames(cds),
         monocle_cluster = as.character(monocle3::clusters(cds, reduction_method = "UMAP")))

if (all(c("X_umap_1", "X_umap_2") %in% colnames(old_coord_df))) {
  old_traj_plot <- function(color_by, title, discrete = TRUE) {
    p <- umap_panel(old_coord_df, color_by, title)
    if (color_by == "outcome_group") p <- p + scale_color_manual(values = OUTCOME_COLORS)
    if (!discrete) p <- p + scale_color_viridis_c()
    p
  }

  old_umap_specs <- list(
    list(ANNOT_COL,          paste0("Original scVI UMAP: ", ANNOT_COL, " (",
                                    COHORT_LABEL, " cohort)"), TRUE,
         "trajectory_by_annot5_original_umap"),
    list("batch",            "Original scVI UMAP: sample (batch)",   TRUE,
         "trajectory_by_sample_original_umap"),
    list("CytoTRACE2_Score", "Original scVI UMAP: CytoTRACE2 score", FALSE,
         "trajectory_by_cytotrace_original_umap"),
    list("pseudotime",       "Original scVI UMAP: pseudotime",       FALSE,
         "trajectory_by_pseudotime_original_umap"),
    list("monocle_cluster",  "Original scVI UMAP: Monocle3 cluster", TRUE,
         "trajectory_by_cluster_original_umap")
  )
  if ("phase" %in% colnames(old_coord_df)) {
    old_umap_specs <- append(old_umap_specs, list(
      list("phase", "Original scVI UMAP: cell-cycle phase", TRUE,
           "trajectory_by_phase_original_umap")
    ))
  }

  for (s in old_umap_specs) {
    save_umap_panel(old_traj_plot(s[[1]], s[[2]], s[[3]]), s[[4]],
                    "trajectory_original_umap")
  }
  message("  Saved ", length(old_umap_specs), " companion plots to trajectory_original_umap/")
} else {
  message("  WARNING: X_umap_1/X_umap_2 missing from colData(cds) -- skipping ",
          "original-UMAP companion plots.")
}

# --- pseudotime distribution by outcome (cell-level visual) ---
pt_df <- as.data.frame(colData(cds)) %>%
  mutate(cell_id = colnames(cds)) %>%
  select(cell_id, batch, outcome_group, pseudotime) %>%
  filter(is.finite(pseudotime))

# `outcome_group` is constant here, so the pooled script's by-outcome density
# would be a single curve. The within-cohort question is how far patients
# spread along their own cohort's trajectory, so split by patient instead.
p_pt_dist <- ggplot(pt_df, aes(x = pseudotime)) +
  geom_density(alpha = 0.5, fill = COHORT_COLOR, color = COHORT_COLOR,
               na.rm = TRUE) +
  labs(title = paste0("Per-cell pseudotime distribution - ", COHORT_LABEL,
                      " cohort"),
       subtitle = "Cell-level visual only; the patient-level summary is patient_level_summary_table.csv",
       x = "Pseudotime", y = "Density") +
  theme_bw(base_size = 12)
save_plot(p_pt_dist, "pseudotime_distribution", "pseudotime",
          width = 7, height = 5.5)

.pt_order <- pt_df %>% group_by(batch) %>%
  summarise(m = median(pseudotime), .groups = "drop") %>% arrange(m) %>% pull(batch)
p_pt_by_patient <- ggplot(pt_df %>% mutate(batch = factor(batch, levels = .pt_order)),
                          aes(x = batch, y = pseudotime)) +
  geom_violin(fill = COHORT_COLOR, alpha = 0.45, color = NA, scale = "width",
              na.rm = TRUE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.35,
               na.rm = TRUE) +
  coord_flip() +
  labs(title = paste0("Pseudotime by patient - ", COHORT_LABEL, " cohort"),
       subtitle = "Cells with finite pseudotime only; patients ordered by median",
       x = "Patient (batch)", y = "Pseudotime") +
  theme_bw(base_size = 12)
save_plot(p_pt_by_patient, "pseudotime_distribution_by_patient", "pseudotime",
          width = 7.5, height = max(4.5, 0.4 * length(.pt_order)))

# --- program scores along continuous pseudotime -------------------------------
# This is the figure that carries the "continuum" claim, so it deliberately
# shows CONTINUOUS pseudotime rather than the operational early/mid/late bins.
# Cells with non-finite pseudotime (outside the root's partition) are absent by
# construction; the count is stated in the subtitle rather than left implicit.
prog_df <- as.data.frame(colData(cds)) %>%
  mutate(cell_id = colnames(cds)) %>%
  select(cell_id, batch, outcome_group,
         any_of(c(ANNOT_COL, ANNOT_COL_LEGACY)), pseudotime,
         score_myoepithelial, score_ductal_notch) %>%
  filter(is.finite(pseudotime))

prog_long <- prog_df %>%
  tidyr::pivot_longer(cols = c(score_myoepithelial, score_ductal_notch),
                      names_to = "program", values_to = "score") %>%
  mutate(program = recode(program,
                          score_myoepithelial = "Myoepithelial-like",
                          score_ductal_notch  = "Ductal / NOTCH-associated"))

PROGRAM_COLORS <- c("Myoepithelial-like" = "#1b6ca8",
                    "Ductal / NOTCH-associated" = "#d1495b")

# Curves are BINNED MEANS, not a LOESS/GAM fit. Two reasons, both practical:
#
#   1. stats::loess does not scale to this cohort -- predLoess() errors outright
#      at ~70,000 points, silently dropping the layer and yielding a blank
#      panel. mgcv (method = "gam") is not installed on this machine, so that
#      fallback is unavailable too.
#   2. A binned mean is exactly what the reader thinks a trend line is: the
#      average program score of the cells at that stage of the trajectory, with
#      a 95% CI on that average. Nothing is smoothed away and nothing depends on
#      a span/knot choice.
#
# Bins hold equal NUMBERS OF CELLS (quantile breaks), not equal pseudotime
# widths, so every point on the curve is estimated from the same amount of data
# and the CI ribbon widens only where cells are genuinely sparse.
bin_program_scores <- function(df, n_bins = 40) {
  n_bins <- max(5, min(n_bins, floor(nrow(df) / (2 * 50))))
  brks <- unique(stats::quantile(df$pseudotime, probs = seq(0, 1, length.out = n_bins + 1),
                                  na.rm = TRUE))
  if (length(brks) < 3) return(NULL)
  df$.bin <- cut(df$pseudotime, breaks = brks, include.lowest = TRUE, labels = FALSE)
  df %>%
    filter(!is.na(.bin), !is.na(score)) %>%
    group_by(program, .bin) %>%
    summarise(pseudotime = stats::median(pseudotime),
              n = dplyr::n(),
              mean_score = mean(score),
              se = stats::sd(score) / sqrt(dplyr::n()),
              .groups = "drop") %>%
    mutate(lo = mean_score - 1.96 * se, hi = mean_score + 1.96 * se)
}

prog_binned <- bin_program_scores(prog_long)

make_program_plot <- function(binned, title, subtitle) {
  ggplot(binned, aes(x = pseudotime, y = mean_score,
                      color = program, fill = program)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = PROGRAM_COLORS) +
    scale_fill_manual(values = PROGRAM_COLORS) +
    labs(title = title, subtitle = subtitle,
         x = "Pseudotime (primary root: CytoTRACE2 cluster mean)",
         y = "Program score (mean of per-gene z-scores)",
         color = "Program", fill = "Program") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

.n_unordered <- sum(!is.finite(colData(cds)$pseudotime))
.sub <- paste0("Mean +/- 95% CI in equal-count pseudotime bins over ",
               format(nrow(prog_df), big.mark = ","),
               " cells with finite pseudotime (",
               format(.n_unordered, big.mark = ","),
               " cells outside the root partition are excluded).\n",
               "Orientation is a chosen root, not proven differentiation.")

if (!is.null(prog_binned)) {
  save_plot(make_program_plot(prog_binned,
              "Prespecified program scores along continuous pseudotime", .sub),
            "program_scores_along_pseudotime", "pseudotime",
            width = 7.5, height = 5.5)

  # The pooled script's by-outcome facet is degenerate here (one level), so it
  # is replaced by the same curves faceted on annot_6 -- does each tumour
  # subtype traverse the cohort's trajectory the same way?
  if (ANNOT_COL %in% colnames(prog_long)) {
    # prog_long is LONG: each cell contributes one row per program, so the
    # cell-count threshold is doubled here.
    .n_rows <- table(factor(as.character(prog_long[[ANNOT_COL]]),
                            levels = intersect(MALIGNANT_LEVELS,
                                    unique(as.character(prog_long[[ANNOT_COL]])))))
    .n_prog <- length(unique(prog_long$program))
    .keep   <- names(.n_rows)[.n_rows >= MIN_CELLS_PER_ANNOT5_FACET * .n_prog]
    .drop   <- setdiff(names(.n_rows), .keep)

    # The two cohorts legitimately end up with different panel counts: Poor has
    # no `Basal cells - Tumor` at all, and its `Epithelial cells - Basal - Tumor`
    # mostly lacks finite pseudotime. Record which classes were panelled so that
    # difference reads as a documented decision, not a missing panel.
    write.csv(data.frame(
                annot_5        = names(.n_rows),
                n_program_rows = as.integer(.n_rows),
                n_cells_approx = as.integer(.n_rows) %/% max(1L, .n_prog),
                min_cells_required = MIN_CELLS_PER_ANNOT5_FACET,
                included       = names(.n_rows) %in% .keep,
                stringsAsFactors = FALSE),
              file.path(OUT_DIR, "pseudotime",
                        sfx("program_scores_facet_inclusion.csv")),
              row.names = FALSE)

    prog_binned_annot <- do.call(rbind, lapply(
      split(prog_long, as.character(prog_long[[ANNOT_COL]])),
      function(d) {
        if (!as.character(d[[ANNOT_COL]][1]) %in% .keep) return(NULL)
        b <- bin_program_scores(d)
        if (is.null(b)) return(NULL)
        b[[ANNOT_COL]] <- as.character(d[[ANNOT_COL]][1])
        b
      }))

    if (!is.null(prog_binned_annot) && nrow(prog_binned_annot) > 0) {
      prog_binned_annot[[ANNOT_COL]] <- factor(prog_binned_annot[[ANNOT_COL]],
                                               levels = MALIGNANT_LEVELS)
      .sub_annot <- if (length(.drop))
        paste0(.sub, " | not panelled (< ", MIN_CELLS_PER_ANNOT5_FACET,
               " cells): ", paste0(.drop, " (n=",
               as.integer(.n_rows[.drop]) %/% max(1L, .n_prog), ")",
               collapse = ", "))
        else .sub
      save_plot(make_program_plot(prog_binned_annot,
                  paste0("Program scores along pseudotime by ", ANNOT_COL,
                         " - ", COHORT_LABEL, " cohort"), .sub_annot) +
                  facet_wrap(as.formula(paste("~", ANNOT_COL))),
                "program_scores_along_pseudotime_by_annot5", "pseudotime",
                width = 12, height = 7)
    }
  }

  write.csv(prog_binned,
            file.path(OUT_DIR, "pseudotime", sfx("program_scores_binned.csv")),
            row.names = FALSE)
  message("  Saved: pseudotime/program_scores_binned.csv (the plotted values)")
} else {
  message("  WARNING: too few ordered cells to bin program scores -- figure skipped")
}

# =============================================================================
# 8. Patient-level table, pseudotime bins, cluster fractions, stats
# =============================================================================

message("=== Step 8: Patient-level summary + stats ===")

cluster_assign <- as.character(monocle3::clusters(cds, reduction_method = "UMAP"))
colData(cds)$monocle_cluster <- cluster_assign

# TERMINOLOGY: these are Leiden clusters in the UMAP embedding. They are NOT
# branches of the principal graph -- the graph's topological branches are a
# different object entirely, and nothing in this script identifies them. The
# variables and outputs are named `cluster` throughout for that reason (they
# were previously, and misleadingly, named `branch`).
colData(cds)$monocle_partition <- as.character(
  monocle3::partitions(cds, reduction_method = "UMAP"))
colData(cds)$connected_to_root <- colData(cds)$monocle_partition == as.character(dominant_partition)

full_df <- as.data.frame(colData(cds)) %>%
  mutate(cell_id = colnames(cds),
         has_finite_pseudotime = is.finite(pseudotime))

# Global pseudotime terciles pooled across ALL malignant cells.
# IMPORTANT: pseudotime is Inf (not merely "large") for cells that are not
# reachable from the chosen root within the same principal-graph partition
# (a direct consequence of learn_graph(..., use_partition = TRUE, ...) --
# each partition gets its own principal tree, and order_cells() only
# propagates finite pseudotime within the root's own partition). Such cells
# have UNDEFINED pseudotime, not "late" pseudotime, so they must be excluded
# from the tercile binning entirely rather than being swept into the top bin
# by cut()'s upper break at +Inf.
finite_pt_idx <- is.finite(full_df$pseudotime)
pt_finite <- full_df$pseudotime[finite_pt_idx]
full_df$pt_bin <- NA_character_
if (length(pt_finite) >= 3) {
  pt_breaks <- stats::quantile(pt_finite, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  pt_breaks[1] <- -Inf; pt_breaks[length(pt_breaks)] <- Inf
  full_df$pt_bin[finite_pt_idx] <- as.character(cut(
    pt_finite, breaks = unique(pt_breaks),
    labels = c("early", "mid", "late")[seq_len(length(unique(pt_breaks)) - 1)],
    include.lowest = TRUE
  ))
}
n_no_pseudotime <- sum(!finite_pt_idx)
if (n_no_pseudotime > 0) {
  message("  NOTE: ", n_no_pseudotime, " / ", nrow(full_df), " cells have ",
          "non-finite pseudotime (unreachable from the root within its ",
          "principal-graph partition) and are excluded from pseudotime-tercile ",
          "and per-patient pseudotime summaries. This is expected behavior of ",
          "learn_graph(use_partition = TRUE) when the UMAP graph fragments into ",
          "multiple partitions (e.g. strong per-patient clustering of malignant ",
          "cells, common in tumor scRNA-seq without batch integration).")
}

# --- pseudotime_bin_summary.csv: per patient fraction in each global tercile ---
pt_bin_summary <- full_df %>%
  filter(!is.na(pt_bin)) %>%
  count(batch, outcome_group, pt_bin, name = "n") %>%
  group_by(batch) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  select(batch, outcome_group, pt_bin, n, frac) %>%
  tidyr::pivot_wider(names_from = pt_bin, values_from = c(n, frac), values_fill = 0)

write.csv(pt_bin_summary,
          file.path(OUT_DIR, "pseudotime", sfx("pseudotime_bin_summary.csv")),
          row.names = FALSE)
message("  Saved: pseudotime/pseudotime_bin_summary.csv")

# --- Monocle cluster fraction by sample / outcome -----------------------------
cluster_counts <- full_df %>%
  count(batch, outcome_group, monocle_cluster, name = "n") %>%
  group_by(batch) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p_cluster_sample <- ggplot(cluster_counts,
                            aes(x = batch, y = frac, fill = monocle_cluster)) +
  geom_col() +
  coord_flip() +
  labs(title = "Monocle cluster fraction by sample", x = "Sample (batch)",
       y = "Fraction of malignant cells", fill = "Monocle cluster") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
save_plot(p_cluster_sample, "cluster_fraction_by_sample", "clusters",
          width = 8, height = max(5, 0.25 * length(unique(cluster_counts$batch))))
# --- gene_modules outputs ---
write.csv(gene_test_res[, c("gene", "status", "p_value", "q_value",
                             intersect("morans_test_statistic", colnames(gene_test_res)),
                             intersect("morans_I", colnames(gene_test_res)))],
          file.path(OUT_DIR, "gene_modules", sfx("graph_test_results.csv")), row.names = FALSE)
message("  Saved: gene_modules/graph_test_results.csv")

# Moran's I over tens of thousands of cells is so overpowered that the q-value
# does essentially no filtering work (the previous run called 9,817 / 10,000
# genes significant at q < 0.05). Significance alone is therefore not a usable
# criterion for "trajectory-associated gene". This companion table ranks by
# EFFECT SIZE (Moran's I) instead, and is the one to read when asking which
# genes actually vary along the trajectory. The module computation above is
# deliberately left on the q-value criterion so this run stays comparable to the
# previous one.
if ("morans_I" %in% colnames(gene_test_res)) {
  moran_tbl <- gene_test_res[order(-gene_test_res$morans_I), ]
  moran_tbl <- moran_tbl[, intersect(c("gene", "status", "morans_I",
                                        "morans_test_statistic", "p_value", "q_value"),
                                      colnames(moran_tbl))]
  moran_tbl$morans_I_rank <- seq_len(nrow(moran_tbl))
  write.csv(utils::head(moran_tbl, 200),
            file.path(OUT_DIR, "gene_modules", sfx("trajectory_genes_top200_by_moransI.csv")),
            row.names = FALSE)
  .mi <- gene_test_res$morans_I
  message("  Saved: gene_modules/trajectory_genes_top200_by_moransI.csv")
  message("    Moran's I distribution: max = ", round(max(.mi, na.rm = TRUE), 3),
          ", n > 0.25 = ", sum(.mi > 0.25, na.rm = TRUE),
          ", n > 0.10 = ", sum(.mi > 0.10, na.rm = TRUE),
          ", n at q < 0.05 = ", length(sig_genes), " / ", nrow(gene_test_res))
} else {
  message("  NOTE: graph_test returned no morans_I column -- effect-size table skipped")
}

if (!is.null(gene_module_df)) {
  write.csv(as.data.frame(gene_module_df),
            file.path(OUT_DIR, "gene_modules", sfx("gene_modules.csv")), row.names = FALSE)
  message("  Saved: gene_modules/gene_modules.csv")

  # module scores aggregated by pseudotime bin (heatmap)
  cell_module_scores <- tryCatch({
    monocle3::aggregate_gene_expression(cds[sig_genes, ], gene_module_df)
  }, error = function(e) { message("  WARNING: aggregate_gene_expression failed: ",
                                    conditionMessage(e)); NULL })

  if (!is.null(cell_module_scores)) {
    mod_df <- as.data.frame(t(as.matrix(cell_module_scores)))
    mod_df$cell_id <- rownames(mod_df)
    mod_df <- merge(mod_df, full_df[, c("cell_id", "pt_bin")], by = "cell_id")
    mod_long <- tidyr::pivot_longer(mod_df, cols = -c(cell_id, pt_bin),
                                     names_to = "module", values_to = "score")
    mod_summary <- mod_long %>%
      filter(!is.na(pt_bin)) %>%
      group_by(module, pt_bin) %>%
      summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

    p_heat <- ggplot(mod_summary, aes(x = pt_bin, y = module, fill = mean_score)) +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(title = "Mean gene-module score by pseudotime bin",
           x = "Pseudotime bin (global terciles)", y = "Module",
           fill = "Mean score") +
      theme_bw(base_size = 11)
    save_plot(p_heat, "module_heatmap_by_pseudotime_bin", "gene_modules",
              width = 6, height = max(4, 0.3 * length(unique(mod_summary$module))))

    p_mod_umap <- traj_plot("pseudotime", "Pseudotime (reference for module UMAP)")
    save_plot(p_mod_umap, "module_score_umap", "gene_modules", width = 7, height = 6.5)
  }
} else {
  message("  gene_module_df is NULL -- skipping gene_modules.csv, module_heatmap, ",
          "and module_score_umap (too few significant genes on this run).")
}

# --- patient_level_summary_table.csv ---
# THIS IS THE SINGLE AUTHORITATIVE PATIENT-LEVEL EXPORT for the whole analysis.
# It is deliberately a SUPERSET of stage 1's patient_level_summary_cytotrace2.csv:
# the CytoTRACE2 and cell-cycle columns are recomputed here from colData(cds) --
# same cells, same global-quantile rule as cytotrace2_malignant_analysis.R -- so
# that one file answers every patient-level question (CytoTRACE2, cell cycle,
# pseudotime, tertile occupancy, cluster occupancy, unreachable-cell fraction)
# without anyone having to join two tables from two output trees. Stage 1's file
# keeps a stage-qualified name so the two can never shadow each other.
cluster_wide <- cluster_counts %>%
  select(batch, monocle_cluster, frac) %>%
  tidyr::pivot_wider(names_from = monocle_cluster, values_from = frac, values_fill = 0,
                      names_prefix = "frac_cluster_")

# Per-patient tertile occupancy, folded into the patient-level table so one file
# answers "where does this patient sit on the shared cohort trajectory?".
# Denominator is that patient's cells WITH finite pseudotime; cells outside the
# root partition are counted separately in frac_no_pseudotime rather than being
# silently absorbed into a bin.
tertile_wide <- full_df %>%
  filter(!is.na(pt_bin)) %>%
  count(batch, pt_bin, name = "n") %>%
  group_by(batch) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  select(batch, pt_bin, frac) %>%
  tidyr::pivot_wider(names_from = pt_bin, values_from = frac, values_fill = 0,
                      names_prefix = "frac_")

# n_cells = TOTAL malignant cells for that patient (the correct denominator
# for the "<10 cells" filter below), NOT just cells with finite pseudotime --
# filtering to finite pseudotime *before* grouping would silently drop entire
# patients whose cells all fall in a different principal-graph partition from
# the root (see note above on learn_graph(use_partition = TRUE)). Instead,
# median/mean pseudotime are computed only over the finite subset per patient
# and are NA for a patient with no finite-pseudotime cells at all -- an
# honest reflection of "pseudotime undefined for this patient's cells", not a
# missing row.
# frac_high_cytotrace_* thresholds are GLOBAL quantiles pooled across malignant
# cells, not per-patient quantiles -- the question is "what share of this
# patient's cells sit in the cohort's top decile/quartile/tercile of CytoTRACE2",
# so the cut points must be shared. The pool is restricted to patients that
# survive the <10-cell filter applied below, which is exactly the rule
# cytotrace2_malignant_analysis.R uses, so the two files agree numerically.
.keep_batches <- full_df %>% count(batch, name = "n") %>% filter(n >= 10) %>% pull(batch)
.cyto_pool <- full_df$CytoTRACE2_Score[full_df$batch %in% .keep_batches]
q_top10 <- quantile(.cyto_pool, probs = 0.90,   na.rm = TRUE)
q_top25 <- quantile(.cyto_pool, probs = 0.75,   na.rm = TRUE)
q_top33 <- quantile(.cyto_pool, probs = 0.6667, na.rm = TRUE)

has_phase <- "phase" %in% colnames(full_df)
if (!has_phase)
  message("  NOTE: no 'phase' column -- frac_G1/frac_S/frac_G2M omitted from the ",
          "patient-level table.")

patient_summary <- full_df %>%
  group_by(batch, outcome_group) %>%
  summarise(n_cells = n(),
            n_cells_with_pseudotime = sum(is.finite(pseudotime)),
            n_cells_no_pseudotime = sum(!is.finite(pseudotime)),
            frac_no_pseudotime = mean(!is.finite(pseudotime)),
            median_pseudotime = if (any(is.finite(pseudotime)))
              median(pseudotime[is.finite(pseudotime)]) else NA_real_,
            mean_pseudotime = if (any(is.finite(pseudotime)))
              mean(pseudotime[is.finite(pseudotime)]) else NA_real_,
            mean_cytotrace = mean(CytoTRACE2_Score, na.rm = TRUE),
            median_cytotrace = median(CytoTRACE2_Score, na.rm = TRUE),
            frac_high_cytotrace_top10 = mean(CytoTRACE2_Score >= q_top10, na.rm = TRUE),
            frac_high_cytotrace_top25 = mean(CytoTRACE2_Score >= q_top25, na.rm = TRUE),
            frac_high_cytotrace_top33 = mean(CytoTRACE2_Score >= q_top33, na.rm = TRUE),
            frac_G1  = if (has_phase) mean(phase == "G1",  na.rm = TRUE) else NA_real_,
            frac_S   = if (has_phase) mean(phase == "S",   na.rm = TRUE) else NA_real_,
            frac_G2M = if (has_phase) mean(phase == "G2M", na.rm = TRUE) else NA_real_,
            mean_score_myoepithelial = mean(score_myoepithelial, na.rm = TRUE),
            mean_score_ductal_notch = mean(score_ductal_notch, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(tertile_wide, by = "batch") %>%
  left_join(cluster_wide, by = "batch")

if (!has_phase)
  patient_summary <- patient_summary %>% select(-frac_G1, -frac_S, -frac_G2M)

n_before <- nrow(patient_summary)
patient_summary_full <- patient_summary
low_n <- patient_summary_full %>% filter(n_cells < 10)
if (nrow(low_n) > 0) {
  message("  WARNING: dropping ", nrow(low_n), " patient(s) with <10 malignant cells: ",
          paste(low_n$batch, collapse = ", "))
}
patient_summary <- patient_summary_full %>% filter(n_cells >= 10)

write.csv(patient_summary,
          file.path(OUT_DIR, "patient_level", sfx("patient_level_summary_table.csv")),
          row.names = FALSE)
message("  Saved: patient_level/patient_level_summary_table.csv (",
        nrow(patient_summary), " patients x ", ncol(patient_summary), " columns, ",
        n_before - nrow(patient_summary), " dropped for <10 cells)")
message("    This is the authoritative patient-level table: CytoTRACE2 + cell cycle ",
        "+ pseudotime + tertile + cluster fractions + frac_no_pseudotime.")
# =============================================================================
# 9. Full per-cell metadata export
# =============================================================================

message("=== Step 9: Per-cell metadata export ===")

# Monocle's own UMAP coordinates, alongside the original scVI ones already
# carried in colData -- needed to show the disconnected population in both
# embeddings and decide whether it was distinct before alignment.
.um <- SingleCellExperiment::reducedDims(cds)[["UMAP"]]
full_df$monocle_umap_1 <- .um[, 1]
full_df$monocle_umap_2 <- .um[, 2]

# QC / doublet / malignancy-profile fields, read from the h5ad obs table only.
message("  Pulling ", length(EXTRA_OBS_COLS), " extra obs columns from the h5ad")
extra_obs <- read_h5ad_obs_subset(H5AD_PATH, full_df$cell_id, EXTRA_OBS_COLS)
.n_ok <- vapply(extra_obs, function(x) sum(!is.na(x)), integer(1))
message("    Non-NA counts: ",
        paste(names(.n_ok), .n_ok, sep = "=", collapse = ", "))
for (cl in colnames(extra_obs)) full_df[[cl]] <- extra_obs[[cl]]

percell_export <- full_df %>%
  select(cell_id, batch, outcome_group,
         any_of(c(ANNOT_COL, ANNOT_COL_LEGACY)),
         monocle_cluster, monocle_partition, connected_to_root,
         pseudotime, has_finite_pseudotime, pt_bin,
         any_of(c("pseudotime_root_myo", "pseudotime_root_top5")),
         any_of("CytoTRACE2_Score"),
         score_myoepithelial, score_ductal_notch,
         any_of(c("phase", "S_score", "G2M_score")),
         any_of(c("X_umap_1", "X_umap_2")),
         monocle_umap_1, monocle_umap_2,
         any_of(c("n_genes", "total_counts", "pct_counts_mt")),
         any_of(EXTRA_OBS_COLS))

write.csv(percell_export,
          file.path(OUT_DIR, "tables", sfx("cds_per_cell_metadata.csv")),
          row.names = FALSE)
message("  Saved: tables/cds_per_cell_metadata.csv (", nrow(percell_export),
        " cells x ", ncol(percell_export), " columns)")

ANNOT5_SUMMARY_PATH <- file.path(OUT_DIR, "tables", sfx("annot5_class_summary.csv"))
# =============================================================================
# PER-CLASS SUMMARY BY annot_5
# =============================================================================
# One row per malignant class. This is what turns "the root is myoepithelial"
# into a reportable ordering along the axis, and it carries three caveats on the
# same row as the numbers they qualify:
#
#   frac_finite_pseudotime -- cells outside the root's partition have UNDEFINED
#     pseudotime and are silently absent from every pseudotime mean. Dropout is
#     not random across classes, so a class mean computed on 63% of its cells is
#     not comparable to one computed on 99%.
#   n_patients / n_patients_ge10 -- a class carried by one patient is a
#     statement about that patient. `Basal cells - Tumor` is 274/276 one patient.
#   mean_n_genes -- sequencing depth is an unadjusted confounder of CytoTRACE2
#     (ANALYSIS_NOTES.md 6.3). A per-class depth difference is the first
#     alternative explanation for a per-class potency difference.
.cd_sum <- as.data.frame(colData(cds))
.cd_sum$.a5 <- factor(as.character(.cd_sum[[ANNOT_COL]]), levels = MALIGNANT_LEVELS)
.cd_sum$.pt <- .cd_sum$pseudotime

.mean0 <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
.med0  <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)

annot5_summary <- do.call(rbind, lapply(levels(.cd_sum$.a5), function(lv) {
  d <- .cd_sum[!is.na(.cd_sum$.a5) & .cd_sum$.a5 == lv, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  fin <- is.finite(d$.pt)
  pat <- table(as.character(d$batch))
  data.frame(
    scope = SCOPE_LABEL,
    annot_5 = lv,
    n_cells = nrow(d),
    frac_of_malignant = nrow(d) / nrow(.cd_sum),
    n_patients = length(pat),
    n_patients_ge10 = sum(pat >= 10),
    top_patient = names(which.max(pat)),
    top_patient_frac = max(pat) / nrow(d),
    mean_cytotrace = .mean0(d$CytoTRACE2_Score),
    median_cytotrace = .med0(d$CytoTRACE2_Score),
    mean_score_myoepithelial = .mean0(d$score_myoepithelial),
    mean_score_ductal_notch = .mean0(d$score_ductal_notch),
    n_finite_pseudotime = sum(fin),
    frac_finite_pseudotime = mean(fin),
    mean_pseudotime = .mean0(d$.pt[fin]),
    median_pseudotime = .med0(d$.pt[fin]),
    frac_G1 = if ("phase" %in% names(d)) mean(d$phase == "G1", na.rm = TRUE) else NA_real_,
    frac_S = if ("phase" %in% names(d)) mean(d$phase == "S", na.rm = TRUE) else NA_real_,
    frac_G2M = if ("phase" %in% names(d)) mean(d$phase == "G2M", na.rm = TRUE) else NA_real_,
    mean_n_genes = if ("n_genes" %in% names(d)) .mean0(d$n_genes) else NA_real_,
    stringsAsFactors = FALSE, row.names = NULL)
}))
annot5_summary <- annot5_summary[order(annot5_summary$median_pseudotime), ]
write.csv(annot5_summary, ANNOT5_SUMMARY_PATH, row.names = FALSE)
message("  Saved: tables/", basename(ANNOT5_SUMMARY_PATH), " (",
        nrow(annot5_summary), " classes)")

# Save final cds object for downstream reuse
saveRDS(cds, file.path(OUT_DIR, "tables", sfx("cds_final.rds")))
message("  Saved: tables/cds_final.rds")

# =============================================================================
# 10. Run configuration / provenance
# =============================================================================

message("=== Step 10: Run configuration ===")

.n_finite <- sum(is.finite(colData(cds)$pseudotime))
.dom_frac <- mean(as.character(monocle3::partitions(cds, reduction_method = "UMAP")) ==
                    as.character(dominant_partition))

run_config <- data.frame(
  key = c("cohort", "outcome_label", "n_tumour_cells", "n_patients",
          "n_hvgs", "n_pca_dims", "alignment_status",
          "n_clusters", "n_partitions", "frac_in_dominant_partition",
          "root_cluster", "root_n_cells",
          "n_cells_finite_pseudotime", "frac_cells_finite_pseudotime",
          "n_graph_test_significant", "n_gene_modules",
          "source_cache", "source_cytotrace_scores", "run_timestamp",
          "annot_col", "n_annot5_classes_present",
          "root_top_annot5", "root_top_annot5_frac", "root_frac_myo_lineage",
          "root_top_patient", "root_top_patient_frac"),
  value = c(COHORT, COHORT_LABEL, as.character(ncol(cds)),
            as.character(length(unique(colData(cds)$batch))),
            as.character(length(malignant_hvgs)), as.character(n_dim),
            alignment_status,
            as.character(length(unique(cluster_assign_all))),
            as.character(n_partitions), as.character(round(.dom_frac, 4)),
            as.character(root_cluster), as.character(length(root_cluster_idx)),
            as.character(.n_finite),
            as.character(round(.n_finite / ncol(cds), 4)),
            as.character(length(sig_genes)),
            as.character(if (is.null(gene_module_df)) 0L else
              length(unique(gene_module_df$module))),
            CACHE_PATH, CYTO_CSV_PATH,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            # The headline answer, lifted onto the one file the driver greps.
            ANNOT_COL,
            as.character(length(unique(as.character(colData(cds)[[ANNOT_COL]])))),
            as.character(.rs_primary$root_top_annot5),
            as.character(round(.rs_primary$root_top_annot5_frac, 4)),
            as.character(round(.rs_primary$root_frac_myo_lineage, 4)),
            as.character(.rs_primary$root_top_patient),
            as.character(round(.rs_primary$root_top_patient_frac, 4))),
  stringsAsFactors = FALSE
)
write.csv(run_config, file.path(OUT_DIR, "tables", sfx("run_config_monocle3.csv")),
          row.names = FALSE)
message("  Saved: tables/", sfx("run_config_monocle3.csv"))
print(run_config)

# --- Terminal health checks, stated in the log so the driver can grep them ---
if (identical(alignment_status, "FAILED_USED_PCA")) {
  message("*** alignment_status = FAILED_USED_PCA -- batch integration failed ",
          "at every alignment_k. The trajectory in this tree is built on an ",
          "UNALIGNED embedding and must NOT be interpreted. ***")
} else {
  message("  align_cds status: ", alignment_status)
}
if (.dom_frac < 0.90)
  message("*** WARNING: the dominant principal-graph partition holds only ",
          round(100 * .dom_frac, 1), "% of cells. The graph is fragmented and ",
          "pseudotime is undefined for the remainder. ***")
if (.n_finite / ncol(cds) < 0.90)
  message("*** WARNING: only ", round(100 * .n_finite / ncol(cds), 1),
          "% of cells have finite pseudotime. ***")

message("=== monocle3_cohort.R complete for the ", COHORT_LABEL, " cohort ===")
message("  Outputs: ", OUT_DIR)
