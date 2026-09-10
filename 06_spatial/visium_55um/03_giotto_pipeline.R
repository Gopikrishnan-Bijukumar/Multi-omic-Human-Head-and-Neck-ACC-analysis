# =============================================================================
# Regular Visium (55 um) — per-sample Giotto pipeline
#
# RETIREMENT NOTE. CytoSPACE and Interaction Changed Genes/Features (ICG/ICF)
# were removed from this pipeline on 2026-08-11, bringing it into line with the
# Visium HD (16 um) pipeline. Unlike the HD side — where CytoSPACE never ran and
# the code was dead — CytoSPACE genuinely RAN here and produced real results.
# Nothing was deleted: the results, the retired code, verbatim pre-strip copies
# of this script, and the full rationale live in
#     archive/cytospace_icf/ARCHIVE_CONTEXT.md
# Read that before considering re-adding any of it. Do not re-add it here: a
# 55 um spot holds ~5-20 cells, which cannot support per-cell-type interaction
# attribution. Cell-type COMPOSITION (Section 8c) is what replaced it.
#
# PARAMETERISATION. This one script runs every 55 um sample. The four original
# per-sample scripts were identical apart from three Section 1 values, which now
# come from config/spatial_visium_samples.csv, keyed by the sample given on the
# command line:
#     module_k            metagene module count, chosen per sample from the
#                         Section 6b silhouette/WSS sweep (5 / 6 / 4 / 4)
#     svg_plot_midpoint   display-only gradient midpoint on the top-SVG plots;
#                         NA reproduces the one sample finalised without it
#     outcome             clinical outcome group
#
# Usage:  Rscript 03_giotto_pipeline.R <SAMPLE>
# =============================================================================

# set working directory
# Resolve this script's own directory so the repo config is found wherever the
# repository lives, and before setwd() moves into the data tree.
.this <- sub("^--file=", "",
             grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
SCRIPT_DIR <- if (is.na(.this)) getwd() else dirname(normalizePath(.this))
REPO_ROOT  <- dirname(dirname(SCRIPT_DIR))

# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) != 1L) stop("usage: Rscript 03_giotto_pipeline.R <SAMPLE>")
SAMPLE_NAME <- .args[[1]]

.cfg <- read.csv(file.path(REPO_ROOT, "config", "spatial_visium_samples.csv"),
                 stringsAsFactors = FALSE)
if (!SAMPLE_NAME %in% .cfg$sample)
  stop(sprintf("sample %s is not in config/spatial_visium_samples.csv", SAMPLE_NAME))
.row <- .cfg[.cfg$sample == SAMPLE_NAME, ]

setwd(file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace"))


# -----------------------------------------------------------------------------
# Paths. BASE_DIR and OUT_DIR are FROZEN and must not be renamed: saveGiotto()
# writes absolute paths INSIDE gobject.RDS (verified — decompressing a saved
# object shows .../giotto_results/<SAMPLE> and its .giotto_scratch shapefile
# paths). The folder is still called "giotto_w_cytospace" after CytoSPACE was
# retired; that name is a frozen historical artifact, left deliberately, not an
# oversight. data/, scripts/ and logs/ move freely — no saved object refers to
# them.
# -----------------------------------------------------------------------------
BASE_DIR    <- file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace")
DATA_DIR    <- file.path(BASE_DIR, "data")
VISIUM_DIR  <- file.path(DATA_DIR, paste0("sample_", SAMPLE_NAME))
SC_REF      <- file.path(DATA_DIR, "sc_ref", "inprogress_3.h5ad")
GENESET_DIR <- file.path(DATA_DIR, "genesets")
OUT_DIR     <- file.path(BASE_DIR, "giotto_results", SAMPLE_NAME)
CKPT_DIR    <- file.path(OUT_DIR, "checkpoints")
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(CKPT_DIR, recursive = TRUE, showWarnings = FALSE)

OUTCOME_MAP    <- setNames(.cfg$outcome, .cfg$sample)
SAMPLE_OUTCOME <- OUTCOME_MAP[[SAMPLE_NAME]]
SC_ANNOT_COL   <- "annot_8_c2l"
N_CORES        <- 16

# -----------------------------------------------------------------------------
# Cross-section constants.
#
# THE RULE: any constant read outside the section that defines it must live
# here, because Section 1 is the only code every execution path traverses.
# This is the single most common bug shape in this codebase. Concretely,
# CT_COLUMN used to be defined in Section 8.6b (a CytoSPACE section) and read by
# Section 9 — deleting 8.6b would have left Section 9 referencing an undefined
# symbol. MODULE_K was hardcoded in ~13 places across Sections 6/6b, where three
# of the four scripts had comments that disagreed with the code.
# -----------------------------------------------------------------------------
CT_COLUMN   <- "dominant_celltype"  # read by Sections 8c, 9
N_PROX_SIMS <- 1000L                # read by Section 9
MODULE_K          <- as.integer(.row$module_k)      # metagene modules, from Section 6b
SVG_PLOT_MIDPOINT <- as.numeric(.row$svg_plot_midpoint)  # display only; NA = omit

# Section 8c spatial blocking. Re-derived for 55 um — this is NOT a units
# conversion from the HD side. Regular Visium spots sit on a 100 um pitch, so
# HD's 25x25-bin / min-100 block is unsatisfiable here (a 400 um block holds at
# most 16 spots, so every block would be dropped and the CSV would be empty).
# Measured over all four samples' real cell_IDs:
#   4x4 (400um): 284-312 blocks, median occupancy 14-15  -> 1/16 quantisation, too coarse
#   6x6 (600um): 138-142 blocks, median occupancy 28-33  <- chosen
#   8x8 (800um):  79- 80 blocks, median occupancy 47-57  -> too few blocks
COMPOSITION_BLOCK_SPOTS <- 6L   # block edge, in spots (600 um at 100 um pitch)
COMPOSITION_MIN_BLOCK   <- 18L  # of a possible 36; drops sparse edge blocks

# One mode flag. Skips Sections 2-8c and regenerates 6e/8/8b/8c from the
# after_s8 checkpoint by calling the SAME functions the full run calls.
RESUME_S8 <- nzchar(Sys.getenv("RESUME_S8"))

# Helper: build save_param list for PNG (600 dpi) or SVG
sparam <- function(name, fmt, ...) {
  base <- list(save_name = name, save_format = fmt, ...)
  if (fmt == "png") base$dpi <- 600
  base
}

# =============================================================================
# SECTION 1 — LOAD LIBRARIES AND SET UP GIOTTO
# =============================================================================
suppressPackageStartupMessages({
  library(Giotto)
  library(GiottoData)
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(viridis)    # magma palette
  library(png)        # readPNG / writePNG for grayscale conversion
  library(jsonlite)   # read scalefactors_json
  library(cluster)
  library(patchwork)
  library(ggrastr)    # rasterised point layers for the large spatial figures
})

# python_path = NULL is correct on this server: Giotto auto-discovers the
# Python-only conda env (miniconda3/envs/giotto_env, with leidenalg/igraph).
# The R packages themselves come from ~/R/x86_64-pc-linux-gnu-library/4.4 via
# ~/.Renviron, so this script must be run with plain `Rscript` — do NOT
# `conda activate giotto_env` first, that would shadow the system R.
giotto_instrs <- createGiottoInstructions(
  save_dir    = OUT_DIR,
  save_plot   = TRUE,
  show_plot   = FALSE,
  return_plot = FALSE,
  python_path = NULL
)

# =============================================================================
# SECTION 1b — INSTRUMENTATION, THREAD PINNING, CHECKPOINTS, RESULT CACHING
# =============================================================================

# Pin BLAS/OpenMP threads. This box has no cgroup limit, so an unpinned BLAS
# will happily oversubscribe every core when several samples run back to back.
Sys.setenv(OMP_NUM_THREADS      = N_CORES,
           OPENBLAS_NUM_THREADS = N_CORES,
           MKL_NUM_THREADS      = N_CORES)
data.table::setDTthreads(N_CORES)

# [MEM] instrumentation. Cheap, and it is the evidence behind every runtime and
# memory claim in the notes. Kept from the HD pipeline; the HD *gating* logic
# (gene-chunked HVF/binSpect) is deliberately NOT ported — at ~4k spots the
# dense matrix is well under 1 GB and none of it is load-bearing here.
.MEM_T0 <- Sys.time()
mem_mark <- function(label) {
  gc_df <- gc(full = FALSE)
  used  <- sum(gc_df[, "used"] * c(8, 8)[seq_len(nrow(gc_df))]) / 1024^2
  rss   <- tryCatch({
    as.numeric(strsplit(readLines("/proc/self/statm", warn = FALSE)[1], " ")[[1]][2]) *
      4096 / 1024^3
  }, error = function(e) NA_real_)
  cat(sprintf("  [MEM] %-52s R=%6.2f GB  RSS=%6.2f GB  t=%6.1f min\n",
              label, used / 1024, rss,
              as.numeric(difftime(Sys.time(), .MEM_T0, units = "mins"))))
  invisible(NULL)
}

# Cache an expensive result to CKPT_DIR keyed by name. Used for the 1000-sim
# proximity result, which is deterministic given the object and takes minutes.
cached <- function(name, expr) {
  f <- file.path(CKPT_DIR, paste0(name, ".rds"))
  if (file.exists(f)) {
    cat("  [CACHE] hit:", name, "\n")
    return(readRDS(f))
  }
  val <- force(expr)
  saveRDS(val, f)
  cat("  [CACHE] stored:", name, "\n")
  val
}

checkpoint_gobject <- function(gobj, name) {
  d <- file.path(CKPT_DIR, name)
  tryCatch({
    saveGiotto(gobject = gobj, foldername = name, dir = CKPT_DIR, overwrite = TRUE)
    cat("  [CKPT] saved:", d, "\n")
  }, error = function(e) {
    cat("  [CKPT] !! save FAILED for ", name, ": ", conditionMessage(e), "\n", sep = "")
    warning("checkpoint save failed for ", name, ": ", conditionMessage(e))
  })
  invisible(gobj)
}

load_checkpoint <- function(name) {
  d <- file.path(CKPT_DIR, name)
  if (!dir.exists(d)) stop("checkpoint not found: ", d,
                           " — run without RESUME_S8 first.")
  loadGiotto(path_to_folder = d)
}

# =============================================================================
# HELPER REGION — canonical cell-type colours, figure savers, and the Section
# 6e/8/8b/8c bodies.
#
# Everything here is defined ABOVE the `if (!RESUME_S8)` block on purpose. The
# resume path calls these SAME functions rather than copies, so an output can
# never drift between a full run and the fast loop. Do not let a "quick"
# resume-path variant appear.
# =============================================================================

# -----------------------------------------------------------------------------
# Canonical cell-type palette, shared with the Visium HD pipeline so the same
# cell type is the same colour in both. Colours are grouped semantically
# (tumour = warm, epithelium = purple, stroma = green/teal, immune = blue/yellow)
# so the map is readable without the legend.
CT_PALETTE <- c(
  # tumour — warm
  ADC___Tumor              = "#C1272D",
  ADMEC___Tumor            = "#E85D75",
  Basal_cells___Tumor      = "#8B3A62",
  Epithelial_cells___Tumor = "#F07C24",
  # normal epithelium
  Epithelial_cells         = "#7B3FA0",
  # stroma — green/teal
  Fibroblast_cells         = "#2E8B57",
  Mural_muscle_cells       = "#9ACD32",
  Endothelial_cells        = "#17A2A2",
  # immune — blue/yellow
  Myeloid_cells            = "#1F6FB4",
  T_cells                  = "#6BAED6",
  B_Plasma_cells           = "#FFD92F"
)
# Pale blue-grey, deliberately absent from both backgrounds this column is drawn
# over (white tissue voids in Section 8, greyscale H&E in Section 12).
CT_NA_COLOR <- "#CFD8DC"
CT_NA_LABEL <- "No majority"

# Build the plotting column: NA becomes a real, labelled factor level rather
# than being fought through ggplot's na.value. The ORIGINAL dominant_celltype
# column keeps real NA — Sections 8c/9 depend on is.na() working there.
#
# NOTE FOR REGULAR VISIUM: this level is used far less than on the HD side.
# n_cell = 20 makes the DWLS purity gate 1/20 = 5%, the shipped rows sum to 1,
# and dominant_celltype is an argmax, so every SOLVED spot gets a call — versus
# ~10% abstention on HD, whose 50% gate leaves most bins with no majority.
# It is NOT always empty here: a spot whose DWLS row is entirely zero has no
# solution and is abstained in Section 8. Measured across the cohort:
# P16 0/4150, P01 0/4021, P09 121/3341 (3.6%),
# P22 0/3328.
ct_plot_factor <- function(x) {
  lv <- c(names(CT_PALETTE), CT_NA_LABEL)
  unknown <- setdiff(unique(x[!is.na(x)]), names(CT_PALETTE))
  if (length(unknown)) {
    warning("CT_PALETTE has no colour for: ", paste(unknown, collapse = ", "),
            " — falling back to grey. Update CT_PALETTE if the reference changed.")
    lv <- c(names(CT_PALETTE), unknown, CT_NA_LABEL)
  }
  factor(ifelse(is.na(x), CT_NA_LABEL, as.character(x)), levels = lv)
}
ct_plot_colors <- function(x) {
  unknown <- setdiff(unique(x[!is.na(x)]), names(CT_PALETTE))
  c(CT_PALETTE,
    if (length(unknown)) setNames(rep("#999999", length(unknown)), unknown),
    setNames(CT_NA_COLOR, CT_NA_LABEL))
}

# -----------------------------------------------------------------------------
# Save a ggplot as PNG (fully raster) + SVG with only the point/tile layers
# rasterised, so the SVG stays small but axes/titles/legends remain real vector
# text that can still be edited in Illustrator.
RASTER_DPI <- 600
save_rasterized_plot <- function(plot_obj, name, dir = OUT_DIR,
                                 width = 8, height = 7, units = "in",
                                 layers = c("Point", "Tile")) {
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot = plot_obj,
                  width = width, height = height, units = units,
                  dpi = RASTER_DPI, bg = "white")
  plot_rast <- ggrastr::rasterise(plot_obj, layers = layers, dpi = RASTER_DPI)
  ggplot2::ggsave(file.path(dir, paste0(name, ".svg")), plot = plot_rast,
                  width = width, height = height, units = units,
                  device = "svg", bg = "white")
}

# -----------------------------------------------------------------------------
# Section 8 figures, as a function so the RESUME_S8=1 path calls the SAME code.
#
# The 11-panel continuous-proportion figure is the same spatCellPlot2D over the
# DWLS enrichment it has always been. At n_cell = 20 those proportions ARE the
# compositional estimate, so there is nothing to replace them with — this differs
# from the HD pipeline, whose n_cell = 2 output is one-hot and needs a separate
# presence-panel rendering.
#
# CHANGED: point_size 1.75 -> 1.2. point_size is a PAGE size (mm on the figure),
# not a data size, so it says nothing about the real spot footprint — at 1.75 on
# these 16x12in / 4-column panels the dot rendered ~33.5px against a 23.2px spot
# pitch, i.e. 1.45x the centre-to-centre spacing. Spots therefore overlapped and
# OCCLUDED each other, so whichever spot was drawn last hid its neighbours'
# proportions (worst in Fibroblast_cells and ADC___Tumor, where high and zero
# spots interleave and the red looked more extensive than it is). Giotto draws
# shape 21 with point_border_stroke = 0.1, so the rendered diameter is
# 0.75 * (point_size * .pt + 0.1 * .stroke / 2) points; solving that for the
# 23.2px pitch gives 1.2. Spots now tile edge-to-edge with no occlusion.
# The value holds for all four samples: their pitches (199-203 full-res px) and
# tissue spans (12704-13769 units) differ by under 3%.
plot_dwls_figures <- function(gobj, dwls_result, cell_type_cols, dominant_ct) {
  for (fmt in c("png", "svg")) {
    tryCatch(
      spatCellPlot2D(gobj, spat_enr_names = "DWLS",
                     cell_annotation_values = cell_type_cols,
                     cow_n_col = 4, coord_fix_ratio = 1, point_size = 1.2,
                     save_param = c(sparam("DWLS_celltype_proportions", fmt),
                                    list(base_width = 16, base_height = 12))),
      error = function(e) {
        cat("  [8] !! spatCellPlot2D FAILED: ", conditionMessage(e), "\n", sep = "")
        warning("spatCellPlot2D failed: ", conditionMessage(e))
      }
    )
  }

  tryCatch({
    # Derive the plotting column here rather than upstream so the resume path is
    # self-sufficient: a checkpoint written before this column existed would
    # otherwise error, the handler would downgrade it to a deferred warning, and
    # the STALE figure from the previous run would be left on disk looking fresh.
    if (!"dominant_celltype_plot" %in% colnames(pDataDT(gobj))) {
      gobj <- addCellMetadata(gobj,
        new_metadata   = data.table(cell_ID                = dwls_result$cell_ID,
                                    dominant_celltype_plot = ct_plot_factor(dominant_ct)),
        by_column      = TRUE,
        column_cell_ID = "cell_ID")
    }
    p <- spatPlot2D(gobj, cell_color = "dominant_celltype_plot",
                    color_as_factor = TRUE,
                    cell_color_code = ct_plot_colors(dominant_ct),
                    point_size = 2.5,
                    save_plot = FALSE, return_plot = TRUE, show_plot = FALSE)
    save_rasterized_plot(p, "dominant_celltype_spatial")
  }, error = function(e) {
    # cat as well as warning(): warnings are deferred to the end of the run, so a
    # figure failure would otherwise be invisible while the job is in flight.
    cat("  [8] !! dominant_celltype spatPlot2D FAILED: ", conditionMessage(e), "\n", sep = "")
    warning("dominant_celltype spatPlot2D failed: ", conditionMessage(e))
  })
}

# -----------------------------------------------------------------------------
# SECTION 8b — CELL-TYPE CALL VALIDATION
#
# Scores every DWLS call against reference markers that were HELD OUT of the
# DWLS signature, and assigns each type a confidence tier. Section 8c joins
# these on and composition_compare.R de-emphasises anything not `solid`, so a
# large apparent composition difference in a badly-called type cannot quietly
# carry the headline.
validate_celltype_calls <- function(gobj, common_genes) {
  cat("--- [8b] Cell-type call validation ---\n")
  tryCatch({
    deg_path <- file.path(DATA_DIR, "sc_ref", "annot_8_c2l_DEGs.csv")
    stopifnot(file.exists(deg_path))
    X_val  <- getExpression(gobj, values = "normalized", output = "matrix")
    md_val <- pDataDT(gobj)
    ct_val <- md_val$dominant_celltype[match(colnames(X_val), md_val$cell_ID)]
    keep   <- !is.na(ct_val)
    ind_ct <- Matrix::sparse.model.matrix(~ 0 + factor(ct_val[keep]))
    colnames(ind_ct) <- levels(factor(ct_val[keep]))
    n_ct_spots <- Matrix::colSums(ind_ct)
    types      <- colnames(ind_ct)

    # (1) held-out DEG markers. `common_genes` is the exact DWLS signature, so
    #     excluding it makes this a genuinely independent check.
    deg <- fread(deg_path)
    deg[, cell_type := gsub("[^A-Za-z0-9_.]", "_", cell_type)]
    deg <- deg[!gene %in% common_genes & gene %in% rownames(X_val)]
    N_VAL_MARKERS <- 50L
    mk_deg <- deg[order(cell_type, -score)][, head(.SD, N_VAL_MARKERS), by = cell_type]

    # (2) canonical panel. PTPRC is deliberately absent (pan-leukocyte, so
    #     myeloid outranking T cells on it is correct biology, not a failure).
    canon_list <- list(
        ADC___Tumor              = "MYB",
        ADMEC___Tumor            = character(0),   # no accepted canonical marker
        Basal_cells___Tumor      = c("KRT5", "TP63"),
        Epithelial_cells___Tumor = c("EPCAM", "KRT8"),
        Epithelial_cells         = c("KRT19", "KRT7"),
        Fibroblast_cells         = c("COL1A1", "DCN", "LUM"),
        Mural_muscle_cells       = c("ACTA2", "MYH11", "TAGLN"),
        Endothelial_cells        = c("PECAM1", "VWF", "CDH5"),
        Myeloid_cells            = c("LYZ", "CD68", "AIF1"),
        T_cells                  = c("CD3E", "CD3D", "CD2"),
        B_Plasma_cells           = c("IGKC", "JCHAIN", "MZB1"))
    mk_canon <- rbindlist(lapply(names(canon_list), function(k)
        data.table(cell_type = k, gene = intersect(canon_list[[k]], rownames(X_val)))))

    score_sets <- function(mk) {
      mk <- mk[cell_type %in% types & gene %in% rownames(X_val)]
      setl <- lapply(types, function(tt) unique(mk[cell_type == tt, gene]))
      names(setl) <- types
      ok <- vapply(setl, length, 1L) > 0
      if (!any(ok)) return(NULL)
      gi <- Matrix::sparseMatrix(
        i = match(unlist(setl[ok]), rownames(X_val)),
        j = rep(seq_len(sum(ok)), times = vapply(setl[ok], length, 1L)),
        x = 1, dims = c(nrow(X_val), sum(ok)),
        dimnames = list(rownames(X_val), types[ok]))
      sums  <- as.matrix(Matrix::t(gi) %*% (X_val[, keep, drop = FALSE] %*% ind_ct))
      denom <- outer(vapply(setl[ok], length, 1L), n_ct_spots, "*")
      m <- sums / denom
      # z-score WITHIN each marker set so rows are comparable to each other
      t(apply(m, 1, function(r) if (sd(r) > 0) (r - mean(r)) / sd(r) else r * 0))
    }
    Z_deg   <- score_sets(mk_deg)
    Z_canon <- score_sets(mk_canon)
    rank_of <- function(Z, tt) if (is.null(Z) || !tt %in% rownames(Z)) NA_integer_
                               else which(colnames(Z)[order(-Z[tt, ])] == tt)

    val <- rbindlist(lapply(types, function(tt) {
      rd <- rank_of(Z_deg, tt); rc <- rank_of(Z_canon, tt)
      # Two-source agreement rule. A type is only downgraded when EVERY available
      # source agrees — but a missing source must not launder a bad call, so when
      # only one source exists (ADMEC___Tumor has no accepted canonical marker)
      # that single source decides on its own.
      rk <- c(rd, rc); rk <- rk[!is.na(rk)]
      tier <- if (!length(rk))        "untested"
              else if (all(rk == 1L)) "solid"
              else if (all(rk >  3L)) "unreliable"
              else                    "ambiguous"
      data.table(cell_type = tt, n_spots = as.integer(n_ct_spots[[tt]]),
                 diag_z_deg   = if (!is.null(Z_deg)   && tt %in% rownames(Z_deg))   Z_deg[tt, tt]   else NA_real_,
                 rank_deg     = rd,
                 diag_z_canon = if (!is.null(Z_canon) && tt %in% rownames(Z_canon)) Z_canon[tt, tt] else NA_real_,
                 rank_canon   = rc,
                 n_markers_deg   = sum(mk_deg$cell_type   == tt),
                 n_markers_canon = sum(mk_canon$cell_type == tt),
                 best_competitor = if (!is.null(Z_deg) && tt %in% rownames(Z_deg))
                                     colnames(Z_deg)[order(-Z_deg[tt, ])][1] else NA_character_,
                 confidence_tier = tier)
    }))[order(factor(confidence_tier,
                     levels = c("solid", "ambiguous", "unreliable", "untested")), -n_spots)]
    if (!is.null(Z_deg)) {
      val[, margin := diag_z_deg - vapply(seq_len(.N), function(i)
            max(Z_deg[cell_type[i], setdiff(colnames(Z_deg), cell_type[i])]), 0.0)]
    }
    fwrite(val, file.path(OUT_DIR, "celltype_validation.csv"))
    print(val[, .(cell_type, n_spots, rank_deg, rank_canon, best_competitor, confidence_tier)])

    bad <- val[confidence_tier != "solid"]
    if (nrow(bad)) {
      cat("\n  *** CELL-TYPE CALL WARNING ***\n")
      for (i in seq_len(nrow(bad)))
        cat(sprintf("  %-26s %-10s n=%6d  DEG rank %s, canonical rank %s (best: %s)\n",
                    bad$cell_type[i], bad$confidence_tier[i], bad$n_spots[i],
                    bad$rank_deg[i], ifelse(is.na(bad$rank_canon[i]), "-", bad$rank_canon[i]),
                    bad$best_competitor[i]))
      cat("  Section 8c/9 results involving these types carry that uncertainty;\n")
      cat("  composition_compare.R de-emphasises them in the Poor-vs-Good figures.\n\n")
    }

    # Concordance heatmap — the diagonal should win every row.
    tryCatch({
      zl <- as.data.table(as.table(Z_deg))
      setnames(zl, c("marker_set", "called_celltype", "z"))
      ph <- ggplot2::ggplot(zl, ggplot2::aes(called_celltype, marker_set, fill = z)) +
        ggplot2::geom_tile() +
        ggplot2::geom_tile(data = zl[as.character(marker_set) == as.character(called_celltype)],
                           colour = "black", linewidth = 0.5, fill = NA) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                      midpoint = 0, name = "z") +
        ggplot2::labs(title = "Cell-type call validation (held-out reference markers)",
                      subtitle = "row = markers OF type; column = spots CALLED that type; boxed = diagonal",
                      x = "called cell type", y = "marker set") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      save_rasterized_plot(ph, "celltype_validation_heatmap", width = 9, height = 7)
    }, error = function(e) warning("validation heatmap: ", conditionMessage(e)))

    rm(X_val); gc(full = TRUE)
    mem_mark("08b cell-type call validation")
    val
  }, error = function(e) {
    cat("  [8b] !! validation FAILED: ", conditionMessage(e), "\n", sep = "")
    warning("Cell-type validation skipped: ", conditionMessage(e))
    NULL
  })
}

# -----------------------------------------------------------------------------
# SECTION 8c — CELL-TYPE COMPOSITION
#
# WHAT IT PRODUCES: the cell-type composition of the tissue, three ways, plus the
# per-block tables composition_compare.R turns into a Poor-vs-Good contrast.
#
# WHAT THE TWO ESTIMATES MEAN — and how this differs from the HD pipeline.
# n_cell in runDWLSDeconv is a PURITY GATE, not a cell count: types below
# 1/n_cell are zeroed and the row renormalised. The HD pipeline uses n_cell = 2,
# so its gate is 50%, at most one type can clear it, and its shipped output is
# one-hot — which is why HD has to recover the discarded continuous stage-1
# matrix separately.
#
# This pipeline uses n_cell = 20, correct for a 55 um spot holding ~5-20 cells.
# The gate is therefore 5%, and the shipped output is ALREADY the continuous
# compositional estimate. Verified on P16: row sums exactly 1.0, 1-5
# types per spot, only 0.38% of entries negative (vs 28% on HD). So:
#   mean_prop_continuous  = column means of the shipped DWLS matrix  <- PRIMARY
#   frac_spots_dominant   = fraction of spots whose argmax type is X <- secondary
# No stage-1 recovery is needed and runDWLSDeconv_keep_stage1() is deliberately
# NOT ported: it pins three unexported Giotto internals to buy a number this
# pipeline already has.
#
# UNCERTAINTY IS SPATIALLY BLOCKED, never per-spot. Neighbouring spots carry
# largely the same tissue's signal, so a binomial interval over ~4k spots is
# meaninglessly tight. Composition is therefore also computed per block and it is
# the between-block spread that carries the error bars. But note what blocking
# does and does not buy: it removes the absurd interval, yet no within-section
# blocking makes blocks truly independent, because tumour nests and stromal bands
# are themselves hundreds of um across. Treat the spread as a HETEROGENEITY
# DESCRIPTOR, not a standard error.

# Percent axis labels without pulling in scales:: — this pipeline has never
# depended on it directly and a figure is not worth a new import.
pct_lab <- function(x) paste0(signif(100 * x, 3), "%")

# Array row/col for every spot. Regular Visium barcodes are random 16-mers
# (GTCACTTCCTTCTAGA-1) and carry NO positional information — the HD pipeline's
# trick of parsing them out of the bin barcode (s_016um_<row>_<col>) cannot work
# here. They come from tissue_positions.csv instead.
#
# Do NOT substitute pixel coordinates: pxl_row_in_fullres is in full-resolution
# image pixels whose um scale differs per slide, so a block size in pixels is not
# comparable across samples. Array indices are on a fixed lattice.
spot_array_coords <- function(ids) {
  tp_path <- file.path(VISIUM_DIR, "spatial", "tissue_positions.csv")
  if (!file.exists(tp_path)) stop("tissue_positions.csv not found: ", tp_path)
  tp <- fread(tp_path)
  setnames(tp, tolower(names(tp)))
  stopifnot(all(c("barcode", "array_row", "array_col") %in% names(tp)))
  ac <- tp[match(ids, barcode), .(arow = array_row, acol = array_col)]
  stopifnot(!anyNA(ac$arow), !anyNA(ac$acol))
  ac
}

summarize_composition <- function(gobj, dwls_result, cell_type_cols, dominant_ct,
                                  ct_tiers = NULL) {
  cat("--- [8c] Cell-type composition ---\n")
  tryCatch({
    md      <- pDataDT(gobj)
    ids     <- dwls_result$cell_ID
    n_total <- length(dominant_ct)
    ct_call <- ifelse(is.na(dominant_ct), CT_NA_LABEL, dominant_ct)

    # ---- (1) majority-call composition, abstentions kept as a real category ---
    comp <- data.table(cell_type = ct_call)[, .(n_spots_dominant = .N), by = cell_type]
    comp[, frac_spots_dominant := n_spots_dominant / n_total]
    cat(sprintf("  %d spots total | %d assigned | %d unassigned (%.2f%%)\n",
                n_total, sum(!is.na(dominant_ct)), sum(is.na(dominant_ct)),
                100 * mean(is.na(dominant_ct))))

    # ---- (2) continuous composition, straight off the shipped DWLS matrix ----
    s1       <- as.matrix(dwls_result[, ..cell_type_cols])   # spots x cell types
    raw_sums <- rowSums(s1, na.rm = TRUE)
    n_neg    <- sum(s1 < 0, na.rm = TRUE)
    cat(sprintf("  DWLS row sums BEFORE clamp: median %.4f, IQR %.4f-%.4f, range %.4f-%.4f\n",
                median(raw_sums), quantile(raw_sums, .25), quantile(raw_sums, .75),
                min(raw_sums), max(raw_sums)))
    cat(sprintf("  QP solver noise: %d of %d entries < 0 (%.3g%%), most negative %.3g\n",
                n_neg, length(s1), 100 * n_neg / length(s1),
                suppressWarnings(min(s1, na.rm = TRUE))))
    # The clamp+renormalise is a near-no-op at n_cell = 20 (rows already sum to
    # 1, ~0.4% of entries are negative at ~1e-15). It is kept anyway, and the
    # pre-clamp distribution above is printed, so any REAL drift away from 1.0
    # stays visible instead of being silently normalised away.
    s1  <- pmax(s1, 0)
    rs  <- rowSums(s1, na.rm = TRUE)
    ok  <- is.finite(rs) & rs > 0
    s1  <- sweep(s1[ok, , drop = FALSE], 1, rs[ok], "/")
    cat(sprintf("  continuous composition over %d of %d spots (%d dropped: zero/NA row)\n",
                sum(ok), length(rs), sum(!ok)))

    cont <- data.table(
      cell_type            = colnames(s1),
      mean_prop_continuous = colMeans(s1),
      sd_prop_continuous   = apply(s1, 2, sd),
      median_prop          = apply(s1, 2, median),
      q25_prop             = apply(s1, 2, quantile, .25),
      q75_prop             = apply(s1, 2, quantile, .75))
    comp <- merge(comp, cont, by = "cell_type", all = TRUE)
    comp[is.na(n_spots_dominant),   n_spots_dominant   := 0L]
    comp[is.na(frac_spots_dominant), frac_spots_dominant := 0]

    # ---- (3) join the Section 8b confidence tiers -----------------------------
    if (!is.null(ct_tiers) && "confidence_tier" %in% colnames(ct_tiers)) {
      comp <- merge(comp, ct_tiers[, .(cell_type, confidence_tier, rank_deg,
                                       rank_canon, best_competitor)],
                    by = "cell_type", all.x = TRUE)
    } else {
      comp[, confidence_tier := NA_character_]
    }
    # The abstention row is not a cell type and has no tier or marker rank.
    comp[cell_type == CT_NA_LABEL, confidence_tier := "not_a_celltype"]

    comp[, `:=`(sample = SAMPLE_NAME, outcome = SAMPLE_OUTCOME)]
    setorder(comp, -mean_prop_continuous)
    fwrite(comp, file.path(OUT_DIR, "composition_summary.csv"))
    print(comp[, .(cell_type, n_spots_dominant,
                   frac = round(frac_spots_dominant, 4),
                   mean_cont = round(mean_prop_continuous, 4), confidence_tier)])

    # ---- (4) composition within each Leiden domain ---------------------------
    if ("leiden_clus" %in% colnames(md)) {
      dl  <- merge(data.table(cell_ID = ids, cell_type = ct_call),
                   md[, .(cell_ID, leiden_clus)], by = "cell_ID")
      byl <- dl[, .N, by = .(leiden_clus, cell_type)]
      byl[, frac := N / sum(N), by = leiden_clus]
      byl[, `:=`(sample = SAMPLE_NAME, outcome = SAMPLE_OUTCOME)]
      setorder(byl, leiden_clus, -frac)
      fwrite(byl, file.path(OUT_DIR, "composition_by_leiden.csv"))
      cat(sprintf("  composition_by_leiden.csv: %d domains x %d types\n",
                  uniqueN(byl$leiden_clus), uniqueN(byl$cell_type)))
    } else {
      byl <- NULL
      cat("  leiden_clus absent — skipping per-domain composition.\n")
    }

    # ---- (5) per-block composition (the unit of uncertainty) -----------------
    blocks <- NULL
    tryCatch({
      ac <- spot_array_coords(ids)
      # THE HEX LATTICE. Regular Visium's array is a staggered hex grid, not a
      # square one: odd rows are offset by half a spot and array_col steps by 2
      # WITHIN a row (verified on this data: 0, 2, 4, 6, 8, 10...). Integer-
      # dividing array_col directly would therefore give blocks half as wide as
      # they are tall. Halve the column index first so blocks are square in
      # physical space.
      bt <- data.table(
        cell_type = ct_call,
        block = paste0(ac$arow %/% COMPOSITION_BLOCK_SPOTS, "_",
                       (ac$acol %/% 2L) %/% COMPOSITION_BLOCK_SPOTS))
      bt[, n_block := .N, by = block]
      n_small <- uniqueN(bt[n_block < COMPOSITION_MIN_BLOCK, block])
      bt <- bt[n_block >= COMPOSITION_MIN_BLOCK]
      if (!nrow(bt)) stop("every block fell below COMPOSITION_MIN_BLOCK (",
                          COMPOSITION_MIN_BLOCK, ") — block size needs re-deriving")
      blocks <- bt[, .(n_spots = .N, n_block = n_block[1]), by = .(block, cell_type)]
      # A type absent from a block is a 0, not a missing row — otherwise the
      # between-block mean is taken only over blocks where the type appears and
      # is biased upward. Complete the grid.
      blocks <- blocks[CJ(block = unique(blocks$block),
                          cell_type = unique(bt$cell_type), unique = TRUE),
                       on = .(block, cell_type)]
      blocks[is.na(n_spots), n_spots := 0L]
      blocks[, n_block := max(n_block, na.rm = TRUE), by = block]
      blocks[, frac := n_spots / n_block]
      # Stamp the geometry into the file. composition_compare.R titles its figure
      # from this rather than hardcoding a size that goes stale the moment
      # COMPOSITION_BLOCK_SPOTS changes.
      blocks[, `:=`(sample = SAMPLE_NAME, outcome = SAMPLE_OUTCOME,
                    block_spots = COMPOSITION_BLOCK_SPOTS,
                    block_um    = 100L * COMPOSITION_BLOCK_SPOTS)]
      setorder(blocks, block, -frac)
      fwrite(blocks, file.path(OUT_DIR, "composition_spatial_blocks.csv"))
      occ <- unique(blocks[, .(block, n_block)])$n_block
      cat(sprintf(paste0("  composition_spatial_blocks.csv: %d blocks of >=%d spots ",
                         "(%d sparse blocks dropped), %dx%d spots = %dx%d um\n"),
                  uniqueN(blocks$block), COMPOSITION_MIN_BLOCK, n_small,
                  COMPOSITION_BLOCK_SPOTS, COMPOSITION_BLOCK_SPOTS,
                  100 * COMPOSITION_BLOCK_SPOTS, 100 * COMPOSITION_BLOCK_SPOTS))
      cat(sprintf("    block occupancy: median %d, range %d-%d (of a possible %d)\n",
                  as.integer(median(occ)), min(occ), max(occ),
                  COMPOSITION_BLOCK_SPOTS^2))
    }, error = function(e) {
      cat("  [8c] !! blocked composition FAILED: ", conditionMessage(e), "\n", sep = "")
      warning("Section 8c blocked composition failed: ", conditionMessage(e))
    })

    # ---- (6) figures ---------------------------------------------------------
    # Stacked bar: both estimates side by side, so the difference between a
    # majority-call composition and a proportional one is visible rather than
    # implicit. Same canonical palette as every other cell-type figure.
    tryCatch({
      pl <- rbind(
        comp[, .(cell_type, value = frac_spots_dominant,
                 estimate = "Majority call\n(fraction of spots)")],
        comp[!is.na(mean_prop_continuous),
             .(cell_type, value = mean_prop_continuous,
               estimate = "DWLS proportion\n(mean over spots)")])
      pl[, cell_type := ct_plot_factor(
            ifelse(cell_type == CT_NA_LABEL, NA_character_, cell_type))]
      pl <- pl[, .(value = sum(value)), by = .(cell_type, estimate)]
      p_comp <- ggplot2::ggplot(pl, ggplot2::aes(x = estimate, y = value,
                                                 fill = cell_type)) +
        ggplot2::geom_col(width = 0.6, colour = "white", linewidth = 0.2) +
        ggplot2::scale_fill_manual(values = ct_plot_colors(dominant_ct),
                                   name = NULL, drop = FALSE) +
        ggplot2::scale_y_continuous(labels = pct_lab, expand = c(0, 0)) +
        ggplot2::labs(
          title    = sprintf("%s (%s outcome) — cell-type composition",
                             SAMPLE_NAME, SAMPLE_OUTCOME),
          subtitle = sprintf(paste0("%s spots at 55um. Left: fraction of spots by dominant ",
                                    "(argmax) type.\nRight: mean SpatialDWLS proportion at ",
                                    "n_cell = 20 (5%% purity gate) — the primary estimate."),
                             format(n_total, big.mark = ",")),
          x = NULL, y = "share of tissue") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
      save_rasterized_plot(p_comp, "composition_barplot", width = 9, height = 7)
    }, error = function(e) {
      cat("  [8c] !! composition barplot FAILED: ", conditionMessage(e), "\n", sep = "")
      warning("composition barplot failed: ", conditionMessage(e))
    })

    # Per-domain heatmap: where in the tissue each cell type concentrates.
    if (!is.null(byl)) tryCatch({
      ph <- ggplot2::ggplot(byl, ggplot2::aes(x = factor(leiden_clus),
                                              y = cell_type, fill = frac)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        ggplot2::scale_fill_viridis_c(option = "magma", direction = -1,
                                      labels = pct_lab, name = "share of\ndomain") +
        ggplot2::labs(title = sprintf("%s — cell-type composition per Leiden domain",
                                      SAMPLE_NAME),
                      subtitle = "column-wise shares; each domain sums to 100%",
                      x = "Leiden cluster", y = NULL) +
        ggplot2::theme_minimal(base_size = 10)
      save_rasterized_plot(ph, "composition_by_leiden_heatmap",
                           width = 10, height = 6)
    }, error = function(e) {
      cat("  [8c] !! per-domain heatmap FAILED: ", conditionMessage(e), "\n", sep = "")
      warning("composition per-domain heatmap failed: ", conditionMessage(e))
    })

    mem_mark("08c cell-type composition")
    invisible(comp)
  }, error = function(e) {
    # cat as well as warning(): warnings are deferred to the end of the run, so a
    # silent failure here would leave the PREVIOUS run's CSVs on disk looking new.
    cat("  [8c] !! composition summary FAILED: ", conditionMessage(e), "\n", sep = "")
    warning("Section 8c composition failed: ", conditionMessage(e))
    NULL
  })
}

# =============================================================================
# SECTION 6e HELPERS — METAGENE FUNCTIONAL ENRICHMENT
#
# Finding spatial co-expression modules says WHERE coordinated expression is,
# not WHAT process it represents. These functions answer the second question.
#
# Two analyses, deliberately:
#   ORA   — hypergeometric overlap of the module's own genes with a pathway.
#           Interpretable and auditable (overlap_genes names every driver), but
#           with small modules it has almost no power.
#   GSEA  — rank ALL genes by correlation with the per-spot metagene score and
#           test the whole transcriptome. This is the PRIMARY analysis: it
#           escapes the small-module problem and asks the better question —
#           what else in the tissue tracks this spatial region?
#
# BOTH ARE NON-FATAL. Missing genesets/ or missing fgsea skips this section with
# an actionable message; no optional renderer should ever abort a long job.
# =============================================================================
read_gmt <- function(path) {
  ln <- readLines(path, warn = FALSE); ln <- ln[nzchar(ln)]
  p  <- strsplit(ln, "\t", fixed = TRUE)
  setNames(lapply(p, function(x) unique(x[-c(1, 2)][nzchar(x[-c(1, 2)])])),
           vapply(p, `[`, "", 1L))
}

# Pearson r between a dense length-n vector y and every ROW of a sparse
# genes x spots matrix, WITHOUT densifying: every term comes from rowSums /
# rowSums-of-squares and ONE sparse matvec. sx/sxx are computed once and reused
# for every module.
sparse_pearson_stats <- function(X) {
  list(n = ncol(X), sx = Matrix::rowSums(X), sxx = Matrix::rowSums(X * X))
}
sparse_pearson <- function(X, y, st) {
  stopifnot(length(y) == st$n, !anyNA(y))
  xy  <- as.numeric(X %*% y)
  sy  <- sum(y); syy <- sum(y * y)
  vx  <- st$sxx - st$sx^2 / st$n
  vy  <- syy    - sy^2     / st$n
  r   <- (xy - st$sx * sy / st$n) / (sqrt(pmax(vx, 0)) * sqrt(vy))
  r[!is.finite(r) | vx <= 0] <- NA_real_
  setNames(r, rownames(X))
}

# Hypergeometric over-representation. This is bit-for-bit what
# clusterProfiler::enricher() computes; using phyper directly avoids the entire
# DOSE/AnnotationDbi/org.Hs.eg.db chain AND the symbol->Entrez mapping step,
# which silently drops 5-10% of genes. Identifiers here are HGNC symbols
# throughout, so no mapping is needed at all.
ora_one <- function(module_genes, sets, universe, min_size = 10L, max_size = 500L) {
  U  <- unique(universe); N <- length(U)
  mg <- intersect(unique(module_genes), U); n <- length(mg)
  if (n == 0L) return(NULL)
  res <- rbindlist(lapply(names(sets), function(s) {
    inU <- intersect(sets[[s]], U); K <- length(inU)
    if (K < min_size || K > max_size) return(NULL)
    ov <- intersect(mg, inU); x <- length(ov)
    data.table(term = s, set_size_in_universe = K, overlap = x,
               overlap_genes = paste(ov, collapse = ";"),
               fold_enrichment = (x / n) / (K / N),
               pvalue = stats::phyper(x - 1L, K, N - K, n, lower.tail = FALSE))
  }))
  if (!nrow(res)) return(NULL)
  res[, n_module := n][]
}

run_metagene_enrichment <- function(gobj) tryCatch({
  ENR_DIR <- file.path(OUT_DIR, "enrichment")
  dir.create(ENR_DIR, recursive = TRUE, showWarnings = FALSE)

  gmts <- list.files(GENESET_DIR, pattern = "\\.gmt$", full.names = TRUE)
  if (!length(gmts)) {
    cat("  [6e] no GMT files in ", GENESET_DIR, " — skipping enrichment.\n", sep = "")
    cat("       Run:  Rscript scripts/fetch_genesets.R   (needs internet)\n")
    return(invisible(NULL))
  }
  # Re-verify the SHA-256 recorded at download time. A silently swapped or
  # truncated GMT would otherwise produce plausible-looking but wrong biology.
  man_path <- file.path(GENESET_DIR, "MANIFEST.json")
  if (file.exists(man_path) && requireNamespace("jsonlite", quietly = TRUE) &&
      requireNamespace("digest", quietly = TRUE)) {
    man <- jsonlite::fromJSON(man_path)
    for (f in gmts) {
      bn <- basename(f)
      if (!is.null(man$files[[bn]]$sha256)) {
        got <- digest::digest(file = f, algo = "sha256")
        if (!identical(got, man$files[[bn]]$sha256))
          stop("GMT checksum mismatch for ", bn, " — genesets/ has changed since download.")
      }
    }
    cat("  [6e] gene sets verified against MANIFEST.json (release ",
        if (is.null(man$release)) "?" else man$release, ")\n", sep = "")
  }
  collections <- setNames(lapply(gmts, read_gmt),
                          sub("\\.v.*$", "", basename(gmts)))
  cat(sprintf("  [6e] %d collections: %s\n", length(collections),
              paste(sprintf("%s (%d sets)", names(collections),
                            vapply(collections, length, 1L)), collapse = ", ")))

  mod_path <- file.path(OUT_DIR, "metagene_module_genes.csv")
  if (!file.exists(mod_path)) {
    cat("  [6e] metagene_module_genes.csv not found — skipping.\n"); return(invisible(NULL))
  }
  mods <- fread(mod_path)
  svg_path <- file.path(OUT_DIR, "SVG_rank.csv")
  X  <- getExpression(gobj, values = "normalized", output = "matrix")
  # UNIVERSE 1 = every gene that survived QC and was testable by binSpect.
  # Answers "is this module enriched relative to everything measurable here?"
  U_qc  <- if (file.exists(svg_path)) fread(svg_path)$feats else rownames(X)
  # UNIVERSE 2 = the SVGs the modules were carved from. Answers the sharper
  # question "given a gene is ALREADY spatially variable, does this module differ
  # from its siblings?" — which controls for the fact that all module genes were
  # pre-selected for spatial variability. Low power by construction; reported as
  # annotation, not inference.
  U_svg <- unique(mods$feat)
  cat(sprintf("  [6e] universes: qc = %d genes | svg = %d genes\n",
              length(U_qc), length(U_svg)))

  # ---------------- ORA ------------------------------------------------------
  ora_all <- rbindlist(lapply(names(collections), function(cn) {
    rbindlist(lapply(c("qc_universe", "svg_universe"), function(un) {
      U <- if (un == "qc_universe") U_qc else U_svg
      rbindlist(lapply(sort(unique(mods$module)), function(k) {
        r <- ora_one(mods[module == k, feat], collections[[cn]], U,
                     min_size = if (un == "qc_universe") 10L else 5L,
                     max_size = if (un == "qc_universe") 500L else Inf)
        if (is.null(r)) return(NULL)
        # BH within (module, collection, universe): the collections are families
        # of wildly different size (50 Hallmark vs ~7,500 GO:BP) and pooling
        # would let GO:BP's mass annihilate Hallmark. Every set is TESTED so the
        # denominator stays honest; only overlap >= 3 is REPORTED.
        r[, `:=`(module = k, collection = cn, universe = un,
                 padj_within = p.adjust(pvalue, "BH"))][]
      }))
    }))
  }), fill = TRUE)
  if (nrow(ora_all)) {
    ora_all[, padj_global := p.adjust(pvalue, "BH"), by = universe]
    for (cn in unique(ora_all$collection)) for (un in unique(ora_all$universe)) {
      out <- ora_all[collection == cn & universe == un & overlap >= 3L][order(module, pvalue)]
      if (nrow(out)) fwrite(out[, .(module, term, n_module, set_size_in_universe, overlap,
                                    overlap_genes, fold_enrichment, pvalue,
                                    padj_within, padj_global)],
                            file.path(ENR_DIR, sprintf("ora_%s_%s.csv", cn, un)))
    }
    cat(sprintf("  [6e] ORA: %d tests, %d reported (overlap>=3), %d significant (padj_within<0.05)\n",
                nrow(ora_all), nrow(ora_all[overlap >= 3L]),
                nrow(ora_all[overlap >= 3L & padj_within < 0.05])))
  }

  # ---------------- correlation-ranked GSEA ----------------------------------
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    cat("  [6e] fgsea not installed — ORA written, GSEA skipped.\n")
    return(invisible(NULL))
  }
  md      <- pDataDT(gobj)
  mg_cols <- grep("^metagene_", colnames(md), value = TRUE)
  if (!length(mg_cols)) { cat("  [6e] no metagene_* columns — GSEA skipped.\n"); return(invisible(NULL)) }
  # loadGiotto warns that expression matrices do not share cell_IDs, so match by
  # cell_ID rather than assuming positional alignment.
  idx <- match(colnames(X), md$cell_ID)
  st  <- sparse_pearson_stats(X)
  cat("  [6e] correlation stats computed once for", nrow(X), "genes x", ncol(X), "spots\n")

  # SEQUENCING-DEPTH CONFOUND — measure it, do not assume it away.
  # A metagene score is the mean normalized expression of its member genes, so a
  # module whose genes are simply broadly expressed will track total counts per
  # spot rather than any biological program. One extra sparse_pearson call gives
  # r(gene, depth) for every gene, which yields the partial correlation
  # controlling for depth:
  #   r_xy.z = (r_xy - r_xz r_yz) / sqrt((1 - r_xz^2)(1 - r_yz^2))
  # GSEA is run on BOTH rankings and clearly labelled, so the depth-driven and
  # depth-independent answers can be compared instead of silently conflated.
  depth   <- as.numeric(md$total_expr[idx])
  r_depth <- sparse_pearson(X, depth, st)
  corr_out <- data.table(gene = rownames(X), r_total_expr = r_depth)
  depth_diag <- data.table(metagene = mg_cols,
                           r_with_total_expr = vapply(mg_cols, function(m)
                             stats::cor(as.numeric(md[[m]][idx]), depth), 0.0))
  cat("  [6e] metagene vs sequencing depth (|r| > 0.5 means the module largely tracks depth):\n")
  for (i in seq_len(nrow(depth_diag)))
    cat(sprintf("       %-12s r = %+.3f%s\n", depth_diag$metagene[i],
                depth_diag$r_with_total_expr[i],
                if (abs(depth_diag$r_with_total_expr[i]) > 0.5) "   <-- DEPTH-DRIVEN" else ""))
  fwrite(depth_diag, file.path(ENR_DIR, "metagene_depth_confound.csv"))

  gsea_all <- list()
  for (mc in mg_cols) {
    k <- as.integer(sub("^metagene_", "", mc))
    y <- as.numeric(md[[mc]][idx])
    if (anyNA(y)) { warning("NA in ", mc, " — skipped"); next }
    r <- sparse_pearson(X, y, st)
    corr_out[[paste0("r_", mc)]] <- r
    r_yz  <- stats::cor(y, depth)
    r_adj <- (r - r_depth * r_yz) / sqrt((1 - r_depth^2) * (1 - r_yz^2))
    r_adj[!is.finite(r_adj)] <- NA_real_
    corr_out[[paste0("rpartial_", mc)]] <- r_adj
    members <- intersect(mods[module == k, feat], names(r))
    for (mode in c("excl_members", "depth_adjusted", "incl_members")) {
      # CIRCULARITY: metagene_k is by construction the mean of its own member
      # genes, so those genes are guaranteed to top their own ranking. Dropping
      # them from the RANKED VECTOR (not from the gene sets — that would change
      # the sets per module and break NES comparability) turns this into the
      # question actually worth asking: what ELSE tracks this region?
      rr <- if (mode == "depth_adjusted") r_adj[!is.na(r_adj)] else r[!is.na(r)]
      if (mode != "incl_members") rr <- rr[!names(rr) %in% members]
      for (cn in names(collections)) {
        fg <- tryCatch(
          suppressWarnings(fgsea::fgsea(pathways = collections[[cn]], stats = rr,
                                        minSize = 10, maxSize = 500,
                                        eps = 0.0, nPermSimple = 10000, nproc = 1)),
          error = function(e) NULL)
        if (is.null(fg) || !nrow(fg)) next
        fg <- as.data.table(fg)
        fg[, leadingEdge := vapply(leadingEdge, paste, "", collapse = ";")]
        gsea_all[[length(gsea_all) + 1L]] <-
          fg[, .(module = k, collection = cn, mode = mode, term = pathway,
                 pval, padj, log2err, ES, NES, size, leadingEdge)]
      }
    }
    cat(sprintf("    module %d: r in [%.3f, %.3f], %d member genes excluded\n",
                k, min(r, na.rm = TRUE), max(r, na.rm = TRUE), length(members)))
  }
  fwrite(corr_out, file.path(ENR_DIR, "metagene_gene_correlations.csv.gz"))
  gsea_dt <- rbindlist(gsea_all, fill = TRUE)
  if (nrow(gsea_dt)) {
    for (cn in unique(gsea_dt$collection)) for (mo in unique(gsea_dt$mode))
      fwrite(gsea_dt[collection == cn & mode == mo][order(module, pval)],
             file.path(ENR_DIR, sprintf("gsea_%s_%s.csv", cn, mo)))
    cat(sprintf("  [6e] GSEA: %d results, %d significant (padj<0.05, excl_members)\n",
                nrow(gsea_dt), nrow(gsea_dt[mode == "excl_members" & padj < 0.05])))
  }

  # ---------------- at-a-glance summary --------------------------------------
  top3 <- function(dt, by_col, val_col) {
    if (!nrow(dt)) return("")
    paste(head(dt[order(get(val_col))][[by_col]], 3), collapse = " | ")
  }
  summ <- rbindlist(lapply(sort(unique(mods$module)), function(k)
    rbindlist(lapply(names(collections), function(cn) {
      o <- if (nrow(ora_all)) ora_all[module == k & collection == cn &
                                      universe == "qc_universe" & overlap >= 3L] else data.table()
      g <- if (nrow(gsea_dt)) gsea_dt[module == k & collection == cn &
                                      mode == "excl_members" & NES > 0] else data.table()
      ga <- if (nrow(gsea_dt)) gsea_dt[module == k & collection == cn &
                                       mode == "depth_adjusted" & NES > 0] else data.table()
      nk <- mods[module == k, .N]
      rdep <- depth_diag[metagene == paste0("metagene_", k), r_with_total_expr]
      data.table(module = k, n_genes = nk, collection = cn,
                 n_genes_annotated = length(intersect(mods[module == k, feat],
                                                      unique(unlist(collections[[cn]])))),
                 r_with_total_expr = if (length(rdep)) round(rdep, 3) else NA_real_,
                 ora_top3  = top3(o, "term", "pvalue"),
                 ora_min_padj  = if (nrow(o)) min(o$padj_within) else NA_real_,
                 gsea_top3 = top3(g, "term", "pval"),
                 gsea_min_padj = if (nrow(g)) min(g$padj) else NA_real_,
                 gsea_depth_adjusted_top3 = top3(ga, "term", "pval"),
                 # A module holding >50% of the module genes is an undifferentiated
                 # bucket; its terms will be generic regardless of the statistics.
                 verdict_note = paste(na.omit(c(
                   if (length(rdep) && abs(rdep) > 0.5)
                     sprintf("DEPTH-DRIVEN (r=%+.2f with total_expr): read gsea_depth_adjusted_top3, not gsea_top3", rdep),
                   if (nk > 0.5 * nrow(mods))
                     "LOW INFORMATION: module holds >50% of all module genes",
                   if (nk < 15L)
                     "UNDERPOWERED for ORA (<15 genes); rely on GSEA")), collapse = "; "))
    }))))
  fwrite(summ, file.path(ENR_DIR, "metagene_module_summary.csv"))

  # ---------------- figures ---------------------------------------------------
  # The GSEA dotplots are NOT drawn here. Section 6e writes the gsea_*.csv files
  # above and stops; scripts/plot_enrichment.R reads them and draws the figures,
  # so the plotting code exists in exactly one place and cannot drift between the
  # four per-sample scripts. Run it after all samples finish; to
  # redraw by hand:  Rscript scripts/plot_enrichment.R

  # ---------------- k-sweep (advisory only) -----------------------------------
  # MODULE_K is chosen per sample in Section 1 from the Section 6b diagnostics.
  # clusterSpatialCorFeats() is only a cutree of a cached hclust, so sweeping k
  # is nearly free. This REPORTS; it never mutates MODULE_K, because doing so
  # would silently invalidate every existing metagene figure and CSV.
  tryCatch({
    sc_path <- file.path(CKPT_DIR, "spat_cor.rds")
    # Match go_bp / go.bp / gobp / c5.go.bp — fetch_genesets.R names the file
    # "go_bp" with an UNDERSCORE, which an over-specific dotted pattern misses.
    i_bp <- grep("go[._]?bp", names(collections), ignore.case = TRUE)[1]
    if (!file.exists(sc_path)) {
      cat("  [6e] k-sweep skipped: no", sc_path, "\n"); return(invisible(NULL))
    }
    if (is.na(i_bp)) {
      cat("  [6e] k-sweep skipped: no GO:BP collection among",
          paste(names(collections), collapse = ", "), "\n"); return(invisible(NULL))
    }
    bp  <- collections[[i_bp]]
    sc0 <- readRDS(sc_path)
    sweep <- rbindlist(lapply(2:8, function(kk) {
      sck <- clusterSpatialCorFeats(sc0, k = kk)
      fc  <- sck[["cor_clusters"]][["spat_clus"]]
      dtk <- data.table(feat = names(fc), module = as.integer(fc))
      rbindlist(lapply(sort(unique(dtk$module)), function(m) {
        gs <- intersect(dtk[module == m, feat], rownames(X))
        if (length(gs) < 3L) return(NULL)
        yy <- as.numeric(Matrix::colMeans(X[gs, , drop = FALSE]))
        rr <- sparse_pearson(X, yy, st); rr <- rr[!is.na(rr)]
        rr <- rr[!names(rr) %in% gs]
        fg <- tryCatch(suppressWarnings(fgsea::fgsea(bp, rr, minSize = 10, maxSize = 500,
                                                     eps = 0.0, nPermSimple = 1000, nproc = 1)),
                       error = function(e) NULL)
        if (is.null(fg) || !nrow(fg)) return(NULL)
        fg <- as.data.table(fg)
        data.table(k = kk, module = m, n_genes = length(gs),
                   n_sig = sum(fg$padj < 0.05, na.rm = TRUE),
                   min_padj = min(fg$padj, na.rm = TRUE),
                   top_terms = paste(head(fg[order(pval)]$pathway, 20), collapse = ";"))
      }))
    }))
    if (nrow(sweep)) {
      fwrite(sweep, file.path(OUT_DIR, "metagene_k_sweep_enrichment.csv"))
      agg <- sweep[, .(n_modules = .N, n_modules_with_sig = sum(n_sig > 0),
                       total_sig = sum(n_sig)), by = k][order(k)]
      cat("  [6e] k-sweep (GO:BP, advisory only — MODULE_K is NOT changed):\n")
      print(agg)
      cat(sprintf("  [6e] RECOMMENDED_K (most modules with >=1 significant term): %d",
                  agg[which.max(n_modules_with_sig), k]))
      cat(sprintf("  |  MODULE_K in use: %d\n", MODULE_K))
    }
  }, error = function(e) warning("k-sweep: ", conditionMessage(e)))

  rm(X); gc(full = TRUE)
  mem_mark("06e metagene functional enrichment")
  invisible(TRUE)
}, error = function(e) {
  cat("  [6e] !! enrichment FAILED: ", conditionMessage(e), "\n", sep = "")
  warning("Metagene enrichment skipped: ", conditionMessage(e))
  invisible(NULL)
})


# =============================================================================
# SECTIONS 2-8c — skipped wholesale when RESUME_S8=1 (see Section 1). The
# brace opened here closes just after Section 8c, where the matching else-
# branch reloads the object from the after_s8 checkpoint and regenerates
# Sections 6e/8/8b/8c from it. A top-level `{` creates no new scope in R, so
# everything below still assigns into the global environment exactly as before.
# =============================================================================
if (!RESUME_S8) {

# =============================================================================
# SECTION 2 — LOAD VISIUM DATA INTO A GIOTTO OBJECT
# =============================================================================
gobj <- createGiottoVisiumObject(
  visium_dir                = VISIUM_DIR,
  expr_data                 = "filter",
  h5_visium_path            = file.path(VISIUM_DIR, "filtered_feature_bc_matrix.h5"),
  h5_tissue_positions_path  = file.path(VISIUM_DIR, "spatial", "tissue_positions.csv"),
  h5_image_png_path         = file.path(VISIUM_DIR, "spatial", "tissue_lowres_image.png"),
  h5_json_scalefactors_path = file.path(VISIUM_DIR, "spatial", "scalefactors_json.json"),
  png_name                  = "tissue_lowres_image.png",
  gene_column_index         = 2,
  instructions              = giotto_instrs
)

gobj <- addCellMetadata(gobj,
  new_metadata   = data.table(cell_ID  = pDataDT(gobj)$cell_ID,
                               sample   = SAMPLE_NAME,
                               outcome  = SAMPLE_OUTCOME),
  by_column      = TRUE,
  column_cell_ID = "cell_ID")

on_tissue  <- pDataDT(gobj)$in_tissue == 1
gobj       <- subsetGiotto(gobj, cell_ids = pDataDT(gobj)$cell_ID[on_tissue])
gene_names <- fDataDT(gobj)$feat_ID
gobj       <- subsetGiotto(gobj, feat_ids = gene_names[!grepl("^MT-|^RPS|^RPL", gene_names)])

# =============================================================================
# SECTION 3 — QUALITY CONTROL, FILTERING, NORMALIZATION
# =============================================================================
gobj <- addStatistics(gobj, expression_values = "raw")

gobj <- filterGiotto(
  gobj,
  expression_threshold   = 5,
  feat_det_in_min_cells  = 10,
  min_det_feats_per_cell = 200,
  expression_values      = "raw",
  verbose                = TRUE
)

gobj <- normalizeGiotto(gobj, scalefactor = 10000, verbose = TRUE)
gobj <- addStatistics(gobj, expression_values = "normalized")

options(bitmapType = "cairo")

# QC spatial plots — both formats
for (metric in c("nr_feats", "total_expr")) {
  for (fmt in c("png", "svg")) {
    spatPlot2D(gobj, cell_color = metric, color_as_factor = FALSE,
               point_size = 2.5, save_plot = TRUE,
               save_param = sparam(paste0("qc_", metric, "_spatial"), fmt))
  }
}

# MYB expression — both formats
for (fmt in c("png", "svg")) {
  tryCatch(
    spatFeatPlot2D(gobj, feats = "MYB", expression_values = "normalized",
                   point_size = 2.5, coord_fix_ratio = 1, gradient_midpoint = 3,
                   save_param = sparam("MYB_expr_norm", fmt)),
    error = function(e) warning("MYB spatFeatPlot2D failed: ", conditionMessage(e))
  )
}

# =============================================================================
# SECTION 4 — HIGHLY VARIABLE GENES, PCA, UMAP, SPATIAL NETWORK
# =============================================================================
gobj <- calculateHVF(
  gobject     = gobj,
  method      = "cov_loess",
  save_plot   = TRUE,
  show_plot   = FALSE,
  return_plot = FALSE,
  save_param  = list(save_name = "hvf_plot", base_width = 8,
                     base_height = 6, units = "in",
                     dpi = 600, save_format = "png"))

# SVG version (second call; HVF result identical, just saves SVG)
calculateHVF(
  gobject     = gobj,
  method      = "cov_loess",
  save_plot   = TRUE,
  show_plot   = FALSE,
  return_plot = FALSE,
  save_param  = list(save_name = "hvf_plot", base_width = 8,
                     base_height = 6, units = "in",
                     save_format = "svg"))

gobj <- runPCA(gobj, feats_to_use = "hvf", expression_values = "normalized")

for (fmt in c("png", "svg")) {
  screePlot(gobj, ncp = 30,
            save_param = sparam("screeplot", fmt))
}

gobj <- runUMAP(gobj, dimensions_to_use = 1:15)
gobj <- createNearestNetwork(gobj, dimensions_to_use = 1:15, k = 15)

gobj <- createSpatialNetwork(gobj, method = "Delaunay", minimum_k = 2,
                              name = "Delaunay_network")

for (fmt in c("png", "svg")) {
  spatPlot2D(gobj, show_network = TRUE, network_color = "red",
             spatial_network_name = "Delaunay_network",
             point_size = 2.5,
             save_param = sparam("spatial_network", fmt))
}

# =============================================================================
# SECTION 5 — LEIDEN CLUSTERING + HMRF SPATIAL DOMAINS
# =============================================================================
gobj <- doLeidenCluster(gobj, resolution = 0.5, n_iterations = 1000,
                        name = "leiden_clus")

for (fmt in c("png", "svg")) {
  spatPlot2D(gobj, cell_color = "leiden_clus", point_size = 2.5,
             save_plot = TRUE,
             save_param = sparam("leiden_spatial", fmt))
  plotUMAP(gobj, cell_color = "leiden_clus",
           save_param = sparam("leiden_umap", fmt))
}

# -- 5 QC: spots per Leiden cluster -------------------------------------------
leiden_tab <- sort(table(pDataDT(gobj)$leiden_clus), decreasing = TRUE)
for (fmt in c("png", "svg")) {
  out_file <- file.path(OUT_DIR, paste0("leiden_spots_per_cluster.", fmt))
  p <- ggplot(data.frame(cluster = names(leiden_tab), n = as.integer(leiden_tab)),
              aes(x = reorder(cluster, -n), y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Leiden cluster", y = "Number of spots",
         title = "Spots per Leiden cluster") +
    theme_classic(base_size = 12)
  ggsave(out_file, p, width = 8, height = 5,
         dpi = if (fmt == "png") 600 else 96)
}

# -- 5b. binSpect for HMRF spatial gene input (COMMENTED OUT) -----------------
# spatial_genes_init <- binSpect(gobj, bin_method = "rank",
#                                 spatial_network_name = "Delaunay_network",
#                                 expression_values    = "normalized")
# hmrf_input_genes <- spatial_genes_init[1:500, ]$feats

# -- 5c. HMRF over k = 3, 4, 5 (COMMENTED OUT) --------------------------------
# HMRF_K_VALS  <- 3:5
# HMRF_BETAS   <- c(0, 5, 10, 15, 20)
#
# for (k_val in HMRF_K_VALS) {
#   cat(sprintf("--- HMRF k=%d ---\n", k_val))
#
#   HMRF_out <- doHMRF(
#     gobj                 = gobj,
#     expression_values    = "normalized",
#     spatial_network_name = "Delaunay_network",
#     spatial_genes        = hmrf_input_genes,
#     k                    = k_val,
#     betas                = c(0, 5, 5),
#     output_folder        = file.path(OUT_DIR, paste0("HMRF_k", k_val))
#   )
#
#   gobj <- addHMRF(
#     gobj         = gobj,
#     HMRFoutput   = HMRF_out,
#     k            = k_val,
#     betas_to_add = HMRF_BETAS,
#     hmrf_name    = "HMRF"
#   )
#
#   hmrf_cols <- paste0("HMRF_k", k_val, "_b.", HMRF_BETAS)
#
#   for (i in seq_along(HMRF_BETAS)) {
#     for (fmt in c("png", "svg")) {
#       spatPlot2D(gobj,
#                  cell_color      = hmrf_cols[i],
#                  point_size      = 2.5,
#                  coord_fix_ratio = 1,
#                  save_plot       = TRUE,
#                  save_param      = sparam(
#                    paste0("HMRF_k", k_val, "_beta", HMRF_BETAS[i], "_spatial"),
#                    fmt))
#     }
#   }
#
#   b10_col <- paste0("HMRF_k", k_val, "_b.10")
#   if (b10_col %in% colnames(pDataDT(gobj))) {
#     for (fmt in c("png", "svg")) {
#       plotUMAP(gobj, cell_color = b10_col,
#                save_param = sparam(paste0("HMRF_k", k_val, "_b10_umap"), fmt))
#     }
#   }
#
#   valid_cols <- intersect(hmrf_cols, colnames(pDataDT(gobj)))
#   hmrf_meta  <- pDataDT(gobj)[, c("cell_ID", valid_cols), with = FALSE]
#   fwrite(hmrf_meta,
#          file.path(OUT_DIR, paste0("HMRF_k", k_val, "_domain_assignments.csv")))
# }

# -- 5d. HMRF cluster validation metrics (COMMENTED OUT) ----------------------
# cat("--- Section 5d: Cluster validation metrics (WSS, Silhouette, Gap) ---\n")
# pca_mat <- as.matrix(getDimReduction(gobj,
#                                       reduction_method = "pca",
#                                       name             = "pca",
#                                       output           = "matrix"))[, 1:15, drop = FALSE]
# pca_dist <- dist(pca_mat, method = "euclidean")
# expr_full  <- getExpression(gobj, values = "normalized", output = "matrix")
# expr_svg_t <- t(as.matrix(expr_full[intersect(hmrf_input_genes, rownames(expr_full)), ]))
# wss_vals     <- numeric(length(HMRF_K_VALS))
# sil_avg_vals <- numeric(length(HMRF_K_VALS))
# for (idx in seq_along(HMRF_K_VALS)) {
#   k_val   <- HMRF_K_VALS[idx]
#   b10_col <- paste0("HMRF_k", k_val, "_b.10")
#   meta_dt <- pDataDT(gobj)
#   if (!b10_col %in% colnames(meta_dt)) {
#     warning(sprintf("Column %s not found — skipping k=%d", b10_col, k_val))
#     wss_vals[idx] <- NA; sil_avg_vals[idx] <- NA; next
#   }
#   clus_vec        <- meta_dt[[b10_col]]
#   names(clus_vec) <- meta_dt$cell_ID
#   clus_int        <- as.integer(clus_vec[rownames(pca_mat)])
#   expr_aligned    <- expr_svg_t[rownames(pca_mat), , drop = FALSE]
#   wss_k <- 0
#   for (cl in unique(clus_int)) {
#     spot_idx <- which(clus_int == cl)
#     if (length(spot_idx) < 2) next
#     cl_mat   <- expr_aligned[spot_idx, , drop = FALSE]
#     centroid <- colMeans(cl_mat)
#     diffs    <- sweep(cl_mat, 2, centroid, "-")
#     wss_k    <- wss_k + sum(diffs^2)
#   }
#   wss_vals[idx] <- wss_k
#   sil_obj           <- cluster::silhouette(clus_int, pca_dist)
#   sil_avg_vals[idx] <- mean(sil_obj[, "sil_width"])
#   cat(sprintf("  k=%2d  WSS=%.2f  Sil=%.4f\n", k_val, wss_k, sil_avg_vals[idx]))
# }
# set.seed(42)
# gap_result <- cluster::clusGap(pca_mat, FUN = kmeans, nstart = 25,
#                                 K.max = max(HMRF_K_VALS), B = 50, verbose = FALSE)
# gap_tab     <- as.data.frame(gap_result$Tab)
# gap_vals    <- gap_tab$gap[HMRF_K_VALS]
# gap_se_vals <- gap_tab$SE.sim[HMRF_K_VALS]
# validation_df <- data.frame(k = HMRF_K_VALS, wss = wss_vals,
#                              silhouette_avg = sil_avg_vals,
#                              gap_stat = gap_vals, gap_se = gap_se_vals)
# fwrite(as.data.table(validation_df),
#        file.path(OUT_DIR, "HMRF_cluster_validation_metrics.csv"))
# for (fmt in c("png", "svg")) {
#   out_file <- file.path(OUT_DIR, paste0("HMRF_validation_wss.", fmt))
#   p_wss <- ggplot(validation_df, aes(x = k, y = wss)) +
#     geom_line(color = "steelblue", linewidth = 0.8) +
#     geom_point(color = "steelblue", size = 2.5) +
#     scale_x_continuous(breaks = HMRF_K_VALS) +
#     labs(x = "Number of domains (k)", y = "WSS",
#          title = "HMRF elbow plot — WSS vs k (beta=10)") +
#     theme_classic(base_size = 12)
#   ggsave(out_file, p_wss, width = 8, height = 5,
#          dpi = if (fmt == "png") 600 else 96)
# }
# best_sil_k <- HMRF_K_VALS[which.max(sil_avg_vals)]
# for (fmt in c("png", "svg")) {
#   out_file <- file.path(OUT_DIR, paste0("HMRF_validation_silhouette.", fmt))
#   p_sil <- ggplot(validation_df, aes(x = k, y = silhouette_avg)) +
#     geom_line(color = "darkorange", linewidth = 0.8) +
#     geom_point(color = "darkorange", size = 2.5) +
#     geom_vline(xintercept = best_sil_k, linetype = "dashed",
#                color = "red", linewidth = 0.7) +
#     annotate("text", x = best_sil_k + 0.4, y = max(sil_avg_vals, na.rm = TRUE),
#              label = paste0("best k=", best_sil_k), hjust = 0, size = 3.5) +
#     scale_x_continuous(breaks = HMRF_K_VALS) +
#     labs(x = "Number of domains (k)", y = "Mean silhouette width",
#          title = "HMRF silhouette score vs k (beta=10)") +
#     theme_classic(base_size = 12)
#   ggsave(out_file, p_sil, width = 8, height = 5,
#          dpi = if (fmt == "png") 600 else 96)
# }
# gap_png_path <- file.path(OUT_DIR, "HMRF_validation_gap.png")
# png(gap_png_path, width = 8, height = 5, units = "in", res = 600)
# plot(gap_result, main = "HMRF gap statistic vs k (kmeans proxy on PCA)",
#      xlab = "Number of clusters (k)", ylab = "Gap statistic")
# dev.off()
# gap_svg_path <- file.path(OUT_DIR, "HMRF_validation_gap.svg")
# svg(gap_svg_path, width = 8, height = 5)
# plot(gap_result, main = "HMRF gap statistic vs k (kmeans proxy on PCA)",
#      xlab = "Number of clusters (k)", ylab = "Gap statistic")
# dev.off()
# cat("--- Section 5d complete.\n")

# =============================================================================
# SECTION 6 — SPATIALLY VARIABLE GENES + SPATIAL CO-EXPRESSION MODULES
# =============================================================================
SVG_rank   <- binSpect(gobj, bin_method = "rank",
                        spatial_network_name = "Delaunay_network",
                        expression_values    = "normalized")
SVG_kmeans <- binSpect(gobj, bin_method = "kmeans",
                        spatial_network_name = "Delaunay_network",
                        expression_values    = "normalized")

fwrite(SVG_rank,   file.path(OUT_DIR, "SVG_rank.csv"))
fwrite(SVG_kmeans, file.path(OUT_DIR, "SVG_kmeans.csv"))

# -- 6 QC: top-5 SVGs spatial feature plots -----------------------------------
top5_svgs <- SVG_rank$feats[1:5]
for (g in top5_svgs) {
  for (fmt in c("png", "svg")) {
    .p_args <- list(gobj, feats = g, expression_values = "normalized",
                    point_size = 2.5, coord_fix_ratio = 1,
                    save_param = sparam(paste0("top_SVG_", g), fmt))
    if (!is.na(SVG_PLOT_MIDPOINT)) .p_args$gradient_midpoint <- SVG_PLOT_MIDPOINT
    do.call(spatFeatPlot2D, .p_args)
  }
}

# -- 6 QC: SVG rank vs adjusted p-value scatter -------------------------------
svg_top200 <- SVG_rank[1:200, ]
for (fmt in c("png", "svg")) {
  out_file <- file.path(OUT_DIR, paste0("SVG_rank_pvalue_scatter.", fmt))
  p <- ggplot(svg_top200, aes(x = seq_len(nrow(svg_top200)),
                               y = -log10(adj.p.value))) +
    geom_point(size = 1.2, color = "steelblue") +
    labs(x = "SVG rank", y = "-log10(adj. p-value)",
         title = "Top 200 spatially variable genes") +
    theme_classic(base_size = 12)
  ggsave(out_file, p, width = 7, height = 5,
         dpi = if (fmt == "png") 600 else 96)
}

# Top 250 SVGs for co-expression module detection
top_svg <- SVG_rank[1:250, ]$feats

spat_cor <- detectSpatialCorFeats(
  gobj,
  method               = "network",
  spatial_network_name = "Delaunay_network",
  subset_feats         = top_svg
)

spat_cor <- clusterSpatialCorFeats(spat_cor, k = MODULE_K)
# Cache for the Section 6e k-sweep, which re-cuts this tree at several k.
saveRDS(spat_cor, file.path(CKPT_DIR, "spat_cor.rds"))

# =============================================================================
# SECTION 6b — METAGENE MODULE k-SELECTION VALIDATION
# Justifies MODULE_K (set in Section 1) for clusterSpatialCorFeats using three complementary metrics
# applied to the gene-gene spatial correlation matrix (1 - |r| distance):
#   1. WSS elbow plot (k=2:12)
#   2. Silhouette width sweep (k=2:12)
#   3. Annotated dendrogram with MODULE_K clusters highlighted via rect.hclust
#   4. Cophenetic correlation coefficient (tree fit quality)
# =============================================================================
cat("--- Section 6b: Metagene module k-selection validation ---\n")

# Reconstruct the correlation matrix exactly as clusterSpatialCorFeats does:
# spat_cor$cor_DT is a long data.table (feat_ID x variable x spat_cor)
cor_DT    <- spat_cor[["cor_DT"]]
cor_DT_dc <- data.table::dcast.data.table(cor_DT,
               formula = feat_ID ~ variable, value.var = "spat_cor")
cor_mat   <- as.matrix(cor_DT_dc[, -1L, with = FALSE])
rownames(cor_mat) <- cor_DT_dc$feat_ID
feat_ord  <- spat_cor[["feat_order"]]
cor_mat   <- cor_mat[feat_ord, feat_ord]
cat(sprintf("  Correlation matrix: %d x %d genes\n", nrow(cor_mat), ncol(cor_mat)))

# Reuse the hclust already computed by clusterSpatialCorFeats (ward.D, 1-r distance)
hc        <- spat_cor[["cor_hclust"]][["spat_clus"]]
gene_dist <- as.dist(1 - cor_mat)   # same distance used internally

# Cophenetic correlation — how faithfully the tree preserves pairwise distances
coph_r <- cor(gene_dist, cophenetic(hc))
cat(sprintf("  Cophenetic correlation coefficient: %.4f\n", coph_r))

# WSS and silhouette for k = 2:12
METAGENE_K_RANGE <- 2:12
wss_mg           <- numeric(length(METAGENE_K_RANGE))
sil_mg           <- numeric(length(METAGENE_K_RANGE))
dist_mat_mg      <- as.matrix(gene_dist)

for (i in seq_along(METAGENE_K_RANGE)) {
  ki     <- METAGENE_K_RANGE[i]
  labels <- cutree(hc, k = ki)

  # WSS in dissimilarity space: sum of squared intra-cluster distances / 2n
  wss_k <- 0
  for (cl in unique(labels)) {
    idx <- which(labels == cl)
    if (length(idx) < 2) next
    wss_k <- wss_k + sum(dist_mat_mg[idx, idx]^2) / (2 * length(idx))
  }
  wss_mg[i] <- wss_k

  sil_obj  <- cluster::silhouette(labels, gene_dist)
  sil_mg[i] <- mean(sil_obj[, "sil_width"])
  cat(sprintf("  k=%2d  WSS=%.2f  Sil=%.4f\n", ki, wss_k, sil_mg[i]))
}

# Save validation table
metagene_val_df <- data.frame(k              = METAGENE_K_RANGE,
                               wss            = wss_mg,
                               silhouette_avg = sil_mg)
fwrite(as.data.table(metagene_val_df),
       file.path(OUT_DIR, "metagene_k_validation.csv"))

best_sil_mg <- METAGENE_K_RANGE[which.max(sil_mg)]
cat(sprintf("  Best silhouette k = %d  |  chosen k = %d\n",
            best_sil_mg, MODULE_K))

# Plot 1: WSS elbow
for (fmt in c("png", "svg")) {
  out_file <- file.path(OUT_DIR, paste0("metagene_k_wss.", fmt))
  p_wss <- ggplot(metagene_val_df, aes(x = k, y = wss)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(aes(color = (k == MODULE_K)), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
    geom_vline(xintercept = MODULE_K, linetype = "dashed", color = "red",
               linewidth = 0.7) +
    annotate("text", x = MODULE_K + 0.25, y = max(wss_mg) * 0.98,
             label = sprintf("k = %d", MODULE_K),
             color = "red", hjust = 0, size = 3.5) +
    scale_x_continuous(breaks = METAGENE_K_RANGE) +
    labs(x = "Number of modules (k)",
         y = "Within-cluster sum of squares (WSS)",
         title = "Metagene co-expression module elbow plot") +
    theme_classic(base_size = 12)
  ggsave(out_file, p_wss, width = 8, height = 5,
         dpi = if (fmt == "png") 600 else 96)
}

# Plot 2: Silhouette sweep
for (fmt in c("png", "svg")) {
  out_file <- file.path(OUT_DIR, paste0("metagene_k_silhouette.", fmt))
  p_sil <- ggplot(metagene_val_df, aes(x = k, y = silhouette_avg)) +
    geom_line(color = "darkorange", linewidth = 0.8) +
    geom_point(aes(color = (k == MODULE_K)), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "darkorange", "TRUE" = "red")) +
    geom_vline(xintercept = MODULE_K, linetype = "dashed", color = "red",
               linewidth = 0.7) +
    annotate("text", x = MODULE_K + 0.25,
             y = sil_mg[METAGENE_K_RANGE == MODULE_K],
             label = sprintf("k=%d (sil=%.3f)", MODULE_K,
                             sil_mg[METAGENE_K_RANGE == MODULE_K]),
             color = "red", hjust = 0, size = 3.5) +
    scale_x_continuous(breaks = METAGENE_K_RANGE) +
    labs(x = "Number of modules (k)", y = "Mean silhouette width",
         title = "Metagene co-expression module silhouette score vs k") +
    theme_classic(base_size = 12)
  ggsave(out_file, p_sil, width = 8, height = 5,
         dpi = if (fmt == "png") 600 else 96)
}

# Plot 3: Dendrogram with MODULE_K clusters highlighted
dend_title <- sprintf(
  "Gene co-expression dendrogram  |  top %d SVGs  |  cophenetic r = %.3f  |  ward.D",
  length(top_svg), coph_r)
dend_colors <- c("firebrick", "steelblue", "darkgreen", "darkorange",
                 "purple", "brown", "darkcyan", "magenta4",
                 "olivedrab", "sienna", "slateblue", "gold3")[seq_len(MODULE_K)]

dend_png <- file.path(OUT_DIR, "metagene_k_dendrogram.png")
png(dend_png, width = 12, height = 6, units = "in", res = 300)
par(mar = c(4, 4, 3, 1))
plot(hc, labels = FALSE, hang = -1, main = dend_title,
     xlab = "Genes", ylab = "Height  (1 - |r|)", sub = "")
rect.hclust(hc, k = MODULE_K, border = dend_colors)
legend("topright", legend = paste0("Module ", seq_len(MODULE_K)),
       fill = dend_colors, bty = "n", cex = 0.85)
dev.off()

dend_svg <- file.path(OUT_DIR, "metagene_k_dendrogram.svg")
svg(dend_svg, width = 12, height = 6)
par(mar = c(4, 4, 3, 1))
plot(hc, labels = FALSE, hang = -1, main = dend_title,
     xlab = "Genes", ylab = "Height  (1 - |r|)", sub = "")
rect.hclust(hc, k = MODULE_K, border = dend_colors)
legend("topright", legend = paste0("Module ", seq_len(MODULE_K)),
       fill = dend_colors, bty = "n", cex = 0.85)
dev.off()

cat(sprintf("--- Section 6b complete. Cophenetic r=%.4f | best sil k=%d | chosen k=%d\n",
            coph_r, best_sil_mg, MODULE_K))

feat_clus_vector <- spat_cor[["cor_clusters"]][["spat_clus"]]
expr_mat         <- getExpression(gobj, values = "normalized", output = "matrix")
feat_clus_dt     <- data.table(feat   = names(feat_clus_vector),
                                module = as.integer(feat_clus_vector))
module_ids       <- sort(unique(feat_clus_dt$module))

meta_list <- lapply(module_ids, function(k) {
  genes_k <- intersect(feat_clus_dt[module == k, feat], rownames(expr_mat))
  scores  <- colMeans(as.matrix(expr_mat[genes_k, , drop = FALSE]))
  rng     <- range(scores)
  if (diff(rng) > 0) (scores - rng[1]) / diff(rng) else scores
})

meta_dt <- as.data.table(setNames(meta_list, paste0("metagene_", module_ids)))
meta_dt[, cell_ID := colnames(expr_mat)]
gobj    <- addCellMetadata(gobj, new_metadata = meta_dt,
                           by_column = TRUE, column_cell_ID = "cell_ID")
fwrite(feat_clus_dt, file.path(OUT_DIR, "metagene_module_genes.csv"))

# -- 6 QC: genes per module barplot -------------------------------------------
module_sizes <- feat_clus_dt[, .N, by = module][order(module)]
for (fmt in c("png", "svg")) {
  out_file <- file.path(OUT_DIR, paste0("module_gene_counts.", fmt))
  p <- ggplot(module_sizes, aes(x = factor(module), y = N)) +
    geom_bar(stat = "identity", fill = "coral") +
    labs(x = "Module", y = "Number of genes",
         title = "Genes per spatial co-expression module") +
    theme_classic(base_size = 12)
  ggsave(out_file, p, width = 6, height = 4,
         dpi = if (fmt == "png") 600 else 96)
}

# Metagene spatial plots — magma colour map, both formats
for (k in module_ids) {
  col_name <- paste0("metagene_", k)
  for (fmt in c("png", "svg")) {
    spatPlot2D(
      gobj,
      cell_color          = col_name,
      color_as_factor     = FALSE,
      cell_color_gradient = viridis::magma(256),
      point_size          = 2.5,
      coord_fix_ratio     = 1,
      save_plot           = TRUE,
      save_param          = sparam(paste0("metagene_module_", k), fmt)
    )
  }
}

# Spatial correlation heatmap — both formats
for (fmt in c("png", "svg")) {
  heatmSpatialCorFeats(
    gobj,
    spatCorObject = spat_cor,
    use_clus_name = "spat_clus",
    save_param    = sparam("spat_cor_heatmap", fmt)
  )
}

# =============================================================================
# SECTION 6e — METAGENE MODULE FUNCTIONAL ENRICHMENT (GSEA + ORA)
# Body is defined in the helper region; the resume path calls the same one.
# =============================================================================
run_metagene_enrichment(gobj)

# =============================================================================
# SECTION 7 — LOAD scRNA REFERENCE AND FIND MARKERS ON annot_7_c2l
# =============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("rhdf5",    quietly = TRUE)) BiocManager::install("rhdf5")
if (!requireNamespace("HDF5Array", quietly = TRUE)) BiocManager::install("HDF5Array")
library(rhdf5); library(HDF5Array)

sce <- zellkonverter::readH5AD(SC_REF, use_hdf5 = TRUE)

stopifnot(SC_ANNOT_COL %in% colnames(colData(sce)))
cell_types_present <- unique(as.character(colData(sce)[[SC_ANNOT_COL]]))
cat("Cell types in", SC_ANNOT_COL, ":", length(cell_types_present), "\n")
print(table(colData(sce)[[SC_ANNOT_COL]]))

counts_mat <- if ("counts" %in% assayNames(sce)) assay(sce, "counts") else assay(sce, 1)

sc_gobj <- createGiottoObject(expression = counts_mat, instructions = giotto_instrs)

sc_meta <- data.table(
  cell_ID   = colnames(counts_mat),
  cell_type = gsub("[^A-Za-z0-9_.]", "_",
                   as.character(colData(sce)[[SC_ANNOT_COL]]))
)
sc_gobj <- addCellMetadata(sc_gobj, new_metadata = sc_meta,
                           by_column = TRUE, column_cell_ID = "cell_ID")
sc_gobj <- normalizeGiotto(sc_gobj, scalefactor = 10000)

markers_gini <- findMarkers_one_vs_all(
  gobject        = sc_gobj,
  method         = "gini",
  expression_values = "normalized",
  cluster_column = "cell_type",
  min_feats      = 1
)
fwrite(markers_gini, file.path(OUT_DIR, "scRNA_markers_gini.csv"))

top_markers      <- markers_gini[, head(.SD, 30), by = "cluster"]
top_marker_genes <- unique(top_markers$feats)
cat("Markers selected:", length(top_marker_genes), "genes across",
    length(unique(top_markers$cluster)), "cell types\n\n")

# Cache the gini markers so the RESUME_S8 path can rebuild the DWLS
# signature for Section 8b without reloading the 8 GB reference h5ad.
saveRDS(markers_gini, file.path(CKPT_DIR, "markers_gini.rds"))

# =============================================================================
# SECTION 7c — scRNA PER-CELL-TYPE EXPRESSION PROFILE
#
# Per-cell-type mean expression, percent-expressed and a tau specificity
# score. This table was originally built for the ICF screen, which has since
# been retired (archive/cytospace_icf/) — it is kept because it is already
# generated, is file.exists()-guarded so it costs nothing on a re-run, and is
# directly useful for interpreting Section 8b validation and Section 8c
# composition: it is how you check whether a cell type that scores badly in
# 8b actually has distinguishing expression in the reference at all.
# =============================================================================
sc_profile_path <- file.path(DATA_DIR, "sc_ref", "sc_celltype_profile.tsv.gz")
if (!file.exists(sc_profile_path)) {
  cat("--- [7c] Building scRNA per-cell-type profile ---\n")
  tryCatch({
    Xsc  <- getExpression(sc_gobj, values = "normalized", output = "matrix")
    ctsc <- pDataDT(sc_gobj)$cell_type
    ind  <- Matrix::sparse.model.matrix(~ 0 + factor(ctsc))
    colnames(ind) <- levels(factor(ctsc))
    n_per <- Matrix::colSums(ind)
    mean_expr <- as.matrix((Xsc %*% ind) %*% Matrix::Diagonal(x = 1 / n_per))
    colnames(mean_expr) <- colnames(ind)
    pct_expr <- as.matrix((methods::as(Xsc > 0, "dMatrix") %*% ind) %*%
                          Matrix::Diagonal(x = 1 / n_per))
    colnames(pct_expr) <- colnames(ind)
    # tau specificity: 0 = ubiquitous, 1 = expressed in exactly one type.
    mx  <- apply(mean_expr, 1, max)
    tau <- ifelse(mx > 0,
                  rowSums(1 - mean_expr / ifelse(mx > 0, mx, 1)) /
                    (ncol(mean_expr) - 1),
                  NA_real_)
    prof <- data.table(gene = rownames(mean_expr), tau_specificity = tau)
    for (cc in colnames(mean_expr)) {
      prof[[paste0("mean_", cc)]] <- mean_expr[, cc]
      prof[[paste0("pct_",  cc)]] <- pct_expr[, cc]
    }
    fwrite(prof, sc_profile_path, sep = "\t")
    cat(sprintf("  [7c] wrote %s (%d genes x %d cell types)\n",
                sc_profile_path, nrow(prof), ncol(mean_expr)))
    rm(Xsc, ind, mean_expr, pct_expr); gc(full = TRUE)
  }, error = function(e) {
    cat("  [7c] !! profile FAILED: ", conditionMessage(e), "\n", sep = "")
    warning("Section 7c scRNA profile failed: ", conditionMessage(e))
  })
} else {
  cat("--- [7c] scRNA profile already exists, skipping ---\n")
}

# =============================================================================
# SECTION 8 — SPATIALDWLS DECONVOLUTION
# =============================================================================
cat("--- [8] SpatialDWLS deconvolution ---\n")

DWLS_matrix <- makeSignMatrixDWLSfromMatrix(
  matrix    = getExpression(sc_gobj, values = "normalized", output = "matrix"),
  cell_type = pDataDT(sc_gobj)$cell_type,
  sign_gene = top_marker_genes
)

common_genes <- intersect(rownames(DWLS_matrix), rownames(gobj))
cat("Common genes (scRNA ∩ Visium):", length(common_genes), "\n")
DWLS_matrix  <- DWLS_matrix[common_genes, ]

if (!requireNamespace("quadprog", quietly = TRUE)) install.packages("quadprog")
if (!requireNamespace("Rfast",    quietly = TRUE)) install.packages("Rfast")
library(quadprog); library(Rfast)

gobj <- runDWLSDeconv(gobject = gobj, sign_matrix = DWLS_matrix,
                      n_cell = 20, cluster_column = "leiden_clus")

dwls_result    <- getSpatialEnrichment(gobj, name = "DWLS", output = "data.table")
fwrite(dwls_result, file.path(OUT_DIR, "DWLS_proportions.csv"))

cell_type_cols <- setdiff(colnames(dwls_result), "cell_ID")

# dominant_celltype is an argmax over the DWLS proportions. At n_cell = 20
# the rows sum to 1 and the 5% gate leaves 1-5 non-zero types per spot, so
# this is always defined — unlike the HD pipeline, where the 50% gate leaves
# ~10% of bins with no majority at all. Measured here: 0 unassigned spots.
.dwls_mat <- as.matrix(dwls_result[, ..cell_type_cols])
.dwls_rs  <- rowSums(.dwls_mat, na.rm = TRUE)
dominant_ct <- cell_type_cols[max.col(.dwls_mat, ties.method = "first")]
# A spot whose DWLS row is ENTIRELY ZERO has no solution at all, and max.col()
# (like which.max()) silently returns column 1 for it — fabricating a confident
# call for the alphabetically-first cell type. Measured: 121 of 3,341 spots in
# P09 (3.6%) hit this and were all being labelled
# Epithelial_cells___Tumor; the other three samples have none. Those spots are
# genuine abstentions and are marked NA, which is what CT_NA_LABEL exists for.
dominant_ct[!is.finite(.dwls_rs) | .dwls_rs <= 0] <- NA_character_
rm(.dwls_mat, .dwls_rs)
gobj <- addCellMetadata(gobj,
  new_metadata   = data.table(cell_ID                = dwls_result$cell_ID,
                              dominant_celltype      = dominant_ct,
                              dominant_celltype_plot = ct_plot_factor(dominant_ct)),
  by_column      = TRUE,
  column_cell_ID = "cell_ID")

plot_dwls_figures(gobj, dwls_result, cell_type_cols, dominant_ct)
mem_mark("08 SpatialDWLS deconvolution")

# =============================================================================
# SECTION 8b — CELL-TYPE CALL VALIDATION
# =============================================================================
CT_VALIDATION_TIER <- validate_celltype_calls(gobj, common_genes)

# =============================================================================
# SECTION 8c — CELL-TYPE COMPOSITION
# =============================================================================
summarize_composition(gobj, dwls_result, cell_type_cols, dominant_ct,
                      CT_VALIDATION_TIER)

# Checkpoint AFTER 8c so RESUME_S8=1 can regenerate 6e/8/8b/8c from here.
checkpoint_gobject(gobj, "after_s8")

} else {
  # ---- RESUME_S8=1 resume path (see Section 1) -----------------------------
  # Regenerates Sections 6e/8/8b/8c from the checkpoint by calling the SAME
  # functions the full run calls, not copies, so an output can never drift
  # between the full run and the fast loop.
  gobj <- load_checkpoint("after_s8")
  cat(sprintf("  [CKPT] resumed: %d spots x %d feats | CT_COLUMN=%s\n",
              nrow(pDataDT(gobj)), nrow(fDataDT(gobj)), CT_COLUMN))
  stopifnot(CT_COLUMN %in% colnames(pDataDT(gobj)))
  mem_mark("08 resumed from checkpoint (Sections 2-8b skipped)")

  # Sections 8/8b/8c all need the DWLS results. Read them back from disk
  # rather than re-solving one quadratic program per spot.
  dwls_csv <- file.path(OUT_DIR, "DWLS_proportions.csv")
  if (file.exists(dwls_csv)) {
    dwls_result    <- fread(dwls_csv)
    cell_type_cols <- setdiff(colnames(dwls_result), "cell_ID")
    md_po          <- pDataDT(gobj)
    dominant_ct    <- md_po$dominant_celltype[match(dwls_result$cell_ID, md_po$cell_ID)]
    plot_dwls_figures(gobj, dwls_result, cell_type_cols, dominant_ct)
  } else {
    warning("RESUME_S8: ", dwls_csv, " not found; Sections 8/8b/8c skipped.")
  }

  # Section 8b needs the DWLS signature gene list, which Section 8 builds.
  # Rebuild it from the cached gini markers rather than reloading the 8 GB h5ad.
  gini_cache <- file.path(CKPT_DIR, "markers_gini.rds")
  if (file.exists(gini_cache) && exists("dwls_result")) {
    mg <- readRDS(gini_cache)
    common_genes <- intersect(unique(mg[, head(.SD, 30), by = "cluster"]$feats),
                              rownames(gobj))
    CT_VALIDATION_TIER <- validate_celltype_calls(gobj, common_genes)
    rm(mg)
  } else {
    warning("RESUME_S8: ", gini_cache, " not found; Section 8b validation skipped.")
    CT_VALIDATION_TIER <- NULL
  }

  if (exists("dwls_result")) {
    summarize_composition(gobj, dwls_result, cell_type_cols, dominant_ct,
                          CT_VALIDATION_TIER)
  }

  run_metagene_enrichment(gobj)
}

# =============================================================================
# SECTION 9 — CELL PROXIMITY ENRICHMENT
# =============================================================================
cat("--- [9] Cell proximity enrichment ---\n")

na_prox_spots <- pDataDT(gobj)$cell_ID[is.na(pDataDT(gobj)[[CT_COLUMN]])]
if (length(na_prox_spots) > 0) {
  cat(sprintf("  Dropping %d spot(s) with NA in %s before proximity analysis.\n",
              length(na_prox_spots), CT_COLUMN))
  gobj_prox <- subsetGiotto(gobj, cell_ids = setdiff(pDataDT(gobj)$cell_ID, na_prox_spots))
} else {
  gobj_prox <- gobj
}

# N_PROX_SIMS is a Section 1 constant, not a literal, because this section runs
# on BOTH execution paths and the resume path must not silently use a different
# simulation count from the full run. Cached: the result is deterministic given
# the object, and re-running 1000 permutations on every RESUME_S8 loop is pure
# waste. Delete checkpoints/prox_sim.rds to force a recompute.
cell_proximities <- cached("prox_sim", cellProximityEnrichment(
  gobject               = gobj_prox,
  cluster_column        = CT_COLUMN,
  spatial_network_name  = "Delaunay_network",
  adjust_method         = "fdr",
  number_of_simulations = N_PROX_SIMS
))
fwrite(cell_proximities$enrichm_res,
       file.path(OUT_DIR, "cell_proximity_enrichment.csv"))

# Save the observed interaction counts table (backs barplot & heatmap)
if (!is.null(cell_proximities$raw_sim_table)) {
  fwrite(as.data.table(cell_proximities$raw_sim_table),
         file.path(OUT_DIR, "cell_proximity_raw_sim_table.csv"))
}
# Save all slots as a single merged CSV for completeness
cp_slots <- names(cell_proximities)
for (sl in cp_slots) {
  obj <- cell_proximities[[sl]]
  if (is.data.frame(obj) || is.data.table(obj)) {
    fwrite(as.data.table(obj),
           file.path(OUT_DIR, paste0("cell_proximity_", sl, ".csv")))
  }
}

for (fmt in c("png", "svg")) {
  tryCatch(
    cellProximityBarplot(gobj_prox, CPscore = cell_proximities, min_sim_ints = 5,
                         save_param = sparam("cell_proximity_barplot", fmt)),
    error = function(e) warning("cellProximityBarplot failed: ", conditionMessage(e))
  )
  tryCatch(
    cellProximityHeatmap(gobj_prox, CPscore = cell_proximities, order_cell_types = TRUE,
                         save_param = sparam("cell_proximity_heatmap", fmt)),
    error = function(e) warning("cellProximityHeatmap failed: ", conditionMessage(e))
  )
}

# =============================================================================
# SECTION 11 — SAVE THE GIOTTO OBJECT
# =============================================================================
cat("--- [11] Saving Giotto object ---\n")
saveGiotto(gobject = gobj, foldername = "giotto_object",
           dir = OUT_DIR, overwrite = TRUE)

# =============================================================================
# SECTION 12 — TISSUE OVERLAY PLOTS (hires greyscale background)
# =============================================================================
cat("--- [12] Tissue overlay plots ---\n")

TP_DIR   <- file.path(OUT_DIR, "tissue_plots")
PT_SIZE  <- 1.75
PT_ALPHA <- 0.9

# Alpha of the greyscale tissue image against the white page: 1 = as before,
# lower = fainter. spatPlot2D has no image-alpha argument (it exposes only
# point_alpha and background_color), so this is baked into the luminance below.
# Compositing a value c over white at alpha a is exactly 1 - a * (1 - c), so this
# constant IS the alpha. Only the background changes; point size/alpha are untouched.
TISSUE_ALPHA <- 0.75
dir.create(TP_DIR, recursive = TRUE, showWarnings = FALSE)

tparam <- function(name, fmt, ...) {
  base <- list(save_name = name, save_format = fmt, save_dir = TP_DIR, ...)
  if (fmt == "png") base$dpi <- 600
  base
}

# Build greyscale hires image and attach
scalef_json_tp  <- jsonlite::fromJSON(
  file.path(VISIUM_DIR, "spatial", "scalefactors_json.json"))
hires_src       <- file.path(VISIUM_DIR, "spatial", "tissue_hires_image.png")
gray_hires_path <- file.path(TP_DIR, "tissue_hires_gray.png")

cat("  Converting hires image to greyscale...\n")
img_px_tp <- png::readPNG(hires_src)
lum <- if (length(dim(img_px_tp)) == 3 && dim(img_px_tp)[3] >= 3) {
  0.2989 * img_px_tp[,,1] + 0.5870 * img_px_tp[,,2] + 0.1140 * img_px_tp[,,3]
} else {
  img_px_tp
}
# Composite the tissue over white at TISSUE_ALPHA (see the constant above).
lum <- 1 - TISSUE_ALPHA * (1 - lum)
gray_rgb      <- array(0, dim = c(nrow(lum), ncol(lum), 3))
gray_rgb[,,1] <- lum
gray_rgb[,,2] <- lum
gray_rgb[,,3] <- lum
png::writePNG(gray_rgb, gray_hires_path)

sf       <- 1 / scalef_json_tp$tissue_hires_scalef
gray_img <- createGiottoLargeImage(
  raster_object = gray_hires_path,
  scale_factor  = sf,
  negative_y    = TRUE,
  name          = "tissue_hires_gray")
gobj <- addGiottoImage(gobj, images = list(gray_img))
bg   <- "tissue_hires_gray"
cat("  Greyscale hires attached. Scale factor:", round(sf, 3),
    "| tissue alpha:", TISSUE_ALPHA, "\n")

# 12a — HMRF domains at beta=10 (COMMENTED OUT — HMRF section disabled)
# cat("  12a: HMRF domains (beta=10, all k)...\n")
# for (k_val in HMRF_K_VALS) {
#   col_b10 <- paste0("HMRF_k", k_val, "_b.10")
#   if (!col_b10 %in% colnames(pDataDT(gobj))) next
#   for (fmt in c("png", "svg")) {
#     tryCatch(
#       spatPlot2D(gobj, show_image = TRUE, image_name = bg,
#                  cell_color = col_b10, color_as_factor = TRUE,
#                  point_size = PT_SIZE, point_alpha = PT_ALPHA,
#                  coord_fix_ratio = 1, save_plot = TRUE,
#                  save_param = tparam(paste0("HMRF_k", k_val, "_b10_on_tissue"), fmt,
#                                      base_width = 8, base_height = 7)),
#       error = function(e)
#         warning(sprintf("spatPlot2D k=%d fmt=%s: %s", k_val, fmt, conditionMessage(e)))
#     )
#   }
#   cat(sprintf("    k=%d done\n", k_val))
# }

# 12b — Beta sweep on greyscale tissue (COMMENTED OUT — HMRF section disabled)
# cat("  12b: Beta sweep (k=3,4,5)...\n")
# for (k_val in c(3, 4, 5)) {
#   for (beta in HMRF_BETAS) {
#     col_name <- paste0("HMRF_k", k_val, "_b.", beta)
#     if (!col_name %in% colnames(pDataDT(gobj))) next
#     for (fmt in c("png", "svg")) {
#       tryCatch(
#         spatPlot2D(gobj, show_image = TRUE, image_name = bg,
#                    cell_color = col_name, color_as_factor = TRUE,
#                    point_size = PT_SIZE, point_alpha = PT_ALPHA,
#                    coord_fix_ratio = 1, save_plot = TRUE,
#                    save_param = tparam(paste0("HMRF_k", k_val, "_b", beta, "_on_tissue"), fmt,
#                                        base_width = 8, base_height = 7)),
#         error = function(e)
#           warning(sprintf("Beta sweep k=%d b=%d fmt=%s: %s", k_val, beta, fmt, conditionMessage(e)))
#       )
#     }
#   }
# }

# 12c — Multi-panel grid (COMMENTED OUT — HMRF section disabled)
# cat("  12c: Multi-panel grid...\n")
# panels <- list()
# for (k_val in HMRF_K_VALS) {
#   col_b10 <- paste0("HMRF_k", k_val, "_b.10")
#   if (!col_b10 %in% colnames(pDataDT(gobj))) next
#   p <- tryCatch(
#     spatPlot2D(gobj, show_image = TRUE, image_name = bg,
#                cell_color = col_b10, color_as_factor = TRUE,
#                point_size = PT_SIZE, point_alpha = PT_ALPHA,
#                coord_fix_ratio = 1, save_plot = FALSE,
#                return_plot = TRUE, show_plot = FALSE,
#                title = paste0("k = ", k_val)),
#     error = function(e) { warning(conditionMessage(e)); NULL }
#   )
#   if (!is.null(p)) panels[[length(panels) + 1]] <- p
# }
# if (length(panels) > 0) {
#   combined <- patchwork::wrap_plots(panels, ncol = 3) +
#     patchwork::plot_annotation(
#       title    = paste0(SAMPLE_NAME, " — HMRF spatial domains (beta=10)"),
#       subtitle = "Each panel shows a different number of domains (k)",
#       theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
#                        plot.subtitle = element_text(size = 10)))
#   ggsave(file.path(TP_DIR, "HMRF_multipanel_tissue_HE.png"),
#          combined, width = 18, height = 12, dpi = 300)
#   ggsave(file.path(TP_DIR, "HMRF_multipanel_tissue_HE.svg"),
#          combined, width = 18, height = 12)
#   cat("  Multi-panel saved.\n")
# }

# 12d — Leiden clusters on tissue
cat("  12d: Leiden clusters...\n")
if ("leiden_clus" %in% colnames(pDataDT(gobj))) {
  for (fmt in c("png", "svg")) {
    tryCatch(
      spatPlot2D(gobj, show_image = TRUE, image_name = bg,
                 cell_color = "leiden_clus", color_as_factor = TRUE,
                 point_size = PT_SIZE, point_alpha = PT_ALPHA,
                 coord_fix_ratio = 1, save_plot = TRUE,
                 save_param = tparam("leiden_on_tissue", fmt,
                                     base_width = 8, base_height = 7)),
      error = function(e) warning("Leiden tissue plot: ", conditionMessage(e))
    )
  }
}

# 12e — Dominant cell type on tissue
cat("  12e: Dominant cell type...\n")
if ("dominant_celltype" %in% colnames(pDataDT(gobj))) {
  # CHANGED: this block used to define its own `bright_palette` and assign it by
  # POSITION over sort(unique(dominant_celltype)). That gave the same cell type a
  # different colour here than in the Section 8/8c figures, and a different
  # colour again in another sample whose set of present types differs — so two
  # tissue overlays could not be compared by eye at all. It now uses the
  # canonical CT_PALETTE via ct_plot_colors(), the same as every other
  # cell-type figure in this pipeline AND in the Visium HD pipeline.
  md_ct <- pDataDT(gobj)
  if (!"dominant_celltype_plot" %in% colnames(md_ct)) {
    gobj <- addCellMetadata(gobj,
      new_metadata   = data.table(cell_ID = md_ct$cell_ID,
                                  dominant_celltype_plot =
                                    ct_plot_factor(md_ct$dominant_celltype)),
      by_column      = TRUE,
      column_cell_ID = "cell_ID")
  }
  for (fmt in c("png", "svg")) {
    tryCatch(
      spatPlot2D(gobj, show_image = TRUE, image_name = bg,
                 cell_color = "dominant_celltype_plot", color_as_factor = TRUE,
                 cell_color_code = ct_plot_colors(pDataDT(gobj)$dominant_celltype),
                 point_size = PT_SIZE, point_alpha = PT_ALPHA,
                 coord_fix_ratio = 1, save_plot = TRUE,
                 save_param = tparam("dominant_celltype_on_tissue", fmt,
                                     base_width = 10, base_height = 8)),
      error = function(e) {
        cat("  12e !! dominant cell type plot FAILED: ", conditionMessage(e), "\n", sep = "")
        warning("Dominant cell type plot: ", conditionMessage(e))
      }
    )
  }
}

# 12f — MYB expression on tissue
cat("  12f: MYB expression...\n")
if ("MYB" %in% fDataDT(gobj)$feat_ID) {
  for (fmt in c("png", "svg")) {
    tryCatch(
      spatFeatPlot2D(gobj, show_image = TRUE, image_name = bg,
                     feats = "MYB", expression_values = "normalized",
                     point_size = PT_SIZE, point_alpha = PT_ALPHA,
                     coord_fix_ratio = 1, gradient_midpoint = 3,
                     save_plot = TRUE,
                     save_param = tparam("MYB_expr_on_tissue", fmt,
                                         base_width = 8, base_height = 7)),
      error = function(e) warning("MYB tissue plot: ", conditionMessage(e))
    )
  }
}

# 12g — Metagene modules on tissue
cat("  12g: Metagene modules...\n")
meta_cols <- grep("^metagene_", colnames(pDataDT(gobj)), value = TRUE)
for (col_name in meta_cols) {
  for (fmt in c("png", "svg")) {
    tryCatch(
      spatPlot2D(gobj, show_image = TRUE, image_name = bg,
                 cell_color = col_name, color_as_factor = FALSE,
                 cell_color_gradient = viridis::magma(256),
                 point_size = PT_SIZE, point_alpha = PT_ALPHA,
                 coord_fix_ratio = 1, save_plot = TRUE,
                 save_param = tparam(paste0(col_name, "_on_tissue"), fmt,
                                     base_width = 8, base_height = 7)),
      error = function(e) warning(col_name, " tissue plot: ", conditionMessage(e))
    )
  }
}

cat("--- Section 12 complete. Tissue plots in:", TP_DIR, "\n")

cat("\n========== DONE with", SAMPLE_NAME, "==========\n")
cat("All outputs in:", OUT_DIR, "\n")
if (RESUME_S8) {
  cat("Mode: RESUME_S8=1 — Sections 2-8b were SKIPPED and the object was loaded\n")
  cat("      from checkpoints/after_s8. Regenerated: 6e, 8 figures, 8b, 8c, 9, 11, 12.\n")
  cat("      Sections 3-6b outputs on disk are from the last FULL run.\n")
} else {
  cat("Mode: full run — Sections 2-12 all executed; checkpoints/after_s8 written.\n")
  cat("      Re-run downstream sections cheaply with:  RESUME_S8=1 Rscript <this script>\n")
}
cat("Cross-sample contrast (after all 4 samples):  Rscript scripts/composition_compare.R\n\n")

# =============================================================================
# END OF SCRIPT
# =============================================================================
