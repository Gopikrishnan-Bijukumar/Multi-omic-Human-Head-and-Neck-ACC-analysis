#!/usr/bin/env Rscript
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# cluster_characterization.R
# -----------------------------------------------------------------------------
# Post-hoc characterisation of the Monocle3 malignant trajectory. Runs AFTER
# monocle3_malignant_trajectory.R has completed and written cds_final.rds.
#
# WHY THIS IS A SEPARATE SCRIPT: everything here consumes the finished
# trajectory and touches nothing upstream of it, so it can be re-run and
# iterated in minutes instead of paying the ~35-45 minute cost of re-deriving
# HVGs, PCA, alignment, UMAP and the principal graph.
#
# It answers three questions:
#
#   1. CLUSTER CROSSWALK. Leiden renumbers clusters on every run, so cluster
#      IDs quoted in earlier reports ("cluster 3", "cluster 19") do not survive
#      re-running. This maps old IDs to new ones by cell-ID overlap, using the
#      snapshot cds_per_cell_metadata_prev_unblocked.csv taken before the
#      patient-blocked-HVG rerun.
#
#   2. FOCAL CLUSTER. The cluster whose per-patient fraction best separates the
#      outcome groups (chosen by monocle3_malignant_trajectory.R and written to
#      tables/focal_cluster_id.txt). Reports its patient composition, its
#      patient-aware pseudobulk markers, and its scores on the prespecified
#      programs -- enough to decide whether it earns a biological label.
#
#   3. DISCONNECTED COMPONENT. The malignant cells that fall outside the
#      root-containing principal-graph partition and therefore have no finite
#      pseudotime. In the previous run these were overwhelmingly one patient
#      (P19). This script does not assume that: it identifies the component
#      structurally, reports which patients actually contribute, and works
#      through the genuine-state / misannotation / doublet / low-quality
#      alternatives with the evidence available.
#
# THE PATIENT IS THE INDEPENDENT UNIT for any outcome-facing comparison. Where
# a comparison is made at the level of individual cells (the within-patient
# marker analysis, which has no patient replication by construction), that is
# stated explicitly in the output and the result is descriptive only.
#
# NOT DONE HERE, DELIBERATELY: no classifier, no predictive accuracy, no
# prognostic claim, no survival model, no claim that one malignant lineage
# differentiates into another, and no copy-number inference (the gene space is
# a 10,000-gene HVG subset, which is not a usable substrate for CNV calling).
# =============================================================================

suppressPackageStartupMessages({
  library(monocle3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(scran)
  library(scuttle)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
  library(ggrastr)
})
set.seed(42)

for (pkg in c("DESeq2", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required package '", pkg, "' is not installed. Both are already ",
         "dependencies of the companion scripts in this directory; install it ",
         "with BiocManager::install('", pkg, "') before re-running.")
}

# =============================================================================
# 0. Constants
# =============================================================================

BASE_DIR <- ACC_DATA_ROOT
SMOKE_TEST <- FALSE

OUT_DIR <- if (SMOKE_TEST)
  file.path(BASE_DIR, "outputs", "smoke_test", "monocle3_analysis") else
  file.path(BASE_DIR, "outputs", "monocle3_analysis")

CHAR_DIR <- file.path(OUT_DIR, "cluster_characterization")
dir.create(CHAR_DIR, recursive = TRUE, showWarnings = FALSE)

CDS_PATH      <- file.path(OUT_DIR, "tables", "cds_final.rds")
PERCELL_PATH  <- file.path(OUT_DIR, "tables", "cds_per_cell_metadata.csv")
FOCAL_ID_PATH <- file.path(OUT_DIR, "tables", "focal_cluster_id.txt")
# Snapshot of the SUPERSEDED annot_6 run's per-cell table, taken 2026-08-26
# immediately before the annot_5 migration overwrote it. This is what makes the
# old -> new cluster crosswalk possible; Leiden renumbers on every run, so
# without it the mapping is unrecoverable. See SUPERSEDED_annot6_RUN.md.
PREV_PERCELL_PATH <- file.path(OUT_DIR, "tables",
                                "cds_per_cell_metadata_prev_annot6.csv")

# Annotation of record. annot_6 is retained as a passenger column so this
# script's crosstab shows both levels side by side.
ANNOT_COL        <- "annot_5"
ANNOT_COL_LEGACY <- "annot_6"

OUTCOME_COLORS <- c("Good" = "forestgreen", "Poor" = "firebrick3")

# Must match PROGRAM_GENES in monocle3_malignant_trajectory.R.
PROGRAM_GENES <- list(
  myoepithelial = c("TP63", "ACTA2", "TAGLN", "CNN1", "MYH11"),
  ductal_notch  = c("KIT", "KRT7", "KRT19", "GABRP", "ELF5",
                    "NOTCH3", "HEY1", "HES1"),
  cell_cycle    = c("MKI67", "TOP2A", "CDK1", "CCNB1", "BIRC5")
)

# Lineage panels for the contamination / misannotation check. If the
# disconnected population is misannotated stroma or immune cells rather than
# tumour, it should score high on one of the non-epithelial panels.
# NOTE: ACTA2/MYH11/TAGLN appear in BOTH the myoepithelial and mural panels --
# that is real biology (myoepithelial cells are contractile), not an error, and
# it means a high mural score alone cannot distinguish a mural cell from a
# myoepithelial-like tumour cell. PTPRC / PECAM1 / DCN carry the discriminating
# weight instead.
LINEAGE_GENES <- list(
  epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"),
  immune           = c("PTPRC", "CD3E", "CD68", "LYZ", "CD14"),
  fibroblast       = c("COL1A1", "COL1A2", "DCN", "LUM", "FBLN1"),
  endothelial      = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  mural            = c("RGS5", "ACTA2", "MYH11", "TAGLN", "NOTCH3")
)

# The previous run's disconnected population was patient P19's Monocle
# cluster 19. Output filenames keep that label because the review checklist
# names those files explicitly; the CONTENT always reports the patient and
# cluster actually found in this run, and a prominent warning fires if they
# differ from the historical ones.
DISCONNECTED_LABEL   <- "P19_cluster19"
EXPECTED_FOCUS_PATIENT <- "P19"

MIN_CELLS_PER_PSEUDOBULK_GROUP <- 10
MIN_PATIENTS_FOR_PSEUDOBULK    <- 3

# =============================================================================
# Helpers
# =============================================================================

save_plot <- function(p, name, width = 8, height = 6, dpi = 600) {
  ggsave(file.path(CHAR_DIR, paste0(name, ".png")), plot = p, width = width,
         height = height, dpi = dpi, units = "in", bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(CHAR_DIR, paste0(name, ".svg")), plot = p, width = width,
         height = height, units = "in", bg = "white", device = svglite::svglite)
  message("  Saved: ", name)
}

write_out <- function(df, name) {
  write.csv(df, file.path(CHAR_DIR, name), row.names = FALSE)
  message("  Saved: ", name, " (", nrow(df), " rows)")
}

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  # Rank-based identity for Cliff's delta: avoids the n1 x n2 outer product,
  # which would be a 1,718 x 2,464 matrix here and far worse on cohort-scale
  # comparisons.
  n1 <- length(x); n2 <- length(y)
  r <- rank(c(x, y))
  u1 <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u1 / (n1 * n2) - 1
}

compare_metric <- function(df, metric, group_col, g1, g2) {
  x <- df[[metric]][df[[group_col]] == g1]
  y <- df[[metric]][df[[group_col]] == g2]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) {
    return(data.frame(metric = metric, n_group1 = length(x), n_group2 = length(y),
                      mean_group1 = NA_real_, mean_group2 = NA_real_,
                      median_group1 = NA_real_, median_group2 = NA_real_,
                      p_value = NA_real_, cliffs_delta = NA_real_,
                      stringsAsFactors = FALSE))
  }
  wt <- tryCatch(stats::wilcox.test(x, y, exact = FALSE), error = function(e) NULL)
  data.frame(
    metric = metric, n_group1 = length(x), n_group2 = length(y),
    mean_group1 = mean(x), mean_group2 = mean(y),
    median_group1 = stats::median(x), median_group2 = stats::median(y),
    p_value = if (!is.null(wt)) wt$p.value else NA_real_,
    cliffs_delta = cliffs_delta(x, y),
    stringsAsFactors = FALSE
  )
}

score_program <- function(logmat, genes, label = "") {
  present <- intersect(genes, rownames(logmat))
  if (length(present) == 0) {
    message("    [", label, "] no genes available -- all NA")
    return(rep(NA_real_, ncol(logmat)))
  }
  sub <- as.matrix(logmat[present, , drop = FALSE])
  mu <- rowMeans(sub); sdv <- apply(sub, 1, stats::sd)
  keep <- is.finite(sdv) & sdv > 0
  if (!any(keep)) return(rep(NA_real_, ncol(logmat)))
  z <- (sub[keep, , drop = FALSE] - mu[keep]) / sdv[keep]
  colMeans(z)
}

# --- patient-aware pseudobulk differential expression -------------------------
# Aggregates raw counts to one profile per (patient, group) and fits DESeq2 with
# the patient as a blocking factor, so the PATIENT is the replicate rather than
# the cell. Per-cell tests across a multi-patient cohort are pseudoreplication
# and would report implausibly small p-values.
#
# Returns NULL (with a logged reason) when the design is not estimable -- too
# few patients contribute cells to both arms. Callers must handle NULL rather
# than assume a table.
pseudobulk_de <- function(counts_mat, meta, group_col, ref_level, alt_level,
                           patient_col = "batch",
                           min_cells = MIN_CELLS_PER_PSEUDOBULK_GROUP,
                           min_patients = MIN_PATIENTS_FOR_PSEUDOBULK) {
  stopifnot(identical(colnames(counts_mat), rownames(meta)))
  meta <- meta[meta[[group_col]] %in% c(ref_level, alt_level), , drop = FALSE]

  tab <- table(meta[[patient_col]], meta[[group_col]])
  eligible <- rownames(tab)[tab[, ref_level] >= min_cells &
                              tab[, alt_level] >= min_cells]
  message("    Patients with >= ", min_cells, " cells in BOTH arms: ",
          length(eligible), " / ", nrow(tab))
  if (length(eligible) < min_patients) {
    message("    NOT ESTIMABLE: need >= ", min_patients,
            " patients for a patient-blocked design, have ", length(eligible),
            ". Skipping pseudobulk DE.")
    return(list(result = NULL, n_patients = length(eligible),
                patients = eligible,
                reason = paste0("not estimable: only ", length(eligible),
                                " patient(s) contribute >= ", min_cells,
                                " cells to both arms")))
  }

  meta <- meta[meta[[patient_col]] %in% eligible, , drop = FALSE]
  sub_counts <- counts_mat[, rownames(meta), drop = FALSE]

  key <- paste(meta[[patient_col]], meta[[group_col]], sep = "__")
  key_f <- factor(key)
  ind <- Matrix::sparseMatrix(i = seq_along(key_f), j = as.integer(key_f),
                              x = 1, dims = c(length(key_f), nlevels(key_f)),
                              dimnames = list(NULL, levels(key_f)))
  pb <- as.matrix(sub_counts %*% ind)

  pb_coldata <- data.frame(
    sample = colnames(pb),
    patient = sub("__.*$", "", colnames(pb)),
    group = sub("^.*__", "", colnames(pb)),
    stringsAsFactors = FALSE
  )
  rownames(pb_coldata) <- pb_coldata$sample
  pb_coldata$group <- factor(pb_coldata$group, levels = c(ref_level, alt_level))
  pb_coldata$patient <- factor(pb_coldata$patient)

  keep_genes <- rowSums(pb) >= 10 & rowSums(pb > 0) >= 2
  pb <- pb[keep_genes, , drop = FALSE]
  message("    Pseudobulk matrix: ", nrow(pb), " genes x ", ncol(pb),
          " samples (", length(eligible), " patients)")

  # DESeq2 requires integer-valued counts but tolerates doubles that are exactly
  # integers. round() rather than storage.mode(): pseudobulk column sums can be
  # large, and an integer coercion that overflows would silently produce NA.
  pb <- round(pb)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pb, colData = pb_coldata,
                                         design = ~ patient + group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("group", alt_level, ref_level))
  out <- as.data.frame(res)
  out$gene <- rownames(out)
  out <- out[order(out$padj, -abs(out$log2FoldChange)), ]
  out <- out[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat",
                 "pvalue", "padj")]
  list(result = out, n_patients = length(eligible), patients = eligible,
       reason = NA_character_)
}

# --- within-group cell-level markers ------------------------------------------
# Used ONLY for the within-patient comparison, where there is exactly one
# patient and therefore no possibility of patient-level replication. Reports
# AUC (the rank-based effect size) alongside the p-value, because with
# thousands of cells the p-value is essentially a function of n and carries no
# information about effect magnitude.
cell_level_markers <- function(logmat, group_vec, g1, g2) {
  sce <- SingleCellExperiment(assays = list(logcounts = logmat))
  fm <- scran::findMarkers(sce, groups = group_vec, test.type = "wilcox",
                            direction = "any", pval.type = "all")
  res <- as.data.frame(fm[[g1]])
  auc_col <- grep("^AUC", colnames(res), value = TRUE)[1]
  out <- data.frame(gene = rownames(res),
                    AUC = res[[auc_col]],
                    p_value = res$p.value,
                    p_adj_BH = res$FDR,
                    stringsAsFactors = FALSE)
  m1 <- Matrix::rowMeans(logmat[, group_vec == g1, drop = FALSE])
  m2 <- Matrix::rowMeans(logmat[, group_vec == g2, drop = FALSE])
  out$mean_logcounts_group1 <- m1[out$gene]
  out$mean_logcounts_group2 <- m2[out$gene]
  out$logFC_group1_vs_group2 <- out$mean_logcounts_group1 - out$mean_logcounts_group2
  out[order(-out$AUC), ]
}

dot_plot <- function(logmat, group_vec, gene_sets, title, subtitle = NULL) {
  # Standard marker dot plot: dot size = fraction of cells with non-zero
  # expression, colour = mean expression z-scored across the groups shown.
  genes <- unlist(gene_sets, use.names = FALSE)
  panel <- rep(names(gene_sets), lengths(gene_sets))
  keep <- genes %in% rownames(logmat)
  genes <- genes[keep]; panel <- panel[keep]

  rows <- lapply(unique(group_vec), function(g) {
    sub <- logmat[genes, group_vec == g, drop = FALSE]
    data.frame(gene = genes, panel = panel, group = g,
               mean_expr = as.numeric(Matrix::rowMeans(sub)),
               frac_expr = as.numeric(Matrix::rowMeans(sub > 0)),
               stringsAsFactors = FALSE)
  })
  # Group by (gene, panel), not gene alone: several markers appear in more than
  # one panel (ACTA2/MYH11/TAGLN/NOTCH3 are both myoepithelial and mural), and
  # grouping on gene alone would z-score across the duplicated copies.
  df <- do.call(rbind, rows) %>%
    group_by(gene, panel) %>%
    mutate(z = if (stats::sd(mean_expr) > 0)
             (mean_expr - mean(mean_expr)) / stats::sd(mean_expr) else 0) %>%
    ungroup() %>%
    mutate(gene = factor(gene, levels = rev(unique(genes))),
           panel = factor(panel, levels = names(gene_sets)))

  ggplot(df, aes(x = group, y = gene, size = frac_expr, color = z)) +
    geom_point() +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                          midpoint = 0) +
    scale_size_continuous(range = c(0.5, 6), limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL,
         size = "Fraction\nexpressing", color = "Mean expr.\n(z across groups)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text.y = element_text(angle = 0))
}

# =============================================================================
# 1. Load
# =============================================================================

message("=== Step 1: Loading trajectory outputs ===")

for (p in c(CDS_PATH, PERCELL_PATH)) {
  if (!file.exists(p))
    stop("Required input not found: ", p,
         "\nRun monocle3_malignant_trajectory.R to completion first.")
}

cds <- readRDS(CDS_PATH)
percell <- read.csv(PERCELL_PATH, stringsAsFactors = FALSE)
message("  cds: ", ncol(cds), " cells x ", nrow(cds), " genes")
message("  per-cell metadata: ", nrow(percell), " rows x ", ncol(percell), " cols")
stopifnot(nrow(percell) == ncol(cds))

rownames(percell) <- percell$cell_id
percell <- percell[colnames(cds), , drop = FALSE]

focal_cluster <- if (file.exists(FOCAL_ID_PATH))
  trimws(readLines(FOCAL_ID_PATH)[1]) else NA_character_
message("  Focal cluster: ", ifelse(is.na(focal_cluster), "<none recorded>", focal_cluster))

counts_mat <- SingleCellExperiment::counts(cds)
sce <- SingleCellExperiment(assays = list(counts = counts_mat))
sce <- scuttle::logNormCounts(sce)
logmat <- SummarizedExperiment::assay(sce, "logcounts")
rm(sce); gc()

# Cell-cycle program is scored here rather than in the main script because it is
# only used for cluster/component characterisation. S_score and G2M_score from
# the upstream Scanpy pipeline are carried in percell and reported alongside it.
message("  Scoring cell-cycle program")
percell$score_cell_cycle <- score_program(logmat, PROGRAM_GENES$cell_cycle, "cell_cycle")

# =============================================================================
# 2. Old -> new cluster crosswalk
# =============================================================================

message("=== Step 2: Cluster crosswalk (previous run -> this run) ===")

if (file.exists(PREV_PERCELL_PATH)) {
  prev <- read.csv(PREV_PERCELL_PATH, stringsAsFactors = FALSE)
  message("  Previous run metadata: ", nrow(prev), " cells")

  joined <- percell %>%
    select(cell_id, new_cluster = monocle_cluster) %>%
    inner_join(prev %>% select(cell_id, old_cluster = monocle_cluster),
               by = "cell_id")
  message("  Cells shared between runs: ", nrow(joined), " / ", nrow(percell))

  n_old <- joined %>% count(old_cluster, name = "n_old")
  n_new <- joined %>% count(new_cluster, name = "n_new")

  crosswalk <- joined %>%
    count(old_cluster, new_cluster, name = "n_shared") %>%
    left_join(n_old, by = "old_cluster") %>%
    left_join(n_new, by = "new_cluster") %>%
    mutate(jaccard = n_shared / (n_old + n_new - n_shared),
           frac_of_old = n_shared / n_old,
           frac_of_new = n_shared / n_new) %>%
    arrange(old_cluster, desc(n_shared))

  write_out(crosswalk, "cluster_crosswalk_old_vs_new.csv")

  best_match <- crosswalk %>%
    group_by(old_cluster) %>%
    slice_max(jaccard, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(jaccard))
  write_out(best_match, "cluster_crosswalk_best_match.csv")

  crosswalk_lines <- c(
    "Monocle cluster crosswalk: previous (pooled-HVG) run -> current (patient-blocked-HVG) run",
    "",
    paste0("Cells shared between runs: ", nrow(joined)),
    paste0("Clusters in previous run: ", nrow(n_old),
           "; clusters in current run: ", nrow(n_new)),
    "",
    "Best match for the clusters named in the review:"
  )
  for (oc in c("3", "19")) {
    row <- best_match[best_match$old_cluster == oc, ]
    if (nrow(row) == 1) {
      crosswalk_lines <- c(crosswalk_lines, sprintf(
        "  old cluster %s (n = %d) -> new cluster %s (n = %d): %d shared cells, Jaccard = %.3f, %.1f%% of the old cluster",
        oc, row$n_old, row$new_cluster, row$n_new, row$n_shared, row$jaccard,
        100 * row$frac_of_old))
    } else {
      crosswalk_lines <- c(crosswalk_lines,
                           paste0("  old cluster ", oc, ": not present in the snapshot"))
    }
  }
  crosswalk_lines <- c(crosswalk_lines, "",
    "A low Jaccard does not mean the biology changed. Leiden splits and merges",
    "clusters freely between runs, so a coherent cell state can be distributed",
    "across several new clusters. Read this table together with the program",
    "scores rather than as a claim about stability.")
  writeLines(crosswalk_lines, file.path(CHAR_DIR, "cluster_crosswalk_summary.txt"))
  message("  Saved: cluster_crosswalk_summary.txt")
} else {
  message("  NOTE: no snapshot at ", PREV_PERCELL_PATH,
          " -- crosswalk to the previous run's cluster numbering is unavailable.")
  writeLines(c(
    "Cluster crosswalk unavailable.",
    "",
    paste0("Expected snapshot not found: ", PREV_PERCELL_PATH),
    "The previous run's per-cell cluster labels were not preserved, so old",
    "cluster IDs (e.g. 'cluster 3', 'cluster 19') cannot be mapped onto this",
    "run's numbering. Clusters here must be identified by their markers and",
    "composition instead."),
    file.path(CHAR_DIR, "cluster_crosswalk_summary.txt"))
}

# =============================================================================
# 3. Focal cluster characterisation
# =============================================================================

message("=== Step 3: Focal cluster characterisation ===")

if (is.na(focal_cluster)) {
  message("  SKIPPED: no focal cluster recorded by the trajectory script.")
} else {
  percell$in_focal <- ifelse(percell$monocle_cluster == focal_cluster,
                             "focal", "other")
  n_focal <- sum(percell$in_focal == "focal")
  message("  Cluster ", focal_cluster, ": ", n_focal, " cells (",
          round(100 * n_focal / nrow(percell), 2), "% of malignant cells)")

  # --- patient composition ---
  focal_comp <- percell %>%
    filter(in_focal == "focal") %>%
    count(batch, outcome_group, name = "n_cells_in_cluster") %>%
    left_join(percell %>% count(batch, name = "n_cells_total"), by = "batch") %>%
    mutate(frac_of_cluster = n_cells_in_cluster / sum(n_cells_in_cluster),
           frac_of_patient = n_cells_in_cluster / n_cells_total) %>%
    arrange(desc(n_cells_in_cluster))
  write_out(focal_comp, "focal_cluster_patient_composition.csv")
  message("  Contributing patients: ", nrow(focal_comp), " / ",
          length(unique(percell$batch)),
          "; largest single-patient share = ",
          round(100 * max(focal_comp$frac_of_cluster), 1), "%")

  # --- program scores, cluster vs rest ---
  prog_metrics <- intersect(
    c("CytoTRACE2_Score", "score_myoepithelial", "score_ductal_notch",
      "score_cell_cycle", "S_score", "G2M_score", "pseudotime"),
    colnames(percell))
  focal_prog <- do.call(rbind, lapply(prog_metrics, function(m) {
    d <- percell
    if (m == "pseudotime") d <- d[is.finite(d$pseudotime), ]
    compare_metric(d, m, "in_focal", "focal", "other")
  }))
  names(focal_prog)[names(focal_prog) == "n_group1"] <- "n_focal"
  names(focal_prog)[names(focal_prog) == "n_group2"] <- "n_other"
  focal_prog$comparison <- paste0("Monocle cluster ", focal_cluster,
                                   " vs all other malignant cells (cell-level, descriptive)")
  write_out(focal_prog, "focal_cluster_program_scores.csv")

  # --- patient-aware pseudobulk markers ---
  message("  Pseudobulk markers (design = ~ patient + in_focal)")
  pb_meta <- percell[, c("batch", "in_focal")]
  rownames(pb_meta) <- percell$cell_id
  focal_pb <- pseudobulk_de(counts_mat, pb_meta, "in_focal",
                             ref_level = "other", alt_level = "focal")

  if (!is.null(focal_pb$result)) {
    write_out(focal_pb$result, "focal_cluster_pseudobulk_markers.csv")
    write_out(data.frame(
      design = "~ patient + in_focal",
      contrast = "in_focal: focal vs other",
      n_patients = focal_pb$n_patients,
      patients = paste(focal_pb$patients, collapse = ";"),
      min_cells_per_group = MIN_CELLS_PER_PSEUDOBULK_GROUP,
      n_genes_tested = nrow(focal_pb$result),
      n_padj_lt_0.05 = sum(focal_pb$result$padj < 0.05, na.rm = TRUE),
      gene_space_note = paste0("10,000-gene HVG subset selected upstream in ",
                                "Scanpy, not the full transcriptome"),
      stringsAsFactors = FALSE
    ), "focal_cluster_pseudobulk_markers_design.csv")

    top_up <- focal_pb$result %>% filter(!is.na(padj), padj < 0.05,
                                          log2FoldChange > 0) %>%
      arrange(desc(log2FoldChange)) %>% head(25)
    top_dn <- focal_pb$result %>% filter(!is.na(padj), padj < 0.05,
                                          log2FoldChange < 0) %>%
      arrange(log2FoldChange) %>% head(25)
    message("    Significant genes (padj < 0.05): ",
            sum(focal_pb$result$padj < 0.05, na.rm = TRUE), " / ",
            nrow(focal_pb$result))
    message("    Top up in cluster ", focal_cluster, ": ",
            paste(head(top_up$gene, 12), collapse = ", "))
    message("    Top down: ", paste(head(top_dn$gene, 12), collapse = ", "))
  } else {
    writeLines(c("Focal-cluster pseudobulk DE not estimable.", "",
                 focal_pb$reason),
               file.path(CHAR_DIR, "focal_cluster_pseudobulk_markers_NOT_ESTIMABLE.txt"))
    top_up <- data.frame(gene = character(0), log2FoldChange = numeric(0))
    top_dn <- data.frame(gene = character(0), log2FoldChange = numeric(0))
  }

  # --- dot plot over the prespecified programs ---
  p_focal_dot <- dot_plot(
    logmat, percell$in_focal, PROGRAM_GENES,
    title = paste0("Prespecified programs: Monocle cluster ", focal_cluster,
                    " vs other malignant cells"),
    subtitle = paste0(n_focal, " cells in cluster ", focal_cluster, " from ",
                       nrow(focal_comp), " patients. Cluster selected by ",
                       "scanning all clusters -- exploratory."))
  save_plot(p_focal_dot, "focal_cluster_program_dotplot", width = 6.5, height = 7)

  # --- per-patient contribution figure ---
  p_focal_comp <- ggplot(focal_comp,
                          aes(x = reorder(batch, frac_of_patient),
                              y = frac_of_patient, fill = outcome_group)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = OUTCOME_COLORS) +
    labs(title = paste0("Fraction of each patient's malignant cells in cluster ",
                         focal_cluster),
         subtitle = paste0(nrow(focal_comp), " of ",
                            length(unique(percell$batch)),
                            " patients contribute cells to this cluster"),
         x = "Patient (batch)", y = "Fraction of that patient's malignant cells",
         fill = "Outcome") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  save_plot(p_focal_comp, "focal_cluster_patient_contribution",
            width = 7, height = max(5, 0.3 * nrow(focal_comp)))

  # --- interpretation note, written from the evidence rather than asserted ---
  myo_row <- focal_prog[focal_prog$metric == "score_myoepithelial", ]
  duc_row <- focal_prog[focal_prog$metric == "score_ductal_notch", ]
  cyc_row <- focal_prog[focal_prog$metric == "score_cell_cycle", ]
  cyt_row <- focal_prog[focal_prog$metric == "CytoTRACE2_Score", ]
  .d <- function(r) if (nrow(r) == 1 && is.finite(r$cliffs_delta)) r$cliffs_delta else NA_real_
  .lab <- function(d) if (is.na(d)) "unknown" else
    if (d > 0.3) "higher" else if (d < -0.3) "lower" else "comparable"

  label_lines <- c(
    paste0("Monocle cluster ", focal_cluster, ": evidence summary"),
    strrep("=", 60), "",
    paste0("Size: ", n_focal, " cells (",
           round(100 * n_focal / nrow(percell), 2),
           "% of malignant cells) from ", nrow(focal_comp), " of ",
           length(unique(percell$batch)), " patients."),
    paste0("Largest single-patient share of the cluster: ",
           round(100 * max(focal_comp$frac_of_cluster), 1), "% (",
           focal_comp$batch[1], ")."),
    "",
    "Program scores versus all other malignant cells (Cliff's delta, cell-level):",
    paste0("  myoepithelial-like    : ", .lab(.d(myo_row)),
           "  (delta = ", round(.d(myo_row), 3), ")"),
    paste0("  ductal/NOTCH          : ", .lab(.d(duc_row)),
           "  (delta = ", round(.d(duc_row), 3), ")"),
    paste0("  cell cycle            : ", .lab(.d(cyc_row)),
           "  (delta = ", round(.d(cyc_row), 3), ")"),
    paste0("  CytoTRACE2            : ", .lab(.d(cyt_row)),
           "  (delta = ", round(.d(cyt_row), 3), ")"),
    "",
    if (nrow(top_up) > 0)
      paste0("Top pseudobulk markers up in this cluster: ",
             paste(head(top_up$gene, 15), collapse = ", "))
    else "Pseudobulk markers: not estimable (see the NOT_ESTIMABLE file).",
    if (nrow(top_dn) > 0)
      paste0("Top pseudobulk markers down in this cluster: ",
             paste(head(top_dn$gene, 15), collapse = ", "))
    else "",
    "",
    "How to read this:",
    "  A biological label is warranted only if the pseudobulk markers and the",
    "  program scores point the same way. If the cluster is myoepithelial-high,",
    "  CytoTRACE2-high and ductal-low, 'early myoepithelial-like state' is",
    "  supportable. If the markers are mixed, or the cluster is dominated by one",
    "  or two patients, report it as an unlabelled cluster and say so.",
    "",
    "  This cluster was selected by scanning every cluster for outcome",
    "  association. Its p-value is therefore not a valid confirmatory test of",
    "  this cluster, and nothing here establishes prognostic utility."
  )
  writeLines(label_lines, file.path(CHAR_DIR, "focal_cluster_interpretation.txt"))
  message("  Saved: focal_cluster_interpretation.txt")
}

# =============================================================================
# 4. Disconnected component
# =============================================================================

message("=== Step 4: Disconnected malignant population ===")

percell$component <- ifelse(percell$has_finite_pseudotime %in% c(TRUE, "TRUE"),
                            "connected", "disconnected")
n_disc <- sum(percell$component == "disconnected")
message("  Cells without finite pseudotime: ", n_disc, " / ", nrow(percell),
        " (", round(100 * n_disc / nrow(percell), 2), "%)")

if (n_disc == 0) {
  writeLines(c(
    "No disconnected malignant population in this run.",
    "",
    "Every malignant cell received a finite pseudotime, i.e. all cells lie in",
    "the root-containing principal-graph partition. The population previously",
    "reported as patient P19 / Monocle cluster 19 did NOT reappear after",
    "patient-blocked HVG selection.",
    "",
    "This is itself the answer to the review question: the disconnection was a",
    "consequence of pooled-HVG feature selection, not a stable property of the",
    "data."),
    file.path(CHAR_DIR, paste0(DISCONNECTED_LABEL, "_disconnection_status.txt")))
  message("  No disconnected cells -- wrote status file and skipped the rest of Step 4.")
} else {

  # --- 4.1 patient composition ---
  disc_comp <- percell %>%
    filter(component == "disconnected") %>%
    count(batch, outcome_group, name = "n_disconnected") %>%
    left_join(percell %>% count(batch, name = "n_cells_total"), by = "batch") %>%
    mutate(frac_of_component = n_disconnected / sum(n_disconnected),
           frac_of_patient = n_disconnected / n_cells_total) %>%
    arrange(desc(n_disconnected))
  write_out(disc_comp, "monocle_cluster19_patient_composition.csv")

  focus_patient <- disc_comp$batch[1]
  message("  Dominant patient: ", focus_patient, " (",
          disc_comp$n_disconnected[1], " cells, ",
          round(100 * disc_comp$frac_of_component[1], 1),
          "% of the component; ",
          round(100 * disc_comp$frac_of_patient[1], 1),
          "% of that patient's malignant cells)")
  message("  Patients contributing: ", nrow(disc_comp))
  if (!identical(focus_patient, EXPECTED_FOCUS_PATIENT))
    message("  *** NOTE: the dominant patient is ", focus_patient,
            ", NOT the historically reported ", EXPECTED_FOCUS_PATIENT,
            ". Output filenames retain the historical label for continuity; ",
            "the contents describe ", focus_patient, ". ***")

  disc_clusters <- percell %>%
    filter(component == "disconnected") %>%
    count(monocle_cluster, monocle_partition, name = "n") %>%
    arrange(desc(n))
  write_out(disc_clusters, "disconnected_component_cluster_composition.csv")

  # --- 4.2 within-patient markers (PRIMARY comparison) ---
  # Controls for patient-specific expression by construction: both arms come
  # from the same patient. There is exactly one patient here, so no pseudobulk
  # design exists and this is necessarily a cell-level, descriptive comparison.
  message("  Within-patient markers: ", focus_patient,
          " disconnected vs connected")
  fp_idx <- which(percell$batch == focus_patient)
  fp_meta <- percell[fp_idx, ]
  n_fp_disc <- sum(fp_meta$component == "disconnected")
  n_fp_conn <- sum(fp_meta$component == "connected")
  message("    ", focus_patient, ": ", n_fp_disc, " disconnected vs ",
          n_fp_conn, " connected malignant cells")

  if (n_fp_disc >= 20 && n_fp_conn >= 20) {
    fp_markers <- cell_level_markers(logmat[, fp_idx, drop = FALSE],
                                      fp_meta$component,
                                      "disconnected", "connected")
    names(fp_markers) <- sub("group1", "disconnected", names(fp_markers))
    names(fp_markers) <- sub("group2", "connected", names(fp_markers))
    fp_markers$patient <- focus_patient
    fp_markers$comparison <- paste0(focus_patient,
      " disconnected vs connected malignant cells -- WITHIN-PATIENT, CELL-LEVEL, DESCRIPTIVE.",
      " One patient affords no replication, so p-values reflect cell counts;",
      " read AUC as the effect size.")
    # The checklist names this file explicitly, so that name is always written.
    # If the dominant patient turns out not to be P19, a correctly-named
    # copy is written alongside it rather than instead of it.
    write_out(fp_markers, "P19_cluster19_vs_other_P19_markers.csv")
    actual_name <- paste0(DISCONNECTED_LABEL, "_vs_other_", focus_patient,
                          "_markers.csv")
    if (actual_name != "P19_cluster19_vs_other_P19_markers.csv")
      write_out(fp_markers, actual_name)
    message("    Top enriched in the disconnected population: ",
            paste(head(fp_markers$gene, 15), collapse = ", "))
  } else {
    message("    SKIPPED: too few cells in one arm for a marker comparison.")
    fp_markers <- NULL
  }

  # --- 4.3 cohort comparison, patient-aware where replication permits ---
  message("  Cohort comparison: disconnected vs all other malignant cells")
  pb_meta2 <- percell[, c("batch", "component")]
  rownames(pb_meta2) <- percell$cell_id
  disc_pb <- pseudobulk_de(counts_mat, pb_meta2, "component",
                            ref_level = "connected", alt_level = "disconnected")
  if (!is.null(disc_pb$result)) {
    disc_pb$result$design_note <- paste0("~ patient + component, ",
                                          disc_pb$n_patients, " patients")
    write_out(disc_pb$result, "cohort_disconnected_vs_rest_pseudobulk.csv")
  } else {
    write_out(data.frame(
      status = "NOT_ESTIMABLE",
      reason = disc_pb$reason,
      n_patients_eligible = disc_pb$n_patients,
      patients_eligible = paste(disc_pb$patients, collapse = ";"),
      min_cells_per_group = MIN_CELLS_PER_PSEUDOBULK_GROUP,
      min_patients_required = MIN_PATIENTS_FOR_PSEUDOBULK,
      interpretation = paste0(
        "The disconnected component is confined to too few patients for a ",
        "patient-blocked cohort design. That concentration is itself the ",
        "finding: this is a patient-restricted population, not a shared ",
        "malignant state. Use the within-patient comparison instead."),
      stringsAsFactors = FALSE
    ), "cohort_disconnected_vs_rest_pseudobulk.csv")
  }

  # --- 4.4 QC, doublet and program summary ---
  message("  QC / doublet / program comparison")
  qc_metrics <- intersect(
    c("total_counts", "n_genes", "n_genes_by_counts", "pct_counts_mt",
      "total_counts_mt", "doublet_score", "dd_score", "scrublet_score",
      "solo_score", "CytoTRACE2_Score", "score_myoepithelial",
      "score_ductal_notch", "score_cell_cycle", "S_score", "G2M_score"),
    colnames(percell))

  # Within the focus patient (the comparison that controls for patient), then
  # cohort-wide for context.
  qc_within <- do.call(rbind, lapply(qc_metrics, function(m)
    compare_metric(fp_meta, m, "component", "disconnected", "connected")))
  qc_within$scope <- paste0("within patient ", focus_patient)

  qc_cohort <- do.call(rbind, lapply(qc_metrics, function(m)
    compare_metric(percell, m, "component", "disconnected", "connected")))
  qc_cohort$scope <- "all patients"

  qc_all <- rbind(qc_within, qc_cohort)
  names(qc_all)[names(qc_all) == "n_group1"] <- "n_disconnected"
  names(qc_all)[names(qc_all) == "n_group2"] <- "n_connected"
  names(qc_all)[names(qc_all) == "mean_group1"] <- "mean_disconnected"
  names(qc_all)[names(qc_all) == "mean_group2"] <- "mean_connected"
  names(qc_all)[names(qc_all) == "median_group1"] <- "median_disconnected"
  names(qc_all)[names(qc_all) == "median_group2"] <- "median_connected"
  qc_all$note <- "cell-level comparison; p-values scale with cell count, read cliffs_delta"
  write_out(qc_all, paste0(DISCONNECTED_LABEL, "_QC_and_program_summary.csv"))

  # Doublet-call rates, which are categorical and so not covered above.
  if ("doublet_consensus_2of3_str" %in% colnames(percell)) {
    dbl_tab <- fp_meta %>%
      count(component, doublet_consensus_2of3_str, name = "n") %>%
      group_by(component) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
    write_out(dbl_tab, paste0(DISCONNECTED_LABEL, "_doublet_consensus_rates.csv"))
  }

  # --- 4.5 malignancy and annotation check ---
  message("  Malignancy / annotation check")
  ann_cols <- intersect(c(ANNOT_COL, ANNOT_COL_LEGACY, "prect_2"), colnames(percell))
  ann_tab <- fp_meta %>%
    count(across(all_of(c("component", ann_cols))), name = "n") %>%
    group_by(component) %>%
    mutate(frac_of_component = n / sum(n)) %>%
    ungroup()

  lineage_scores <- lapply(names(LINEAGE_GENES), function(nm) {
    s <- score_program(logmat[, fp_idx, drop = FALSE], LINEAGE_GENES[[nm]], nm)
    d <- data.frame(component = fp_meta$component, score = s,
                    stringsAsFactors = FALSE)
    r <- compare_metric(d, "score", "component", "disconnected", "connected")
    r$metric <- paste0("lineage_score_", nm)
    r
  })
  lineage_df <- do.call(rbind, lineage_scores)
  names(lineage_df)[names(lineage_df) == "n_group1"] <- "n_disconnected"
  names(lineage_df)[names(lineage_df) == "n_group2"] <- "n_connected"
  names(lineage_df)[names(lineage_df) == "mean_group1"] <- "mean_disconnected"
  names(lineage_df)[names(lineage_df) == "mean_group2"] <- "mean_connected"
  names(lineage_df)[names(lineage_df) == "median_group1"] <- "median_disconnected"
  names(lineage_df)[names(lineage_df) == "median_group2"] <- "median_connected"
  lineage_df$scope <- paste0("within patient ", focus_patient)
  lineage_df$caveat <- paste0(
    "prect_2 is a 1:1 deterministic recode of the annot_5 malignant flag (every '- Tumor' cell is ",
    "'Cancer profile'), so it is a relabelling of the annotation and NOT ",
    "independent malignancy evidence. No CNV inference exists for this project.")

  write.csv(ann_tab, file.path(CHAR_DIR,
    paste0(DISCONNECTED_LABEL, "_annotation_crosstab.csv")), row.names = FALSE)
  message("  Saved: ", DISCONNECTED_LABEL, "_annotation_crosstab.csv")
  write_out(lineage_df, "P19_malignancy_and_annotation_check.csv")

  p_lineage_dot <- dot_plot(
    logmat[, fp_idx, drop = FALSE], fp_meta$component, LINEAGE_GENES,
    title = paste0("Lineage markers: ", focus_patient,
                    " disconnected vs connected malignant cells"),
    subtitle = paste0("Contamination / misannotation check. ACTA2, MYH11, TAGLN ",
                       "and NOTCH3 are shared between the mural and\n",
                       "myoepithelial panels -- PTPRC, PECAM1 and DCN carry the ",
                       "discriminating weight."))
  save_plot(p_lineage_dot, paste0(DISCONNECTED_LABEL, "_lineage_dotplot"),
            width = 6.5, height = 8)

  p_prog_dot <- dot_plot(
    logmat[, fp_idx, drop = FALSE], fp_meta$component, PROGRAM_GENES,
    title = paste0("Prespecified programs: ", focus_patient,
                    " disconnected vs connected malignant cells"))
  save_plot(p_prog_dot, paste0(DISCONNECTED_LABEL, "_program_dotplot"),
            width = 6.5, height = 7)

  # --- 4.6 integration check: both embeddings ---
  # If the population is already separate in the original scVI embedding, it was
  # distinct before this script's alignment ever ran; if it is only separate in
  # the Monocle embedding, integration or graph construction isolated it.
  message("  Integration check: plotting both embeddings")
  fp_plot <- fp_meta
  emb_specs <- list(
    list("X_umap_1", "X_umap_2", "Original scVI UMAP (pre-alignment)"),
    list("monocle_umap_1", "monocle_umap_2", "Monocle3 UMAP (patient-aligned)")
  )
  plots <- list()
  for (es in emb_specs) {
    if (!all(c(es[[1]], es[[2]]) %in% colnames(fp_plot))) next
    p <- ggplot(fp_plot, aes(x = .data[[es[[1]]]], y = .data[[es[[2]]]],
                              color = component)) +
      geom_point(size = 0.5, stroke = 0, alpha = 0.7, na.rm = TRUE) +
      scale_color_manual(values = c(connected = "grey65",
                                     disconnected = "#d1495b")) +
      labs(title = es[[3]], x = NULL, y = NULL, color = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
    plots[[length(plots) + 1]] <- ggrastr::rasterise(p, dpi = 600)
  }
  if (length(plots) > 0) {
    p_emb <- if (requireNamespace("patchwork", quietly = TRUE) && length(plots) > 1)
      patchwork::wrap_plots(plots, nrow = 1) +
        patchwork::plot_annotation(
          title = paste0(focus_patient,
                          ": disconnected vs connected malignant cells"),
          subtitle = paste0(n_fp_disc, " disconnected, ", n_fp_conn,
                             " connected. Separation in the LEFT panel means the ",
                             "population was already distinct before alignment."))
    else plots[[1]]
    save_plot(p_emb, paste0(DISCONNECTED_LABEL, "_connected_vs_disconnected_umap"),
              width = 11, height = 5.8)
  }

  # Cohort-wide view, so the component can be located in the full embedding.
  p_cohort_emb <- ggplot(percell %>% arrange(component == "disconnected"),
                          aes(x = monocle_umap_1, y = monocle_umap_2,
                              color = component)) +
    geom_point(size = 0.35, stroke = 0, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = c(connected = "grey80", disconnected = "#d1495b")) +
    labs(title = "Disconnected malignant cells in the Monocle3 embedding",
         subtitle = paste0(n_disc, " of ", nrow(percell),
                            " malignant cells lie outside the root partition and ",
                            "have no finite pseudotime"),
         x = "Monocle UMAP 1", y = "Monocle UMAP 2", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  save_plot(ggrastr::rasterise(p_cohort_emb, dpi = 600),
            "disconnected_component_cohort_umap", width = 7, height = 6.5)

  # --- 4.7 disconnection status vs the previous run ---
  status <- c(
    paste0("Disconnected malignant population: status after patient-blocked HVG selection"),
    strrep("=", 78), "",
    paste0("Run: patient-blocked HVGs (scran::modelGeneVar(block = batch)), ",
           "full Monocle3 re-derivation."),
    "",
    "FINDING",
    paste0("  ", n_disc, " of ", nrow(percell), " malignant cells (",
           round(100 * n_disc / nrow(percell), 2),
           "%) still fall outside the root-containing partition and have no",
           " finite pseudotime."),
    paste0("  Dominant patient: ", focus_patient, " with ",
           disc_comp$n_disconnected[1], " cells (",
           round(100 * disc_comp$frac_of_component[1], 1),
           "% of the component, ",
           round(100 * disc_comp$frac_of_patient[1], 1),
           "% of that patient's malignant cells)."),
    paste0("  Patients contributing at all: ", nrow(disc_comp), " (",
           paste(sprintf("%s=%d", disc_comp$batch, disc_comp$n_disconnected),
                 collapse = ", "), ")."),
    paste0("  Monocle clusters involved: ",
           paste(sprintf("%s (n=%d)", disc_clusters$monocle_cluster,
                         disc_clusters$n), collapse = ", "), "."),
    ""
  )

  if (exists("crosswalk") && is.data.frame(crosswalk)) {
    old19 <- percell %>%
      inner_join(prev %>% select(cell_id, old_cluster = monocle_cluster),
                 by = "cell_id") %>%
      filter(old_cluster == "19")
    if (nrow(old19) > 0) {
      still_disc <- sum(old19$component == "disconnected")
      status <- c(status,
        "PERSISTENCE VERSUS THE PREVIOUS RUN",
        paste0("  The previous run's cluster 19 held ", nrow(old19),
               " cells that are still present in this run."),
        paste0("  Of those, ", still_disc, " (",
               round(100 * still_disc / nrow(old19), 1),
               "%) remain disconnected here."),
        paste0("  Verdict: ", if (still_disc / nrow(old19) > 0.8)
                 "the population PERSISTS -- blocked HVG selection did not dissolve it."
               else if (still_disc / nrow(old19) < 0.2)
                 "the population is LARGELY RESOLVED -- it was substantially an artefact of pooled HVG selection."
               else
                 "the population is PARTIALLY resolved -- inspect the composition table before drawing conclusions."),
        "")
    }
  } else {
    status <- c(status,
      "PERSISTENCE VERSUS THE PREVIOUS RUN",
      "  Not assessable: no snapshot of the previous run's cluster labels.",
      "")
  }

  status <- c(status,
    "HOW THIS IS HANDLED DOWNSTREAM",
    "  These cells are reported separately and excluded from the pseudotime",
    "  tertile summaries; per-patient exclusion is visible as frac_no_pseudotime",
    "  in patient_level_summary_table.csv. They are NOT assigned an artificial",
    "  pseudotime and NOT forced into the root partition.",
    "",
    paste0("  Note that ", focus_patient, " loses ",
           round(100 * disc_comp$frac_of_patient[1], 1),
           "% of its malignant cells from every pseudotime summary as a result."),
    "",
    "MALIGNANCY EVIDENCE",
    "  No copy-number inference (inferCNV / copyKAT / numbat) exists for this",
    "  project, and the gene space here is a 10,000-gene HVG subset that is not",
    "  a usable substrate for CNV calling. The only malignancy label available,",
    "  prect_2, is a 1:1 deterministic recode of the annot_5 malignant flag -- every '- Tumor' cell is",
    "  'Cancer profile' -- so it restates the annotation and cannot independently",
    "  confirm malignancy. See the annotation crosstab and lineage marker check",
    "  for what the expression data alone can say.",
    "",
    "WHAT TO CONCLUDE",
    "  Read this together with:",
    paste0("    - ", DISCONNECTED_LABEL, "_QC_and_program_summary.csv (quality/doublets)"),
    paste0("    - P19_malignancy_and_annotation_check.csv (contamination)"),
    paste0("    - ", DISCONNECTED_LABEL, "_connected_vs_disconnected_umap (integration)"),
    "  A population that is QC-normal, epithelial-marker-positive, immune- and",
    "  stroma-negative, and already separate in the pre-alignment embedding is a",
    "  genuine patient-restricted malignant state and should be reported as a",
    "  separate component. One with elevated doublet scores, low gene counts or",
    "  mixed lineage markers is a technical artefact and should be excluded."
  )
  writeLines(status, file.path(CHAR_DIR,
    paste0(DISCONNECTED_LABEL, "_disconnection_status.txt")))
  message("  Saved: ", DISCONNECTED_LABEL, "_disconnection_status.txt")
}

message("=== cluster_characterization.R complete -> ", CHAR_DIR, " ===")
