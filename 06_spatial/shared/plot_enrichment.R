# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# plot_enrichment.R — figures for the Section 6e metagene GSEA
#
#   Rscript plot_enrichment.R                 # all samples, all collections
#   PE_SAMPLE=P16 Rscript plot_enrichment.R
#
# Reads the CSVs written by Section 6e,
#   giotto_results/<SAMPLE>/enrichment/gsea_<collection>_<mode>.csv
# and draws the dotplots (no Giotto, fgsea or gobject needed).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace")
RES_DIR  <- Sys.getenv("PE_RES_DIR", file.path(BASE_DIR, "giotto_results"))

SAMPLES <- c("P16", "P01", "P09", "P22")
if (nzchar(Sys.getenv("PE_SAMPLE"))) SAMPLES <- Sys.getenv("PE_SAMPLE")

COLLECTIONS <- c(hallmark = "Hallmark", go_bp = "GO: Biological Process",
                 go_mf = "GO: Molecular Function", reactome = "Reactome")

# `excl_members`: a module's own member genes are dropped from the ranked
# vector. Use "depth_adjusted" for the sequencing-depth-controlled version.
MODE  <- Sys.getenv("PE_MODE", "excl_members")
N_TOP <- 5L        # top terms per module, by p-value
SIG   <- 0.05

# Diverging NES scale, same blue/red the rest of the project uses.
COL_LO <- "#2166AC"; COL_HI <- "#B2182B"

# -----------------------------------------------------------------------------
# ACRONYMS + pretty_term() (default wrap width 42), from this folder.
.script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "enrichment_labels.R"))

# -----------------------------------------------------------------------------
plot_one <- function(sample_id, coll_key) {
  enr_dir <- file.path(RES_DIR, sample_id, "enrichment")
  f <- file.path(enr_dir, sprintf("gsea_%s_%s.csv", coll_key, MODE))
  if (!file.exists(f)) {
    cat(sprintf("  -- %-12s %-9s : no %s\n", sample_id, coll_key, basename(f)))
    return(invisible(NULL))
  }
  g <- fread(f)
  if (!nrow(g)) {
    cat(sprintf("  -- %-12s %-9s : empty\n", sample_id, coll_key))
    return(invisible(NULL))
  }

  tp <- g[order(module, pval)][, head(.SD, N_TOP), by = module]
  tp[, label := pretty_term(term, coll_key)]
  tp[, significant := padj < SIG]

  # Order terms by the module where each is strongest, then by NES within that
  # module. This is what gives a dotplot its readable diagonal; ordering by NES
  # alone interleaves modules and the structure disappears.
  anchor <- tp[order(pval), .SD[1], by = label][, .(label, a_mod = module, a_nes = NES)]
  tp <- merge(tp, anchor, by = "label")
  lv <- unique(anchor[order(-a_mod, a_nes), label])
  tp[, label := factor(label, levels = lv)]

  n_terms  <- uniqueN(tp$label)
  n_mod    <- uniqueN(tp$module)
  # Symmetric limits so white is NES = 0 in every panel.
  nes_lim  <- max(abs(tp$NES), na.rm = TRUE)
  n_sig    <- tp[significant == TRUE, .N]

  # When every point clears padj < SIG the shape legend has one level and says
  # nothing the subtitle has not already said; drop it and keep the space.
  any_ns <- tp[significant == FALSE, .N] > 0

  p <- ggplot(tp, aes(x = factor(module), y = label)) +
    geom_point(aes(size = size, fill = NES, shape = significant),
               colour = "grey25", stroke = 0.45) +
    scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 1),
      breaks = c(TRUE, FALSE),
      labels = c(sprintf("padj < %.2f", SIG), "not significant"),
      name = NULL, guide = if (any_ns) "legend" else "none") +
    scale_fill_gradient2(low = COL_LO, mid = "white", high = COL_HI,
                         midpoint = 0, limits = c(-nes_lim, nes_lim),
                         name = "NES") +
    scale_size_continuous(range = c(1.6, 7), name = "gene set size") +
    guides(
      shape = if (any_ns) guide_legend(
                order = 1, override.aes = list(size = 3.4, fill = "grey70"))
              else "none",
      fill  = guide_colourbar(order = 2, barheight = grid::unit(62, "pt"),
                              barwidth = grid::unit(9, "pt")),
      size  = guide_legend(order = 3, override.aes = list(shape = 21, fill = "grey70"))) +
    labs(
      x = "metagene module", y = NULL,
      title = sprintf("Correlation-ranked GSEA · %s", COLLECTIONS[[coll_key]]),
      subtitle = sprintf(
        "%s  —  top %d terms per module (k = %d), %d of %d at padj < %.2f\nmodule member genes excluded from the ranking",
        sample_id, N_TOP, n_mod, n_sig, nrow(tp), SIG)) +
    theme_bw(base_size = 11) +
    theme(
      panel.border       = element_rect(colour = "grey35", linewidth = 0.5),
      panel.grid.major   = element_line(colour = "grey88", linetype = "dotted",
                                        linewidth = 0.35),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 9, lineheight = 0.92),
      axis.text.x        = element_text(size = 10),
      axis.title.x       = element_text(size = 10, margin = margin(t = 6)),
      axis.ticks         = element_line(colour = "grey45", linewidth = 0.35),
      plot.title         = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle      = element_text(size = 9, colour = "grey30",
                                        lineheight = 1.15,
                                        margin = margin(b = 8)),
      plot.title.position = "plot",
      legend.title       = element_text(size = 9),
      legend.text        = element_text(size = 8.5),
      legend.key.size    = grid::unit(11, "pt"),
      plot.margin        = margin(8, 10, 6, 6))

  # Canvas sized from its contents: width is the term labels, the dot panel and
  # the legend (the label column varies most between collections).
  extra_lines  <- sum(grepl("\n", levels(tp$label)))
  max_lab_char <- max(nchar(unlist(strsplit(levels(tp$label), "\n", fixed = TRUE))))
  lab_col <- 0.065 * max_lab_char            # ~9pt characters
  panel_w <- max(2.4, 0.52 * n_mod)
  w <- lab_col + panel_w + 1.55              # + legend column
  h <- 1.9 + 0.30 * n_terms + 0.10 * extra_lines

  name <- file.path(enr_dir, sprintf("gsea_dotplot_%s", coll_key))
  tryCatch({
    ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
    ggsave(paste0(name, ".svg"), p, width = w, height = h, device = "svg",
           bg = "white")
    cat(sprintf("  ok %-12s %-9s : %2d terms x %d modules, %2d sig  (%.1f x %.1f in)\n",
                sample_id, coll_key, n_terms, n_mod, n_sig, w, h))
  }, error = function(e) {
    # report a failed save immediately
    cat(sprintf("  !! %s %s FAILED: %s\n", sample_id, coll_key,
                conditionMessage(e)))
    warning(sample_id, " ", coll_key, ": ", conditionMessage(e))
  })
  invisible(NULL)
}

# -----------------------------------------------------------------------------
cat(sprintf("--- GSEA dotplots (mode = %s, top %d per module) ---\n", MODE, N_TOP))
for (s in SAMPLES) for (cn in names(COLLECTIONS)) plot_one(s, cn)
cat("\n========== enrichment figures DONE ==========\n\n")
