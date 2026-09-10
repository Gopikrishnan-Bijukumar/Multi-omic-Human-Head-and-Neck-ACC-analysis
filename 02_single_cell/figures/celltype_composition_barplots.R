# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# processed_celltype_composition_barplots.R — scRNA cell-type composition
#
# Three stacked-bar figures from inprogress_3.h5ad:
#   1. one bar per sample (batch), annot_8_c2l labels, good-outcome samples left
#   2. one bar per outcome group, annot_8_c2l labels, all cells pooled
#   3. one bar per outcome group, annot_5 collapsed to tumour / stromal / immune
#
#   Rscript named_scripts/processed_celltype_composition_barplots.R
#
# The palette and the plot styling are lifted from the spatial figure
#   spatial/reg_visium/giotto_w_cytospace/scripts/composition_compare.R
# so a cell type carries the SAME colour in the scRNA and Visium figures.
#
# The h5ad is opened read-only and only obs/ is read — the count matrix is
# never touched, so this runs in seconds on the 8 GB file and cannot alter it.
#
# NOTE on the reader: every obs/ dataset in this h5ad is LZF-compressed (h5py's
# own filter). hdf5r and rhdf5 cannot decode LZF — it needs an external HDF5
# filter plugin, which is not installed on this box — so the obs columns
# are pulled through reticulate + h5py. That is the only Python in the file;
# all counting, plotting and writing is R.
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(data.table)
  library(ggplot2)
})

H5AD    <- file.path(ACC_DATA_ROOT, "outputs/h5ad_files/inprogress_3.h5ad")
OUT_DIR <- file.path(ACC_DATA_ROOT, "outputs/images")

# The env the project's notebooks use ("spatial" kernel); it has h5py.
PY_BIN <- Sys.getenv("COMPOSITION_PYTHON",
                     "python")

CELLTYPE_KEY <- "annot_8_c2l"
BATCH_KEY    <- "batch"
OUTCOME_KEY  <- "clinical_outcome"

# Figure 3 only. annot_5 is the annotation of record; annot_7_c2l is its coarse
# form and is read purely to cross-check the hand-written compartment map.
ANNOT5_KEY   <- "annot_5"
COARSE_KEY   <- "annot_7_c2l"

OUTCOME_LVL <- c("Good", "Poor")   # Good first, so it sits on the left

# The project cell-type palette. Copied verbatim from CT_PALETTE in
# giotto_w_cytospace/scripts/composition_compare.R, which is the source of
# truth — do not recolour here, edit it there and copy across.
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

# -----------------------------------------------------------------------------
# FIGURE 3 constants — the three-compartment collapse of annot_5.
#
# Fixed segment order rather than the abundance sort used for Figures 1-2, so
# the three bands stay in the same place whatever the data does.
COMPARTMENT_LVL <- c("Tumour", "Stromal / structural", "Immune")

# One colour per compartment, taken from the families already in CT_PALETTE
# (tumour warm, stroma green, immune blue) so Figure 3 reads as part of the set.
COMPARTMENT_PALETTE <- c(
  "Tumour"               = "#C1272D",
  "Stromal / structural" = "#2E8B57",
  "Immune"               = "#1F6FB4"
)

# Every annot_5 label -> its compartment. Keyed on the raw label, NOT the
# palette_key() form. Two judgement calls are baked in here:
#
#  * "Actively Dividing cells" carries no tumour tag at annot_5, but annot_6
#    labels the same cells "ADC - Tumor" and annot_7_c2l puts them in
#    "Tumor cells". They are tumour.
#  * The non-malignant epithelium (Basal/EA/Ductal/Glandular/Secretory/Ciliated)
#    is neither tumour nor immune. annot_7_c2l keeps it as its own class; here it
#    joins the stromal/structural compartment, which the subtitle states.
COMPARTMENT_MAP <- c(
  # tumour
  "Myoepithelial cells - Tumor"                   = "Tumour",
  "Epithelial cells - Tumor"                      = "Tumour",
  "Actively Dividing cells"                       = "Tumour",
  "Actively Dividing Myoepithelial cells - Tumor" = "Tumour",
  "Epithelial cells - Basal - Tumor"              = "Tumour",
  "Basal cells - Tumor"                           = "Tumour",
  # stroma proper
  "Fibroblast cells"                              = "Stromal / structural",
  "Fibroblast cells - CAF"                        = "Stromal / structural",
  "Fibroblast cells - Inflammatory"               = "Stromal / structural",
  "Endothelial cells"                             = "Stromal / structural",
  "Mural cells"                                   = "Stromal / structural",
  "Muscle cells - Tongue"                         = "Stromal / structural",
  # non-malignant epithelium, grouped with the structural compartment
  "Epithelial cells - Basal"                      = "Stromal / structural",
  "Epithelial cells - EA"                         = "Stromal / structural",
  "Epithelial cells - Ductal"                     = "Stromal / structural",
  "Epithelial cells - Glandular"                  = "Stromal / structural",
  "Epithelial cells - Secretory"                  = "Stromal / structural",
  "Epithelial cells - Ciliated"                   = "Stromal / structural",
  # immune
  "Macrophage - M2"                               = "Immune",
  "Macrophage - AP/TAM"                           = "Immune",
  "Mast cells"                                    = "Immune",
  "CD4+ T cells"                                  = "Immune",
  "CD8+ T cells"                                  = "Immune",
  "Memory B cells"                                = "Immune",
  "B cells"                                       = "Immune",
  "Plasma cells"                                  = "Immune"
)

# The same collapse expressed on annot_7_c2l. Used only to verify COMPARTMENT_MAP
# — see the cross-check in the Figure 3 section.
COARSE_MAP <- c(
  "Tumor cells"        = "Tumour",
  "Fibroblast cells"   = "Stromal / structural",
  "Mural/muscle cells" = "Stromal / structural",
  "Endothelial cells"  = "Stromal / structural",
  "Epithelial cells"   = "Stromal / structural",
  "Myeloid cells"      = "Immune",
  "T cells"            = "Immune",
  "B / Plasma cells"   = "Immune"
)

# Print a percentage inside a stacked segment only when the segment is big
# enough to hold the text; below this the labels collide into noise.
LABEL_MIN_FRAC <- 0.05

# -----------------------------------------------------------------------------
# Load. Mode 'r' is read-only — never 'r+' or 'a' — and only obs/ is addressed,
# so X is not so much as opened. h5py stores a categorical as categories plus
# 0-based codes; both come back raw and the factor is rebuilt on the R side.
use_python(PY_BIN, required = TRUE)
py_run_string("
import h5py

def _read_obs_categoricals(path, keys):
    out = {}
    with h5py.File(path, 'r') as f:
        for k in keys:
            if k not in f['obs']:
                raise KeyError(\"obs has no column '%s'\" % k)
            g = f['obs'][k]
            cats = [c.decode() if isinstance(c, bytes) else str(c)
                    for c in g['categories'][:]]
            codes = g['codes'][:]
            if (codes < 0).any():
                raise ValueError(\"obs['%s'] has unassigned (-1) codes\" % k)
            out[k] = {'categories': cats, 'codes': codes.astype('int32')}
    return out
")

read_obs_factor <- function(raw, key) {
  cats <- unlist(raw[[key]]$categories, use.names = FALSE)
  factor(cats[raw[[key]]$codes + 1L], levels = cats)
}

# Cell-type label -> CT_PALETTE key: " - " marks a tumour subtype (___), and
# every remaining space or slash becomes a single underscore.
palette_key <- function(x) gsub("[ /]", "_", gsub(" - ", "___", x))

cat("--- Reading obs from the h5ad (read-only) ---\n")
raw <- py$`_read_obs_categoricals`(H5AD, list(BATCH_KEY, OUTCOME_KEY, CELLTYPE_KEY))

obs <- data.table(
  sample    = read_obs_factor(raw, BATCH_KEY),
  outcome   = read_obs_factor(raw, OUTCOME_KEY),
  cell_type = read_obs_factor(raw, CELLTYPE_KEY)
)
rm(raw)

stopifnot(all(levels(obs$outcome) %in% OUTCOME_LVL))
obs[, outcome := factor(as.character(outcome), levels = OUTCOME_LVL)]
obs[, cell_type := factor(palette_key(as.character(cell_type)),
                          levels = palette_key(levels(cell_type)))]

unknown <- setdiff(levels(obs$cell_type), names(CT_PALETTE))
if (length(unknown)) stop("no colour for cell type(s): ",
                          paste(unknown, collapse = ", "),
                          " — update CT_PALETTE.")

cat(sprintf("  %d cells | %d samples | %d cell types (%s)\n",
            nrow(obs), uniqueN(obs$sample), uniqueN(obs$cell_type), CELLTYPE_KEY))

# One outcome per sample; fail loudly if that ever stops holding.
smp <- obs[, .(outcome = unique(outcome)), by = sample]
if (anyDuplicated(smp$sample)) stop("sample(s) map to more than one ", OUTCOME_KEY)
# Good samples first, then Poor, each keeping the h5ad's category order.
setorder(smp, outcome, sample)
obs[, sample := factor(as.character(sample), levels = as.character(smp$sample))]
cat(sprintf("  %s\n", paste(sprintf("%s n=%d", OUTCOME_LVL,
                                    tabulate(smp$outcome)), collapse = "  |  ")))

# -----------------------------------------------------------------------------
# Counts -> per-bar proportions. Bars are normalised because cell yield varies
# 3.5x across samples (1,052-7,863), so raw counts would mostly show library
# size rather than composition.
# Absent sample x cell-type combinations are simply dropped: an empty segment
# draws nothing either way, and every bar is normalised over the cells it has.
composition <- function(dt) {
  out <- dt[, .N, by = .(sample, outcome, cell_type)]
  out[, proportion := N / sum(N), by = sample]
  stopifnot(all(abs(out[, sum(proportion), by = sample]$V1 - 1) < 1e-9))
  out[]
}

per_sample <- composition(obs)
stopifnot(sum(per_sample$N) == nrow(obs))

# Pool every cell of a group into a single bar by making the group the "sample".
pooled <- composition(copy(obs)[, sample := outcome])

# Stack the abundant types nearest the axis; legend follows the same order.
# Colours are pinned by name in scale_fill_manual, so ordering never moves them.
ct_order <- per_sample[, .(m = mean(proportion)), by = cell_type][order(-m), cell_type]
for (d in list(per_sample, pooled)) d[, cell_type := factor(cell_type, levels = ct_order)]
setorder(per_sample, sample, cell_type)
setorder(pooled, sample, cell_type)

cat("\n  pooled composition (%):\n")
print(dcast(pooled, cell_type ~ sample, value.var = "proportion", fill = 0)[
        , lapply(.SD, function(x) if (is.numeric(x)) round(100 * x, 2) else x)])

# -----------------------------------------------------------------------------
# FIGURE BUILDER — one stacked bar per sample, faceted by outcome so the two
# groups sit in labelled panels with a real gap between them. Styling is the
# same as composition_compare.R::stacked_bar so the scRNA and spatial figures
# read as one pair.
stacked_bar <- function(dt, title, subtitle, label_min, label_fmt,
                        x_angle = 0, y_expand_top = 0.02,
                        palette = CT_PALETTE,
                        legend_labels = function(x)
                          gsub("_", " ", gsub("___", " / ", x))) {
  ggplot(dt, aes(x = sample, y = proportion, fill = cell_type)) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
    # Blank the small labels rather than filtering the rows: position_stack
    # recomputes the stack from whatever data it is given, so passing a filtered
    # subset here silently shifts every label off its own segment.
    geom_text(
      aes(label = fifelse(proportion >= label_min,
                          sprintf(label_fmt, 100 * proportion), "")),
      position = position_stack(vjust = 0.5),
      size = 3.1, colour = "white", fontface = "bold") +
    facet_grid(~ outcome, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = palette, name = NULL, labels = legend_labels) +
    scale_y_continuous(labels = function(x) paste0(signif(100 * x, 3), "%"),
                       expand = expansion(mult = c(0, y_expand_top))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "share of cells") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", colour = "grey80"),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      strip.text         = element_text(face = "bold", size = 12),
      axis.text.x        = element_text(size = 10, angle = x_angle,
                                        hjust = if (x_angle) 1 else 0.5),
      plot.title         = element_text(face = "bold"),
      legend.key.size    = unit(11, "pt"),
      legend.text        = element_text(size = 9.5))
}

save_fig <- function(p, name, width, height) {
  tryCatch({
    ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
           width = width, height = height, dpi = 600, bg = "white")
    ggsave(file.path(OUT_DIR, paste0(name, ".svg")), p,
           width = width, height = height, device = "svg", bg = "white")
    cat(sprintf("  wrote %s.{png,svg}\n", file.path(OUT_DIR, name)))
  }, error = function(e) {
    # cat as well as warning(): a deferred warning would let a stale figure from
    # a previous run sit on disk looking freshly generated.
    cat("  !! ", name, " FAILED: ", conditionMessage(e), "\n", sep = "")
    warning(name, " failed: ", conditionMessage(e))
  })
}

cat("\n--- Figures ---\n")

# FIGURE 1 — every sample. 24 bars, so the sample labels are turned 45 degrees.
save_fig(
  stacked_bar(per_sample,
              title    = "Cell-type composition by sample",
              subtitle = sprintf("scRNA-seq, %s labels, normalised per sample",
                                 CELLTYPE_KEY),
              label_min = LABEL_MIN_FRAC, label_fmt = "%.0f%%",
              x_angle = 45),
  "processed_celltype_composition_by_batch", width = 16, height = 6.5)

# FIGURE 2 — the two outcome groups, every cell pooled. A lower label threshold
# than Figure 1: with only two bars the small compartments have room to be read.
save_fig(
  stacked_bar(pooled,
              title    = "Cell-type composition by clinical outcome",
              subtitle = sprintf("scRNA-seq, %s labels, all cells pooled per group",
                                 CELLTYPE_KEY),
              label_min = 0.03, label_fmt = "%.1f%%"),
  "processed_celltype_composition_by_outcome", width = 6, height = 6)

# -----------------------------------------------------------------------------
# FIGURE 3 — tumour / stromal / immune.
#
# The 26 annot_5 labels collapsed to three compartments via COMPARTMENT_MAP. The
# map is written by hand (annot_5's tumour tagging is incomplete — see the
# comment on COMPARTMENT_MAP), so before it is used it is checked against
# annot_7_c2l, the coarse annotation that already carries these classes. The two
# must agree for every cell; if an upstream relabel ever breaks that, the script
# stops here rather than quietly drawing a wrong figure.
cat("\n--- Compartments (", ANNOT5_KEY, " -> tumour / stromal / immune) ---\n", sep = "")

# Same file and same row order as the first read, so these bind straight on.
raw <- py$`_read_obs_categoricals`(H5AD, list(ANNOT5_KEY, COARSE_KEY))
obs[, annot_5 := read_obs_factor(raw, ANNOT5_KEY)]
obs[, coarse  := read_obs_factor(raw, COARSE_KEY)]
rm(raw)

unmapped <- setdiff(levels(obs$annot_5), names(COMPARTMENT_MAP))
if (length(unmapped)) stop("no compartment for ", ANNOT5_KEY, " label(s): ",
                           paste(unmapped, collapse = ", "),
                           " — update COMPARTMENT_MAP.")

obs[, compartment := factor(COMPARTMENT_MAP[as.character(annot_5)],
                            levels = COMPARTMENT_LVL)]
stopifnot(!anyNA(obs$compartment))

# Cross-check against the coarse annotation.
unmapped <- setdiff(levels(obs$coarse), names(COARSE_MAP))
if (length(unmapped)) stop("no compartment for ", COARSE_KEY, " label(s): ",
                           paste(unmapped, collapse = ", "),
                           " — update COARSE_MAP.")
coarse_cmp <- factor(COARSE_MAP[as.character(obs$coarse)], levels = COMPARTMENT_LVL)
n_disagree <- sum(coarse_cmp != obs$compartment)
if (n_disagree) {
  print(table(annot_5_derived = obs$compartment, coarse_derived = coarse_cmp))
  stop(sprintf("COMPARTMENT_MAP disagrees with %s for %d of %d cells — see the table above",
               COARSE_KEY, n_disagree, nrow(obs)))
}
cat(sprintf("  cross-check OK: agrees with %s for all %d cells\n", COARSE_KEY, nrow(obs)))

n_cmp <- tabulate(obs$compartment, nbins = length(COMPARTMENT_LVL))
cat(sprintf("  %s\n", paste(sprintf("%s n=%d (%.1f%%)", COMPARTMENT_LVL,
                                    n_cmp, 100 * n_cmp / nrow(obs)),
                            collapse = "  |  ")))

# Pool per outcome group, reusing composition(): the compartment takes the place
# of cell_type, the outcome takes the place of sample.
pooled_cmp <- composition(copy(obs)[, `:=`(sample = outcome, cell_type = compartment)])
stopifnot(sum(pooled_cmp$N) == nrow(obs))
pooled_cmp[, cell_type := factor(as.character(cell_type), levels = COMPARTMENT_LVL)]
setorder(pooled_cmp, sample, cell_type)

cat("\n  pooled compartment composition (%):\n")
print(dcast(pooled_cmp, cell_type ~ sample, value.var = "proportion", fill = 0)[
        , lapply(.SD, function(x) if (is.numeric(x)) round(100 * x, 2) else x)])

# The compartment names are already display-ready, so the legend labels pass
# through unchanged instead of going through the palette_key() un-mangling.
save_fig(
  stacked_bar(pooled_cmp,
              title    = "Compartment composition by clinical outcome",
              subtitle = sprintf(paste0("scRNA-seq, %s collapsed to tumour / stromal / immune\n",
                                        "stromal includes non-malignant epithelium"),
                                 ANNOT5_KEY),
              label_min = 0.03, label_fmt = "%.1f%%",
              palette = COMPARTMENT_PALETTE, legend_labels = identity),
  "processed_compartment_composition_by_outcome", width = 6, height = 6)

cat("\n========== composition bar plots DONE ==========\n")
cat("Outputs in:", OUT_DIR, "\n\n")
