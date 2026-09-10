# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# composition_compare.R — cell-type composition, poor vs good outcome
#                         Regular Visium (55 um), 4-sample cohort
#
# One stacked bar per sample: the SpatialDWLS cell-type proportions, good-outcome
# samples on the left, poor-outcome on the right.
#
#   Rscript scripts/composition_compare.R
#
# Reads the per-sample Section 8c output. Giotto-free: data.table + ggplot2 only,
# no gobject load, no h5ad. Runs in seconds.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace")
# Overridable so the script can be pointed at a scratch copy for testing without
# touching giotto_results/. Unset in normal use.
RES_DIR  <- Sys.getenv("COMPOSITION_RES_DIR", file.path(BASE_DIR, "giotto_results"))
OUT_DIR  <- file.path(RES_DIR, "composition_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Good first so the good-outcome samples sit on the left of the figure.
SAMPLES     <- c(P01 = "Good", P09 = "Good",
                 P16 = "Poor", P22 = "Poor")
OUTCOME_LVL <- c("Good", "Poor")
stopifnot(all(SAMPLES %in% OUTCOME_LVL))

# The project cell-type palette, identical to CT_PALETTE in the per-sample
# scripts, so these bars match every other cell-type figure in the project.
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

# Print a percentage inside a stacked segment only when the segment is big
# enough to hold the text; below this the labels collide into noise.
LABEL_MIN_FRAC <- 0.05

# -----------------------------------------------------------------------------
# Load. Fail loudly and specifically: the usual reason a file is missing is that
# a sample has not been re-run since Section 8c was added, and that is worth
# saying outright rather than dying inside a merge 40 lines later.
cat("--- Reading per-sample composition ---\n")
comp <- rbindlist(lapply(names(SAMPLES), function(s) {
  f <- file.path(RES_DIR, s, "composition_summary.csv")
  if (!file.exists(f)) {
    stop(sprintf(paste0(
      "\n  Missing: %s\n",
      "  Sample %s has no Section 8c output yet. Run the pipeline for it first:\n",
      "    Rscript ../visium_55um/03_giotto_pipeline.R %s\n"), f, s, s))
  }
  dt <- fread(f)
  dt[, sample := s][, outcome := SAMPLES[[s]]]
  dt[, .(sample, outcome, cell_type, mean_prop_continuous)]
}), use.names = TRUE)

# `No majority` is not a cell type — it is the spots where DWLS returned all
# zeros (121 in P09) and so it carries no proportion. Dropping it and
# renormalising the real cell types below is what keeps the four bars directly
# comparable to one another.
n_drop <- comp[is.na(mean_prop_continuous), .N]
if (n_drop) {
  cat(sprintf("  dropped %d non-cell-type row(s): %s\n", n_drop,
              paste(unique(comp[is.na(mean_prop_continuous),
                                paste0(cell_type, " (", sample, ")")]),
                    collapse = ", ")))
  comp <- comp[!is.na(mean_prop_continuous)]
}

comp[, sample  := factor(sample,  levels = names(SAMPLES))]
comp[, outcome := factor(outcome, levels = OUTCOME_LVL)]
stopifnot(uniqueN(comp$sample) == length(SAMPLES))

unknown <- setdiff(unique(comp$cell_type), names(CT_PALETTE))
if (length(unknown)) stop("no colour for cell type(s): ",
                          paste(unknown, collapse = ", "),
                          " — update CT_PALETTE.")

# -----------------------------------------------------------------------------
# Normalise each bar to exactly 1. mean_prop_continuous already sums to 1.0 per
# sample (n_cell = 20 makes the shipped DWLS output a proportional composition),
# so this only absorbs rounding and the dropped non-cell-type row above. No other
# transform is applied — the bars are the DWLS proportions.
comp[, proportion := mean_prop_continuous / sum(mean_prop_continuous), by = sample]
stopifnot(all(abs(comp[, sum(proportion), by = sample]$V1 - 1) < 1e-9))

# Stack the abundant types nearest the axis; legend follows the same order.
ct_order <- comp[, .(m = mean(proportion)), by = cell_type][order(-m), cell_type]
comp[, cell_type := factor(cell_type, levels = ct_order)]
setorder(comp, sample, cell_type)

cat(sprintf("  %d samples x %d cell types\n",
            uniqueN(comp$sample), uniqueN(comp$cell_type)))
print(dcast(comp, cell_type ~ sample, value.var = "proportion")[
        , lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])

fwrite(comp[, .(sample, outcome, cell_type, proportion)],
       file.path(OUT_DIR, "composition_stacked_bar.csv"))
cat(sprintf("\n  wrote %s\n", file.path(OUT_DIR, "composition_stacked_bar.csv")))

# -----------------------------------------------------------------------------
# FIGURE BUILDER — one stacked bar per sample, faceted by outcome so the two
# groups sit in labelled panels with a real gap between them. Both figures below
# go through this, so their styling cannot drift apart.
stacked_bar <- function(dt, title, subtitle, label_min, label_fmt,
                        y_expand_top = 0.02) {
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
    scale_fill_manual(values = CT_PALETTE, name = NULL,
                      labels = function(x) gsub("_", " ", gsub("___", " / ", x))) +
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
      axis.text.x        = element_text(size = 10),
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
    cat(sprintf("  wrote %s.{png,svg}\n", name))
  }, error = function(e) {
    # cat as well as warning(): a deferred warning would let a stale figure from
    # a previous run sit on disk looking freshly generated.
    cat("  !! ", name, " FAILED: ", conditionMessage(e), "\n", sep = "")
    warning(name, " failed: ", conditionMessage(e))
  })
}

cat("\n--- Figures ---\n")

# FIGURE 1 — all cell types.
save_fig(
  stacked_bar(comp,
              title    = "Cell-type composition by outcome",
              subtitle = "SpatialDWLS proportions, normalised per sample",
              label_min = LABEL_MIN_FRAC, label_fmt = "%.0f%%"),
  "composition_stacked_bar", width = 8.5, height = 6)

# -----------------------------------------------------------------------------
# FIGURE 2 — the immune compartment only.
#
# Same three colours as Figure 1, and deliberately NOT renormalised to 100%:
# `proportion` stays the share of ALL cells, so a bar height here means the same
# thing it means in Figure 1 and the figure answers "how much immune infiltrate
# is in this tumour". Renormalising would make P16 (1.4% immune in total)
# look like P22 (11.7%).
IMMUNE_TYPES <- c("Myeloid_cells", "T_cells", "B_Plasma_cells")
IMMUNE_LABEL_MIN <- 0.01   # 1% of all cells; Figure 1's 5% would blank nearly all
stopifnot(all(IMMUNE_TYPES %in% levels(comp$cell_type)))

imm <- comp[cell_type %in% IMMUNE_TYPES]
imm[, cell_type := droplevels(cell_type)]
imm_tot <- imm[, .(total = sum(proportion)), by = .(sample, outcome)]

cat(sprintf("  immune total per sample: %s\n",
            paste(sprintf("%s %.2f%%", imm_tot$sample, 100 * imm_tot$total),
                  collapse = "  |  ")))

p_imm <- stacked_bar(
    imm,
    title    = "Immune-cell composition by outcome",
    subtitle = "SpatialDWLS proportions, as a share of ALL cells (not renormalised)",
    label_min = IMMUNE_LABEL_MIN, label_fmt = "%.1f%%",
    y_expand_top = 0.12) +
  # The total is the number this figure exists to communicate; the stacked
  # segments alone make it hard to add up by eye.
  geom_text(data = imm_tot, inherit.aes = FALSE,
            aes(x = sample, y = total, label = sprintf("%.1f%%", 100 * total)),
            vjust = -0.6, size = 3.6, fontface = "bold", colour = "grey20")

save_fig(p_imm, "composition_immune_stacked_bar", width = 7.5, height = 5.2)

cat("\n========== composition comparison DONE ==========\n")
cat("Outputs in:", OUT_DIR, "\n\n")
