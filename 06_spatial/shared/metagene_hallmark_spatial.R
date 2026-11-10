# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# metagene_hallmark_spatial.R — metagene territories on the tissue, with each
# territory's top enriched programs beside it in the same color. One figure per
# gene-set collection: hallmark, GO biological process, Reactome.
#
#   Rscript metagene_hallmark_spatial.R                    # 4 samples x 3 collections
#   MHS_SAMPLE=P16 Rscript metagene_hallmark_spatial.R
#   MHS_COLL=go_bp,reactome Rscript metagene_hallmark_spatial.R
#   MHS_MODE=excl_members  Rscript metagene_hallmark_spatial.R
#
# The figure shows the greyscale H&E with every metagene territory in its own
# color, and that territory's top enriched terms in the same color on the right.
#
# Reads the saved gobject and the Section 6e CSVs only (k is read from the
# object). The gobject is loaded once per sample and all three figures are drawn
# from it.
#
# INPUTS  giotto_results/<S>/giotto_object/                      (loadGiotto)
#         giotto_results/<S>/enrichment/gsea_<collection>_{excl_members,
#                                                          depth_adjusted}.csv
#         giotto_results/<S>/enrichment/metagene_depth_confound.csv
#         data/sample_<S>/spatial/tissue_hires_image.png + scalefactors_json.json
# OUTPUTS giotto_results/<S>/tissue_plots/metagene_domains_<collection>.{png,svg}
#         giotto_results/<S>/enrichment/metagene_domain_top_<collection>.csv
#
# In the SVG, all text is live <text> (see the ggsave call at the bottom) and only
# the spot lattice is rasterized.
# =============================================================================

suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(png)
  library(svglite)
  library(ggrastr)
})

BASE_DIR <- file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace")
DATA_DIR <- file.path(BASE_DIR, "data")
RES_DIR  <- file.path(BASE_DIR, "giotto_results")

.script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "enrichment_labels.R"))

SAMPLES <- c("P16", "P01", "P09", "P22")
if (nzchar(Sys.getenv("MHS_SAMPLE"))) SAMPLES <- Sys.getenv("MHS_SAMPLE")

# Collections to draw. Keys match the Section 6e CSV names
# (gsea_<key>_<mode>.csv, as in plot_enrichment.R); values are the display
# string used in the title, footnote and empty-block note.
COLLECTIONS <- c(hallmark = "hallmark",
                 go_bp    = "GO biological-process",
                 reactome = "Reactome")
if (nzchar(Sys.getenv("MHS_COLL"))) {
  want <- trimws(strsplit(Sys.getenv("MHS_COLL"), ",", fixed = TRUE)[[1]])
  bad  <- setdiff(want, names(COLLECTIONS))
  if (length(bad)) stop("MHS_COLL: unknown collection(s) ", paste(bad, collapse = ", "),
                        " — known: ", paste(names(COLLECTIONS), collapse = ", "))
  COLLECTIONS <- COLLECTIONS[want]
}

# Same map as ../visium_55um/03_giotto_pipeline.R Section 1.
OUTCOME_MAP <- c(P16 = "Poor", P01 = "Good",
                 P09 = "Good", P22 = "Poor")

N_TOP     <- 3L      # terms shown per metagene
SIG       <- 0.05
DEPTH_CUT <- 0.5     # |r| with total counts above which a module is depth-driven

# Spot footprint as a fraction of the measured spot pitch: 0.88 draws a
# near-touching lattice. The page-space point_size is solved per figure (see the
# tissue panel).
SPOT_FILL    <- 0.88
# Opacity of the greyscale H&E over the white page, the same constant as
# Section 12 of ../visium_55um/03_giotto_pipeline.R.
TISSUE_ALPHA <- 0.75

# ggplot's geom_text `size` is in millimetres; svglite writes `size * .pt` as the
# SVG font-size, and Illustrator reads that number as points. So 9 pt type is
# 9/.pt, not 9 — verified off the rendered SVG (`font-size: 9.00px`).
PT9   <- 9  / .pt    # 3.1631 — terms, NES/padj, spot counts
PT11  <- 11 / .pt    # module header
PT8   <- 8  / .pt    # italic notes and the footnote

# -----------------------------------------------------------------------------
# PALETTE
#
# Fixed by module index, so a module keeps its color between the PNG, the SVG
# and the CSV.
#
# Checked for all pairs: worst all-pairs CVD dE 10.9 (deutan), worst
# normal-vision dE 22.1. Territories also carry numbered tags.
#
# SPOT_PAL is the mark color. TEXT_PAL is the same hue darkened for >= 4.5:1
# contrast on white (4.52-5.71:1). The header swatch carries the exact SPOT_PAL
# value.
SPOT_PAL <- c("#da751a", "#064af6", "#01aaa0", "#b41e72", "#9f7ce6", "#5a6700")
TEXT_PAL <- c("#b95d00", "#1552fe", "#09847c", "#ce3b87", "#8663ca", "#6e7c21")

INK_DIM  <- "#5f5e5b"   # NES / padj columns, italic notes
INK_MID  <- "#52514e"   # subtitle, footnote
INK_HI   <- "#0b0b0b"   # title

# =============================================================================
# HELPERS
# =============================================================================

# Greyscale tissue background composited over the white page with Section 12's
# formula and constant (../visium_55um/03_giotto_pipeline.R):
#
#     lum  = 0.2989R + 0.5870G + 0.1140B
#     grey = 1 - TISSUE_ALPHA * (1 - lum)
#
# Returns the greyscale matrix and its extent; the background is drawn with
# annotation_raster (Giotto's large-image path renormalizes a single-channel
# image to its own intensity range). The extent formula is the one Giotto uses
# for negative_y = TRUE (0, ncol*sf, -nrow*sf, 0). Nothing is written to disk.
build_grey_image <- function(sample_id) {
  vis_dir <- file.path(DATA_DIR, paste0("sample_", sample_id), "spatial")
  src     <- file.path(vis_dir, "tissue_hires_image.png")

  px  <- png::readPNG(src)
  lum <- if (length(dim(px)) == 3 && dim(px)[3] >= 3) {
    0.2989 * px[, , 1] + 0.5870 * px[, , 2] + 0.1140 * px[, , 3]
  } else px
  grey <- 1 - TISSUE_ALPHA * (1 - lum)

  sf <- 1 / jsonlite::fromJSON(
    file.path(vis_dir, "scalefactors_json.json"))$tissue_hires_scalef
  list(raster = grDevices::as.raster(grey),
       xmin = 0, xmax = ncol(grey) * sf,
       ymin = -nrow(grey) * sf, ymax = 0)
}

# Where to place a territory's numbered tag: bin the domain's spots on a coarse
# grid, take the densest bin, and snap to the nearest real spot, so the tag sits
# on tissue belonging to the domain.
densest_spot <- function(x, y, bin) {
  bx <- floor(x / bin); by <- floor(y / bin)
  key <- paste(bx, by, sep = "_")
  tab <- table(key)
  hit <- names(tab)[which.max(tab)]
  sel <- key == hit
  cx  <- mean(x[sel]); cy <- mean(y[sel])
  i   <- which.min((x - cx)^2 + (y - cy)^2)
  c(x[i], y[i])
}

# Per-module GSEA source. Section 6e already writes the verdict into
# metagene_module_summary.csv ("DEPTH-DRIVEN ... read gsea_depth_adjusted_top3");
# this applies it. A module whose score correlates |r| > 0.5 with total counts
# has an unadjusted ranking that partly reports sequencing depth, so it reads
# the depth-controlled file and the figure says so on the header line.
# MHS_MODE forces one file for every module.
pick_modes <- function(enr_dir, modules) {
  forced <- Sys.getenv("MHS_MODE")
  conf_f <- file.path(enr_dir, "metagene_depth_confound.csv")
  r <- setNames(rep(NA_real_, length(modules)), modules)
  if (file.exists(conf_f)) {
    cf <- fread(conf_f)
    cf[, k := as.integer(sub("^metagene_", "", metagene))]
    r[as.character(cf$k)] <- cf$r_with_total_expr
  }
  mode <- ifelse(!is.na(r) & abs(r) > DEPTH_CUT, "depth_adjusted", "excl_members")
  if (nzchar(forced)) mode[] <- forced
  data.table(module = modules, r_with_total_expr = as.numeric(r),
             mode_used = as.character(mode))
}

# =============================================================================
# THE TERMS PANEL
#
# Laid out on a fixed grid of line slots (one slot = LINE_IN inches), so figures
# with different module counts have the same line spacing.
# =============================================================================
LINE_IN    <- 0.205   # one line slot, inches
HEAD_IN    <- 1.35    # title + subtitle band, inches
FOOT_WRAP  <- 78      # footnote wrap column, characters

# Slots reserved at the bottom for the wrapped footnote, from the number of lines
# it wraps to. A footnote line costs ~0.67 slots at the rendered slot height;
# 0.75 per line plus 0.6 for the gap above.
foot_slots <- function(footnote) {
  0.75 * length(strsplit(footnote, "\n", fixed = TRUE)[[1]]) + 0.6
}

# Statistics sit in a left-aligned column starting at STAT_X, just after the
# term text. All three are panel-width fractions. At 9 pt Liberation Sans the
# widest wrapped term line (TERM_WRAP = 34) ends at x = 0.488 and the widest
# statistic is 1.806 in, so STAT_X = 0.58 leaves a 0.44 in gutter and keeps the
# longest statistic inside the panel, for all twelve sample x collection blocks.
STAT_X    <- 0.58
BLOCK_R   <- 0.96     # right edge of the header's spot count, = widest stat end
TERM_WRAP <- 34       # characters, passed to pretty_term()

# Tissue : terms width ratio in the patchwork split (also used for the spot-size
# derivation).
TISSUE_SHARE <- 1.7

terms_panel <- function(blocks, n_slots, sample_id, footnote, FOOT_SLOTS,
                        coll_name) {
  rows <- list(); slot <- 0

  push <- function(x, lab, col, size, face = "plain", hjust = 0, dy = 0) {
    slot_now <- slot + dy
    rows[[length(rows) + 1]] <<- data.table(
      x = x, y = -slot_now, label = lab, colour = col,
      size = size, face = face, hjust = hjust)
  }

  swatches <- list()

  for (b in blocks) {
    slot <- slot + 1
    swatches[[length(swatches) + 1]] <-
      data.table(x = 0.018, y = -slot, fill = b$fill)
    # Header at 11 pt, above the 9 pt statistics.
    push(0.055, sprintf("METAGENE %d", b$module), b$ink, PT11, "bold")
    push(BLOCK_R, sprintf("%s spots · %.1f%%",
                          format(b$n_spots, big.mark = ","), b$pct),
         INK_DIM, PT9, hjust = 1)

    if (!nrow(b$terms)) {
      slot <- slot + 1
      push(0.075, sprintf("no positively enriched %s term", coll_name),
           INK_DIM, PT9, "italic")
    } else {
      for (i in seq_len(nrow(b$terms))) {
        lines <- strsplit(b$terms$label[i], "\n", fixed = TRUE)[[1]]
        for (j in seq_along(lines)) {
          slot <- slot + 1
          push(0.075, lines[j], b$ink, PT9)
          if (j == 1L) {
            # "n.s." is not decoration. A module can have fewer than three
            # significant terms (P01 metagene_2 has one under
            # hallmark), and without the marker the block reads as three
            # findings.
            push(STAT_X, sprintf("NES %+.2f   padj %.1e%s",
                                 b$terms$NES[i], b$terms$padj[i],
                                 if (b$terms$padj[i] >= SIG) "   n.s." else ""),
                 INK_DIM, PT9, hjust = 0)
          }
        }
      }
    }
    if (b$depth_adj) {
      slot <- slot + 1
      push(0.075, "ranking controlled for sequencing depth", INK_DIM, PT8,
           "italic")
    }
    slot <- slot + 0.9   # gap before the next module
  }

  txt <- rbindlist(rows)
  swt <- rbindlist(swatches)

  # Footnote pinned to the bottom slots, not to the content, so it lands in the
  # same place in all four figures. Top-anchored (vjust = 1) so it grows down
  # into the reserved band however many lines it wraps to.
  foot <- data.table(x = 0, y = -(n_slots - FOOT_SLOTS + 0.4), label = footnote,
                     colour = INK_MID, size = PT8, face = "italic", hjust = 0)

  ggplot() +
    geom_point(data = swt, aes(x = x, y = y), shape = 22, size = 2.9,
               fill = swt$fill, colour = swt$fill, stroke = 0) +
    geom_text(data = txt, aes(x = x, y = y, label = label, hjust = hjust),
              colour = txt$colour, size = txt$size, fontface = txt$face,
              lineheight = 0.95) +
    geom_text(data = foot, aes(x = x, y = y, label = label, hjust = hjust),
              colour = foot$colour, size = foot$size, fontface = foot$face,
              vjust = 1, lineheight = 1.15) +
    scale_x_continuous(limits = c(-0.01, 1.02), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-n_slots, 0.4), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.margin      = margin(6, 6, 4, 10))
}

# =============================================================================
# ONE SAMPLE — the parts that do not depend on the collection
#
# loadGiotto, the z-score argmax, the nearest-neighbor pitch measurement and the
# H&E composite run once per sample. Returns a context that draw_collection()
# consumes.
# =============================================================================
sample_context <- function(sample_id) {
  out_dir <- file.path(RES_DIR, sample_id)
  enr_dir <- file.path(out_dir, "enrichment")
  tp_dir  <- file.path(out_dir, "tissue_plots")
  dir.create(tp_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("--- %s ---\n", sample_id))

  gobj <- loadGiotto(file.path(out_dir, "giotto_object"), verbose = FALSE)
  md   <- pDataDT(gobj)

  mg_cols <- grep("^metagene_[0-9]+$", colnames(md), value = TRUE)
  mg_cols <- mg_cols[order(as.integer(sub("^metagene_", "", mg_cols)))]
  if (!length(mg_cols)) {
    cat("  !! no metagene_* columns in the saved object — skipped\n"); return(invisible(NULL))
  }
  modules <- as.integer(sub("^metagene_", "", mg_cols))
  K <- length(modules)
  if (K > length(SPOT_PAL)) {
    stop(sample_id, ": ", K, " modules but the validated palette has ",
         length(SPOT_PAL), " slots — extend and re-validate it, do not cycle.")
  }

  # --- domain assignment: argmax of z-scored scores --------------------------
  # Raw metagene scores are means over gene sets of different size and different
  # expression level, so their magnitudes are not comparable and a raw argmax
  # would just pick whichever module happens to sit highest everywhere.
  scores <- as.matrix(md[, ..mg_cols])
  z <- apply(scores, 2, function(v) {
    s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  })
  dom <- modules[max.col(z, ties.method = "first")]

  n_spots <- as.integer(table(factor(dom, levels = modules)))
  cat(sprintf("  %d spots, %d modules; sizes: %s\n", nrow(md), K,
              paste(sprintf("m%d=%d", modules, n_spots), collapse = " ")))

  # --- spatial layout --------------------------------------------------------
  sl <- getSpatialLocations(gobj, output = "data.table")
  sl <- merge(sl, data.table(cell_ID = md$cell_ID, dom = dom), by = "cell_ID")
  n_cells <- nrow(md)

  # The gobject is not needed below; drop it before drawing the figures.
  rm(gobj, md, scores, z); invisible(gc(FALSE))

  # Spot pitch, measured: coordinates are full-res pixels and the um-per-pixel
  # scale differs per slide. Only the query points are subsampled; neighbors use
  # all spots.
  set.seed(1)
  q     <- sl[sample(.N, min(.N, 400L))]
  pitch <- stats::median(vapply(seq_len(nrow(q)), function(i) {
    dd <- (sl$sdimx - q$sdimx[i])^2 + (sl$sdimy - q$sdimy[i])^2
    sqrt(min(dd[dd > 0]))
  }, numeric(1)))

  # Crop to the spots rather than to the slide (the capture area is larger than the
  # section).
  padx <- 0.03 * diff(range(sl$sdimx)); pady <- 0.03 * diff(range(sl$sdimy))
  xlim <- range(sl$sdimx) + c(-padx, padx)
  ylim <- range(sl$sdimy) + c(-pady, pady)

  dk <- build_grey_image(sample_id)

  # Numbered tags. Placed once: they mark territories, and the territories are
  # the same in all three figures — a tag that moved between collections would
  # imply the domain map had moved with it.
  tags <- rbindlist(lapply(modules, function(k) {
    s <- sl[dom == k]
    if (!nrow(s)) return(NULL)
    p <- densest_spot(s$sdimx, s$sdimy, bin = 8 * pitch)
    data.table(module = k, x = p[1], y = p[2],
               fill = SPOT_PAL[k], ink = TEXT_PAL[k])
  }))

  sl[, dom_f := factor(dom, levels = modules)]

  list(sample_id = sample_id, enr_dir = enr_dir, tp_dir = tp_dir,
       modules = modules, K = K, n_spots = n_spots, n_cells = n_cells,
       sl = sl, tags = tags, pitch = pitch, xlim = xlim, ylim = ylim, dk = dk,
       modes = pick_modes(enr_dir, modules))
}

# =============================================================================
# ONE FIGURE — one sample, one collection
# =============================================================================
draw_collection <- function(ctx, coll_key) {
  coll_name <- COLLECTIONS[[coll_key]]
  enr_dir   <- ctx$enr_dir
  modules   <- ctx$modules
  n_spots   <- ctx$n_spots
  n_cells   <- ctx$n_cells
  K         <- ctx$K
  modes     <- ctx$modes

  # --- terms -----------------------------------------------------------------
  cache <- list()
  read_mode <- function(m) {
    if (is.null(cache[[m]])) {
      f <- file.path(enr_dir, sprintf("gsea_%s_%s.csv", coll_key, m))
      if (!file.exists(f))
        cat(sprintf("  -- %s: no %s\n", coll_key, basename(f)))
      cache[[m]] <<- if (file.exists(f)) fread(f) else data.table()
    }
    cache[[m]]
  }

  blocks <- list(); audit <- list()
  for (i in seq_along(modules)) {
    k  <- modules[i]
    mm <- modes$mode_used[i]
    g  <- read_mode(mm)
    # NES > 0 only: the CSVs carry depleted terms too, and "what is this
    # section doing" means enriched, not depleted. For P09 module 3 the
    # two strongest hallmark terms by p-value are both negative-NES.
    tp <- if (nrow(g)) g[module == k & NES > 0][order(pval)][seq_len(min(N_TOP, .N))]
          else data.table()
    if (nrow(tp)) tp[, label := pretty_term(term, coll_key, width = TERM_WRAP)]

    blocks[[length(blocks) + 1]] <- list(
      module    = k,
      fill      = SPOT_PAL[k],
      ink       = TEXT_PAL[k],
      n_spots   = n_spots[i],
      pct       = 100 * n_spots[i] / n_cells,
      depth_adj = identical(mm, "depth_adjusted"),
      terms     = tp)

    if (nrow(tp)) {
      audit[[length(audit) + 1]] <- data.table(
        sample = ctx$sample_id, module = k, n_spots = n_spots[i],
        pct_spots = round(100 * n_spots[i] / n_cells, 2),
        r_with_total_expr = modes$r_with_total_expr[i],
        mode_used = mm, rank = seq_len(nrow(tp)), term = tp$term,
        NES = tp$NES, padj = tp$padj, significant = tp$padj < SIG)
    }
  }

  if (length(audit)) {
    fwrite(rbindlist(audit),
           file.path(enr_dir, sprintf("metagene_domain_top_%s.csv", coll_key)))
  }

  # --- footnote --------------------------------------------------------------
  # Built here because the canvas height below reserves room for it.
  n_depth <- sum(vapply(blocks, function(b) b$depth_adj, logical(1)))
  n_ns    <- sum(vapply(blocks, function(b)
                        if (nrow(b$terms)) sum(b$terms$padj >= SIG) else 0L,
                        numeric(1)))
  footnote <- paste0(
    "Top ", N_TOP, " MSigDB ", coll_name, " terms per metagene, ranked by GSEA ",
    "p-value over genes correlated with the metagene score; NES > 0 only. ",
    "Module member genes excluded from the ranking",
    if (n_depth) paste0(", except for the ", n_depth,
                        " module(s) marked above, which use the ",
                        "depth-controlled ranking (|r| > ", DEPTH_CUT,
                        " with total counts).")
    else ".",
    if (n_ns) sprintf(" %d term(s) marked n.s. did not reach padj < %.2f.",
                      n_ns, SIG) else "")
  # Hard-wrapped, since the panel clips at its own edge.
  footnote   <- paste(strwrap(footnote, width = FOOT_WRAP), collapse = "\n")
  FOOT_SLOTS <- foot_slots(footnote)

  # --- canvas geometry -------------------------------------------------------
  # Settled before either panel is drawn, because the tissue panel derives its
  # spot size from the physical height it will end up with.
  need_lines <- 0
  for (b in blocks) {
    n_lab <- if (nrow(b$terms))
      sum(vapply(strsplit(b$terms$label, "\n", fixed = TRUE), length, 1L)) else 1L
    need_lines <- need_lines + 1 + n_lab + as.integer(b$depth_adj) + 0.9
  }
  H <- max(8.6, HEAD_IN + LINE_IN * (need_lines + FOOT_SLOTS))
  W <- 14.2
  n_slots <- (H - HEAD_IN) / LINE_IN
  PANEL_H_IN <- H - HEAD_IN - 0.25        # title band + outer margins

  # --- tissue panel ----------------------------------------------------------
  # Raster, spot table, pitch, tags and limits all come from the context: they
  # are properties of the section, not of the collection. Only the point size is
  # recomputed, because it depends on the panel height, which depends on how
  # many lines this collection's terms wrapped to.
  dk   <- ctx$dk
  sl   <- ctx$sl
  tags <- ctx$tags
  xlim <- ctx$xlim
  ylim <- ctx$ylim

  # Spot size from the rendered geometry: inverting
  #   diameter_in = 0.75 * size * .pt / 72.27   (shape 21, stroke 0)
  # for SPOT_FILL of the measured pitch gives a near-touching lattice at this
  # panel's geometry.
  # Tissue panel width follows TISSUE_SHARE: W minus the outer plot margins, times
  # the tissue's share of the patchwork split, minus the tissue panel's own
  # left/right margins.
  panel_w  <- (W - (12 + 12) / 72) * TISSUE_SHARE / (TISSUE_SHARE + 1) - (2 + 4) / 72
  panel_h  <- PANEL_H_IN
  scale_in <- min(panel_w / diff(xlim), panel_h / diff(ylim))
  pt_size  <- scale_in * ctx$pitch * SPOT_FILL / (0.75 * .pt / 72.27)
  cat(sprintf("  %-8s pitch %.0f units, scale %.3g in/unit -> point_size %.2f\n",
              coll_key, ctx$pitch, scale_in, pt_size))

  # The numbered tags (white disc, ring and numeral in the module's colors) are
  # the required relief for the contrast rule — against the grey tissue several
  # fills drop below 3:1, and the rule is discharged by visible direct labels —
  # and they keep the territory-to-block link working in greyscale print. Placed
  # in sample_context(), because they mark territories, not terms.
  #
  # ONLY the spot lattice is rasterized. Four thousand vector circles is what
  # made the old SVG a 2.5 MB document that fights back when you open it, and
  # nothing is gained by keeping them editable — you are never going to retype a
  # spot. Everything above and below this layer stays vector: the tag discs are
  # still real <circle>, and every numeral, term and statistic is still <text>.
  p_tis <- ggplot() +
    annotation_raster(dk$raster, xmin = dk$xmin, xmax = dk$xmax,
                      ymin = dk$ymin, ymax = dk$ymax, interpolate = TRUE) +
    ggrastr::rasterise(
      geom_point(data = sl, aes(x = sdimx, y = sdimy, fill = dom_f),
                 shape = 21, size = pt_size, stroke = 0, alpha = 1),
      dpi = 600) +
    scale_fill_manual(values = setNames(SPOT_PAL[modules],
                                        as.character(modules)),
                      guide = "none") +
    geom_point(data = tags, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 21, size = 6.6, fill = "white", colour = tags$fill,
               stroke = 1.5) +
    geom_text(data = tags, aes(x = x, y = y, label = module),
              inherit.aes = FALSE, colour = tags$ink, size = 3.4,
              fontface = "bold") +
    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void() +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          legend.position  = "none",
          plot.margin      = margin(4, 2, 4, 4))

  # --- compose ---------------------------------------------------------------
  p_txt <- terms_panel(blocks, n_slots, ctx$sample_id, footnote, FOOT_SLOTS,
                       coll_name)

  fig <- (p_tis | p_txt) +
    plot_layout(widths = c(TISSUE_SHARE, 1)) +
    plot_annotation(
      title = sprintf("%s — metagene territories and their %s programs",
                      ctx$sample_id, coll_name),
      subtitle = sprintf(
        "%s outcome  ·  %d co-expression modules  ·  each spot colored by its dominant metagene (argmax of z-scored module score)",
        OUTCOME_MAP[[ctx$sample_id]], K),
      theme = theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.title      = element_text(colour = INK_HI, face = "bold",
                                       size = 15, margin = margin(b = 3)),
        plot.subtitle   = element_text(colour = INK_MID, size = 9.5,
                                       margin = margin(b = 6)),
        plot.margin     = margin(10, 12, 8, 12)))

  stem <- file.path(ctx$tp_dir, sprintf("metagene_domains_%s", coll_key))
  ggsave(paste0(stem, ".png"), fig, width = W, height = H, dpi = 600, bg = "white")

  # fix_text_size = FALSE: svglite's default (TRUE) adds textLength/lengthAdjust to
  # every string, which distorts glyphs when the text is edited in Illustrator.
  # ggsave passes ... to the device.
  ggsave(paste0(stem, ".svg"), fig, width = W, height = H, bg = "white",
         device = svglite::svglite, fix_text_size = FALSE)
  cat(sprintf("  ok  %.1f x %.1f in  ->  %s.{png,svg}\n", W, H, basename(stem)))
  invisible(NULL)
}

# =============================================================================
cat(sprintf("--- metagene territories + programs (%s) ---\n",
            paste(names(COLLECTIONS), collapse = ", ")))
for (s in SAMPLES) {
  tryCatch({
    ctx <- sample_context(s)
    if (!is.null(ctx)) {
      # Per-collection tryCatch as well as the per-sample one: a collection that
      # fails should cost its own figure, not the two that would have been drawn
      # after it from a context that is already built and paid for.
      for (cn in names(COLLECTIONS)) {
        tryCatch(draw_collection(ctx, cn), error = function(e) {
          cat(sprintf("  !! %s %s FAILED: %s\n", s, cn, conditionMessage(e)))
          warning(s, " ", cn, ": ", conditionMessage(e))
        })
      }
    }
    rm(ctx); invisible(gc(FALSE))
  }, error = function(e) {
    # report a failed save immediately
    cat(sprintf("  !! %s FAILED: %s\n", s, conditionMessage(e)))
    warning(s, ": ", conditionMessage(e))
  })
}
cat("\n========== DONE ==========\n\n")
