# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# biomarker_spatial_plots.R
# Loads saved Giotto objects for 4 ACC Visium samples and plots expression of
# good-outcome and poor-outcome biomarkers overlaid on tissue.
# Outputs: biomarkers_on_spatial_tissue/{SAMPLE}/{good|poor}_outcome/{GENE}_expr_on_tissue.png
# All PNGs at 600 dpi; all text at font size 14.

setwd(file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace"))
options(bitmapType = "cairo")

# =============================================================================
# CONSTANTS
# =============================================================================
BASE_DIR <- file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace")
OUT_DIR  <- file.path(BASE_DIR, "biomarkers_on_spatial_tissue")

SAMPLES <- c("P16", "P01", "P09", "P22")
OUTCOME_MAP <- c(P16 = "Poor", P01 = "Good",
                 P09 = "Good", P22 = "Poor")

good_outcome_biomarkers <- c("KRT77", "SCEL", "COL17A1", "TSBP1", "STATH",
                              "PRH2", "MPPED1", "LTF", "HOXA3", "KRT13")
poor_outcome_biomarkers <- c("NCAN", "VAX1", "RBFOX3", "GBX2", "BPIFB4",
                              "PLP1", "LIX1", "POU6F2", "BARX1", "LEFTY2")

FONT_SIZE  <- 14
PLOT_W     <- 8   # inches
PLOT_H     <- 7   # inches
POINT_SIZE <- 1.75
POINT_ALPHA <- 0.9

# Alpha of the greyscale tissue image against the white page: 1 = as before,
# lower = fainter. spatFeatPlot2D has no image-alpha argument, so this is baked
# into the luminance below; compositing a value c over white at alpha a is exactly
# 1 - a * (1 - c), so this constant IS the alpha. Kept in step with the per-sample
# pipeline scripts (Section 12) so both sets of on-tissue figures match.
TISSUE_ALPHA <- 0.75

# =============================================================================
# LIBRARIES
# =============================================================================
suppressPackageStartupMessages({
  library(Giotto)
  library(ggplot2)
  library(png)
  library(jsonlite)
  library(data.table)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

font_theme <- theme(
  text         = element_text(size = FONT_SIZE),
  axis.text    = element_text(size = FONT_SIZE),
  axis.title   = element_text(size = FONT_SIZE),
  legend.text  = element_text(size = FONT_SIZE),
  legend.title = element_text(size = FONT_SIZE),
  plot.title   = element_text(size = FONT_SIZE)
)

# =============================================================================
# MAIN LOOP
# =============================================================================
for (sample in SAMPLES) {
  outcome <- OUTCOME_MAP[[sample]]
  cat(sprintf("\n===== %s  [%s outcome] =====\n", sample, outcome))

  sample_out   <- file.path(OUT_DIR, sample)
  good_out_dir <- file.path(sample_out, "good_outcome")
  poor_out_dir <- file.path(sample_out, "poor_outcome")
  dir.create(good_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(poor_out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Load saved Giotto object
  # ---------------------------------------------------------------------------
  giotto_obj_dir <- file.path(BASE_DIR, "giotto_results", sample, "giotto_object")

  instrs <- createGiottoInstructions(
    save_dir    = sample_out,
    save_plot   = FALSE,
    show_plot   = FALSE,
    return_plot = TRUE,
    python_path = NULL
  )

  cat("  Loading Giotto object...\n")
  gobj <- loadGiotto(path_to_folder = giotto_obj_dir)
  gobj <- replaceGiottoInstructions(gobj, instructions = instrs)

  # ---------------------------------------------------------------------------
  # Build greyscale hires tissue image and attach (in-memory only)
  # ---------------------------------------------------------------------------
  # FIXED: the sample dirs moved under data/ in the Aug 2026 restructure, so this
  # resolved to a non-existent path and jsonlite then tried to parse the path string
  # itself as JSON ("lexical error: invalid char in json text"). Matches the
  # DATA_DIR/VISIUM_DIR convention used by the four per-sample scripts.
  visium_dir  <- file.path(BASE_DIR, "data", paste0("sample_", sample))
  scalef_json <- jsonlite::fromJSON(
    file.path(visium_dir, "spatial", "scalefactors_json.json"))
  hires_src   <- file.path(visium_dir, "spatial", "tissue_hires_image.png")
  gray_path   <- file.path(sample_out, "tissue_hires_gray.png")

  cat("  Building greyscale tissue image...\n")
  img_px <- png::readPNG(hires_src)
  lum <- if (length(dim(img_px)) == 3 && dim(img_px)[3] >= 3) {
    0.2989 * img_px[,,1] + 0.5870 * img_px[,,2] + 0.1140 * img_px[,,3]
  } else {
    img_px
  }
  # Composite the tissue over white at TISSUE_ALPHA (see the constant above).
  lum <- 1 - TISSUE_ALPHA * (1 - lum)
  gray_rgb      <- array(0, dim = c(nrow(lum), ncol(lum), 3))
  gray_rgb[,,1] <- lum
  gray_rgb[,,2] <- lum
  gray_rgb[,,3] <- lum
  png::writePNG(gray_rgb, gray_path)

  sf       <- 1 / scalef_json$tissue_hires_scalef
  gray_img <- createGiottoLargeImage(
    raster_object = gray_path,
    scale_factor  = sf,
    negative_y    = TRUE,
    name          = "tissue_hires_gray"
  )
  gobj <- addGiottoImage(gobj, images = list(gray_img))
  cat(sprintf("  Greyscale image attached (scale factor: %.3f)\n", sf))

  available_feats <- fDataDT(gobj)$feat_ID

  # ---------------------------------------------------------------------------
  # Plot helper
  # ---------------------------------------------------------------------------
  plot_gene <- function(gene, out_dir) {
    if (!gene %in% available_feats) {
      cat(sprintf("  [SKIP] %s — not found in expression matrix\n", gene))
      return(invisible(NULL))
    }

    p <- tryCatch(
      spatFeatPlot2D(
        gobj,
        show_image        = TRUE,
        image_name        = "tissue_hires_gray",
        feats             = gene,
        expression_values   = "normalized",
        cell_color_gradient = c("grey15", "#FF4500", "#FFD700"),
        gradient_style      = "sequential",
        point_size          = POINT_SIZE,
        point_alpha         = POINT_ALPHA,
        coord_fix_ratio     = 1,
        save_plot           = FALSE,
        show_plot           = FALSE
      ),
      error = function(e) {
        warning(sprintf("  [ERROR] %s: %s\n", gene, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(p)) return(invisible(NULL))

    # spatFeatPlot2D may return a list when called with a single gene
    if (is.list(p) && !inherits(p, "gg")) p <- p[[1]]

    p <- p + font_theme

    out_file <- file.path(out_dir, paste0(gene, "_expr_on_tissue.png"))
    ggsave(out_file, plot = p, dpi = 600,
           width = PLOT_W, height = PLOT_H, units = "in")
    cat(sprintf("  Saved: %s/%s\n", basename(out_dir), basename(out_file)))
  }

  # ---------------------------------------------------------------------------
  # Good-outcome biomarkers
  # ---------------------------------------------------------------------------
  cat("  --- Good-outcome biomarkers ---\n")
  for (gene in good_outcome_biomarkers) plot_gene(gene, good_out_dir)

  # ---------------------------------------------------------------------------
  # Poor-outcome biomarkers
  # ---------------------------------------------------------------------------
  cat("  --- Poor-outcome biomarkers ---\n")
  for (gene in poor_outcome_biomarkers) plot_gene(gene, poor_out_dir)

  cat(sprintf("  Finished %s\n", sample))
}

cat("\n========== DONE ==========\n")
cat("All outputs in:", OUT_DIR, "\n\n")
