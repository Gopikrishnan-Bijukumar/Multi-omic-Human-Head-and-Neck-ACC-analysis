#!/usr/bin/env Rscript
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# =============================================================================
# fetch_genesets.R — one-time setup for Section 6e metagene enrichment
#
#     Rscript scripts/fetch_genesets.R
#
# Plain Rscript. Do NOT `conda activate giotto_env` first — on this server the R
# packages live in ~/R/x86_64-pc-linux-gnu-library/4.4 (picked up automatically
# via ~/.Renviron) and giotto_env is Python-only. fgsea and fastmatch are already
# installed as binaries here, so unlike the Anvil/HD version this script does not
# need a compiler and does not build anything from source.
#
# This is the ONLY network-dependent step in the entire pipeline. Running it once
# makes every downstream run fully offline, SHA-pinned and reproducible.
#
# WHAT IT DOES
#   1. checks the two required packages are present (fastmatch + fgsea)
#   2. downloads MSigDB symbol GMTs (Hallmark, GO:BP, GO:MF, Reactome)
#   3. VERIFIES the downloads actually contain real gene sets — not just that a
#      file arrived — including a round-trip fgsea test
#   4. writes genesets/MANIFEST.json with a SHA-256 per file, which Section 6e
#      re-checks on every run
#
# WHY NOT msigdbr: its gene set data was moved into a separate `msigdbdf`
# package hosted on r-universe, and that package is NOT currently resolvable.
# The msigdbr path is a live failure, not a theoretical risk.
#
# WHY NOT clusterProfiler: ~40 packages including the AnnotationDbi/SQLite chain,
# plus a mandatory symbol->Entrez mapping that silently drops 5-10% of symbols.
# Its enricher() is phyper, which Section 6e calls directly in ~10 lines.
# =============================================================================

suppressPackageStartupMessages({library(utils); library(tools)})

BASE_DIR    <- file.path(ACC_DATA_ROOT, "spatial/reg_visium/giotto_w_cytospace")
GENESET_DIR <- file.path(BASE_DIR, "data", "genesets")
dir.create(GENESET_DIR, recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

# -----------------------------------------------------------------------------
# 1. packages
# -----------------------------------------------------------------------------
cat("=== [1/4] R packages ===\n")
giotto_before <- tryCatch(as.character(packageVersion("GiottoClass")), error = function(e) NA)
for (pkg in c("fastmatch", "fgsea", "jsonlite", "digest")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-10s already installed (%s)\n", pkg, packageVersion(pkg)))
  } else {
    cat(sprintf("  %-10s installing...\n", pkg))
    # update = FALSE IS MANDATORY. An unpinned Bioconductor update could pull a
    # newer Matrix/S4Vectors underneath the working Giotto 4.2.2 install and
    # break the whole pipeline for a plotting dependency.
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
giotto_after <- tryCatch(as.character(packageVersion("GiottoClass")), error = function(e) NA)
if (!is.na(giotto_before) && !identical(giotto_before, giotto_after))
  stop("GiottoClass version changed during install (", giotto_before, " -> ", giotto_after,
       "). This pipeline runs on Giotto 4.2.2 — investigate before running anything.")
cat(sprintf("  GiottoClass still %s (unchanged)\n", giotto_after))
stopifnot(requireNamespace("fgsea", quietly = TRUE))

# -----------------------------------------------------------------------------
# 2. gene sets — fallback ladder, stop at the first source that VERIFIES
# -----------------------------------------------------------------------------
cat("\n=== [2/4] gene set collections ===\n")

msigdb_files <- function(rel) c(
  hallmark      = sprintf("h.all.v%s.symbols.gmt",           rel),
  go_bp         = sprintf("c5.go.bp.v%s.symbols.gmt",        rel),
  go_mf         = sprintf("c5.go.mf.v%s.symbols.gmt",        rel),
  reactome      = sprintf("c2.cp.reactome.v%s.symbols.gmt",  rel))

SOURCES <- list(
  list(name = "MSigDB 2026.1.Hs", release = "2026.1.Hs",
       base = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/",
       files = msigdb_files("2026.1.Hs")),
  list(name = "MSigDB 2024.1.Hs", release = "2024.1.Hs",
       base = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/",
       files = msigdb_files("2024.1.Hs")),
  # Genuinely INDEPENDENT source (not just an older Broad snapshot). No Hallmark;
  # we proceed with 3 collections and say so loudly.
  list(name = "Bader Lab (current)", release = "baderlab-current",
       base = "https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/",
       files = c(go_bp    = "GO/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_symbol.gmt",
                 reactome = "Pathways/Human_Reactome_symbol.gmt"))
)

read_gmt <- function(path) {
  ln <- readLines(path, warn = FALSE); ln <- ln[nzchar(ln)]
  p  <- strsplit(ln, "\t", fixed = TRUE)
  setNames(lapply(p, function(x) unique(x[-c(1, 2)][nzchar(x[-c(1, 2)])])),
           vapply(p, `[`, "", 1L))
}

# The verification that matters: not "did a file arrive" but "does it contain
# real, usable, human gene sets". Canaries are checked against the ACTUAL release
# rather than hardcoded across releases, because MSigDB renames sets between
# them — GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION was renamed out of 2026.1.Hs.
verify_gmt <- function(path, label) {
  if (!file.exists(path) || file.size(path) < 10000) return("file missing or too small")
  sets <- tryCatch(read_gmt(path), error = function(e) NULL)
  if (is.null(sets))                 return("unparseable")
  if (length(sets) < 40L)            return(sprintf("only %d sets", length(sets)))
  if (anyDuplicated(names(sets)))    return("duplicate set names")
  syms <- unique(unlist(sets, use.names = FALSE))
  if (length(syms) < 500L)           return(sprintf("only %d unique symbols", length(syms)))
  # human symbols look like HGNC: mostly uppercase alphanumeric
  if (mean(grepl("^[A-Z0-9][A-Z0-9._-]*$", head(syms, 2000))) < 0.9)
    return("symbols do not look like HGNC identifiers")
  NULL
}

fetched <- NULL
for (src in SOURCES) {
  cat(sprintf("\n  trying: %s\n", src$name))
  ok <- TRUE; got <- list()
  for (nm in names(src$files)) {
    url  <- paste0(src$base, src$files[[nm]])
    dest <- file.path(GENESET_DIR, sprintf("%s.v%s.symbols.gmt", nm, src$release))
    r <- tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE),
                  error = function(e) 1L, warning = function(w) 1L)
    bad <- if (!identical(as.integer(r), 0L)) "download failed" else verify_gmt(dest, nm)
    if (!is.null(bad)) {
      cat(sprintf("    %-10s FAILED (%s)\n", nm, bad)); ok <- FALSE
      if (file.exists(dest)) unlink(dest)
      break
    }
    sets <- read_gmt(dest)
    cat(sprintf("    %-10s OK  %5d sets, %6d unique symbols\n",
                nm, length(sets), length(unique(unlist(sets, use.names = FALSE)))))
    got[[nm]] <- list(path = dest, url = url, n_sets = length(sets))
  }
  if (ok && length(got)) { fetched <- list(src = src, files = got); break }
}

if (is.null(fetched)) {
  stop("\n  Could not retrieve any gene set collection.\n",
       "  All sources failed. If this machine is behind a proxy, set:\n",
       "    export https_proxy=http://<host>:<port>\n",
       "    export http_proxy=http://<host>:<port>\n",
       "  and re-run. Verify DNS with:  getent hosts data.broadinstitute.org\n")
}
cat(sprintf("\n  using: %s (%d collections)\n", fetched$src$name, length(fetched$files)))
if (is.null(fetched$files$hallmark))
  cat("  NOTE: this source has no Hallmark collection — proceeding with",
      length(fetched$files), "collections.\n")

# -----------------------------------------------------------------------------
# 3. functional round-trip test — proves the parser, the sets AND fgsea all work
#    together before a single real result is computed
# -----------------------------------------------------------------------------
cat("\n=== [3/4] round-trip functional test ===\n")
tst  <- read_gmt(fetched$files[[1]]$path)
# pick a real set of usable size from the actual downloaded file
cand <- names(tst)[vapply(tst, length, 1L) >= 20L & vapply(tst, length, 1L) <= 200L]
if (!length(cand)) stop("no gene set of testable size found — collection looks wrong")
canary <- cand[1]
members <- tst[[canary]]
universe <- unique(c(unlist(tst, use.names = FALSE)))
universe <- head(unique(c(members, setdiff(universe, members))), 8000)
# synthetic ranking: canary members at the very top
stats <- setNames(c(seq(3, 1, length.out = length(members)),
                    rnorm(length(universe) - length(members), 0, 0.5)),
                  c(members, setdiff(universe, members)))
fg <- fgsea::fgsea(pathways = tst, stats = stats, minSize = 10, maxSize = 500,
                   eps = 0.0, nPermSimple = 1000, nproc = 1)
row <- as.data.frame(fg)[fg$pathway == canary, ]
cat(sprintf("  canary set: %s (%d genes)\n", canary, length(members)))
cat(sprintf("  fgsea NES = %.2f, padj = %.3g\n", row$NES, row$padj))
if (!(is.finite(row$NES) && row$NES > 0 && row$padj < 0.01))
  stop("fgsea round-trip FAILED — a set planted at the top of the ranking was not recovered.")
# and the same overlap must be significant by phyper
ph <- stats::phyper(length(members) - 1L, length(members),
                    length(universe) - length(members), length(members),
                    lower.tail = FALSE)
cat(sprintf("  phyper on the same overlap: p = %.3g\n", ph))
if (!(ph < 1e-10)) stop("phyper round-trip FAILED")
cat("  both engines verified.\n")

# -----------------------------------------------------------------------------
# 4. manifest
# -----------------------------------------------------------------------------
cat("\n=== [4/4] manifest ===\n")
man <- list(
  source     = fetched$src$name,
  release    = fetched$src$release,
  downloaded = format(Sys.time(), tz = "UTC", usetz = TRUE),
  host       = Sys.info()[["nodename"]],
  r_version  = as.character(getRversion()),
  fgsea      = as.character(packageVersion("fgsea")),
  files      = setNames(lapply(names(fetched$files), function(nm) {
    f <- fetched$files[[nm]]
    list(url = f$url, sha256 = digest::digest(file = f$path, algo = "sha256"),
         bytes = as.numeric(file.size(f$path)), n_sets = f$n_sets)
  }), vapply(fetched$files, function(f) basename(f$path), ""))
)
jsonlite::write_json(man, file.path(GENESET_DIR, "MANIFEST.json"),
                     auto_unbox = TRUE, pretty = TRUE)
writeLines(c(
  sprintf("Gene sets for regular-Visium metagene enrichment (Section 6e)"),
  sprintf("source      : %s", man$source),
  sprintf("release     : %s", man$release),
  sprintf("downloaded  : %s on %s", man$downloaded, man$host),
  sprintf("fgsea       : %s   R: %s", man$fgsea, man$r_version),
  "",
  "Files (SHA-256 re-verified by Section 6e on every run):",
  vapply(names(man$files), function(f)
    sprintf("  %-40s %7d sets  %s", f, man$files[[f]]$n_sets, man$files[[f]]$sha256), "")
), file.path(GENESET_DIR, "PROVENANCE.txt"))

cat("  wrote MANIFEST.json + PROVENANCE.txt\n")
cat(sprintf("\nDONE. %d collections in %s\n", length(fetched$files), GENESET_DIR))
cat("Section 6e will now run offline on every pipeline invocation.\n")
