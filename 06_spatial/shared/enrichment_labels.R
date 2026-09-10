# =============================================================================
# enrichment_labels.R — MSigDB term-label prettifying, shared by every script
# that puts a gene-set name on a figure.
#
#   source(file.path(BASE_DIR, "scripts", "enrichment_labels.R"))
#
# Side-effect free on purpose: it defines ACRONYMS and pretty_term() and does
# nothing else, so any script can source it at any point without surprises.
#
# This used to live inside plot_enrichment.R. It was lifted out when a second
# figure script needed the same labels — for exactly the reason plot_enrichment.R
# gives at its own top: two copies of shared plotting code drift apart silently.
# =============================================================================

# Term labels. MSigDB names are SHOUTED_WITH_UNDERSCORES and prefixed by their
# collection, which wastes the entire left margin on a string repeated in every
# row. Strip the prefix, Title Case the rest, then restore the acronyms that
# Title Case would otherwise mangle into "Tnfa" / "Dna".
ACRONYMS <- c("DNA", "RNA", "MRNA", "TRNA", "RRNA", "NADH", "ATP", "ADP", "GTP",
              "TNFA", "TNF", "NFKB", "MYC", "KRAS", "TGF", "TGFB", "WNT", "MTOR",
              "IL2", "IL6", "STAT3", "STAT5", "JAK", "PI3K", "AKT", "MTORC1",
              "EMT", "ECM", "MHC", "IFN", "UV", "ROS", "G2M", "E2F", "P53",
              "APC", "CDC", "ER", "GPCR", "SRP", "NK", "TCA", "II", "III", "IV")

# `width` is the wrap column. plot_enrichment.R's dotplot axis takes 42; a
# narrower side panel passes a smaller number. It stays an argument with the
# original default so the dotplots are byte-identical after the extraction.
pretty_term <- function(x, collection, width = 42) {
  x <- sub(paste0("^", toupper(gsub("_", "", collection)), "_"), "", x)
  x <- sub("^(HALLMARK|GOBP|GOMF|GOCC|REACTOME|KEGG|WP|BIOCARTA)_", "", x)
  x <- gsub("_", " ", x)
  # Title Case each word, then put the acronyms back.
  x <- vapply(strsplit(tolower(x), " ", fixed = TRUE), function(w) {
    w <- ifelse(nchar(w) > 0,
                paste0(toupper(substring(w, 1, 1)), substring(w, 2)), w)
    hit <- toupper(w) %in% ACRONYMS
    w[hit] <- toupper(w[hit])
    paste(w, collapse = " ")
  }, character(1))
  # Wrap rather than truncate: a term cut at 60 chars loses the words that
  # distinguish it from its neighbours.
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}
