# Stage 35 reference deconvolution using the PUBLISHED DWLS package (CRAN 0.1.0).
# Reads a signature matrix and a bulk matrix written by 35_dwls_crosscheck.py and writes
# one row of compartment proportions per bulk sample. No project data is modified.
suppressMessages(library(DWLS))
args <- commandArgs(trailingOnly = TRUE)
sig_f <- args[1]; bulk_f <- args[2]; out_f <- args[3]

S <- as.matrix(read.csv(sig_f, row.names = 1, check.names = FALSE))
B <- as.matrix(read.csv(bulk_f, row.names = 1, check.names = FALSE))
common <- intersect(rownames(S), rownames(B))
S <- S[common, , drop = FALSE]; B <- B[common, , drop = FALSE]
cat("reference DWLS:", length(common), "genes x", ncol(S), "compartments,",
    ncol(B), "samples\n")

res <- matrix(NA_real_, nrow = ncol(B), ncol = ncol(S),
              dimnames = list(colnames(B), colnames(S)))
for (i in seq_len(ncol(B))) {
  w <- tryCatch(DWLS::solveDampenedWLS(S, B[, i]),
                error = function(e) { cat("  sample", colnames(B)[i], "failed:",
                                          conditionMessage(e), "\n"); rep(NA_real_, ncol(S)) })
  w[w < 0] <- 0
  if (sum(w, na.rm = TRUE) > 0) w <- w / sum(w, na.rm = TRUE)
  res[i, ] <- w
}
write.csv(res, out_f)
cat("wrote", out_f, "\n")
