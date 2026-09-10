# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

library(WGCNA)
library(DESeq2)
library(gplots)

allowWGCNAThreads()

# Output location for all plots, tables, and saved R objects
output_dir <- file.path(ACC_DATA_ROOT, "bulk_rna_data/outputs/images/clinical_outcome/wgcna/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# SECTION 1: Load Data
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 1: Loading Data\n")
cat("=============================================================================\n\n")

# Load raw counts (samples as rows, genes as columns)
counts_raw <- read.csv(file.path(ACC_DATA_ROOT, "bulk_rna_data/outputs/csv_outputs/counts_for_WGCNA.csv"), row.names = 1)
cat("Raw counts dimensions:", nrow(counts_raw), "samples x", ncol(counts_raw), "genes\n")

# Load metadata
metadata <- read.csv(file.path(ACC_DATA_ROOT, "bulk_rna_data/outputs/csv_outputs/bulk_meta_data.csv"), row.names = 1)
cat("Metadata dimensions:", nrow(metadata), "samples x", ncol(metadata), "traits\n\n")

# Verify sample order matches
if (!all(rownames(counts_raw) == rownames(metadata))) {
  metadata <- metadata[rownames(counts_raw), ]
  cat("Reordered metadata to match counts\n")
}

# -----------------------------------------------------------------------------
# SECTION 2: Pre-filtering Before Normalization
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 2: Pre-filtering Genes\n")
cat("=============================================================================\n\n")

# Transpose for filtering (genes as rows temporarily)
counts_t <- t(counts_raw)

# Filter: keep genes with >= 10 counts in at least 50% of samples
min_samples <- ceiling(ncol(counts_t) * 0.5)
keep_genes <- rowSums(counts_t >= 10) >= min_samples
counts_filtered <- counts_t[keep_genes, ]

cat("Genes before filtering:", nrow(counts_t), "\n")
cat("Genes after filtering (>=10 counts in >=50% samples):", nrow(counts_filtered), "\n\n")

# -----------------------------------------------------------------------------
# SECTION 3: VST Normalization with DESeq2
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 3: VST Normalization\n")
cat("=============================================================================\n\n")

# Ensure counts are integers
counts_int <- apply(counts_filtered, c(1, 2), as.integer)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData = metadata,
  design = ~ 1
)

cat("Created DESeq2 object\n")

# Apply VST normalization (blind = TRUE for unbiased transformation)
vsd <- vst(dds, blind = TRUE)
cat("Applied VST normalization\n")

# Extract normalized data and transpose for WGCNA (samples as rows, genes as columns)
datExpr <- as.data.frame(t(assay(vsd)))

cat("Normalized expression matrix:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n\n")

# -----------------------------------------------------------------------------
# SECTION 4: Variance Filtering
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 4: Variance Filtering\n")
cat("=============================================================================\n\n")

# Calculate variance for each gene
gene_variance <- apply(datExpr, 2, var)

# Keep top 75% most variable genes (remove bottom 25%)
variance_threshold <- quantile(gene_variance, 0.25)
high_var_genes <- gene_variance > variance_threshold
datExpr <- datExpr[, high_var_genes]

cat("Removed bottom 25% low-variance genes\n")
cat("Genes remaining:", ncol(datExpr), "\n\n")

# -----------------------------------------------------------------------------
# SECTION 5: Quality Control - Check for Bad Genes/Samples
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 5: Quality Control\n")
cat("=============================================================================\n\n")

gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing", sum(!gsg$goodGenes), "problematic genes\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("Removing", sum(!gsg$goodSamples), "problematic samples\n")
  }
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

cat("After QC:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n\n")

# -----------------------------------------------------------------------------
# SECTION 6: Prepare Trait Data
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 6: Prepare Trait Data\n")
cat("=============================================================================\n\n")

# Match samples to expression data
traitData <- metadata[rownames(datExpr), ]

# Convert categorical traits to numeric
datTraits <- data.frame(
  ClinicalOutcome = ifelse(traitData$ClinicalOutcome == "5Plus", 1, 0),
  DEGOutcome = ifelse(traitData$DEGOutcome == "5Plus", 1, 0),
  ModelOutcome1 = ifelse(traitData$ModelOutcome1 == "Average", 1, 0),
  ModelOutcome2 = ifelse(traitData$ModelOutcome2 == "Average", 1, 0),
  row.names = rownames(traitData)
)

cat("Numeric trait matrix created\n")
cat("Trait encoding:\n")
cat("  ClinicalOutcome/DEGOutcome: 5Plus=1, 2Minus=0\n")
cat("  ModelOutcome1/2: Average=1, Poor=0\n\n")

# -----------------------------------------------------------------------------
# SECTION 6B: Initial Sample Clustering (All 20 Samples, Before Outlier Removal)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 6B: Initial Sample Clustering (All 20 Samples)\n")
cat("=============================================================================\n\n")

sampleTree_all20 <- hclust(dist(datExpr), method = "average")

svg(paste0(output_dir, "00_sample_clustering_pre_removal.svg"), width = 12, height = 9)
par(cex = 0.8)
par(mar = c(2, 4, 2, 0))
plot(sampleTree_all20,
     main = "Sample Clustering (All 20 Samples, Before Outlier Removal)",
     sub = "",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.2,
     cex.main = 1.5)
dev.off()

cat("Saved: 00_sample_clustering_pre_removal.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 7: Remove Outlier Sample
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 7: Remove Outlier Sample\n")
cat("=============================================================================\n\n")

# Confirmed via visual inspection of 00_sample_clustering_pre_removal.svg (all 20 samples):
# P12 is a definite outlier; no second outlier was found.
outlier_samples <- c("P12")

cat("Removing outlier sample(s):", paste(outlier_samples, collapse = ", "), "(identified from initial clustering)\n")

# Remove from expression data
datExpr <- datExpr[!(rownames(datExpr) %in% outlier_samples), ]

# Remove from trait data
datTraits <- datTraits[!(rownames(datTraits) %in% outlier_samples), ]

cat("Samples remaining:", nrow(datExpr), "\n")
cat("Sample names:", paste(rownames(datExpr), collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# SECTION 8: Sample Clustering (After Outlier Removal)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 8: Sample Clustering\n")
cat("=============================================================================\n\n")

# Hierarchical clustering of samples
sampleTree <- hclust(dist(datExpr), method = "average")

# Plot sample dendrogram
svg(paste0(output_dir, "01_sample_clustering.svg"), width = 12, height = 9)
par(cex = 0.8)
par(mar = c(2, 4, 2, 0))
plot(sampleTree,
     main = "Sample Clustering (After Outlier Removal)",
     sub = "",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.2,
     cex.main = 1.5)
dev.off()

cat("Saved: 01_sample_clustering.svg\n")

# Plot sample dendrogram with trait heatmap
svg(paste0(output_dir, "02_sample_dendrogram_traits.svg"), width = 12, height = 9)
traitColors <- numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels = colnames(datTraits),
  main = "Sample Dendrogram with Trait Heatmap"
)
dev.off()

cat("Saved: 02_sample_dendrogram_traits.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 9: Soft Threshold Power Selection
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 9: Soft Threshold Power Selection\n")
cat("=============================================================================\n\n")

# Test range of powers
powers <- c(1:10, seq(12, 20, by = 2))

# Calculate scale-free topology fit
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

# Plot results
svg(paste0(output_dir, "03_soft_threshold_selection.svg"), width = 10, height = 5)
par(mfrow = c(1, 2))

# Scale-free topology fit
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R²)",
     type = "n",
     main = "Scale Independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     col = "red",
     cex = 0.9)
abline(h = 0.85, col = "red", lty = 2)
abline(h = 0.80, col = "blue", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     col = "red",
     cex = 0.9)

dev.off()

cat("Saved: 03_soft_threshold_selection.svg\n")

# Select power - use estimated or default to 14 if R² doesn't reach 0.85
# softPower <- sft$powerEstimate
softPower <- 14
# if (is.na(softPower)) {
#   softPower <- 14  # Use 14 based on our analysis (R² plateaus here)
#   cat("NOTE: R² did not reach 0.85. Using power = 14 (where R² plateaus)\n")
# } else {
#   cat("Estimated soft-thresholding power:", softPower, "\n")
# }

cat("Using soft power:", softPower, "\n\n")

# -----------------------------------------------------------------------------
# SECTION 10: Network Construction and Module Detection
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 10: Network Construction & Module Detection\n")
cat("=============================================================================\n\n")

# IMPORTANT: Override cor function to use WGCNA's version
cor <- WGCNA::cor

net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 50,
  reassignThreshold = 0,
  mergeCutHeight = 0.40,
  deepSplit = 1,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = paste0(output_dir, "TOM"),
  verbose = 3
)

# Restore base R cor function after blockwiseModules
cor <- stats::cor

# Module summary
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
nModules <- length(unique(moduleColors)) - 1  # Exclude grey

cat("\nModules detected:", nModules, "(excluding grey/unassigned)\n")
cat("\nModule sizes:\n")
print(table(moduleColors))

# Plot gene dendrogram with modules
svg(paste0(output_dir, "04_gene_dendrogram_modules.svg"), width = 12, height = 9)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram with Module Colors"
)
dev.off()

cat("\nSaved: 04_gene_dendrogram_modules.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 11: Module-Trait Relationships
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 11: Module-Trait Relationships\n")
cat("=============================================================================\n\n")

# Calculate module eigengenes
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Correlate modules with traits
nSamples <- nrow(datExpr)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Create annotation matrix
textMatrix <- paste0(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1),
  ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

# Plot heatmap
svg(paste0(output_dir, "05_module_trait_relationships.svg"), width = 10, height = 8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships\n(correlation, p-value)"
)
dev.off()

cat("Saved: 05_module_trait_relationships.svg\n")

# Report significant associations
cat("\nSignificant module-trait associations (p < 0.05):\n")
for (trait in colnames(datTraits)) {
  sig_modules <- names(MEs)[moduleTraitPvalue[, trait] < 0.05]
  if (length(sig_modules) > 0) {
    cat("  ", trait, ":", paste(sig_modules, collapse = ", "), "\n")
  }
}
cat("\n")

# -----------------------------------------------------------------------------
# SECTION 12: Gene Significance and Module Membership (ClinicalOutcome — main)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 12: Gene Significance & Module Membership\n")
cat("=============================================================================\n\n")

# Calculate module membership (kME) for all genes
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

colnames(geneModuleMembership) <- paste0("MM.", gsub("ME", "", names(MEs)))
colnames(MMPvalue) <- paste0("p.MM.", gsub("ME", "", names(MEs)))

# Calculate gene significance for ClinicalOutcome
geneTraitSig <- as.data.frame(cor(datExpr, datTraits$ClinicalOutcome, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSig), nSamples))

colnames(geneTraitSig) <- "GS.ClinicalOutcome"
colnames(GSPvalue) <- "p.GS.ClinicalOutcome"

cat("Calculated gene significance for ClinicalOutcome\n")
cat("Calculated module membership for all modules\n\n")

# -----------------------------------------------------------------------------
# SECTION 13: Create Gene Information Table (ClinicalOutcome — main)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 13: Export Gene Information\n")
cat("=============================================================================\n\n")

# Compile gene info
geneInfo <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  GS.ClinicalOutcome = geneTraitSig$GS.ClinicalOutcome,
  p.GS.ClinicalOutcome = GSPvalue$p.GS.ClinicalOutcome
)

# Add module membership columns
geneInfo <- cbind(geneInfo, geneModuleMembership, MMPvalue)

# Sort by gene significance
geneInfo <- geneInfo[order(-abs(geneInfo$GS.ClinicalOutcome)), ]

# Save full table
write.csv(geneInfo, paste0(output_dir, "06_gene_info_complete.csv"), row.names = FALSE)
cat("Saved: 06_gene_info_complete.csv\n")

# -----------------------------------------------------------------------------
# SECTION 14: Extract Hub Genes (ClinicalOutcome — main)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 14: Hub Gene Identification\n")
cat("=============================================================================\n\n")

# Get hub genes for each module (top 20 by module membership)
hub_genes_all <- data.frame()

for (module in unique(moduleColors)) {
  if (module == "grey") next

  # Get genes in this module
  module_genes <- geneInfo[geneInfo$Module == module, ]

  # Get the module membership column for this module
  mm_col <- paste0("MM.", module)

  if (mm_col %in% colnames(module_genes) && nrow(module_genes) > 0) {
    # Sort by absolute module membership
    module_genes <- module_genes[order(-abs(module_genes[, mm_col])), ]

    # Take top 20 (or fewer if module is smaller)
    n_top <- min(20, nrow(module_genes))
    top_hubs <- module_genes[1:n_top, ]

    # Create standardized output
    hub_df <- data.frame(
      Gene = top_hubs$Gene,
      Module = top_hubs$Module,
      Rank = 1:n_top,
      GS.ClinicalOutcome = top_hubs$GS.ClinicalOutcome,
      p.GS.ClinicalOutcome = top_hubs$p.GS.ClinicalOutcome,
      ModuleMembership = top_hubs[, mm_col]
    )

    hub_genes_all <- rbind(hub_genes_all, hub_df)
  }
}

write.csv(hub_genes_all, paste0(output_dir, "07_hub_genes_top20.csv"), row.names = FALSE)
cat("Saved: 07_hub_genes_top20.csv\n")
cat("Extracted top 20 hub genes per module\n\n")

# -----------------------------------------------------------------------------
# SECTION 15: Module Eigengene Visualization
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 15: Module Eigengene Plots\n")
cat("=============================================================================\n\n")

# Eigengene barplot by sample
svg(paste0(output_dir, "08_module_eigengenes_barplot.svg"), width = 14, height = 10)

nMods <- ncol(MEs)
nRows <- ceiling(sqrt(nMods))
nCols <- ceiling(nMods / nRows)
par(mfrow = c(nRows, nCols), mar = c(4, 3, 2, 1))

for (i in 1:ncol(MEs)) {
  module_name <- gsub("ME", "", names(MEs)[i])
  barplot(
    MEs[, i],
    names.arg = rownames(datExpr),
    las = 2,
    main = paste0("ME", module_name),
    col = module_name,
    ylab = "Eigengene",
    cex.names = 0.7
  )
}

dev.off()
cat("Saved: 08_module_eigengenes_barplot.svg\n")

# Eigengene heatmap with ClinicalOutcome
svg(paste0(output_dir, "09_module_eigengene_heatmap.svg"), width = 10, height = 8)
plotMEpairs(MEs, y = datTraits$ClinicalOutcome)
dev.off()
cat("Saved: 09_module_eigengene_heatmap.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 16: Save R Objects
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 16: Save R Objects\n")
cat("=============================================================================\n\n")

save(
  datExpr,
  datTraits,
  metadata,
  net,
  MEs,
  moduleColors,
  moduleLabels,
  moduleTraitCor,
  moduleTraitPvalue,
  geneModuleMembership,
  geneTraitSig,
  geneInfo,
  softPower,
  file = paste0(output_dir, "WGCNA_results.RData")
)

cat("Saved: WGCNA_results.RData\n")
cat("Load with: load('WGCNA_results.RData')\n\n")

# -----------------------------------------------------------------------------
# SECTION 17: Module Membership vs Gene Significance Scatter Plots
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 17: Module Membership vs Gene Significance Plots\n")
cat("=============================================================================\n\n")

# Key modules = modules significantly correlated with ClinicalOutcome (p < 0.05), grey excluded
key_modules <- rownames(moduleTraitPvalue)[moduleTraitPvalue[, "ClinicalOutcome"] < 0.05]
key_modules <- gsub("ME", "", key_modules)
key_modules <- setdiff(key_modules, "grey")

key_modules <- key_modules[key_modules %in% unique(moduleColors)]

# Create scatter plots for key modules
svg(paste0(output_dir, "10_MM_vs_GS_scatter.svg"), width = 15, height = 10)
n_km <- length(key_modules)
par(mfrow = c(ceiling(n_km / 4), min(4, max(1, n_km))))  # adaptive grid, up to 4 cols

for (module in key_modules) {
  column <- match(module, substring(names(geneModuleMembership), 4))
  if (is.na(column)) next

  moduleGenes <- moduleColors == module

  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSig[moduleGenes, 1]),
    xlab = paste("Module Membership (", module, ")", sep = ""),
    ylab = "Gene Significance (ClinicalOutcome)",
    main = paste(module, "module"),
    cex.main = 1.2,
    cex.lab = 1.1,
    cex.axis = 1.1,
    col = module
  )
}

dev.off()
cat("Saved: 10_MM_vs_GS_scatter.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 18: TOM-based Network Heatmap
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 18: TOM Network Heatmap\n")
cat("=============================================================================\n\n")

# ADJ (unsigned abs-correlation adjacency) feeds intramodularConnectivity() in Section 20 below;
# kept distinct from the signed 'adjacency' matrix computed next for the TOM plots
cor <- WGCNA::cor
ADJ <- abs(cor(datExpr, use = "p"))^softPower
cor <- stats::cor

cat("Recalculating TOM matrix (this may take a few minutes)...\n")

cor <- WGCNA::cor

# Calculate adjacency matrix
adjacency <- adjacency(datExpr, power = softPower, type = "signed")

# Calculate TOM
TOM <- TOMsimilarity(adjacency, TOMType = "signed")

cor <- stats::cor

# Assign gene names
colnames(TOM) <- colnames(datExpr)
rownames(TOM) <- colnames(datExpr)

cat("TOM matrix dimensions:", dim(TOM), "\n")

# Select top hub genes for visualization (top 400 genes by connectivity)
nSelect <- min(20, ncol(datExpr))

# Calculate connectivity from adjacency
connectivity <- colSums(adjacency) - 1

# Select most connected genes
selectGenes <- order(-connectivity)[1:nSelect]

# Select TOM for these genes
selectTOM <- TOM[selectGenes, selectGenes]
selectColors <- moduleColors[selectGenes]

# Create dissimilarity TOM
dissTOM <- 1 - selectTOM

# Cluster genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot TOM heatmap
svg(paste0(output_dir, "11_TOM_network_heatmap.svg"), width = 12, height = 12)

TOMplot(
  dissTOM,
  geneTree,
  selectColors,
  main = "Network Heatmap (Top 20 Hub Genes)"
)

dev.off()
cat("Saved: 11_TOM_network_heatmap.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 18: TOM-based Network Heatmap (Publication Quality)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 18: TOM Network Heatmap (Top Hub Genes)\n")
cat("=============================================================================\n\n")

# Option A: Top 20 hub genes overall (by connectivity)
nSelect <- 20

# Calculate connectivity from adjacency
connectivity <- colSums(adjacency) - 1

# Select most connected genes
selectGenes <- order(-connectivity)[1:nSelect]

# Select TOM for these genes
selectTOM <- TOM[selectGenes, selectGenes]
selectColors <- moduleColors[selectGenes]
selectGeneNames <- colnames(datExpr)[selectGenes]

# Create dissimilarity TOM
dissTOM <- 1 - selectTOM

# Assign gene names for labeling
colnames(dissTOM) <- selectGeneNames
rownames(dissTOM) <- selectGeneNames

# Cluster genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot TOM heatmap
svg(paste0(output_dir, "11_TOM_network_heatmap_top20.svg"), width = 10, height = 10)

TOMplot(
  dissTOM,
  geneTree,
  selectColors,
  main = "Co-expression Network (Top 20 Hub Genes)"
)

dev.off()
cat("Saved: 11_TOM_network_heatmap_top20.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 18B: TOM Heatmap for Key Module Hub Genes (WITH GENE NAMES)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 18B: TOM Network Heatmap (Key Module Hubs)\n")
cat("=============================================================================\n\n")

# Get top 10 hub genes from turquoise and blue modules
key_modules_for_plot <- c("turquoise", "blue")
n_per_module <- 10

hub_indices <- c()

for (mod in key_modules_for_plot) {
  mod_genes <- which(moduleColors == mod)
  mod_connectivity <- connectivity[mod_genes]
  top_in_mod <- mod_genes[order(-mod_connectivity)[1:n_per_module]]
  hub_indices <- c(hub_indices, top_in_mod)
}

# Subset TOM
selectTOM <- TOM[hub_indices, hub_indices]
selectColors <- moduleColors[hub_indices]
selectGeneNames <- colnames(datExpr)[hub_indices]

# TOM dissimilarity
dissTOM <- 1 - selectTOM
colnames(dissTOM) <- selectGeneNames
rownames(dissTOM) <- selectGeneNames

# Cluster genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot
svg(
  paste0(output_dir, "11_TOM_network_heatmap_key_modules_labeled.svg"),
  width = 12,
  height = 12
)

TOMplot(
  dissTOM,
  geneTree,
  selectColors,
  geneLabels = selectGeneNames,   # <<< THIS IS THE KEY LINE
  main = "Co-expression Network\n(Top 10 Hub Genes from Turquoise & Blue Modules)"
)

dev.off()

cat("Saved: 11_TOM_network_heatmap_key_modules_labeled.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 19: Module Eigengene Network
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 19: Module Eigengene Network\n")
cat("=============================================================================\n\n")

# Calculate eigengene adjacency
MET <- orderMEs(cbind(MEs, datTraits))

svg(paste0(output_dir, "12_eigengene_network.svg"), width = 12, height = 10)
par(mfrow = c(1, 2))

# Dendrogram
plotEigengeneNetworks(
  MET,
  "Eigengene Dendrogram",
  marDendro = c(0, 4, 2, 0),
  plotHeatmaps = FALSE
)

# Heatmap
par(cex = 0.8)
plotEigengeneNetworks(
  MET,
  "Eigengene Adjacency Heatmap",
  marHeatmap = c(6, 6, 2, 2),
  plotDendrograms = FALSE,
  xLabelsAngle = 90
)

dev.off()
cat("Saved: 12_eigengene_network.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 20: Intramodular Connectivity
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 20: Intramodular Connectivity\n")
cat("=============================================================================\n\n")

# Calculate intramodular connectivity
kIN <- intramodularConnectivity(ADJ, moduleColors)

# Add to gene info
geneInfo$kWithin <- kIN$kWithin[match(geneInfo$Gene, colnames(datExpr))]
geneInfo$kOut <- kIN$kOut[match(geneInfo$Gene, colnames(datExpr))]
geneInfo$kTotal <- kIN$kTotal[match(geneInfo$Gene, colnames(datExpr))]

# Save updated gene info
write.csv(geneInfo, paste0(output_dir, "06_gene_info_complete.csv"), row.names = FALSE)
cat("Updated 06_gene_info_complete.csv with connectivity measures\n\n")

# -----------------------------------------------------------------------------
# SECTION 21: Export Network for Cytoscape
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 21: Export Network for Cytoscape\n")
cat("=============================================================================\n\n")

# Function to export module network
exportModuleNetwork <- function(module_name, geneInfo, TOM, moduleColors, threshold = 0.1, max_genes = 50) {

  # Get genes in module
  modGenes <- which(moduleColors == module_name)

  # Sort by module membership
  mm_col <- paste0("MM.", module_name)
  if (!mm_col %in% colnames(geneInfo)) return(NULL)

  module_geneInfo <- geneInfo[geneInfo$Module == module_name, ]
  module_geneInfo <- module_geneInfo[order(-abs(module_geneInfo[, mm_col])), ]

  # Get top genes
  topGenes <- head(module_geneInfo$Gene, max_genes)
  topGeneIdx <- match(topGenes, colnames(datExpr))
  topGeneIdx <- topGeneIdx[!is.na(topGeneIdx)]

  if (length(topGeneIdx) < 2) return(NULL)

  # Get TOM subset
  modTOM <- TOM[topGeneIdx, topGeneIdx]
  colnames(modTOM) <- colnames(datExpr)[topGeneIdx]
  rownames(modTOM) <- colnames(datExpr)[topGeneIdx]

  # Create edge list
  edges <- data.frame()
  for (i in 1:(length(topGeneIdx)-1)) {
    for (j in (i+1):length(topGeneIdx)) {
      weight <- modTOM[i, j]
      if (weight > threshold) {
        edges <- rbind(edges, data.frame(
          Source = rownames(modTOM)[i],
          Target = rownames(modTOM)[j],
          Weight = weight,
          Module = module_name
        ))
      }
    }
  }

  # Create node attributes
  nodes <- data.frame(
    Gene = topGenes,
    Module = module_name,
    GS.ClinicalOutcome = geneInfo$GS.ClinicalOutcome[match(topGenes, geneInfo$Gene)],
    MM = geneInfo[match(topGenes, geneInfo$Gene), mm_col],
    kWithin = geneInfo$kWithin[match(topGenes, geneInfo$Gene)]
  )

  return(list(edges = edges, nodes = nodes))
}

# Export networks for key modules (ClinicalOutcome-significant, from Section 17)
for (mod in key_modules) {
  if (mod %in% unique(moduleColors)) {
    result <- exportModuleNetwork(mod, geneInfo, TOM, moduleColors, threshold = 0.15, max_genes = 50)
    if (!is.null(result)) {
      write.csv(result$edges, paste0(output_dir, "13_cytoscape_", mod, "_edges.csv"), row.names = FALSE)
      write.csv(result$nodes, paste0(output_dir, "13_cytoscape_", mod, "_nodes.csv"), row.names = FALSE)
      cat("Saved: 13_cytoscape_", mod, "_edges.csv and nodes.csv\n")
    }
  }
}
cat("\n")

# -----------------------------------------------------------------------------
# SECTION 22: Module Size Distribution
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 22: Module Size Distribution\n")
cat("=============================================================================\n\n")

svg(paste0(output_dir, "14_module_sizes.svg"), width = 12, height = 6)
par(mar = c(8, 4, 4, 2))

module_table <- sort(table(moduleColors), decreasing = TRUE)

ymax <- max(module_table)
pad  <- ceiling(0.08 * ymax)   # extra headroom so top label is visible

bp <- barplot(
  module_table,
  col = names(module_table),
  las = 2,
  main = "Number of Genes per Module",
  ylab = "Number of Genes",
  cex.names = 0.8,
  ylim = c(0, ymax + pad)      # ensures turquoise + label fits within y-axis
)

# Add count labels at correct x positions (bar midpoints)
text(
  x = bp,
  y = as.numeric(module_table) + ceiling(0.02 * ymax),
  labels = as.numeric(module_table),
  cex = 0.7,
  pos = 3
)

dev.off()
cat("Saved: 14_module_sizes.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 23: Module Eigengene Expression by Outcome (ClinicalOutcome)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 23: Module Eigengene by Outcome (ClinicalOutcome)\n")
cat("=============================================================================\n\n")

# Create boxplots of eigengene expression by ClinicalOutcome
svg(paste0(output_dir, "15_eigengene_by_outcome.svg"), width = 14, height = 10)

# Get modules significant for ClinicalOutcome (p < 0.05)
sig_modules <- names(MEs)[moduleTraitPvalue[, "ClinicalOutcome"] < 0.05]
sig_modules <- sig_modules[1:min(8, length(sig_modules))]

par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))

for (me in sig_modules) {
  module_name <- gsub("ME", "", me)

  # Create grouping factor based on ClinicalOutcome (1 = 5Plus, 0 = 2Minus; see Section 6)
  # Levels ordered explicitly (not left to alphabetical sort) so col=c(green,red) below
  # maps to (5Plus, 2Minus) -- assumes 5Plus is the favorable/longer-survival category,
  # consistent with this project's green=good/red=poor outcome color convention
  outcome <- factor(ifelse(datTraits$ClinicalOutcome == 1, "5Plus", "2Minus"),
                     levels = c("5Plus", "2Minus"))

  boxplot(
    MEs[, me] ~ outcome,
    main = paste0(module_name, " Module"),
    ylab = "Eigengene Expression",
    # col = c("#E41A1C", "#4DAF4A"),
    col = c("#4DAF4A","#E41A1C"),
    xlab = "ClinicalOutcome"
  )

  # Add points
  stripchart(
    MEs[, me] ~ outcome,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 19,
    col = "black",
    cex = 0.8
  )

  # Add p-value for ClinicalOutcome
  pval <- moduleTraitPvalue[me, "ClinicalOutcome"]
  mtext(paste0("p = ", signif(pval, 2)), side = 3, line = 0, cex = 0.7)
}

dev.off()
cat("Saved: 15_eigengene_by_outcome.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 24: Hub Gene Summary Table (ClinicalOutcome)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 24: Hub Gene Summary (ClinicalOutcome)\n")
cat("=============================================================================\n\n")

# Create comprehensive hub gene table for significant modules
hub_summary <- data.frame()

sig_module_colors <- gsub("ME", "", sig_modules)

for (module in sig_module_colors) {
  if (!module %in% unique(moduleColors)) next

  mm_col <- paste0("MM.", module)
  if (!mm_col %in% colnames(geneInfo)) next

  module_genes <- geneInfo[geneInfo$Module == module, ]
  module_genes <- module_genes[order(-abs(module_genes[, mm_col])), ]

  top_genes <- head(module_genes, 20)

  hub_entry <- data.frame(
    Module = module,
    Rank = 1:nrow(top_genes),
    Gene = top_genes$Gene,
    GS.ClinicalOutcome = round(top_genes$GS.ClinicalOutcome, 4),
    p.GS.ClinicalOutcome = signif(top_genes$p.GS.ClinicalOutcome, 3),
    ModuleMembership = round(top_genes[, mm_col], 4),
    kWithin = round(top_genes$kWithin, 2)
  )

  hub_summary <- rbind(hub_summary, hub_entry)
}

write.csv(hub_summary, paste0(output_dir, "16_hub_genes_summary_clinicaloutcome.csv"), row.names = FALSE)
cat("Saved: 16_hub_genes_summary_clinicaloutcome.csv\n\n")

# -----------------------------------------------------------------------------
# SECTION 25: Gene Significance Volcano Plot (ClinicalOutcome)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 25: Gene Significance Volcano Plot (ClinicalOutcome)\n")
cat("=============================================================================\n\n")

svg(paste0(output_dir, "17_GS_volcano_plot.svg"), width = 12, height = 8)

# Prepare data
gs_values <- geneInfo$GS.ClinicalOutcome
pvalues <- geneInfo$p.GS.ClinicalOutcome
log_pval <- -log10(pvalues + 1e-10)

# Plot
plot(
  gs_values, log_pval,
  pch = 19,
  col = adjustcolor(geneInfo$Module, alpha.f = 0.5),
  xlab = "Gene Significance (ClinicalOutcome)",
  ylab = "-log10(p-value)",
  main = "Gene Significance Volcano Plot (ClinicalOutcome)\n(Colored by Module)",
  cex = 0.6
)

# Add significance lines
abline(h = -log10(0.05), col = "red", lty = 2, lwd = 1.5)
abline(v = c(-0.3, 0.3), col = "gray", lty = 2)

# Label top genes
top_genes <- geneInfo[order(-abs(geneInfo$GS.ClinicalOutcome)), ][1:10, ]
text(
  top_genes$GS.ClinicalOutcome,
  -log10(top_genes$p.GS.ClinicalOutcome + 1e-10),
  labels = substr(top_genes$Gene, 1, 12),
  cex = 0.6,
  pos = 3
)

# Add legend
legend(
  "topright",
  legend = c("p < 0.05", "|GS| > 0.3"),
  lty = 2,
  col = c("red", "gray"),
  cex = 0.8
)

dev.off()
cat("Saved: 17_GS_volcano_plot.svg\n\n")

# -----------------------------------------------------------------------------
# SECTION 26: Module Correlation Heatmap
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 26: Module Correlation Heatmap\n")
cat("=============================================================================\n\n")

svg(paste0(output_dir, "18_module_correlation_heatmap.svg"), width = 10, height = 10)

# Calculate module-module correlations
ME_cor <- cor(MEs, use = "pairwise.complete.obs")

# Create heatmap
labeledHeatmap(
  Matrix = ME_cor,
  xLabels = names(MEs),
  yLabels = names(MEs),
  xSymbols = gsub("ME", "", names(MEs)),
  ySymbols = gsub("ME", "", names(MEs)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  setStdMargins = TRUE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Module-Module Correlation"
)

dev.off()
cat("Saved: 18_module_correlation_heatmap.svg\n\n")

# Significant-only version: restrict both axes to modules whose correlation
# with ClinicalOutcome (the project's primary outcome) is p < 0.05, so
# significant modules are shown correlating only with other significant
# modules. Grey (unassigned genes) is always excluded.
sig_ME <- rownames(moduleTraitPvalue)[moduleTraitPvalue[, "ClinicalOutcome"] < 0.05]
sig_ME <- setdiff(sig_ME, "MEgrey")
sig_ME <- sig_ME[sig_ME %in% colnames(ME_cor)]

if (length(sig_ME) >= 2) {
  ME_cor_sig <- ME_cor[sig_ME, sig_ME]

  svg(paste0(output_dir, "18b_module_correlation_heatmap_significant.svg"), width = 10, height = 10)

  labeledHeatmap(
    Matrix = ME_cor_sig,
    xLabels = sig_ME,
    yLabels = sig_ME,
    xSymbols = gsub("ME", "", sig_ME),
    ySymbols = gsub("ME", "", sig_ME),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    setStdMargins = TRUE,
    cex.text = 0.6,
    zlim = c(-1, 1),
    main = "Module-Module Correlation\n(ClinicalOutcome-significant modules, p < 0.05)"
  )

  dev.off()
  cat("Saved: 18b_module_correlation_heatmap_significant.svg\n\n")
} else {
  cat("Skipped 18b_module_correlation_heatmap_significant.svg: fewer than 2 ClinicalOutcome-significant modules\n\n")
}

# -----------------------------------------------------------------------------
# SECTION 27: Save All R Objects (Final)
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("SECTION 27: Save All R Objects (Final)\n")
cat("=============================================================================\n\n")

save(
  datExpr,
  datTraits,
  metadata,
  net,
  TOM,
  MEs,
  moduleColors,
  moduleLabels,
  moduleTraitCor,
  moduleTraitPvalue,
  geneModuleMembership,
  geneTraitSig,
  geneInfo,
  softPower,
  ADJ,
  kIN,
  file = paste0(output_dir, "WGCNA_results_complete.RData")
)

cat("Saved: WGCNA_results_complete.RData (with TOM and connectivity)\n\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=============================================================================\n\n")

cat("All analyses (gene significance, module membership, hub genes, key module selection) performed with respect to ClinicalOutcome only\n\n")

cat("OUTPUT FILES:\n")
cat("  EIGENGENE VISUALIZATION:\n")
cat("    08_module_eigengenes_barplot.svg    - Eigengene expression by sample\n")
cat("    09_module_eigengene_heatmap.svg     - Eigengene relationships\n")
cat("\n  NETWORK VISUALIZATIONS:\n")
cat("    10_MM_vs_GS_scatter.svg             - Module membership vs gene significance\n")
cat("    11_TOM_network_heatmap.svg          - TOM-based network heatmap\n")
cat("    12_eigengene_network.svg            - Module eigengene network\n")
cat("    13_cytoscape_*_edges/nodes.csv      - Cytoscape export files\n")
cat("    14_module_sizes.svg                 - Module size distribution\n")
cat("    15_eigengene_by_outcome.svg         - Eigengene boxplots by ClinicalOutcome\n")
cat("    17_GS_volcano_plot.svg              - Gene significance volcano plot\n")
cat("    18_module_correlation_heatmap.svg   - Module-module correlations (all modules)\n")
cat("    18b_module_correlation_heatmap_significant.svg - Module-module correlations (ClinicalOutcome-significant modules only)\n")
cat("\n  DATA FILES:\n")
cat("    06_gene_info_complete.csv                      - Gene info with connectivity (ClinicalOutcome)\n")
cat("    16_hub_genes_summary_clinicaloutcome.csv       - Hub gene summary table\n")
cat("    WGCNA_results_complete.RData                   - All R objects\n")

cat("\nNEXT STEPS:\n")
cat("  1. Use 13_cytoscape_*.csv files in Cytoscape for network visualization\n")
cat("  2. Perform GO/KEGG enrichment on module genes\n")
cat("  3. Validate hub genes experimentally\n")

cat("\n=============================================================================\n")
