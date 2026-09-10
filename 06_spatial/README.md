# Spatial transcriptomics

Supplementary Methods, Single cell section 7 (7a–7e). Spatial transcriptomics at two resolutions:
Visium (55 µm²) and Visium HD (2 µm², binned to 16 µm²). All samples were profiled using the
Human Transcriptome Probe Set v2.0 (10x Genomics) and processed independently with Giotto in R.

- **`visium_55um/`** — the standard Visium cohort: 4 samples, 2 good and 2 poor outcome.
- **`visium_hd_16um/`** — the HD cohort: 2 patients, one poor and one good outcome. See the note
  in that folder about where the HD pipeline was executed.
- **`shared/`** — enrichment and comparison stages used by both resolutions.

## Pipeline stages (7a–7e)

**7a Pre-processing and normalisation.** Filtered feature barcode matrices and spatial
coordinates are loaded with `createGiottoVisiumObject` and restricted to spots/bins marked as on
tissue. Mitochondrial and ribosomal genes are removed before quality control filtering. Bins are
retained if at least 50 or 200 genes were detected at 16 µm and 55 µm resolution respectively,
and genes if detected in at least 50 bins in HD and 10 bins at standard resolution. Counts are
normalised and log-transformed with `normalizeGiotto`.

**7b Dimensionality reduction, clustering and spatial network.** Highly variable genes via
`calculateHVF`; PCA on those genes; the first 15 PCs used for UMAP and shared nearest neighbour
graph generation (k = 15); Leiden clustering at resolution 0.5; a Delaunay triangulation spatial
network with a minimum k of 2 used for all downstream spatial analyses.

**7c Spatially variable genes and co-expression modules.** `binSpect` (rank-based binarization)
on the Delaunay network with BH FDR correction. The 250 most significant genes feed
`detectSpatialCorFeats`, partitioned into modules by Ward-linkage hierarchical clustering on
gene–gene correlation distance. Module number k is chosen per sample by a k sweep from 2 to 12
using mean silhouette width and within-cluster sum of squares — the resulting per-sample values
live in `config/spatial_visium_samples.csv`.

**7d Functional enrichment.** Metagene modules tested against MSigDB Hallmark, BP, MF and
Reactome collections (2026.1.Hs). Over-representation uses the exact hypergeometric test against
the QC-surviving background with BH correction. Gene-set enrichment uses `fgsea` with
`minSize = 10`, `maxSize = 500` and `nPermSimple = 10000`, on genes pre-ranked by Pearson
correlation with each module's metagene score. To reduce the chance of a module's own genes
topping its ranking, the primary ranking excludes member genes using `excl_members`.

**7e Deconvolution and immune infiltration.** Cell type composition estimated using spatial
dampened weighted least squares (SpatialDWLS) with the single-cell data as reference. Per-Leiden
domain and per-spatial-block compositions are summarised per sample, using blocks of
400 × 400 µm² for HD and a 600 µm × 600 µm hexagonal lattice for standard Visium. Cross-sample
comparisons between outcome groups are descriptive only; no statistical tests were computed.
Immune infiltration is summarised as the summed proportion of myeloid, T-cell and B/plasma cells
(not normalised).

## A note on CytoSPACE

The pipeline directory in the source project is named `giotto_w_cytospace`, and CytoSPACE and
Interaction Changed Features were part of an earlier version. Both were retired: a 16 µm bin
holds roughly 1–3 cells and a 55 µm spot roughly 5–20, and neither supports per-cell-type
interaction attribution. Cell-type composition (7e) replaced it. No CytoSPACE code is included
here.
