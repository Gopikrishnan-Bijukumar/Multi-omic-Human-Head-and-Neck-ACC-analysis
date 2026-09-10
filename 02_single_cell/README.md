# Single-cell RNA-seq

Supplementary Methods, Single cell sections 1–3. Cohort: 24 samples (19 overlapping with the bulk
cohort) spanning 9 distinct head and neck regions; 13 good and 11 poor outcome.

## Run order

| Script | What it does | Key parameters |
|---|---|---|
| `01_preprocessing_concat.ipynb` | Loads CellRanger output per sample and concatenates. Filtering is done independently per sample | >= 200 detected genes per cell; genes expressed in >= 3 cells; < 10% mitochondrial counts |
| `02_doubletdetection_scrublet.ipynb` | DoubletDetection and Scrublet doublet calls | — |
| `03_solo_doublet_calling.ipynb` | SOLO doublet calls on the DoubletDetection- and Scrublet-completed dataset | — |
| `04_doublet_consensus_removal.ipynb` | Removes cells flagged by 2 of the 3 algorithms | consensus 2/3 |
| `05_top95_count_marking.ipynb` | Marks cells above the 95th percentile of total counts as potential doublets or outliers | 95th percentile |
| `06_final_filtering.ipynb` | Applies the final filters before annotation | — |
| `07_normalization_hvg_scvi.ipynb` | Normalisation, log1p, HVG selection, cell cycle scoring and regression, and scVI batch correction | 10,000 HVGs (`seurat_v3`); S and G2M scores regressed out; scVI 30-PC latent space, 300 epochs, validation every 10 epochs |
| `08_first_pass_annotation.ipynb` | Neighbourhood graph, UMAP and Leiden clustering on the `X_scVI` latent representation, then first-pass manual annotation | Leiden resolution 0.6 |
| `09_marking_cells_for_removal.ipynb` | Marks remaining low-quality clusters | — |
| `10_subcluster_and_annotation.ipynb` | Subclustering and final manual cell-type annotation across the 24 samples | — |
| `11_pseudobulk_pydeseq2.ipynb` | Sample-level pseudobulk differential expression: expression summed across all cells per sample (`adata.obs['batch']`), PyDESeq2 with the design based on clinical outcome | padj < 0.05, abs(log2FC) > 0.5 |
| `12_pseudobulk_gsea.ipynb` | GSEA on the pseudobulk differential expression result | — |

`plotting_utils.py` is a small UMAP plotting helper imported by several notebooks.

## `figures/`

UMAP and PCA views by clinical outcome, a rasterised UMAP for publication, marker gene and
differentially expressed gene heatmaps, and cell-type composition barplots.
