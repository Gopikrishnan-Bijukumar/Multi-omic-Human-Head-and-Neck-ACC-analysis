# Bulk RNA-seq

Supplementary Methods, Bulk RNA-seq sections 1–5. Cohort: 20 tumour samples from 10 anatomical
sites, stratified by survival duration (poor < 2 years, favourable > 5 years; 10 per group).

## Run order

| Script | What it does | Key parameters |
|---|---|---|
| `01_counts_to_h5ad.ipynb` | Converts raw counts to H5AD and annotates with sample metadata | — |
| `02_pydeseq2_differential_expression.ipynb` | PyDESeq2 differential expression between outcome cohorts, design factors from the clinical metadata, contrast via `DeseqStats` | padj < 0.05, abs(log2FC) > 0.5 |
| `03_merge_deg_into_h5ad.ipynb` | Attaches the DE results back onto the H5AD `var` | — |
| `04_gsea_prerank.ipynb` | GSEApy `prerank` on the PyDESeq2 `stat` statistic against locally stored MSigDB 2025 collections (KEGG, Reactome, BP, CC, MF, Oncogenic Signatures, Hallmark) | 1000 permutations, seed 6, min gene set size 15 |
| `05_gsea_plots.ipynb` | Filters on Normalized Enrichment Score and plots. Positive NES denotes enrichment in the good-outcome cohort, negative in the poor-outcome cohort | — |
| `06_wgcna.R` | Weighted gene co-expression network analysis: VST normalisation via the DESeq2 dependency, `blockwiseModules`, module–trait Pearson correlation of eigengenes, kME, hub genes, intramodular connectivity, igraph network layout | genes with counts >= 10 in >= 50% of samples; bottom 25% by variance excluded; soft power 14, signed weighted network; hub genes = top 20 per module by kME; TOM >= 0.225 (blue) and >= 0.185 (turquoise); Fruchterman-Reingold, 5000 iterations |
| `07_ora_modules.ipynb` | Over-representation analysis of a co-expression module's genes via `gseapy.enrich` against the same local MSigDB collections. Set `MODULE` to select the module | p <= 0.05; background = the 16,639 genes surviving WGCNA QC |
| `08_ora_plots.ipynb` | Plots significant enriched terms | — |

One outlier (a good-outcome sample) was identified by clustering and `goodSamplesGenes`
inspection and removed before network construction, leaving 19 samples for WGCNA.

## `figures/`

Presentation figures drawing on the outputs above: volcano plots, differentially expressed gene
heatmaps, an ACC marker gene heatmap, PCA with covariate highlighting, and a cohort overview.
