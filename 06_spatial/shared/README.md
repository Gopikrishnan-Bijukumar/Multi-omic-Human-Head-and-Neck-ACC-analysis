# Shared spatial stages

Used by the Visium pipelines after the per-sample Giotto run.

| Script | Methods element |
|---|---|
| `fetch_genesets.R` | Retrieves and caches the MSigDB Hallmark, BP, MF and Reactome collections (2026.1.Hs) used for enrichment |
| `metagene_hallmark_spatial.R` | Metagene module scoring and the functional enrichment of 7d: hypergeometric over-representation against the QC-surviving background, and `fgsea` on genes ranked by correlation with each module's metagene score |
| `plot_enrichment.R` | Enrichment figures |
| `enrichment_labels.R` | Term label tidying for those figures |
| `composition_compare.R` | Cross-sample composition comparison of 7e, including the immune infiltration summary (summed myeloid, T-cell and B/plasma proportions). Descriptive only; no statistical tests |
| `biomarker_spatial_plots.R` | Expression of good- and poor-outcome biomarkers overlaid on tissue |

## Run order

`fetch_genesets.R` is run once (it is the only network step). The per-sample Giotto pipeline is
then run once per sample, followed by the shared stages across all four:

```bash
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16

Rscript fetch_genesets.R

for s in P16 P01 P09 P22; do
  Rscript ../visium_55um/03_giotto_pipeline.R "$s"
done

Rscript plot_enrichment.R              # set PE_SAMPLE=<sample> to plot just one
Rscript metagene_hallmark_spatial.R    # set MHS_SAMPLE=<sample> to plot just one
Rscript composition_compare.R          # cross-sample; needs all four finished
Rscript biomarker_spatial_plots.R
```

The sample list and each sample's parameters come from `config/spatial_visium_samples.csv`.
Section 6e of the per-sample pipeline writes the enrichment CSVs that `plot_enrichment.R` and
`metagene_hallmark_spatial.R` read.

> [!NOTE]
> Run with plain `Rscript` so the system R library is used; Giotto finds its Python environment
> (leidenalg) by itself. Pin the BLAS/OpenMP thread counts as above.
