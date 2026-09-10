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

The per-sample Giotto pipeline is one parameterised script; run it once per sample, then the
shared stages once across all four:

```bash
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16

for s in P16 P01 P09 P22; do
  Rscript ../visium_55um/03_giotto_pipeline.R "$s"
done

Rscript plot_enrichment.R        # set PE_SAMPLE=<sample> to plot just one
Rscript composition_compare.R    # cross-sample; needs all four finished
```

The sample list and each sample's parameters come from `config/spatial_visium_samples.csv`.
Section 6e writes the enrichment CSVs that `plot_enrichment.R` draws from, so it must follow the
per-sample runs; it is cheap and safe to re-run on its own.

> [!NOTE]
> Use plain `Rscript`. The R packages come from `~/R/x86_64-pc-linux-gnu-library/4.4` via
> `~/.Renviron`. Do **not** `conda activate giotto_env` first — that environment is Python-only
> (Giotto finds it by itself for leidenalg) and activating it shadows the system R and its
> library. Pin the BLAS/OpenMP thread counts as above: an unpinned BLAS will oversubscribe
> every core.
