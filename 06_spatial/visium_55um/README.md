# Visium, 55 µm

4 samples, 2 good and 2 poor outcome.

| Script | What it does |
|---|---|
| `01_preprocessing.ipynb` | Pre-processing and QC for the 55 µm cohort |
| `02_normalization_clustering.ipynb` | Normalization and clustering |
| `03_giotto_pipeline.R` | The full per-sample Giotto pipeline (methods 7a–7e) |

## Running the pipeline

```bash
Rscript 03_giotto_pipeline.R P01
```

One invocation processes one sample. Three per-sample values come from
`config/spatial_visium_samples.csv`:

| Column | Meaning |
|---|---|
| `module_k` | Metagene module count, chosen per sample from the Section 6b silhouette/WSS sweep (5 / 6 / 4 / 4) |
| `svg_plot_midpoint` | Display-only gradient midpoint on the top-SVG plots (`NA` = none) |
| `outcome` | Clinical outcome group |
