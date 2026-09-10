# Visium, 55 µm

4 samples, 2 good and 2 poor outcome.

| Script | What it does |
|---|---|
| `01_preprocessing.ipynb` | Pre-processing exploration for the 55 µm cohort |
| `02_normalization_clustering.ipynb` | Normalisation and clustering exploration |
| `03_giotto_pipeline.R` | The full per-sample Giotto pipeline (methods 7a–7e) |

## Running the pipeline

```bash
Rscript 03_giotto_pipeline.R P01
```

One invocation processes one sample. The four original per-sample scripts were identical apart
from three Section 1 values, which now come from `config/spatial_visium_samples.csv`:

| Column | Meaning |
|---|---|
| `module_k` | Metagene module count, chosen per sample from the Section 6b silhouette/WSS sweep (5 / 6 / 4 / 4) |
| `svg_plot_midpoint` | Display-only gradient midpoint on the top-SVG plots; `NA` reproduces the one sample finalised without it |
| `outcome` | Clinical outcome group |

Substituting a sample's values back into this script reproduces the corresponding original
per-sample script exactly.
