# Visium HD, 16 µm

2 patients, one poor and one good outcome. 2 µm² bins aggregated to 16 µm².

| Script | What it does |
|---|---|
| `convert_tissue_positions.py` | Converts the 10x `tissue_positions.parquet` into the CSV form `createGiottoVisiumObject` expects. Run before the Giotto pipeline |

The HD analysis follows the same steps as `../visium_55um/03_giotto_pipeline.R` (methods 7a–7e),
with the 16 µm QC thresholds (>= 50 genes per bin, genes in >= 50 bins) and memory adaptations
for HD scale: `normalizeGiotto(scale_feats = FALSE, scale_cells = FALSE)` on both the spatial and
single-cell objects, gene-chunked `calculateHVF` and `binSpect` with global FDR applied after
chunking, and a streaming exporter for the counts table.
