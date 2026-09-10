# Visium HD, 16 µm

2 patients, one poor and one good outcome. 2 µm² bins aggregated to 16 µm².

| Script | What it does |
|---|---|
| `convert_tissue_positions.py` | Converts the 10x tissue position table into the CSV form `createGiottoVisiumObject` expects. Must be run before the Giotto pipeline |
| `derive_ncells_per_spot.py` | Derives the number of cells per bin |

## Where the HD Giotto pipeline is

**The R pipeline that produced the published HD results is not in this repository.**

Giotto densifies sparse matrices, and at HD scale that is fatal: after QC the two samples are
155,638 × 16,168 and 138,203 × 13,565 bins × genes, which is 0.68 GB and 0.27 GB as sparse
matrices but 18.75 GiB and 13.97 GiB dense — a 28× and 52× inflation. Straightforward runs were
killed at 220 GB of memory inside `binSpect` and `addStatistics`.

The pipeline was therefore ported to an HPC cluster and modified to avoid densification:
`normalizeGiotto(scale_feats = FALSE, scale_cells = FALSE)` on both the spatial and single-cell
objects, gene-chunked `calculateHVF` and `binSpect` implementations with global FDR applied after
chunking, a streaming exporter for the counts table, thread pinning, and gobject checkpointing.
Chunked implementations were verified against the unchunked ones at `max|diff| = 0` for both
`sd`/`gini` statistics and both `binSpect` bin methods.

That modified pipeline lives on the HPC allocation where it ran, and the local copies in the
source project were the pre-port versions that never completed — their output directories are
empty. Rather than publish code that did not generate the results, this folder documents the
situation. The analysis steps are identical to `../visium_55um/03_giotto_pipeline.R` (methods
7a–7e) apart from the memory adaptations above and the 16 µm QC thresholds (>= 50 genes per bin,
genes in >= 50 bins).
