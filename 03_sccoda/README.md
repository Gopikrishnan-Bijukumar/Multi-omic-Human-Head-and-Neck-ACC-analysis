# scCODA compositional analysis

Supplementary Methods, Single cell section 4. Differential cell-type abundance between the two
clinical outcome groups, using scCODA 0.1.9.

| Script | What it does |
|---|---|
| `01_sccoda_compositional_analysis.ipynb` | Aggregates single-cell annotations into a sample-by-cell-type count matrix and fits the scCODA model with outcome as a categorical covariate and good outcome as the reference level |
| `02_reference_sensitivity.py` | Repeats the analysis using endothelial cells, M2 macrophages and fibroblasts as alternative reference cell types, testing the robustness of the mural-cell-based reference |

Model fitting uses Hamiltonian Monte Carlo with 100,000 posterior samples after 10,000 burn-in
iterations, with random seeds fixed. Effects are assessed by posterior inclusion probability
against scCODA's data-adaptive threshold targeting an estimated false discovery rate < 0.05;
model coefficients, log2 fold changes, posterior standard deviations and 94% highest density
intervals are extracted from the posterior.

Cell-type proportions are additionally compared descriptively with a two-sided Mann-Whitney U
test on sample-level proportions and displayed as box-whisker plots.
