# CellPhoneDB cell–cell communication

Supplementary Methods, Single cell section 5 (5a–5d). CellPhoneDB 5.0.1, statistical analysis
method with interaction scoring enabled. Expression is taken from the raw count layer of the
single-cell object, normalised to 10,000 transcripts per cell and log1p transformed. Statistical
testing uses 1,500 permutations at p < 0.05.

## Run order

The nine notebooks run in order; each ends in an assertion gate that the next one depends on.

| Notebook | Methods element |
|---|---|
| `00_audit_inputs_and_sources.ipynb` | Input provenance audit before anything is written |
| `01_build_cpdb_inputs.ipynb` | Builds the per-patient and pooled CellPhoneDB inputs, including cell-type retention rules |
| `02_run_cellphonedb.ipynb` | Runs CellPhoneDB globally across all samples and independently per patient (5a) |
| `03_extract_patient_endpoints.ipynb` | Defines directed endpoints: a unique ligand–receptor interaction and sender–receiver cell type pair, evaluable for a patient only when both cell types have >= 10 cells and CellPhoneDB returns a finite, significant interaction score (5a) |
| `04_patient_level_inference.ipynb` | Patient-level differential interaction analysis (5b) |
| `05_pooled_exploratory_analysis.ipynb` | Exploratory pooled-cell analysis (5c) |
| `06_sensitivity_analysis.ipynb` | Down-sampling sensitivity across seeds (5c) |
| `07_figures_and_tables.ipynb` | Figures and tables |
| `08_final_verification.ipynb` | End-to-end verification of every reported quantity |

`utils.py` holds the shared implementation; `config.yaml` (validated against
`config_schema.json`) holds every threshold. `tests/` covers the statistical helpers.
`cpdb_interaction_coverage.py` reports interaction coverage.

## Statistical procedure (5b)

Per-patient interaction scores are compared between outcome groups, testing endpoints only when
evaluable in at least six patients per group. Two-sided Mann-Whitney U tests are applied using a
pre-specified procedure: exact tests for tie-free distributions and the asymptotic test with tie
correction when ties are present. Effect sizes are rank-biserial correlations, positive values
denoting higher scores in the good-outcome cohort, reported with median differences. Uncertainty
uses 10,000 stratified bootstrap re-samples, patients re-sampled within outcome groups, giving
95 percentile confidence intervals. Benjamini-Hochberg correction is applied to all evaluable
endpoints, with FDR < 0.05 considered significant. Fisher's exact tests separately identify
endpoints whose evaluability differed between outcome groups, independently B-H adjusted to flag
selection-sensitive comparisons.

## Pooled analysis (5c) and consistency (5d)

The pooled analysis pools all samples within each outcome group separately, retaining cell types
represented by more than 10 cells in each group and down-sampling to a common depth of 42,214
cells per group. The primary analysis uses seed 0 and is repeated with seeds 1–4. A paired
Wilcoxon signed rank test compares the pooled cohorts across shared cell type pairs of each
ligand–receptor interaction, with FDR computed per pooled cohort at q < 0.05. Because the unit of
replication is cell interactions rather than patients, these results are exploratory.

Per-patient consistency (5d) is descriptive only and not statistically tested, given the small
patient cohort size.
