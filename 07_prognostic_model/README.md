# Prognostic model — three-level hierarchical Bayesian shrinkage

Supplementary Methods, section 8.

The bulk RNA-seq dataset (n = 20) is the internal discovery cohort, categorized as poor (< 2 year)
and good (> 5 year) outcome, and is used for model fitting. An independent 54-patient Danish bulk
RNA-seq cohort is the survival validation cohort; a second cohort (N = 54, MD Anderson CCR2020)
provides a predefined molecular subtype label for external non-survival validation; a third
external cohort (Frerich 2018) tests whether the panel distinguishes the high-risk subset. The
single-cell dataset (N = 24, 95,441 cells) supports candidate gene nomination and cell
composition deconvolution.

**No expression or outcome data from external cohorts contributed to gene selection or model
fitting at any stage.**

## The result

A locked four-gene panel — **DSCAM, ODC1, NCAPG, CCNB2** — scored as the unweighted mean of
sign-weighted rank-inverse-normal expression. All four coefficients are positive, so higher means
worse prognosis. The hierarchical Bayesian shrinkage procedure selected the genes and their
directions; it does **not** weight the deployed score, and no posterior quantity enters any
reported C-index or hazard ratio.

## Start here

| File | What it gives you |
|---|---|
| `release/model_spec.json` | The model specification, with the frozen discovery reference and locked threshold |
| `release/score_acc4.py` | The scoring function (cohort-batch and frozen single-sample) |
| `release/results_dictionary.csv` | Every reported number, with its value, uncertainty, source file, generating script and evidence tier |
| `release/cohorts.csv` | Cohorts, endpoints, event counts and roles |
| `release/MANIFEST.json` | sha256 fingerprints of every raw and harmonized input |
| `release/environment.txt`, `release/requirements-frozen.txt` | The environment used to generate the results |
| `release/build_release.py`, `release/build_results_dictionary.py` | Rebuild the release files from the results on disk |

## Core pipeline

| Script | Methods element |
|---|---|
| `00_dk_survival.py` | DK overall-survival table used by the figure scripts |
| `03_data.py` | Data assembly, rank-inverse-normal transform, evidence prior weights |
| `modules.py` | The 52 biological modules (curated programs, WGCNA modules, single-cell DE, evidence sets) forming the middle level of the hierarchy |
| `mlbayes.py` | The model: three-level shrinkage (source → module → gene) fitted by ADVI |
| `common.py` | Input fingerprinting |
| `20_pseudobulk.py` | Single-cell pseudobulk aggregation per cell type and patient, used for deconvolution |
| `11_lock_and_validate.py` | Priors tuned on the discovery cohort alone by leave-one-out CV, refit on all 20 patients, then applied to the validation cohort |
| `13_compact_panel.py` | Panel size selection by the parsimony rule — leave-one-out across three seeds plus repeated five-fold CV, taking the smallest k within 0.02 AUC of the best-performing size (k = 4) |
| `14_stability_check.py` | Gene stability, bootstrap confidence intervals, cut-point sensitivity |
| `16_panel8_lock.py` | 8-gene variant and its leave-one-gene-out robustness test; also writes the locked 4-gene DK score used by later scripts |

Candidate genes were nominated by integrating bulk differential expression (PyDESeq2), single-cell
pseudobulk differential expression in differentially enriched cell types, and WGCNA hub genes from
outcome-associated modules — all computed on internal cohorts — and retained if supported by at
least two of the three sources.

## Validation and diagnostics

| Script | Methods element |
|---|---|
| `23_ccr2020_external.py` | External subtype validation: ACC-I vs ACC-II in the CCR2020 cohort |
| `40_frerich_external.py` | Cross-platform replication of the subtype result in Frerich 2018 |
| `18_brayer_benchmark.py`, `18b_brayer_fairness.py`, `19_recover_subgroup.py`, `30_brayer_published_group.py` | Benchmarking against a published aggressive ACC classifier (Brayer et al. 2023) |
| `24_delong.py` | DeLong's method for correlated ROC curves, comparing the panel against that published classifier in the same patients |
| `28_matched_nulls.py` | 10,000 abundance- and variance-matched random gene panels and the 2,000-permutation null distribution |
| `26_cox_diagnostics.py` | Harrell's C-index and per standard deviation Cox hazard ratios at 5 years and over full follow-up; proportional-hazards checks; bootstrap optimism correction (Harrell's method, b = 1000) |
| `27_threshold.py` | Decision threshold derived solely from the discovery cohort; frozen single-sample scoring |
| `25_nested_cv.py`, `25b_nested_fixedk.py`, `25c_fidelity.py` | Nested cross-validation of the whole procedure, with outcome-derived inputs rebuilt inside each fold |
| `31b_internal_oof.py`, `31e_internal_stability_full.py` | Internal-cohort out-of-fold predictions and run-to-run distribution of the nested-CV estimate |
| `33_vi_diagnostics.py`, `33b_comparator.py` | Variational inference diagnostics and a penalized-regression comparator |
| `34_axis_stability.py` | Stability of the underlying expression axis |
| `21_signature_and_gate.py` | Deconvolution signature matrix and leave-one-patient-out validation gate |
| `17_subtype_discovery.py` | Data-defined DK subgroup discovery |
| `29_joint_contribution.py` | Joint contribution of panel genes |
| `36_unified_metrics.py` | Unified metric tables: AUC and C-index wherever each is computable |

## Deconvolution

| Script | Methods element |
|---|---|
| `22_dk_composition.py` | Bulk cell type composition in the validation cohort |
| `35_dwls_crosscheck.py`, `35r_dwls_reference.R` | Cross-check of the DWLS implementation against the published DWLS package, including leave-one-patient-out synthetic mixtures |
| `32_deconv_figures.py` | Deconvolution figures |

## Figures

`15_final_figures.py`, `31f_stability_figure.py`, `36c_metric_figures.py`,
`52_manuscript_figures.py`.

### Regenerating the manuscript tables and figures

From this directory. Nothing below refits the model — every stage reads the locked scores:

```bash
python 36_unified_metrics.py      # unified metrics
python 36c_metric_figures.py      # figures M1-M3, S1
python 32_deconv_figures.py       # deconvolution figures M4, S2, S3
python 31f_stability_figure.py    # internal-stability figure S4
python 52_manuscript_figures.py   # Figure 1 and Figure S1
```

Run them in that order: `36c_` reads what `36_` writes. Resulting numbers are indexed in
`release/results_dictionary.csv`.
