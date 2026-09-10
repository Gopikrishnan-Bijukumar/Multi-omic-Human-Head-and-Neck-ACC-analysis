# Prognostic model — three-level hierarchical Bayesian shrinkage

Supplementary Methods, section 8.

The bulk RNA-seq dataset (n = 20) is the internal discovery cohort, categorised as poor (< 2 year)
and good (> 5 year) outcome, and is used for model fitting. An independent 54-patient Danish bulk
RNA-seq cohort is the survival validation cohort; a second cohort (N = 54, MD Anderson CCR2020)
provides a predefined molecular subtype label for external non-survival validation; a third
external cohort (Frerich 2018) checks whether the panel distinguishes the high-risk subset. The
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
| `release/model_spec.json` | The single definitive model specification, with the frozen discovery reference and locked threshold |
| `release/results_dictionary.csv` | Every manuscript-facing number, each with its value, uncertainty, source file, generating script and evidence tier |
| `release/DEPRECATED.md` | The do-not-quote map: every superseded artefact and the wrong number it would produce |
| `release/score_acc4.py` | The one scoring function (cohort-batch and frozen single-sample), verified against the stored scores |
| `release/MANIFEST.json` | sha256 fingerprints of every raw and harmonised input |
| `release/environment.txt`, `release/requirements-frozen.txt` | The exact environment behind every published number |

Script filenames are kept exactly as they were when the results were generated, because
`results_dictionary.csv` cites each number's generating script by name. The numbering reflects
the order stages were run, not a clean pipeline order.

## Core pipeline

| Script | Methods element |
|---|---|
| `03_data.py` | Data assembly, rank-inverse-normal transform, evidence prior weights |
| `modules.py` | The 52 biological modules (curated programmes, WGCNA modules, single-cell DE, evidence sets) forming the middle level of the hierarchy |
| `mlbayes.py` | The model: three-level shrinkage (source → module → gene) fitted by ADVI |
| `bayes.py`, `common.py` | Supporting inference and input-fingerprinting utilities |
| `20_pseudobulk.py` | Single-cell pseudobulk differential expression in differentially enriched cell types, one of the three nomination sources |
| `11_lock_and_validate.py` | **The pipeline that produces the final result.** Priors tuned on the discovery cohort alone by leave-one-out CV, refit on all 20 patients, then applied once to the validation cohort |
| `13_compact_panel.py` | Panel size selection by the parsimony rule — leave-one-out across three seeds plus repeated five-fold CV, taking the smallest k within 0.02 AUC of the best-performing size. This gave k = 4 |
| `14_stability_check.py` | Gene stability, bootstrap confidence intervals, cut-point sensitivity |
| `16_panel8_lock.py` | The 8-gene variant, retained as a reported negative result: it failed its own pre-specified leave-one-gene-out robustness gate and adds nothing over 4 genes |

Candidate genes were nominated by integrating bulk differential expression (PyDESeq2), single-cell
pseudobulk differential expression in differentially enriched cell types, and WGCNA hub genes from
outcome-associated modules — all computed on internal cohorts — and retained if supported by at
least two of the three sources.

## Validation and diagnostics

| Script | Methods element |
|---|---|
| `23_ccr2020_external.py` | External subtype validation: ACC-I vs ACC-II in the CCR2020 cohort |
| `40_frerich_external.py` | Cross-platform replication of the subtype result in Frerich 2018 |
| `18_brayer_benchmark.py`, `18b_brayer_fairness.py`, `30_brayer_published_group.py` | Benchmarking against a published aggressive ACC classifier |
| `24_delong.py` | DeLong's method for correlated ROC curves, comparing the panel against that published classifier in the same patients |
| `28_matched_nulls.py` | 10,000 abundance- and variance-matched random gene panels and the 2,000-permutation null distribution |
| `26_cox_diagnostics.py` | Harrell's C-index and per standard deviation Cox hazard ratios at 5 years and over full follow-up; proportional-hazards checks |
| `27_threshold.py` | Decision threshold derived solely from the discovery cohort |
| `25_nested_cv.py`, `25b_nested_fixedk.py`, `25c_fidelity.py` | Nested cross-validation of the whole procedure |
| `31b_internal_oof.py`, `31d/31e_internal_stability*.py` | Internal performance re-monitored for information leakage. Several quantities had been computed once from the 20-patient data rather than within each fold; correcting this lowered the internal AUC without affecting external validation |
| `33_vi_diagnostics.py`, `33b_comparator.py` | Variational inference convergence diagnostics |
| `34_axis_stability.py` | Stability of the underlying expression axis |
| `21_signature_and_gate.py` | Reference single-sample scoring, mapping each patient independently into the discovery cohort's expression pattern rather than relative to the validation batch |
| `17_subtype_discovery.py`, `19_recover_subgroup.py` | Subtype and high-risk subgroup recovery |
| `29_joint_contribution.py` | Joint contribution of panel genes |
| `31_unified_metrics.py`, `36_unified_metrics_v2.py` | Unified metric tables (v2 supersedes v1 for M1–M3 and S1; see `release/DEPRECATED.md`) |
| `00_baseline_dk.py`, `01_explore.py`, `02_feasibility.py`, `04_cv.py`, `05_biology.py`, `06_mlcv.py`, `07_denoise_combine.py`, `08_procedure_cv.py`, `09_final.py`, `10_endpoint_matched.py` | Early diagnostics kept for the record. Several are deliberately exploratory because they learn from validation-cohort survival; `release/DEPRECATED.md` says which |

Bootstrap optimism correction (Harrell's method, b = 1000) was applied to the multivariable model.

## Deconvolution

| Script | Methods element |
|---|---|
| `22_dk_composition.py` | Bulk cell type composition in the validation cohort |
| `35_dwls_crosscheck.py`, `35r_dwls_reference.R` | DWLS-based deconvolution against the internal 24-sample single-cell dataset, by leave-one-patient-out cross-validation against known single-cell composition before application to the validation cohort |
| `32_deconv_figures.py` | Deconvolution figures |

## Figures

`12_compare_figure.py`, `15_final_figures.py`, `31c_metric_figures.py`,
`31f_stability_figure.py`, `36c_metric_figures_v2.py`, `52_manuscript_figures_v13.py`.

### Regenerating the manuscript tables and figures

From this directory. Nothing below refits the model — every stage reads the locked scores:

```bash
python 36_unified_metrics_v2.py   # unified metrics (v2, sourced from 31e)
python 36c_metric_figures_v2.py   # manuscript figures M1-M3, S1
python 32_deconv_figures.py       # deconvolution figures M4, S2, S3
python 31f_stability_figure.py    # internal-stability figure S4
```

Run them in that order: `36c_` reads what `36_` writes. Resulting numbers are indexed in
`release/results_dictionary.csv`; `release/DEPRECATED.md` says which older outputs they
supersede.

The manuscript prose builders from the source project (`38_`–`53_build_manuscript_v*.py`) are not
included: they generate document text and PDFs rather than analysis results.
