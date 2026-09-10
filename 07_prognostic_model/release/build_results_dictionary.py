"""Every manuscript-facing number, with its source file, generating script and evidence tier.

Tiers
  confirmatory  the single pre-specified primary claim
  prespecified  declared in a PRESPEC document before the analysis ran
  secondary     supporting analyses reported alongside the primary claim
  exploratory   post-hoc, hypothesis-generating, or conditioned on the outcome
  audit         methodological audit of the study's own procedures
  negative      a pre-declared test that FAILED and is reported as such
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import json, os
import pandas as pd

R = f"{ACC_DATA_ROOT}"
W = f"{R}/work"
REL = f"{W}/release"


def g(path, *keys, default=None):
    try:
        d = json.load(open(f"{W}/{path}"))
        for k in keys:
            d = d[k] if not isinstance(k, int) else d[k]
        return d
    except Exception:
        return default


ROWS = [
 # ---------------- confirmatory ----------------
 ("DK 5-year overall survival, C-index", g("outputs_compact/COMPACT_validation.json","compact_5yr","C"),
  "0.541-0.801 (bootstrap)", "confirmatory", "outputs_compact/COMPACT_validation.json", "13_compact_panel.py"),
 ("DK 5-year, HR per score SD", g("outputs_compact/COMPACT_validation.json","compact_5yr","hr_per_sd"),
  "1.103-2.745, p=0.0171", "confirmatory", "outputs_metrics_v2/UNIFIED_metrics.json", "36_unified_metrics_v2.py"),
 ("DK time-dependent AUC at 60 months", 0.68831, "0.526-0.832", "confirmatory",
  "outputs_metrics_v2/UNIFIED_metrics.json", "36_unified_metrics_v2.py"),
 ("DK matched random-panel null, empirical p", g("outputs_nulls/NULLS_results.json","dk","p", default=0.0033),
  "32 of 10,000 panels >= observed", "confirmatory", "outputs_nulls/NULLS_results.json", "28_matched_nulls.py"),
 # ---------------- pre-specified external ----------------
 ("CCR2020 ACC-I vs ACC-II subtype AUC", g("outputs_ccr2020/CCR2020_results.json","test1","auc"),
  "0.773-0.968 (bootstrap)", "prespecified", "outputs_ccr2020/CCR2020_results.json", "23_ccr2020_external.py"),
 ("CCR2020 subtype, Holm-adjusted p", 6.137e-06, "-", "prespecified",
  "outputs_ccr2020/CCR2020_results.json", "23_ccr2020_external.py"),
 ("CCR2020 matched random-panel p", 0.013499, "134 of 10,000", "prespecified",
  "outputs_nulls/NULLS_results.json", "28_matched_nulls.py  [supersedes 0.0145 in CCR2020_results.json]"),
 ("CCR2020 time-to-death Spearman rho (decedents only)", -0.32895, "p=0.0575", "negative",
  "outputs_ccr2020/CCR2020_results.json", "23_ccr2020_external.py  [FAILED its criterion; exploratory]"),
 # ---------------- secondary ----------------
 ("DK full-follow-up C-index", g("outputs_compact/COMPACT_validation.json","compact_full","C"),
  "HR/SD 1.527, p=0.0167", "secondary", "outputs_compact/COMPACT_validation.json", "13_compact_panel.py"),
 ("DK adjusted HR (stage + solid histology)", 1.53197, "1.006-2.332, p=0.0466; 6.7 events/parameter",
  "secondary", "outputs_coxdx/COXDX_results.json", "26_cox_diagnostics.py"),
 ("Proportional-hazards test, score", 0.779432, "spline comparison p=0.774", "secondary",
  "outputs_coxdx/COXDX_results.json", "26_cox_diagnostics.py"),
 ("Frozen single-sample score, DK C", 0.6893095768374164,
  "HR/SD 1.936 (1.186-3.162), p=0.0082", "exploratory",
  "outputs_threshold/THRESHOLD_results.json", "27_threshold.py  [POST HOC: designed after DK unblinded]"),
 ("Discovery-locked threshold, DK HR", 3.3829751372814765, "1.433-7.985, log-rank p=0.0033",
  "exploratory", "outputs_threshold/THRESHOLD_results.json", "27_threshold.py  [POST HOC]"),
 # ---------------- composition ----------------
 ("DK three-compartment composition score, C", 0.65590, "HR/SD 1.630, p=0.0391", "prespecified",
  "outputs_deconv/DK_COMPOSITION_results.json", "22_dk_composition.py"),
 ("Composition vs expression score, Spearman", 0.477, "p=2.67e-4; combining changes C by -0.006",
  "secondary", "outputs_deconv/DK_COMPOSITION_results.json", "22_dk_composition.py"),
 ("Deconvolution gate G3 (per-patient composition r)", 0.685, "required 0.70 - FAILED; "
  "rank-based statements only", "negative", "outputs_deconv/GATE_results.json", "21_signature_and_gate.py"),
 # ---------------- internal audit ----------------
 ("Internal leakage-free discovery AUC", g("outputs_metrics/INTERNAL_stability_full.json","summary","leakage_free","mean"),
  "2.5-97.5 pct 0.7475-0.8152 over 14 draws; published 0.810 lies inside", "audit",
  "outputs_metrics/INTERNAL_stability_full.json", "31e_internal_stability_full.py"),
 ("Measured information leak (paired)", g("outputs_metrics/INTERNAL_stability_full.json","summary","leak_paired","mean"),
  "+/- 0.0247", "audit", "outputs_metrics/INTERNAL_stability_full.json", "31e_internal_stability_full.py"),
 ("Originally reported internal AUC (WITHDRAWN)", 0.906, "do not quote", "audit",
  "outputs_nested/NESTED_results.json", "25_nested_cv.py  [WITHDRAWN]"),
 ("Locked-panel recovery, leakage-free 25-fold", "DSCAM 5/25, ODC1 5/25, NCAPG 3/25, CCNB2 0/25",
  "no gene above 5/25", "audit", "outputs_nested/NESTED_results.json", "25b_nested_fixedk.py"),
 # ---------------- new: review response ----------------
 ("ADVI run-to-run: top-4 Jaccard across 20 seeds", g("outputs_vi/VI_DIAGNOSTICS.json","multiseed","topk_jaccard","top4","mean"),
  "beta Pearson 0.948 but |beta| rank Spearman 0.241", "audit",
  "outputs_vi/VI_DIAGNOSTICS.json", "33_vi_diagnostics.py"),
 ("ADVI run-to-run: DK C of each seed's own top-4", g("outputs_vi/VI_DIAGNOSTICS.json","multiseed","dk_c_seed_own_panel","mean"),
  "SD 0.043; locked panel 0.671", "audit", "outputs_vi/VI_DIAGNOSTICS.json", "33_vi_diagnostics.py"),
 ("Hyperparameter sensitivity: DK C over 12 configs x 5 seeds", g("outputs_vi/VI_DIAGNOSTICS.json","sensitivity","dk_c_own_panel","mean"),
  "SD 0.038, range 0.522-0.689; 60/60 distinct panels", "audit",
  "outputs_vi/VI_DIAGNOSTICS.json", "33_vi_diagnostics.py"),
 ("Penalized comparator selected by discovery LOO, DK C", g("outputs_vi/COMPARATOR.json","selected","dk_c_dense"),
  "L1 C=0.3, LOO AUC 0.830; top-4 form 0.584; locked panel 0.671", "audit",
  "outputs_vi/COMPARATOR.json", "33b_comparator.py"),
 ("Axis stability: proliferation-programme enrichment in fold panels", g("outputs_axis/AXIS_stability.json","programme_concentration","proliferation_union","enrichment"),
  "p<1e-4; 16/25 folds vs best single gene 5/25", "audit",
  "outputs_axis/AXIS_stability.json", "34_axis_stability.py"),
 ("Axis stability: fold-panel DK C", g("outputs_axis/AXIS_stability.json","external_discrimination","fold_panel_dk_c","mean"),
  "SD 0.050; random panels 0.509; NOT significantly above random (p=0.094)", "audit",
  "outputs_axis/AXIS_stability.json", "34_axis_stability.py"),
 ("DWLS cross-check: composite score rank agreement", g("outputs_dwls/DWLS_CROSSCHECK.json","composition_survival","score_spearman_local_vs_reference"),
  "published CRAN DWLS 0.1.0 vs local; reference failed on 22/54 samples", "audit",
  "outputs_dwls/DWLS_CROSSCHECK.json", "35_dwls_crosscheck.py"),
 ("DWLS cross-check: synthetic per-compartment recovery", g("outputs_dwls/DWLS_CROSSCHECK.json","synthetic_recovery","local","median_per_compartment_spearman"),
  "reference 0.738 on the same mixtures", "audit",
  "outputs_dwls/DWLS_CROSSCHECK.json", "35_dwls_crosscheck.py"),
 # ---------------- comparison ----------------
 ("Brayer 14-gene classifier, DK C", 0.54009, "our 0.671; NOT a superiority claim - see section 9",
  "secondary", "outputs_metrics_v2/UNIFIED_metrics.json", "18_brayer_benchmark.py"),
 ("8-gene panel leave-one-gene-out gate", "FAILED (p_max 0.333 vs <0.05)", "delta-C vs 4-gene "
  "-0.001 (-0.066, +0.057)", "negative", "outputs_panel8/PANEL8_validation.json", "16_panel8_lock.py"),
 ("Data-derived DK subgroup, best Holm p", 0.23926, "raw 0.0342 - FAILED after correction",
  "negative", "outputs_subtype/SUBTYPE_results.json", "17_subtype_discovery.py"),
]

df = pd.DataFrame(ROWS, columns=["quantity", "value", "uncertainty_or_note", "tier",
                                 "source_file", "generating_script"])
df.to_csv(f"{REL}/results_dictionary.csv", index=False)
print(df.to_string(index=False, max_colwidth=52))
print(f"\n{len(df)} entries -> {REL}/results_dictionary.csv")
print("\ntier counts:"); print(df.tier.value_counts().to_string())
