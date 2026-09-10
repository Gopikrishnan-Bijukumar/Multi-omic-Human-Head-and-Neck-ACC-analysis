"""Build the frozen manuscript release (review item #5).

Creates, from artefacts already on disk, a single self-consistent bundle:
  MANIFEST.json           sha256 of every raw input, including the five external-cohort
                          files that carried no fingerprint anywhere in the project
  model_spec.json         the one model specification, with the frozen discovery reference
  cohorts.csv             cohort manifest with the ROLE of each cohort
  results_dictionary.csv  every manuscript number -> value, source file, script, evidence tier
  environment.txt         the interpreter and package versions actually used

Verifies as it goes: score_acc4.py must reproduce the stored DK score and the reported
frozen-reference C-index, or the build fails.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)
# The numbered build scripts and common.py live one level up, in 07_prognostic_model/.
PARENT = str(_Path(__file__).resolve().parent.parent)

import json, os, subprocess, sys
import numpy as np, pandas as pd
from scipy import stats
from sksurv.metrics import concordance_index_censored

R = f"{ACC_DATA_ROOT}"
REL = HERE
sys.path.insert(0, HERE); sys.path.insert(0, REL); sys.path.insert(0, PARENT)
from common import input_fingerprint
import importlib.util
spec = importlib.util.spec_from_file_location("d3", f"{PARENT}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

PANEL = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]

# ------------------------------------------------------------------ 1. model spec
compact = json.load(open(f"{ACC_DATA_ROOT}/model_outputs_compact/COMPACT_model.json"))
thr = json.load(open(f"{ACC_DATA_ROOT}/model_outputs_threshold/THRESHOLD_results.json"))
D = d3.load("all")
tr_raw, dk_raw = D["tr_raw"], D["dk_raw"]
lock = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta = pd.Series(np.array(lock["beta"]), index=lock["genes"])

spec_out = {
    "name": "ACC four-gene prognostic signature",
    "genes": PANEL,
    "signs": {g: float(np.sign(beta[g])) for g in PANEL},
    "discovery_coefficients": {g: float(compact["weights"][g]) for g in PANEL},
    "scheme": "sign",
    "scoring": ("mean of sign-weighted rank-inverse-normal expression; all four discovery "
                "coefficients are positive, so the deployed score is an unweighted mean. The "
                "hierarchical Bayesian procedure selected genes and directions and does not "
                "weight the score."),
    "input": "log2 depth-normalised expression, as in outputs/03_harmonized/*_norm_log2.csv",
    "frozen_reference": {
        "description": "discovery-cohort raw log2 values for the four genes; the frozen "
                       "single-sample transform maps a new patient into these ranks",
        "n_samples": int(tr_raw.shape[0]),
        "samples": list(tr_raw.index),
        "values": {g: [float(v) for v in tr_raw[g].values] for g in PANEL}},
    "frozen_threshold": {
        "value": float(thr["locked_threshold"]["threshold"]),
        "derivation": "Youden point of the discovery-cohort frozen score; uses no validation data",
        "status": "post hoc - designed after the Danish cohort was unblinded; externally "
                  "evaluated, not prospectively validated",
        "dk_performance": {"hr": thr["locked_threshold"]["survival"]["hr"],
                           "ci": thr["locked_threshold"]["survival"]["ci"],
                           "logrank_p": thr["locked_threshold"]["survival"]["logrank_p"]}},
    "performance": {
        "dk_batch_C_5yr": thr["procedures"]["a_batch_cohort_rint"]["C_5yr"],
        "dk_frozen_C_5yr": thr["procedures"]["c_frozen_jse_reference"]["C_5yr"]},
    "provenance": {"selected_by": "work/13_compact_panel.py from work/11_lock_and_validate.py",
                   "model_code": "work/mlbayes.py (MultiLevelModel, mean-field ADVI in PyTorch)",
                   "locked_artefact": "work/outputs_compact/COMPACT_model.json"}}
json.dump(spec_out, open(f"{REL}/model_spec.json", "w"), indent=2)
print("model_spec.json written")

# ------------------------------------------------------------------ 2. verify the scorer
import score_acc4
T, E = D["T"], D["E"]
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
s_batch = score_acc4.score_batch(dk_raw[PANEL])
stored = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel4_score.npy")
d_batch = float(np.max(np.abs(s_batch - stored)))
c_batch = float(concordance_index_censored(E5.astype(bool), T5, s_batch)[0])
s_frozen = score_acc4.score_frozen(dk_raw[PANEL], tr_raw[PANEL])
c_frozen = float(concordance_index_censored(E5.astype(bool), T5, s_frozen)[0])
checks = {
    "batch_score_matches_dk_panel4_score_npy": {"max_abs_diff": d_batch, "pass": d_batch < 1e-10},
    "batch_C_5yr": {"got": c_batch, "expect": 0.6714922048997772, "pass": abs(c_batch - 0.6714922048997772) < 1e-9},
    "frozen_C_5yr": {"got": c_frozen, "expect": thr["procedures"]["c_frozen_jse_reference"]["C_5yr"],
                     "pass": abs(c_frozen - thr["procedures"]["c_frozen_jse_reference"]["C_5yr"]) < 1e-9}}
for k, v in checks.items():
    print(f"  [{'OK ' if v['pass'] else 'FAIL'}] {k}: {v.get('got', v.get('max_abs_diff'))}")
assert all(v["pass"] for v in checks.values()), "release scorer does not reproduce locked results"

# ------------------------------------------------------------------ 3. manifest
INPUTS = {
  "internal_bulk_rnaseq":      f"{R}/../bulk_rna_data/h5ad_files/bulk_counts.h5ad",
  "single_cell":               f"{R}/../outputs/h5ad_files/inprogress_3.h5ad",
  "wgcna_module_table":        f"{R}/../bulk_rna_data/outputs/images/clinical_outcome/wgcna/06_gene_info_complete.csv",
  "pydeseq2_significant":      f"{R}/../bulk_rna_data/outputs/pydeseq2_out/clinical_outcome/DEG_significant.csv",
  "pydeseq2_all":              f"{R}/../bulk_rna_data/outputs/pydeseq2_out/clinical_outcome/DEG_results.csv",
  "clinical_metadata":         f"{R}/../clinical_info.csv",
  "dk_expression":             f"{R}/external_dk/DK_All_Genes_Data.csv",
  "dk_survival":               f"{R}/external_dk/DK_sample_data_table2.csv",
  "ccr2020_expression":        f"{ACC_DATA_ROOT}/external_ccr2020/ACC_RNAseq.xlsx",
  "ccr2020_clinical":          f"{ACC_DATA_ROOT}/external_ccr2020/CCR2020_Clinical.xlsx",
  "brayer_tableS2":            f"{ACC_DATA_ROOT}/external_brayer_supp/TableS2_Poor_Survival_DE.csv",
  "brayer_tableS1":            f"{ACC_DATA_ROOT}/external_brayer_supp/TableS1_NoMYB_DE.csv",
  "brayer_table5_screenshot":  f"{R}/table_5_ss.png",
  "harmonized_internal":       f"{R}/outputs/03_harmonized/train_norm_log2.csv",
  "harmonized_dk":             f"{R}/outputs/03_harmonized/dk_norm_log2.csv",
  "evidence_table":            f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv",
  "sample_labels":             f"{R}/outputs/01_candidate_table/sample_labels.csv"}
man = {"generated_utc": pd.Timestamp.utcnow().isoformat(),
       "note": ("Fingerprints for every raw and harmonised input. Before this release the five "
                "external-cohort inputs (DK x2, CCR2020 x2, Brayer) carried no hash anywhere in "
                "the project, although common.py already implemented the machinery."),
       "inputs": {}}
for k, p in INPUTS.items():
    p = os.path.normpath(p)
    man["inputs"][k] = input_fingerprint(p) if os.path.exists(p) else {"path": p, "error": "missing"}
    print(f"  hashed {k}")
json.dump(man, open(f"{REL}/MANIFEST.json", "w"), indent=2)
print("MANIFEST.json written")

# ------------------------------------------------------------------ 4. cohorts
pd.DataFrame([
 {"cohort": "internal (discovery)", "n": 20, "endpoint": "Poor vs Good outcome (binary)",
  "events": "10 Poor / 10 Good", "role": "discovery - model fitting and gene selection",
  "survival_times_available": "no - `years` holds only the 2/5 group label",
  "used_for_fitting": "yes"},
 {"cohort": "Danish (DK)", "n": 54, "endpoint": "overall survival, censored at 60 months",
  "events": "21 deaths by 60 months, 0 censored before 60 months",
  "role": "confirmatory external validation - the single primary claim",
  "survival_times_available": "yes", "used_for_fitting": "no - never entered any fit"},
 {"cohort": "CCR2020 (MD Anderson)", "n": 54, "endpoint": "ACC-I vs ACC-II molecular subtype",
  "events": "20 ACC-I / 34 ACC-II", "role": "pre-specified external biological validation - "
  "SUBTYPE, not survival",
  "survival_times_available": "no - censoring times missing for all 20 living patients",
  "used_for_fitting": "no"}]).to_csv(f"{REL}/cohorts.csv", index=False)
print("cohorts.csv written")
