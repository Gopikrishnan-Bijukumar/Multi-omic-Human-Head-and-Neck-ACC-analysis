"""Nested cross-validation at a fixed panel (k=4, scheme=sign).

Adds two fixed-panel arms to `25_nested_cv.py`:

  as_built_fixed      k=4, scheme=sign, all priors/modules/RINT from the full 20 patients
  leakage_free_fixed  k=4, scheme=sign, everything rebuilt in-fold

The machinery is reused from 25_nested_cv.py: the file is executed up to its driver section,
which defines the data, the in-fold rebuilders and the model wrappers.

Merges into model_outputs_nested/NESTED_results.json.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold

R = f"{ACC_DATA_ROOT}"
SRC = f"{HERE}/25_nested_cv.py"

# reuse everything in 25_nested_cv.py above its driver section
_src = open(SRC).read().split('res = {"n_patients"')[0]
exec(compile(_src, SRC, "exec"))

FIXED_K = 4
FIXED_SCHEME = "sign"


def run_fixed(arm, n_repeats=5):
    print("\n" + "=" * 74)
    print(f"ARM: {arm}_fixed   (k={FIXED_K}, scheme={FIXED_SCHEME}, no inner selection)")
    print("=" * 74)
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=7)
    oof = {r: np.full(len(y), np.nan) for r in range(n_repeats)}
    picked = []
    for f, (tr_idx, te_idx) in enumerate(rskf.split(np.zeros(len(y)), y)):
        rep = f // 5
        Ztr, Zte_fn, M, sid, nsrc, s_ev, d_ev = fold_artifacts(tr_idx, arm)
        b = fit_beta(Ztr, y[tr_idx], M, sid, nsrc, s_ev, d_ev, seed=100 + f)
        idx_g = np.argsort(-np.abs(b))[:FIXED_K]
        oof[rep][te_idx] = panel_score(Zte_fn(te_idx), idx_g, b[idx_g], FIXED_SCHEME)
        picked.append([genes[i] for i in idx_g])
        if (f + 1) % 5 == 0:
            print(f"  {f+1}/{n_repeats*5} folds")
    aucs = [roc_auc_score(y, oof[r]) for r in range(n_repeats)]
    flat = pd.Series([g for p in picked for g in p]).value_counts()
    return {"arm": f"{arm}_fixed", "k": FIXED_K, "scheme": FIXED_SCHEME,
            "auc_mean": float(np.mean(aucs)), "auc_sd": float(np.std(aucs)),
            "auc_per_repeat": [float(a) for a in aucs],
            "gene_selection_frequency": {k: int(v) for k, v in flat.head(15).items()},
            "locked_panel_recovery": {g: int(flat.get(g, 0)) for g in
                                      ["DSCAM", "ODC1", "NCAPG", "CCNB2"]}}


a_fix = run_fixed("as_built")
b_fix = run_fixed("leakage_free")

res = json.load(open(f"{OUT}/NESTED_results.json"))
res["arm_as_built_fixed_k4"] = a_fix
res["arm_leakage_free_fixed_k4"] = b_fix
res["optimism_decomposition"] = {
    "as_built_fixed_k4": a_fix["auc_mean"],
    "as_built_nested_k_selection": res["arm_as_built"]["auc_mean"],
    "leakage_free_nested_k_selection": res["arm_leakage_free"]["auc_mean"],
    "leakage_free_fixed_k4": b_fix["auc_mean"],
    "component_k_and_scheme_selection": float(a_fix["auc_mean"] - res["arm_as_built"]["auc_mean"]),
    "component_evidence_prior_module_rint_leak": float(
        res["arm_as_built"]["auc_mean"] - res["arm_leakage_free"]["auc_mean"]),
}

print("\n" + "=" * 74)
print("OPTIMISM DECOMPOSITION")
print("=" * 74)
d = res["optimism_decomposition"]
print(f"  as-built, fixed k=4                     {d['as_built_fixed_k4']:.3f}")
print(f"  as-built, k chosen in inner loop        {d['as_built_nested_k_selection']:.3f}")
print(f"  leakage-free, k chosen in inner loop    {d['leakage_free_nested_k_selection']:.3f}")
print(f"  leakage-free, fixed k=4 (locked panel)  {d['leakage_free_fixed_k4']:.3f}")
print()
print(f"  attributable to k/scheme selection      {d['component_k_and_scheme_selection']:+.3f}")
print(f"  attributable to evidence/prior leakage  {d['component_evidence_prior_module_rint_leak']:+.3f}")
print("\nlocked-panel recovery under leakage-free refitting (out of 25 folds):")
for g, n in b_fix["locked_panel_recovery"].items():
    print(f"  {g:<8s} selected in {n:2d}/25 folds")

json.dump(res, open(f"{OUT}/NESTED_results.json", "w"), indent=1)
print(f"\nupdated {OUT}/NESTED_results.json")
