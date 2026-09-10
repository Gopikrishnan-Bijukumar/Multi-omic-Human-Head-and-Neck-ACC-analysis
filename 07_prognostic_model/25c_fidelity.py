"""
Stage 7c - exact-design fidelity check.

`25b` reported the as-built fixed-k=4 arm at 0.856 against a published 0.906. That gap is not
necessarily a harness bug: `13_compact_panel.py` averages a leave-one-out component (20 folds x 3
seeds) with a 5x5 cross-validation component, whereas 25/25b use 5x5 only.

This script settles it by replicating the ORIGINAL resampling design exactly - LOO x 3 seeds plus
RepeatedStratifiedKFold(5, 5), scored as 0.5*mean(LOO) + 0.5*mean(CV5), at fixed k=4/sign - and
runs the leakage-free arm through the same design so the two are directly comparable.

If the as-built arm reproduces ~0.906, the harness is faithful and the leakage estimate stands.

Reads only; merges into work/outputs_nested/NESTED_results.json.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold

R = f"{ACC_DATA_ROOT}"
SRC = f"{HERE}/25_nested_cv.py"
_src = open(SRC).read().split('res = {"n_patients"')[0]
exec(compile(_src, SRC, "exec"))

FIXED_K, FIXED_SCHEME = 4, "sign"


def run_original_design(arm):
    """Exactly 13_compact_panel.py's folds: LOO x 3 seeds, then RepeatedStratifiedKFold(5,5)."""
    print("\n" + "=" * 74)
    print(f"ARM: {arm}  (original design: LOO x3 seeds + 5x5 CV, k=4/sign)")
    print("=" * 74)
    folds = []
    for seed in range(3):
        for i in range(len(y)):
            folds.append((np.array([j for j in range(len(y)) if j != i]), np.array([i]),
                          seed, "loo"))
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=7)
    for f, (tr, te) in enumerate(rskf.split(np.zeros(len(y)), y)):
        folds.append((tr, te, 100 + f, "cv5"))

    pred_loo = np.full((3, len(y)), np.nan)
    pred_cv = []
    for n, (tr_idx, te_idx, seed, kind) in enumerate(folds):
        Ztr, Zte_fn, M, sid, nsrc, s_ev, d_ev = fold_artifacts(tr_idx, arm)
        b = fit_beta(Ztr, y[tr_idx], M, sid, nsrc, s_ev, d_ev, seed=seed)
        idx_g = np.argsort(-np.abs(b))[:FIXED_K]
        sc = panel_score(Zte_fn(te_idx), idx_g, b[idx_g], FIXED_SCHEME)
        if kind == "loo":
            pred_loo[seed, te_idx[0]] = sc[0]
        else:
            pred_cv.append((te_idx, sc))
        if (n + 1) % 25 == 0:
            print(f"  {n+1}/{len(folds)} fits")

    auc_loo = [roc_auc_score(y, pred_loo[s]) for s in range(3)]
    auc_cv = []
    for rep in range(5):
        p = np.full(len(y), np.nan)
        for te_idx, sc in pred_cv[rep * 5:(rep + 1) * 5]:
            p[te_idx] = sc
        auc_cv.append(roc_auc_score(y, p))
    disc = 0.5 * np.mean(auc_loo) + 0.5 * np.mean(auc_cv)
    print(f"  LOO  AUC = {np.mean(auc_loo):.3f} (sd {np.std(auc_loo):.3f})")
    print(f"  CV5  AUC = {np.mean(auc_cv):.3f} (sd {np.std(auc_cv):.3f})")
    print(f"  disc_auc = {disc:.3f}   <- comparable to the published table")
    return {"loo_auc": float(np.mean(auc_loo)), "loo_sd": float(np.std(auc_loo)),
            "cv5_auc": float(np.mean(auc_cv)), "cv5_sd": float(np.std(auc_cv)),
            "disc_auc": float(disc)}


a = run_original_design("as_built")
b = run_original_design("leakage_free")

res = json.load(open(f"{OUT}/NESTED_results.json"))
res["original_design_fidelity"] = {
    "published_disc_auc_k4_sign": 0.906,
    "published_loo": 0.907, "published_cv5": 0.906,
    "as_built": a, "leakage_free": b,
    "fidelity_gap_vs_published": float(a["disc_auc"] - 0.906),
    "leak_under_original_design": float(a["disc_auc"] - b["disc_auc"]),
}

print("\n" + "=" * 74)
print("FIDELITY")
print("=" * 74)
print(f"  published            0.906  (LOO 0.907 / CV5 0.906)")
print(f"  as-built reproduced  {a['disc_auc']:.3f}  (LOO {a['loo_auc']:.3f} / CV5 {a['cv5_auc']:.3f})")
print(f"  gap                  {a['disc_auc'] - 0.906:+.3f}")
print(f"  leakage-free         {b['disc_auc']:.3f}")
print(f"  leak                 {a['disc_auc'] - b['disc_auc']:+.3f}")

json.dump(res, open(f"{OUT}/NESTED_results.json", "w"), indent=1)
print(f"\nupdated {OUT}/NESTED_results.json")
