"""Out-of-fold predictions for the internal cohort.

Runs the resampling design of `25c_fidelity.py` at fixed k=4/sign - leave-one-out x 3 seeds plus
RepeatedStratifiedKFold(5, 5), scored as 0.5*mean(LOO) + 0.5*mean(CV5) - for both arms, and saves
the per-patient prediction vectors used for the internal-cohort ROC curves.

Writes model_outputs_metrics/internal_oof_predictions.npz.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold

R = f"{ACC_DATA_ROOT}"
SRC = f"{HERE}/25_nested_cv.py"
METRICS_OUT = f"{ACC_DATA_ROOT}/model_outputs_metrics"

# reuse 25_nested_cv.py above its driver section: the data (`y`, `genes`),
# the in-fold rebuilders (`fold_artifacts`) and the model wrappers (`fit_beta`, `panel_score`)
_src = open(SRC).read().split('res = {"n_patients"')[0]
exec(compile(_src, SRC, "exec"))

FIXED_K, FIXED_SCHEME = 4, "sign"


def run_original_design(arm):
    """13_compact_panel.py's folds: LOO x 3 seeds, then RepeatedStratifiedKFold(5,5)."""
    print("\n" + "=" * 74)
    print(f"ARM: {arm}  (LOO x3 seeds + 5x5 CV, k=4/sign)")
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
    pred_cv5 = np.full((5, len(y)), np.nan)
    cv_buf = []
    for n, (tr_idx, te_idx, seed, kind) in enumerate(folds):
        Ztr, Zte_fn, M, sid, nsrc, s_ev, d_ev = fold_artifacts(tr_idx, arm)
        b = fit_beta(Ztr, y[tr_idx], M, sid, nsrc, s_ev, d_ev, seed=seed)
        idx_g = np.argsort(-np.abs(b))[:FIXED_K]
        sc = panel_score(Zte_fn(te_idx), idx_g, b[idx_g], FIXED_SCHEME)
        if kind == "loo":
            pred_loo[seed, te_idx[0]] = sc[0]
        else:
            cv_buf.append((te_idx, sc))
        if (n + 1) % 25 == 0:
            print(f"  {n+1}/{len(folds)} fits")

    for rep in range(5):
        for te_idx, sc in cv_buf[rep * 5:(rep + 1) * 5]:
            pred_cv5[rep, te_idx] = sc

    auc_loo = [roc_auc_score(y, pred_loo[s]) for s in range(3)]
    auc_cv = [roc_auc_score(y, pred_cv5[r]) for r in range(5)]
    disc = 0.5 * np.mean(auc_loo) + 0.5 * np.mean(auc_cv)
    print(f"  LOO  AUC = {np.mean(auc_loo):.3f} (sd {np.std(auc_loo):.3f})")
    print(f"  CV5  AUC = {np.mean(auc_cv):.3f} (sd {np.std(auc_cv):.3f})")
    print(f"  disc_auc = {disc:.3f}")
    return pred_loo, pred_cv5, float(disc), float(np.mean(auc_loo)), float(np.mean(auc_cv))


out = {"y": y.astype(int)}
summary = {}
for arm in ["as_built", "leakage_free"]:
    ploo, pcv, disc, a_loo, a_cv = run_original_design(arm)
    out[f"oof_{arm}_loo"] = ploo
    out[f"oof_{arm}_cv5"] = pcv
    summary[arm] = {"disc_auc": disc, "loo_auc": a_loo, "cv5_auc": a_cv}

np.savez(f"{METRICS_OUT}/internal_oof_predictions.npz", **out)
for arm in ["as_built", "leakage_free"]:
    print(f"  {arm:<14s} disc_auc {summary[arm]['disc_auc']:.3f}")
print(f"\nsaved -> {METRICS_OUT}/internal_oof_predictions.npz")
