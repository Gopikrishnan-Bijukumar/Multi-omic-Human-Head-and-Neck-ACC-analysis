"""STAGE 7b - recover the internal cohort's out-of-fold predictions.

`25_nested_cv.py`, `25b_nested_fixedk.py` and `25c_fidelity.py` all build per-patient out-of-fold
prediction vectors and then discard them, keeping only summary AUCs. That is enough for a table
but not for a ROC curve, so the internal cohort could not appear alongside DK and CCR2020 in the
unified-metrics figures.

This script re-runs the ORIGINAL resampling design at fixed k=4/sign - exactly `25c_fidelity.py`'s
`run_original_design`, i.e. leave-one-out x 3 seeds plus RepeatedStratifiedKFold(5, 5), scored as
0.5*mean(LOO) + 0.5*mean(CV5) - for both arms, and saves the prediction vectors.

It does NOT edit 25b or 25c. The harness is borrowed with the same exec-the-top idiom those two
scripts already use, so no script that produced a published number is touched.

Regression check: the leakage-free arm must reproduce disc_auc = 0.810 +/- 0.04 and the as-built
arm 0.881 +/- 0.04 (work/outputs_nested/log_25c.txt). ADVI is stochastic, hence the tolerance.

Reads only; writes work/outputs_metrics/internal_oof_predictions.npz.
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

# reuse the harness verbatim: everything above the driver line defines the data (`y`, `genes`),
# the in-fold rebuilders (`fold_artifacts`) and the model wrappers (`fit_beta`, `panel_score`)
_src = open(SRC).read().split('res = {"n_patients"')[0]
exec(compile(_src, SRC, "exec"))

FIXED_K, FIXED_SCHEME = 4, "sign"


def run_original_design(arm):
    """13_compact_panel.py's folds: LOO x 3 seeds, then RepeatedStratifiedKFold(5,5)."""
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

print("\n" + "=" * 74)
print("REGRESSION CHECK vs work/outputs_nested/log_25c.txt")
print("=" * 74)
fail = []
for arm, want in [("as_built", 0.881), ("leakage_free", 0.810)]:
    got = summary[arm]["disc_auc"]
    ok = abs(got - want) <= 0.04
    if not ok:
        fail.append((arm, got, want))
    print(f"  [{'OK ' if ok else 'FAIL'}] {arm:<14s} disc_auc got {got:.3f}  expect {want:.3f} "
          f"(+/-0.04, ADVI is stochastic)")
print(f"  measured leak = {summary['as_built']['disc_auc'] - summary['leakage_free']['disc_auc']:+.3f}"
      f"   (published +0.071)")
if fail:
    print("\n*** REGRESSION FAILURE - do not use these predictions ***")
else:
    print(f"\nsaved -> {METRICS_OUT}/internal_oof_predictions.npz")
