"""STAGE 7e - how stable is the internal-cohort leakage-free estimate across runs?

`31b_internal_oof.py` re-ran the original design (LOO x3 seeds + 5x5 CV, fixed k=4/sign) and got
disc_auc = 0.767 for the leakage-free arm against the published 0.810 (work/outputs_nested/
log_25c.txt), failing its +/-0.04 regression tolerance.

The decomposition is informative and is why this script exists:

                     published (25c)     re-run (31b)
  LOO  component     0.843 (sd 0.005)    0.757 (sd 0.037)
  CV5  component     0.776 (sd 0.036)    0.778 (sd 0.051)
  disc_auc           0.810               0.767

The CV5 component reproduces to 0.002. The whole discrepancy sits in the LOO component, which is
an AUC over 20 single-patient predictions, each from a separately fitted ADVI model - so a small
perturbation in any one fit can flip ranked pairs. The published LOO sd of 0.005 across three
seeds looks implausibly tight against a re-run sd of 0.037.

This script settles whether 0.810 vs 0.767 is ordinary run-to-run noise by repeating the
leakage-free arm N_REPEATS times end to end, reporting the spread of each component separately.
No number is adjusted on the basis of the result; the outcome is reported either way.

Reads only; writes work/outputs_metrics/INTERNAL_stability.json.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json

import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold

R = f"{ACC_DATA_ROOT}"
SRC = f"{HERE}/25_nested_cv.py"
OUT_M = f"{ACC_DATA_ROOT}/model_outputs_metrics"

_src = open(SRC).read().split('res = {"n_patients"')[0]
exec(compile(_src, SRC, "exec"))

FIXED_K, FIXED_SCHEME = 4, "sign"
N_REPEATS = 3
ARM = "leakage_free"


def one_run(offset):
    """One full pass of the original design. `offset` shifts every ADVI seed."""
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
    buf = []
    for n, (tr_idx, te_idx, seed, kind) in enumerate(folds):
        Ztr, Zte_fn, M, sid, nsrc, s_ev, d_ev = fold_artifacts(tr_idx, ARM)
        b = fit_beta(Ztr, y[tr_idx], M, sid, nsrc, s_ev, d_ev, seed=seed + offset)
        idx_g = np.argsort(-np.abs(b))[:FIXED_K]
        sc = panel_score(Zte_fn(te_idx), idx_g, b[idx_g], FIXED_SCHEME)
        if kind == "loo":
            pred_loo[seed, te_idx[0]] = sc[0]
        else:
            buf.append((te_idx, sc))
        if (n + 1) % 25 == 0:
            print(f"    {n+1}/{len(folds)} fits")
    for rep in range(5):
        for te_idx, sc in buf[rep * 5:(rep + 1) * 5]:
            pred_cv5[rep, te_idx] = sc
    al = [roc_auc_score(y, pred_loo[s]) for s in range(3)]
    ac = [roc_auc_score(y, pred_cv5[r]) for r in range(5)]
    return (float(np.mean(al)), float(np.std(al)), float(np.mean(ac)), float(np.std(ac)),
            float(.5 * np.mean(al) + .5 * np.mean(ac)), pred_loo, pred_cv5)


runs, keep = [], {}
for k in range(N_REPEATS):
    off = 1000 * (k + 1)
    print(f"\n=== repeat {k+1}/{N_REPEATS}  (seed offset {off}) ===")
    torch.manual_seed(off)
    np.random.seed(off)
    lm, ls, cm, cs, disc, ploo, pcv = one_run(off)
    print(f"  LOO {lm:.3f} (sd {ls:.3f})   CV5 {cm:.3f} (sd {cs:.3f})   disc_auc {disc:.3f}")
    runs.append({"offset": off, "loo_auc": lm, "loo_sd": ls, "cv5_auc": cm, "cv5_sd": cs,
                 "disc_auc": disc})
    keep[f"run{k}_loo"] = ploo
    keep[f"run{k}_cv5"] = pcv

d = np.array([r["disc_auc"] for r in runs])
l = np.array([r["loo_auc"] for r in runs])
c = np.array([r["cv5_auc"] for r in runs])

print("\n" + "=" * 74)
print(f"leakage-free, {N_REPEATS} independent repeats of the original design")
print("=" * 74)
print(f"  LOO component  {l.mean():.3f} +/- {l.std():.3f}   (published 0.843, 31b 0.757)")
print(f"  CV5 component  {c.mean():.3f} +/- {c.std():.3f}   (published 0.776, 31b 0.778)")
print(f"  disc_auc       {d.mean():.3f} +/- {d.std():.3f}   (published 0.810, 31b 0.767)")
lo, hi = float(d.min()), float(d.max())
covers = lo - 0.005 <= 0.810 <= hi + 0.005
print(f"  observed range [{lo:.3f}, {hi:.3f}] "
      f"{'covers' if covers else 'DOES NOT cover'} the published 0.810")

np.savez(f"{OUT_M}/internal_oof_repeats.npz", y=y.astype(int), **keep)
json.dump({
    "arm": ARM, "design": "LOO x3 seeds + 5x5 CV, fixed k=4/sign", "n_repeats": N_REPEATS,
    "published_25c": {"loo": 0.843, "cv5": 0.776, "disc_auc": 0.810},
    "run_31b": {"loo": 0.757, "cv5": 0.778, "disc_auc": 0.767},
    "repeats": runs,
    "summary": {"loo_mean": float(l.mean()), "loo_sd": float(l.std()),
                "cv5_mean": float(c.mean()), "cv5_sd": float(c.std()),
                "disc_mean": float(d.mean()), "disc_sd": float(d.std()),
                "disc_min": lo, "disc_max": hi,
                "range_covers_published_0810": bool(covers)},
    "interpretation": (
        "The CV5 component is reproducible; the LOO component is not. LOO AUC here is computed "
        "over 20 single-patient out-of-fold scores, each from a separately fitted ADVI model, so "
        "it is far more sensitive to fit-level noise than the 5x5 CV. Any internal-cohort figure "
        "should therefore be driven by the CV5 component, and the published 0.810 should carry "
        "an explicit run-to-run uncertainty rather than being quoted as a point estimate."),
}, open(f"{OUT_M}/INTERNAL_stability.json", "w"), indent=2)
print(f"\nsaved -> {OUT_M}/INTERNAL_stability.json")
