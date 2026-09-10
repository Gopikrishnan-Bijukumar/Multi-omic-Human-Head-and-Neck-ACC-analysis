"""STAGE 7f - a reproducible internal-cohort figure, and a reproducible leak estimate.

`31d_internal_stability.py` showed that the published leakage-free 0.810 (work/outputs_nested/
log_25c.txt) is a single stochastic draw that does not reproduce: 3 repeats gave 0.780 +/- 0.014,
range [0.761, 0.795]. That settled *whether* there is a problem but left two things unfinished.

1. The replacement estimate rests on 4 draws (31b's one + 31d's three). Four draws pin a mean to
   about +/-0.007 and cannot support a percentile interval at all. A number that goes into a
   manuscript in place of a withdrawn one has to be better resolved than the one it replaces.

2. The published **leak** is +0.071 = 0.881 (as_built) - 0.810 (leakage_free). If the leakage-free
   term moves, so does the leak - and 0.881 is itself a single draw. Quoting a corrected 0.78
   against an uncorrected 0.881 would silently restate the leak as +0.10 with no evidence behind
   either end. The leak has to be re-estimated the same way, and it has to be estimated **paired**:
   both arms run under the same seed offset, leak taken within-repeat, so the shared fit-level noise
   cancels instead of being counted twice.

So this script runs N_REPEATS paired repeats of the original design (LOO x3 seeds + 5x5 CV, fixed
k=4/sign) over BOTH arms, and reports for each arm and for the paired leak: mean, sd, standard
error, and a 2.5-97.5 percentile interval over draws.

The three leakage-free draws already in INTERNAL_stability.json (offsets 1000/2000/3000) and 31b's
unoffset draw are pooled into the arm-level leakage-free summary - they are independent draws from
the same design, so discarding them would throw away evidence. Offsets here start at 4000 to
guarantee no seed collision. The paired leak uses only this script's repeats, because pairing
requires both arms from one offset.

No number is adjusted on the basis of the result. Writes work/outputs_metrics/
INTERNAL_stability_full.json (rewritten after every repeat, so a partial run is still usable) and
internal_stability_full_preds.npz. Reads only otherwise.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import os

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
N_REPEATS = 10
OFFSET0, OFFSET_STEP = 4000, 1000
ARMS = ["as_built", "leakage_free"]

# independent draws already on disk, folded into the arm-level summaries (not into the paired leak)
PRIOR_DRAWS = {"as_built": [0.881], "leakage_free": [0.767, 0.795, 0.785, 0.761]}
PUBLISHED = {"as_built": 0.881, "leakage_free": 0.810, "leak": 0.071}


def one_run(arm, offset):
    """One full pass of the original design for one arm. `offset` shifts every ADVI seed."""
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
        Ztr, Zte_fn, M, sid, nsrc, s_ev, d_ev = fold_artifacts(tr_idx, arm)
        b = fit_beta(Ztr, y[tr_idx], M, sid, nsrc, s_ev, d_ev, seed=seed + offset)
        idx_g = np.argsort(-np.abs(b))[:FIXED_K]
        sc = panel_score(Zte_fn(te_idx), idx_g, b[idx_g], FIXED_SCHEME)
        if kind == "loo":
            pred_loo[seed, te_idx[0]] = sc[0]
        else:
            buf.append((te_idx, sc))
    for rep in range(5):
        for te_idx, sc in buf[rep * 5:(rep + 1) * 5]:
            pred_cv5[rep, te_idx] = sc

    al = [roc_auc_score(y, pred_loo[s]) for s in range(3)]
    ac = [roc_auc_score(y, pred_cv5[r]) for r in range(5)]
    return {"loo_auc": float(np.mean(al)), "loo_sd": float(np.std(al)),
            "cv5_auc": float(np.mean(ac)), "cv5_sd": float(np.std(ac)),
            "disc_auc": float(.5 * np.mean(al) + .5 * np.mean(ac))}, pred_loo, pred_cv5


def describe(v, label):
    """mean / sd / se / percentile interval over independent draws."""
    v = np.asarray(v, float)
    lo, hi = (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) >= 5 \
        else (float(v.min()), float(v.max()))
    d = {"label": label, "n_draws": int(len(v)), "mean": float(v.mean()),
         "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
         "se": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0,
         "min": float(v.min()), "max": float(v.max()),
         "pct_lo": lo, "pct_hi": hi,
         "interval_kind": "2.5-97.5 percentile" if len(v) >= 5 else "observed range"}
    return d


repeats, keep = [], {}
for k in range(N_REPEATS):
    off = OFFSET0 + OFFSET_STEP * k
    torch.manual_seed(off)
    np.random.seed(off)
    row = {"offset": off}
    for arm in ARMS:
        m, ploo, pcv = one_run(arm, off)
        row[arm] = m
        keep[f"r{k}_{arm}_loo"] = ploo
        keep[f"r{k}_{arm}_cv5"] = pcv
    row["leak"] = row["as_built"]["disc_auc"] - row["leakage_free"]["disc_auc"]
    repeats.append(row)
    print(f"[{k+1:2d}/{N_REPEATS}] off={off}  as_built {row['as_built']['disc_auc']:.3f}   "
          f"leakage_free {row['leakage_free']['disc_auc']:.3f}   leak {row['leak']:+.3f}",
          flush=True)

    # rewrite after every repeat so an interrupted run still leaves a usable artifact
    arm_draws = {a: [r[a]["disc_auc"] for r in repeats] + PRIOR_DRAWS[a] for a in ARMS}
    leak_paired = [r["leak"] for r in repeats]
    summ = {a: describe(arm_draws[a], f"{a} disc_auc, pooled independent draws") for a in ARMS}
    summ["leak_paired"] = describe(leak_paired, "as_built - leakage_free, paired within repeat")
    for a in ARMS:
        summ[a]["covers_published"] = bool(summ[a]["min"] <= PUBLISHED[a] <= summ[a]["max"])
        summ[a]["loo_mean"] = float(np.mean([r[a]["loo_auc"] for r in repeats]))
        summ[a]["cv5_mean"] = float(np.mean([r[a]["cv5_auc"] for r in repeats]))
    summ["leak_paired"]["covers_published"] = bool(
        summ["leak_paired"]["min"] <= PUBLISHED["leak"] <= summ["leak_paired"]["max"])

    json.dump({
        "design": "LOO x3 seeds + 5x5 CV, fixed k=4/sign; both arms paired per seed offset",
        "n_repeats_paired": len(repeats), "offsets": [r["offset"] for r in repeats],
        "prior_independent_draws_pooled": PRIOR_DRAWS,
        "published_25c": PUBLISHED, "repeats": repeats, "summary": summ,
    }, open(f"{OUT_M}/INTERNAL_stability_full.json", "w"), indent=2)

np.savez(f"{OUT_M}/internal_stability_full_preds.npz", y=y.astype(int), **keep)

print("\n" + "=" * 78)
print(f"internal cohort, original design, {N_REPEATS} paired repeats")
print("=" * 78)
for key in ["as_built", "leakage_free", "leak_paired"]:
    s = summ[key]
    print(f"  {key:<14s} {s['mean']:+.3f} +/- {s['sd']:.3f} (se {s['se']:.3f})  "
          f"[{s['pct_lo']:.3f}, {s['pct_hi']:.3f}] {s['interval_kind']}, n={s['n_draws']}  "
          f"published {PUBLISHED[key.replace('_paired','')]:.3f} "
          f"{'covered' if s['covers_published'] else 'NOT covered'}")
print(f"\nsaved -> {OUT_M}/INTERNAL_stability_full.json")
