"""Internal-cohort nested-CV estimate and the as-built minus leakage-free difference, over repeats.

Runs N_REPEATS paired repeats of the design (LOO x3 seeds + 5x5 CV, fixed k=4/sign) for both
arms, and reports for each arm and for the paired difference: mean, sd, standard error, and a
2.5-97.5 percentile interval over draws. Both arms in a repeat use the same seed offset, so the
difference is taken within repeat.

Arm-level summaries also pool earlier independent draws of the same design (PRIOR_DRAWS; seed
offsets 0-3000). Offsets here start at 4000. The paired difference uses only this script's repeats.

Writes model_outputs_metrics/INTERNAL_stability_full.json (rewritten after every repeat) and
internal_stability_full_preds.npz.
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

# earlier independent draws of the same design, pooled into the arm-level summaries only
PRIOR_DRAWS = {"as_built": [0.881], "leakage_free": [0.767, 0.795, 0.785, 0.761]}


def one_run(arm, offset):
    """One full pass of the design for one arm. `offset` shifts every ADVI seed."""
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
        summ[a]["loo_mean"] = float(np.mean([r[a]["loo_auc"] for r in repeats]))
        summ[a]["cv5_mean"] = float(np.mean([r[a]["cv5_auc"] for r in repeats]))

    json.dump({
        "design": "LOO x3 seeds + 5x5 CV, fixed k=4/sign; both arms paired per seed offset",
        "n_repeats_paired": len(repeats), "offsets": [r["offset"] for r in repeats],
        "prior_independent_draws_pooled": PRIOR_DRAWS,
        "repeats": repeats, "summary": summ,
    }, open(f"{OUT_M}/INTERNAL_stability_full.json", "w"), indent=2)

np.savez(f"{OUT_M}/internal_stability_full_preds.npz", y=y.astype(int), **keep)

print("\n" + "=" * 78)
print(f"internal cohort, {N_REPEATS} paired repeats")
print("=" * 78)
for key in ["as_built", "leakage_free", "leak_paired"]:
    s = summ[key]
    print(f"  {key:<14s} {s['mean']:+.3f} +/- {s['sd']:.3f} (se {s['se']:.3f})  "
          f"[{s['pct_lo']:.3f}, {s['pct_hi']:.3f}] {s['interval_kind']}, n={s['n_draws']}")
print(f"\nsaved -> {OUT_M}/INTERNAL_stability_full.json")
