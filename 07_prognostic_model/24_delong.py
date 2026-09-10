"""
Stage 6 - is the DK AUC gap between our panel and Brayer's statistically significant?

The four signatures are scored on the SAME 44 DK patients, so the AUCs are correlated and an
unpaired comparison would be wrong. This uses DeLong's test for two correlated ROC curves
(DeLong, DeLong & Clarke-Pearson, Biometrics 1988), plus a paired bootstrap on delta-AUC as a
distribution-free check.

Endpoint is the discovery definition transferred unchanged to DK (13_compact_panel.py:205):
    early = (T <= 24 months) AND died      late = (T > 60 months)
Patients between 24 and 60 months are excluded, exactly as in the original analysis.

Nothing is refitted. Scores are the ones already locked on disk; brayer49 is recomputed with the
published Table 5 coefficients and the same JSE-only orientation used in 18_brayer_benchmark.py.
Reads only; writes to work/outputs_delong/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.fonttype"] = "none"

ROOT = f"{ACC_DATA_ROOT}"
sys.path.insert(0, HERE)
import importlib
d3 = importlib.import_module("03_data")

OUT = f"{ACC_DATA_ROOT}/model_outputs_delong"
RNG = np.random.default_rng(20260815)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


# ---------------- DeLong ----------------
def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(N, dtype=float)
    out[J] = T
    return out


def delong_cov(scores, y):
    """AUCs and their covariance matrix for k classifiers on the same samples.

    scores : (k, n) array;  y : (n,) 0/1 with 1 = positive class.
    """
    pos = scores[:, y == 1]
    neg = scores[:, y == 0]
    k, m = pos.shape
    n = neg.shape[1]
    tx = np.array([_midrank(p) for p in pos])
    ty = np.array([_midrank(q) for q in neg])
    tz = np.array([_midrank(np.concatenate([p, q])) for p, q in zip(pos, neg)])
    aucs = (tz[:, :m].sum(1) / m - (m + 1) / 2.0) / n
    v01 = (tz[:, :m] - tx) / n                 # structural components, positives
    v10 = 1.0 - (tz[:, m:] - ty) / m           # structural components, negatives
    s01 = np.cov(v01)
    s10 = np.cov(v10)
    cov = s01 / m + s10 / n
    return aucs, np.atleast_2d(cov)


def delong_test(s_a, s_b, y):
    """Two-sided p for AUC_a - AUC_b on paired samples."""
    aucs, cov = delong_cov(np.vstack([s_a, s_b]), y)
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var <= 0:
        return float(aucs[0]), float(aucs[1]), np.nan, 1.0
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    return float(aucs[0]), float(aucs[1]), float(z), float(2 * stats.norm.sf(abs(z)))


# ---------------- data ----------------
print("=" * 72)
print("Stage 6 - DeLong test, DK binary endpoint")
print("=" * 72)

D = d3.load("all")          # must match 18_brayer_benchmark.py:77 - the full gene universe,
                            # otherwise most of their 49 genes are missing from the lookup
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]

early = (T <= 24) & (E == 1)
late = T > 60
mask = early | late
lab = early[mask].astype(int)
print(f"DK patients scored : {len(T)}")
print(f"early (<=24mo, died): {int(early.sum())}   late (>60mo): {int(late.sum())}   "
      f"excluded (24-60mo): {int((~mask).sum())}")

s4 = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel4_score.npy")
s8 = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel8_score.npy")
sb14 = np.load(f"{ACC_DATA_ROOT}/model_outputs_brayer/dk_brayer14_score.npy")

# brayer49 was not cached - recompute exactly as 18_brayer_benchmark.py did
coef = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_brayer/brayer_table5_coefs.csv")
Zd_all = d3.rint(dk_raw).values
Zt_all = d3.rint(tr_raw).values
gi_all = {g: i for i, g in enumerate(genes_all)}
sub = coef[coef.gene.map(lambda g: g in gi_all)]
ii = [gi_all[g] for g in sub.gene]
ww = sub.coefficient.values
sb49 = Zd_all[:, ii] @ ww
sgn = 1.0 if roc_auc_score(y, Zt_all[:, ii] @ ww) >= 0.5 else -1.0
sb49 = sgn * sb49
print(f"brayer49 recomputed from {len(sub)}/49 available genes (orientation {sgn:+.0f}, JSE-only)")

CAND = {"ours_4gene": s4, "ours_8gene": s8, "brayer14": sb14, "brayer49": sb49}
for k, v in CAND.items():
    print(f"  {k:<12s} AUC = {roc_auc_score(lab, v[mask]):.3f}")

# ---------------- paired tests vs ours_4gene ----------------
print("\n" + "=" * 72)
print("DeLong, paired on the same 44 patients   (reference = ours_4gene)")
print("=" * 72)

res = {"endpoint": "DK early (<=24mo & died) vs late (>60mo)",
       "n_early": int(early.sum()), "n_late": int(late.sum()),
       "n_excluded_24_60mo": int((~mask).sum()),
       "method": "DeLong 1988, two correlated ROC curves, paired",
       "auc": {k: float(roc_auc_score(lab, v[mask])) for k, v in CAND.items()},
       "comparisons": {}}


def paired_boot(a, b, y_, n=10000):
    idx = np.arange(len(y_))
    d = []
    for _ in range(n):
        bi = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y_[bi])) < 2:
            continue
        d.append(roc_auc_score(y_[bi], a[bi]) - roc_auc_score(y_[bi], b[bi]))
    d = np.array(d)
    return float(d.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), \
        float(2 * min((d <= 0).mean(), (d >= 0).mean()))


ref = s4[mask]
for name in ["brayer14", "brayer49", "ours_8gene"]:
    other = CAND[name][mask]
    a1, a2, z, p = delong_test(ref, other, lab)
    bm, blo, bhi, bp = paired_boot(ref, other, lab)
    res["comparisons"][f"ours_4gene_vs_{name}"] = {
        "auc_ours": a1, "auc_other": a2, "delta": a1 - a2,
        "delong_z": z, "delong_p": p,
        "boot_delta_mean": bm, "boot_ci": [blo, bhi], "boot_p": bp}
    print(f"ours_4gene ({a1:.3f}) vs {name} ({a2:.3f}) : "
          f"dAUC = {a1 - a2:+.3f}  DeLong z = {z:+.2f}  p = {p:.4f}   "
          f"| bootstrap dAUC {bm:+.3f} [{blo:+.3f}, {bhi:+.3f}] p = {bp:.4f}")

# 8-gene as reference too, since it is the better of ours on this endpoint
print("\n(reference = ours_8gene)")
ref8 = s8[mask]
for name in ["brayer14", "brayer49"]:
    other = CAND[name][mask]
    a1, a2, z, p = delong_test(ref8, other, lab)
    bm, blo, bhi, bp = paired_boot(ref8, other, lab)
    res["comparisons"][f"ours_8gene_vs_{name}"] = {
        "auc_ours": a1, "auc_other": a2, "delta": a1 - a2,
        "delong_z": z, "delong_p": p,
        "boot_delta_mean": bm, "boot_ci": [blo, bhi], "boot_p": bp}
    print(f"ours_8gene ({a1:.3f}) vs {name} ({a2:.3f}) : "
          f"dAUC = {a1 - a2:+.3f}  DeLong z = {z:+.2f}  p = {p:.4f}   "
          f"| bootstrap dAUC {bm:+.3f} [{blo:+.3f}, {bhi:+.3f}] p = {bp:.4f}")

# each signature against chance
print("\nagainst chance (AUC = 0.5):")
res["vs_chance"] = {}
for k, v in CAND.items():
    a, cov = delong_cov(v[mask][None, :], lab)
    se = float(np.sqrt(cov[0, 0]))
    z = (a[0] - 0.5) / se
    p = float(2 * stats.norm.sf(abs(z)))
    res["vs_chance"][k] = {"auc": float(a[0]), "se": se, "z": float(z), "p": p,
                           "ci": [float(a[0] - 1.96 * se), float(a[0] + 1.96 * se)]}
    print(f"  {k:<12s} AUC = {a[0]:.3f} [{a[0]-1.96*se:.3f}, {a[0]+1.96*se:.3f}]  p = {p:.4f}")

# ---------------- figure ----------------
C = {"ours_4gene": "#B4436C", "ours_8gene": "#7D3C6B", "brayer14": "#3C6E8F", "brayer49": "#6E8FA8"}
fig, ax = plt.subplots(1, 2, figsize=(9.5, 4))

for k, v in CAND.items():
    s = v[mask]
    o = np.argsort(-s)
    tp = np.concatenate([[0], np.cumsum(lab[o]) / lab.sum()])
    fp = np.concatenate([[0], np.cumsum(1 - lab[o]) / (1 - lab).sum()])
    ax[0].plot(fp, tp, lw=2 if k.startswith("ours") else 1.4,
               ls="-" if k.startswith("ours") else "--", color=C[k],
               label=f"{k}  {roc_auc_score(lab, s):.3f}")
ax[0].plot([0, 1], [0, 1], ls=":", c="0.6", lw=1)
ax[0].set_xlabel("false positive rate"); ax[0].set_ylabel("true positive rate")
ax[0].set_title(f"DK external, early vs late ({int(early.sum())} vs {int(late.sum())})", fontsize=10)
ax[0].legend(frameon=False, fontsize=8, loc="lower right")

names = ["brayer14", "brayer49", "ours_8gene"]
d = [res["comparisons"][f"ours_4gene_vs_{n}"] for n in names]
yy = np.arange(len(names))
ax[1].axvline(0, c="0.6", lw=1, ls=":")
for i, (n, r) in enumerate(zip(names, d)):
    ax[1].plot([r["boot_ci"][0], r["boot_ci"][1]], [i, i], c="0.35", lw=1.6)
    ax[1].plot(r["delta"], i, "o", ms=7, color=C[n])
    ax[1].text(r["boot_ci"][1] + 0.02, i, f"DeLong p = {r['delong_p']:.3f}", va="center", fontsize=8)
ax[1].set_yticks(yy); ax[1].set_yticklabels([f"vs {n}" for n in names])
ax[1].set_xlabel("$\\Delta$AUC (ours_4gene $-$ comparator)")
ax[1].set_title("paired difference, 10,000-resample CI", fontsize=10)
ax[1].set_xlim(-0.35, 0.75)
for a_ in ax:
    a_.spines[["top", "right"]].set_visible(False)
save(fig, "DELONG_dk_roc_and_delta")

json.dump(res, open(f"{OUT}/DELONG_results.json", "w"), indent=1)
print(f"\nwrote {OUT}/DELONG_results.json")
