"""
Stage 10 - null panels matched on abundance, variance and sign structure.

WHY
---
The review accepted the finite-permutation p-value form but objected to how the null panels were
drawn. In `13_compact_panel.py` the null is

    panel_score(Zd, rng.choice(len(genes), K, replace=False), rng.choice([-1., 1.], K), "sign")

which differs from the real panel in two ways that both make the null easier to beat:

  * random +/-1 signs, whereas the real panel is all-positive and therefore coherent;
  * no matching on expression abundance or variance, so a low p-value could partly reflect having
    picked well-measured, variable genes rather than these four genes specifically.

This script re-runs both nulls with the sign structure preserved and the panels matched on
abundance and variance, at B = 10,000, and reports the unmatched version alongside so the effect
of matching is visible. The counting form (b+1)/(B+1) is retained, and b and B are reported.

Applied to both external results: the DK 5-year C-index and the CCR2020 ACC-I subtype AUC.

Reads only; writes to work/outputs_nulls/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import sys
import importlib.util
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sksurv.metrics import concordance_index_censored
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_nulls"
EXT = f"{ACC_DATA_ROOT}/external_ccr2020"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3)

RNG = np.random.default_rng(20260815)
B = 10000
PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


def matched_bins(mean_v, var_v, n_bins=10):
    """Joint abundance x variance bin id for every gene."""
    mb = np.digitize(mean_v, np.quantile(mean_v, np.linspace(0, 1, n_bins + 1)[1:-1]))
    vb = np.digitize(var_v, np.quantile(var_v, np.linspace(0, 1, n_bins + 1)[1:-1]))
    return mb * n_bins + vb


def draw_matched(bins, panel_idx, pool_by_bin, rng):
    """One null panel: same size, one gene drawn from each real gene's abundance/variance bin."""
    out = []
    for gidx in panel_idx:
        pool = pool_by_bin[bins[gidx]]
        pick = pool[rng.integers(len(pool))]
        while pick in out:
            pick = pool[rng.integers(len(pool))]
        out.append(pick)
    return np.array(out)


def pval(null, obs, higher_is_better=True):
    b = int(np.sum(null >= obs)) if higher_is_better else int(np.sum(null <= obs))
    return {"b": b, "B": int(len(null)), "p": float((b + 1) / (len(null) + 1)),
            "null_mean": float(np.mean(null)), "null_p95": float(np.percentile(null, 95))}


print("=" * 74)
print("Stage 10 - matched null panels, B = 10,000")
print("=" * 74)

res = {"B": B, "panel": PANEL4,
       "null_design": {
           "sign_structure": "all-positive, matching the real panel (previously random +/-1)",
           "matching": "joint decile bins of mean expression and variance within the cohort",
           "p_value_form": "(b+1)/(B+1), b and B reported"}}

# ================= DK: 5-year C-index =================
D = d3.load("all")
genes_all, T, E = D["genes"], D["T"], D["E"]
dk_raw = D["dk_raw"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (D["tr_raw"].mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
gi = {g: i for i, g in enumerate(genes)}
Xd = dk_raw[genes]
Zd = d3.rint(Xd).values
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
pidx = np.array([gi[g] for g in PANEL4])

obs_c = concordance_index_censored(E5.astype(bool), T5, Zd[:, pidx].mean(1))[0]
print(f"\nDK universe: {len(genes)} genes "
      f"(filter: DK mean>1, DK detected in >60%, JSE mean>1)")
print(f"observed 4-gene C-index (5-year OS) = {obs_c:.4f}")

mean_v, var_v = Xd.values.mean(0), Xd.values.var(0)
bins = matched_bins(mean_v, var_v)
pool_by_bin = {b_: np.where(bins == b_)[0] for b_ in np.unique(bins)}
print("  panel gene bins: " + ", ".join(
    f"{g}(n={len(pool_by_bin[bins[gi[g]]])})" for g in PANEL4))

null_m, null_u = np.empty(B), np.empty(B)
for i in range(B):
    idx = draw_matched(bins, pidx, pool_by_bin, RNG)
    null_m[i] = concordance_index_censored(E5.astype(bool), T5, Zd[:, idx].mean(1))[0]
    idx = RNG.choice(len(genes), 4, replace=False)
    null_u[i] = concordance_index_censored(E5.astype(bool), T5, Zd[:, idx].mean(1))[0]
    if (i + 1) % 2500 == 0:
        print(f"  {i+1}/{B}")

dk = {"observed_C": float(obs_c), "n_universe": len(genes),
      "matched": pval(null_m, obs_c), "unmatched_same_signs": pval(null_u, obs_c)}
res["dk_c_index_5yr"] = dk
print(f"  matched   null: mean {dk['matched']['null_mean']:.4f}  "
      f"b={dk['matched']['b']}/{B}  p={dk['matched']['p']:.4f}")
print(f"  unmatched null: mean {dk['unmatched_same_signs']['null_mean']:.4f}  "
      f"b={dk['unmatched_same_signs']['b']}/{B}  p={dk['unmatched_same_signs']['p']:.4f}")

# ================= CCR2020: ACC-I subtype AUC =================
print("\n" + "-" * 74)
print("CCR2020 subtype AUC")
expr = pd.read_excel(f"{EXT}/ACC_RNAseq.xlsx", sheet_name="ACC_Mitani_data_RPKM")
expr = expr[expr["gene_name"].map(lambda g: isinstance(g, str))].copy()
expr["gene_name"] = expr["gene_name"].astype(str)
expr = expr.groupby("gene_name").sum()
clin = pd.read_excel(f"{EXT}/CCR2020_Clinical.xlsx", sheet_name="Annotations").set_index("TID")
samples = [s for s in expr.columns if s in clin.index]
expr = expr[samples]
y_acc1 = (clin.loc[samples, "ACC Subtype"].astype(int) == 1).astype(int).values

# Universe: detected (RPKM > 0) in at least 25% of samples. A stricter RPKM > 1 filter would
# exclude DSCAM itself (mean RPKM 0.42, above 1 in only 11% of samples), which would make the
# null incomparable to the real panel. Abundance matching, not filtering, is what controls for
# expression level here - so a low-abundance panel gene is matched against low-abundance nulls.
u2 = (expr.values > 0).mean(1) >= 0.25
g2 = expr.index[u2].tolist()
missing = [g for g in PANEL4 if g not in set(g2)]
assert not missing, f"panel genes excluded by the universe filter: {missing}"
gi2 = {g: i for i, g in enumerate(g2)}
X2 = expr.loc[g2].T
Z2 = d3.rint(X2).values
p2 = np.array([gi2[g] for g in PANEL4])
obs_auc = roc_auc_score(y_acc1, Z2[:, p2].mean(1))
print(f"universe: {len(g2)} genes; observed 4-gene AUC = {obs_auc:.4f} "
      f"({int(y_acc1.sum())} ACC-I vs {int((1-y_acc1).sum())} ACC-II)")

m2, v2 = X2.values.mean(0), X2.values.var(0)
bins2 = matched_bins(m2, v2)
pool2 = {b_: np.where(bins2 == b_)[0] for b_ in np.unique(bins2)}
nm2, nu2 = np.empty(B), np.empty(B)
for i in range(B):
    idx = draw_matched(bins2, p2, pool2, RNG)
    nm2[i] = roc_auc_score(y_acc1, Z2[:, idx].mean(1))
    idx = RNG.choice(len(g2), 4, replace=False)
    nu2[i] = roc_auc_score(y_acc1, Z2[:, idx].mean(1))
    if (i + 1) % 2500 == 0:
        print(f"  {i+1}/{B}")

cc = {"observed_AUC": float(obs_auc), "n_universe": len(g2),
      "matched": pval(nm2, obs_auc), "unmatched_same_signs": pval(nu2, obs_auc)}
res["ccr2020_subtype_auc"] = cc
print(f"  matched   null: mean {cc['matched']['null_mean']:.4f}  "
      f"b={cc['matched']['b']}/{B}  p={cc['matched']['p']:.4f}")
print(f"  unmatched null: mean {cc['unmatched_same_signs']['null_mean']:.4f}  "
      f"b={cc['unmatched_same_signs']['b']}/{B}  p={cc['unmatched_same_signs']['p']:.4f}")

# ================= figure =================
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
for a, (nm, nu, obs, lab, p_) in zip(ax, [
        (null_m, null_u, obs_c, "DK C-index (5-year OS)", dk["matched"]["p"]),
        (nm2, nu2, obs_auc, "CCR2020 ACC-I AUC", cc["matched"]["p"])]):
    a.hist(nu, bins=60, color="0.85", label="unmatched null", edgecolor="none")
    a.hist(nm, bins=60, color="#3C6E8F", alpha=.75, label="matched null", edgecolor="none")
    a.axvline(obs, color="#B4436C", lw=2, label=f"observed {obs:.3f}")
    a.set_xlabel(lab); a.set_ylabel("null panels")
    a.set_title(f"matched p = {p_:.4f}", fontsize=9)
    a.legend(frameon=False, fontsize=8)
    a.spines[["top", "right"]].set_visible(False)
save(fig, "NULLS_matched_vs_unmatched")

json.dump(res, open(f"{OUT}/NULLS_results.json", "w"), indent=1)
print(f"\nwrote {OUT}/NULLS_results.json")
