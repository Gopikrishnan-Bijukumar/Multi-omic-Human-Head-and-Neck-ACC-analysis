"""
Stage 11 - do the four genes contribute independently?

WHY
---
The consolidated document claimed "all four genes contribute" on the strength of four separate
single-gene AUCs (DSCAM 0.851, ODC1 0.854, NCAPG 0.776, CCNB2 0.832 in CCR2020). Individual
performance says nothing about *independent* contribution: NCAPG and CCNB2 are both G2/M genes and
are expected to be correlated, so four good marginal AUCs are consistent with one shared axis.

This replaces that claim with three things that can actually support it:

  1. the inter-gene correlation structure in both external cohorts;
  2. joint models - logistic for CCR2020 subtype, Cox for DK survival - with adjusted per-gene
     effects and confidence intervals;
  3. drop-one performance change with bootstrap CIs, plus conditional permutation importance.

Whatever survives is what gets written. Reads only; writes to work/outputs_joint/.
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
from lifelines import CoxPHFitter
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_joint"
EXT = f"{ACC_DATA_ROOT}/external_ccr2020"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3)

RNG = np.random.default_rng(20260815)
PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]
NBOOT = 2000


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


print("=" * 74)
print("Stage 11 - independent contribution of the four genes")
print("=" * 74)

res = {"panel": PANEL4, "n_boot": NBOOT,
       "why": "individual AUCs cannot establish independent contribution under collinearity"}

# ---------------- load both cohorts ----------------
D = d3.load("all")
T, E = D["T"], D["E"]
Zd_all = d3.rint(D["dk_raw"]).values
gid = {g: i for i, g in enumerate(D["genes"])}
Zdk = Zd_all[:, [gid[g] for g in PANEL4]]
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

expr = pd.read_excel(f"{EXT}/ACC_RNAseq.xlsx", sheet_name="ACC_Mitani_data_RPKM")
expr = expr[expr["gene_name"].map(lambda g: isinstance(g, str))].copy()
expr["gene_name"] = expr["gene_name"].astype(str)
expr = expr.groupby("gene_name").sum()
clin = pd.read_excel(f"{EXT}/CCR2020_Clinical.xlsx", sheet_name="Annotations").set_index("TID")
sm_ = [s for s in expr.columns if s in clin.index]
Zcc = d3.rint(expr.loc[PANEL4, sm_].T).values
y_acc1 = (clin.loc[sm_, "ACC Subtype"].astype(int) == 1).astype(int).values
print(f"DK n={len(T)} (events5y={int(E5.sum())});  CCR2020 n={len(sm_)} "
      f"({int(y_acc1.sum())} ACC-I)")

# ================= 1. correlation structure =================
print("\n" + "=" * 74)
print("1. inter-gene correlation")
print("=" * 74)
cor = {}
for nm, Z in [("DK", Zdk), ("CCR2020", Zcc)]:
    C = np.corrcoef(Z.T)
    cor[nm] = pd.DataFrame(C, index=PANEL4, columns=PANEL4).round(3).to_dict()
    off = C[np.triu_indices(4, 1)]
    print(f"\n{nm}  (mean |r| off-diagonal = {np.abs(off).mean():.3f}, max = {np.abs(off).max():.3f})")
    print(pd.DataFrame(C, index=PANEL4, columns=PANEL4).round(3).to_string())
    cor[f"{nm}_mean_abs_offdiag"] = float(np.abs(off).mean())
    cor[f"{nm}_max_abs_offdiag"] = float(np.abs(off).max())
res["correlation"] = cor
print(f"\nNCAPG-CCNB2 (both G2/M): DK r = {np.corrcoef(Zdk.T)[2,3]:.3f}, "
      f"CCR2020 r = {np.corrcoef(Zcc.T)[2,3]:.3f}")

# ================= 2. joint models =================
print("\n" + "=" * 74)
print("2. joint models - adjusted per-gene effects")
print("=" * 74)

X = sm.add_constant(Zcc)
lg = sm.Logit(y_acc1, X).fit(disp=0)
ci = lg.conf_int()
joint_cc = {}
print("\nCCR2020, joint logistic (outcome = ACC-I):")
for i, g in enumerate(PANEL4, start=1):
    joint_cc[g] = {"beta": float(lg.params[i]),
                   "or": float(np.exp(lg.params[i])),
                   "ci": [float(np.exp(ci[i][0])), float(np.exp(ci[i][1]))],
                   "p": float(lg.pvalues[i])}
    v = joint_cc[g]
    print(f"  {g:<8s} OR {v['or']:8.3f}  [{v['ci'][0]:.3f}-{v['ci'][1]:.3f}]  p = {v['p']:.4f}")
n_sig_cc = sum(v["p"] < 0.05 for v in joint_cc.values())
print(f"  -> {n_sig_cc}/4 genes retain p < 0.05 when adjusted for the other three")

dfd = pd.DataFrame(Zdk, columns=PANEL4)
dfd["time"], dfd["event"] = T5, E5
cph = CoxPHFitter().fit(dfd, "time", "event")
joint_dk = {}
print("\nDK, joint Cox (5-year OS):")
for g in PANEL4:
    s_ = cph.summary.loc[g]
    joint_dk[g] = {"hr": float(s_["exp(coef)"]),
                   "ci": [float(s_["exp(coef) lower 95%"]), float(s_["exp(coef) upper 95%"])],
                   "p": float(s_["p"])}
    v = joint_dk[g]
    print(f"  {g:<8s} HR {v['hr']:6.3f}  [{v['ci'][0]:.3f}-{v['ci'][1]:.3f}]  p = {v['p']:.4f}")
n_sig_dk = sum(v["p"] < 0.05 for v in joint_dk.values())
print(f"  -> {n_sig_dk}/4 genes retain p < 0.05 when adjusted for the other three")
res["joint_models"] = {"ccr2020_logistic": joint_cc, "n_significant_ccr2020": n_sig_cc,
                       "dk_cox": joint_dk, "n_significant_dk": n_sig_dk}

# ================= 3. drop-one and conditional permutation =================
print("\n" + "=" * 74)
print("3. drop-one performance change (bootstrap CI) and permutation importance")
print("=" * 74)


def auc_of(Z, cols):
    return roc_auc_score(y_acc1, Z[:, cols].mean(1))


def c_of(Z, cols):
    return concordance_index_censored(E5.astype(bool), T5, Z[:, cols].mean(1))[0]


full_cc = auc_of(Zcc, list(range(4)))
full_dk = c_of(Zdk, list(range(4)))
print(f"full 4-gene: CCR2020 AUC = {full_cc:.4f},  DK C = {full_dk:.4f}")

drop = {}
for j, g in enumerate(PANEL4):
    rest = [i for i in range(4) if i != j]
    d_cc = auc_of(Zcc, rest) - full_cc
    d_dk = c_of(Zdk, rest) - full_dk
    bc, bd = [], []
    for _ in range(NBOOT):
        b1 = RNG.choice(len(y_acc1), len(y_acc1), replace=True)
        if len(np.unique(y_acc1[b1])) > 1:
            bc.append(roc_auc_score(y_acc1[b1], Zcc[b1][:, rest].mean(1))
                      - roc_auc_score(y_acc1[b1], Zcc[b1][:, :].mean(1)))
        b2 = RNG.choice(len(T5), len(T5), replace=True)
        if E5[b2].sum() > 3:
            try:
                bd.append(concordance_index_censored(E5[b2].astype(bool), T5[b2],
                                                     Zdk[b2][:, rest].mean(1))[0]
                          - concordance_index_censored(E5[b2].astype(bool), T5[b2],
                                                       Zdk[b2].mean(1))[0])
            except Exception:
                pass
    # conditional permutation: shuffle this gene, keep the others intact
    perm_cc, perm_dk = [], []
    for _ in range(500):
        Zp = Zcc.copy(); Zp[:, j] = RNG.permutation(Zp[:, j])
        perm_cc.append(auc_of(Zp, list(range(4))))
        Zq = Zdk.copy(); Zq[:, j] = RNG.permutation(Zq[:, j])
        perm_dk.append(c_of(Zq, list(range(4))))
    drop[g] = {
        "ccr2020_delta_auc": float(d_cc),
        "ccr2020_delta_ci": [float(np.percentile(bc, 2.5)), float(np.percentile(bc, 97.5))],
        "dk_delta_C": float(d_dk),
        "dk_delta_ci": [float(np.percentile(bd, 2.5)), float(np.percentile(bd, 97.5))],
        "ccr2020_perm_drop": float(full_cc - np.mean(perm_cc)),
        "dk_perm_drop": float(full_dk - np.mean(perm_dk)),
    }
    v = drop[g]
    print(f"  drop {g:<8s} CCR2020 dAUC {v['ccr2020_delta_auc']:+.4f} "
          f"[{v['ccr2020_delta_ci'][0]:+.3f},{v['ccr2020_delta_ci'][1]:+.3f}]   "
          f"DK dC {v['dk_delta_C']:+.4f} "
          f"[{v['dk_delta_ci'][0]:+.3f},{v['dk_delta_ci'][1]:+.3f}]")
res["drop_one"] = drop
print("\n  conditional permutation (mean performance lost when that gene alone is shuffled):")
for g in PANEL4:
    print(f"    {g:<8s} CCR2020 {drop[g]['ccr2020_perm_drop']:+.4f}   "
          f"DK {drop[g]['dk_perm_drop']:+.4f}")

n_ci_excl_cc = sum(drop[g]["ccr2020_delta_ci"][1] < 0 for g in PANEL4)
n_ci_excl_dk = sum(drop[g]["dk_delta_ci"][1] < 0 for g in PANEL4)
res["summary"] = {
    "genes_with_significant_adjusted_effect_ccr2020": n_sig_cc,
    "genes_with_significant_adjusted_effect_dk": n_sig_dk,
    "genes_whose_removal_significantly_hurts_ccr2020": n_ci_excl_cc,
    "genes_whose_removal_significantly_hurts_dk": n_ci_excl_dk,
    "statement": (
        "Individual AUCs do not establish independent contribution. Report the adjusted joint "
        "effects and drop-one intervals; where the panel behaves as a single correlated axis, say "
        "so rather than claiming four independent contributors.")}
print(f"\n  genes whose removal significantly hurts (CI excludes 0): "
      f"CCR2020 {n_ci_excl_cc}/4, DK {n_ci_excl_dk}/4")

# ================= figure =================
fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.6))
im = ax[0].imshow(np.corrcoef(Zcc.T), vmin=-1, vmax=1, cmap="RdBu_r")
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(PANEL4, rotation=45, ha="right")
ax[0].set_yticks(range(4)); ax[0].set_yticklabels(PANEL4)
for i in range(4):
    for j in range(4):
        ax[0].text(j, i, f"{np.corrcoef(Zcc.T)[i,j]:.2f}", ha="center", va="center", fontsize=7)
ax[0].set_title("gene correlation, CCR2020", fontsize=9)
plt.colorbar(im, ax=ax[0], fraction=.046)

yy = np.arange(4)
for a, key, cikey, lab in [(ax[1], "ccr2020_delta_auc", "ccr2020_delta_ci", "$\\Delta$AUC, CCR2020"),
                           (ax[2], "dk_delta_C", "dk_delta_ci", "$\\Delta$C, DK")]:
    a.axvline(0, ls=":", c="0.6", lw=1)
    for i, g in enumerate(PANEL4):
        lo, hi = drop[g][cikey]
        a.plot([lo, hi], [i, i], c="0.35", lw=1.6)
        a.plot(drop[g][key], i, "o", ms=7, color="#B4436C")
    a.set_yticks(yy); a.set_yticklabels([f"drop {g}" for g in PANEL4])
    a.set_xlabel(lab); a.set_title("effect of removing each gene", fontsize=9)
    a.spines[["top", "right"]].set_visible(False)
save(fig, "JOINT_contribution")

json.dump(res, open(f"{OUT}/JOINT_results.json", "w"), indent=1)
print(f"\nwrote {OUT}/JOINT_results.json")
