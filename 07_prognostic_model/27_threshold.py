"""
Stage 9 - deployability: single-sample scoring and a transferable threshold.

WHY
---
Two linked criticisms from the review:

  (#3) The rank-inverse-normal transform is cohort-relative. A patient's score depends on who else
       is in the batch, a single new patient cannot be scored at all, and the document's claim that
       the score is "invariant to platform, library size, and distributional shape" is overclaimed.
  (#6) The DK median split is computed inside DK, so it does not validate a clinical decision
       threshold - it only visualises the continuous score.

This script replaces both caveats with a concrete, testable procedure: freeze reference quantiles
from the discovery cohort and score each patient independently against them.

  (a) batch cohort-RINT     - what was reported; every DK patient influences every other
  (b) leave-one-out DK ref  - each DK patient scored against the other 53; measures how much a
                              score depends on cohort composition, without cross-platform effects
  (c) frozen JSE reference  - each DK patient mapped into the JSE rank space one at a time. This
                              is the deployable procedure, and it is what makes a JSE-derived
                              threshold meaningful in DK.

Sensitivity/specificity/PPV/NPV are then reported at the JSE-locked Youden cut-point under (c).
DK survival is used only to evaluate; the threshold is derived from JSE alone.

Reads only; writes to work/outputs_threshold/.
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
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sksurv.metrics import concordance_index_censored
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_threshold"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3)

PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


def rint_batch(X):
    r = stats.rankdata(X, axis=0)
    n = X.shape[0]
    return stats.norm.ppf((r - 0.375) / (n + 0.25))


def rint_vs_reference(x_new, X_ref):
    """Map one or more samples into the reference cohort's rank space, one sample at a time.

    Each row of `x_new` is transformed using only `X_ref` - never using the other rows. This is
    the operation a deployed assay performs on a single incoming patient.
    """
    ref = np.sort(X_ref, axis=0)
    n = ref.shape[0]
    out = np.empty(x_new.shape, dtype=float)
    for j in range(x_new.shape[1]):
        pos = np.searchsorted(ref[:, j], x_new[:, j], side="left")
        out[:, j] = stats.norm.ppf((pos + 0.5) / (n + 1.0))
    return out


def clopper_pearson(k, n, alpha=0.05):
    lo = stats.beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


print("=" * 74)
print("Stage 9 - deployability: single-sample scoring and threshold transfer")
print("=" * 74)

D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]
gi = {g: i for i, g in enumerate(genes_all)}
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta = dict(zip(LOCK["genes"], LOCK["beta"]))
sgn = np.sign([beta[g] for g in PANEL4])
Xd = dk_raw[PANEL4].values.astype(float)
Xt = tr_raw[PANEL4].values.astype(float)
print(f"DK n = {len(T)}, JSE n = {len(y)}, panel = {PANEL4}, signs = {sgn.astype(int).tolist()}")

res = {"panel": PANEL4, "n_dk": int(len(T)), "n_jse": int(len(y))}


# ================= the three scoring procedures =================
s_batch = (rint_batch(Xd) * sgn).mean(1)

s_loo = np.empty(len(T))
for i in range(len(T)):
    others = np.delete(np.arange(len(T)), i)
    s_loo[i] = (rint_vs_reference(Xd[[i]], Xd[others]) * sgn).mean()

s_frozen = (rint_vs_reference(Xd, Xt) * sgn).mean(1)

# JSE scores under its own batch transform, for the threshold
s_jse = (rint_batch(Xt) * sgn).mean(1)

ref = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel4_score.npy")
print(f"\nbatch cohort-RINT reproduces the locked score: "
      f"max |diff| = {np.abs(s_batch - ref).max():.2e}")


def surv_summary(sc, tag):
    c5 = concordance_index_censored(E5.astype(bool), T5, sc)[0]
    cf = concordance_index_censored(E.astype(bool), T, sc)[0]
    m = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "s": stats.zscore(sc)}),
                          "time", "event")
    s_ = m.summary.iloc[0]
    out = {"C_5yr": float(c5), "C_full": float(cf), "hr_per_sd": float(s_["exp(coef)"]),
           "ci": [float(s_["exp(coef) lower 95%"]), float(s_["exp(coef) upper 95%"])],
           "p": float(s_["p"])}
    print(f"  {tag:<26s} C(5y)={out['C_5yr']:.3f}  C(full)={out['C_full']:.3f}  "
          f"HR/SD={out['hr_per_sd']:.3f} [{out['ci'][0]:.3f}-{out['ci'][1]:.3f}] p={out['p']:.4g}")
    return out


print("\n" + "=" * 74)
print("discrimination under each scoring procedure")
print("=" * 74)
res["procedures"] = {
    "a_batch_cohort_rint": surv_summary(s_batch, "(a) batch cohort-RINT"),
    "b_loo_dk_reference": surv_summary(s_loo, "(b) leave-one-out DK ref"),
    "c_frozen_jse_reference": surv_summary(s_frozen, "(c) frozen JSE reference"),
}

res["agreement"] = {
    "batch_vs_loo_spearman": float(stats.spearmanr(s_batch, s_loo).correlation),
    "batch_vs_frozen_spearman": float(stats.spearmanr(s_batch, s_frozen).correlation),
    "loo_vs_frozen_spearman": float(stats.spearmanr(s_loo, s_frozen).correlation),
}
print(f"\nrank agreement: batch vs LOO rho = {res['agreement']['batch_vs_loo_spearman']:.4f}, "
      f"batch vs frozen rho = {res['agreement']['batch_vs_frozen_spearman']:.4f}")
print("  -> (b) near 1.0 means a single patient's score barely depends on cohort composition;")
print("     (c) is the number a deployed assay would actually achieve.")


# ================= JSE-locked threshold =================
print("\n" + "=" * 74)
print("JSE-locked Youden cut-point applied under the frozen-reference transform")
print("=" * 74)

cuts = np.unique(s_jse)
J = [((s_jse[y == 1] > c).mean() - (s_jse[y == 0] > c).mean()) for c in cuts]
thr = float(cuts[int(np.argmax(J))])
print(f"threshold derived on JSE only: {thr:+.5f}  (Youden J = {np.max(J):.3f})")

hi = s_frozen > thr
print(f"DK patients flagged high-risk: {int(hi.sum())} of {len(hi)} ({hi.mean()*100:.0f}%)")

# operating characteristics against the 5-year outcome, restricted to patients with
# a determinate outcome by 60 months (died by 60mo, or followed beyond 60mo)
died5 = (T <= 60) & (E == 1)
alive5 = T > 60
det = died5 | alive5
tp = int((hi & died5).sum()); fp = int((hi & alive5).sum())
fn = int((~hi & died5).sum()); tn = int((~hi & alive5).sum())
sens = tp / max(tp + fn, 1); spec = tn / max(tn + fp, 1)
ppv = tp / max(tp + fp, 1); npv = tn / max(tn + fn, 1)
oc = {
    "threshold": thr, "n_determinate": int(det.sum()),
    "n_died_by_60mo": int(died5.sum()), "n_alive_beyond_60mo": int(alive5.sum()),
    "TP": tp, "FP": fp, "FN": fn, "TN": tn,
    "sensitivity": sens, "sensitivity_ci": clopper_pearson(tp, tp + fn),
    "specificity": spec, "specificity_ci": clopper_pearson(tn, tn + fp),
    "PPV": ppv, "PPV_ci": clopper_pearson(tp, tp + fp),
    "NPV": npv, "NPV_ci": clopper_pearson(tn, tn + fn),
}
print(f"  determinate outcome at 60 months: n = {oc['n_determinate']} "
      f"({oc['n_died_by_60mo']} died, {oc['n_alive_beyond_60mo']} alive)")
print(f"  TP={tp} FP={fp} FN={fn} TN={tn}")
for k in ["sensitivity", "specificity", "PPV", "NPV"]:
    lo, hi_ci = oc[f"{k}_ci"]
    print(f"  {k:<12s} {oc[k]:.3f}  [{lo:.3f}-{hi_ci:.3f}]")
res["locked_threshold"] = oc

lr = logrank_test(T5[hi], T5[~hi], E5[hi], E5[~hi])
hrs = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "g": hi.astype(int)}),
                        "time", "event")
ci = np.exp(hrs.confidence_intervals_.iloc[0].values)
res["locked_threshold"]["survival"] = {
    "hr": float(np.exp(hrs.params_.iloc[0])), "ci": [float(ci[0]), float(ci[1])],
    "logrank_p": float(lr.p_value)}
print(f"  survival at this locked cut: HR {np.exp(hrs.params_.iloc[0]):.2f} "
      f"[{ci[0]:.2f}-{ci[1]:.2f}], log-rank p = {lr.p_value:.4f}")

# median split, for contrast only
med = s_batch > np.median(s_batch)
lr_med = logrank_test(T5[med], T5[~med], E5[med], E5[~med])
res["median_split_for_contrast"] = {
    "note": "DK-internal cut; a descriptive visualisation, not a validated decision rule",
    "n_high": int(med.sum()), "logrank_p": float(lr_med.p_value)}
print(f"  (contrast) DK-internal median split: n_high = {int(med.sum())}, "
      f"log-rank p = {lr_med.p_value:.4f} - descriptive only")

res["statement"] = (
    "Under a frozen JSE reference each DK sample is transformed independently of every other, so "
    "a single new patient can be scored and a discovery-derived threshold is transferable. The "
    "batch cohort-RINT figures reported elsewhere are not invariant to cohort composition; the "
    "frozen-reference C-index is the value a deployed assay would achieve.")

# ================= figure =================
C1, C2, C3 = "#B4436C", "#3C6E8F", "#7D8F3C"
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))

ax[0].scatter(s_batch, s_frozen, s=26, color=C1, alpha=.85, edgecolor="none")
lims = [min(s_batch.min(), s_frozen.min()), max(s_batch.max(), s_frozen.max())]
ax[0].plot(lims, lims, ls=":", c="0.6", lw=1)
ax[0].set_xlabel("batch cohort-RINT score"); ax[0].set_ylabel("frozen JSE-reference score")
ax[0].set_title(f"scoring procedures, rho = {res['agreement']['batch_vs_frozen_spearman']:.3f}",
                fontsize=9)

names = ["(a) batch\ncohort-RINT", "(b) LOO\nDK ref", "(c) frozen\nJSE ref"]
vals = [res["procedures"][k]["C_5yr"] for k in
        ["a_batch_cohort_rint", "b_loo_dk_reference", "c_frozen_jse_reference"]]
ax[1].bar(names, vals, color=[C2, C3, C1], width=.6)
ax[1].axhline(0.5, ls=":", c="0.6", lw=1)
ax[1].set_ylim(0.4, 0.8); ax[1].set_ylabel("C-index, 5-year OS")
ax[1].set_title("discrimination by procedure", fontsize=9)
for i, v in enumerate(vals):
    ax[1].text(i, v + 0.008, f"{v:.3f}", ha="center", fontsize=8)

from lifelines import KaplanMeierFitter
km = KaplanMeierFitter()
for m, lab, c in [(hi, f"high risk (n={int(hi.sum())})", C1),
                  (~hi, f"low risk (n={int((~hi).sum())})", C2)]:
    km.fit(T5[m], E5[m], label=lab)
    km.plot_survival_function(ax=ax[2], ci_show=False, color=c, lw=2)
ax[2].set_title(f"JSE-locked threshold, log-rank p = {lr.p_value:.3f}", fontsize=9)
ax[2].set_xlabel("months"); ax[2].set_ylabel("overall survival")
ax[2].legend(frameon=False, fontsize=8)
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
save(fig, "THRESHOLD_deployability")

pd.DataFrame({"dk_sample": dk_raw.index, "score_batch": s_batch, "score_loo": s_loo,
              "score_frozen_jse": s_frozen, "high_risk_locked": hi.astype(int),
              "T": T, "E": E}).to_csv(f"{OUT}/dk_scores_by_procedure.csv", index=False)
json.dump(res, open(f"{OUT}/THRESHOLD_results.json", "w"), indent=1)
print(f"\nwrote {OUT}/THRESHOLD_results.json")
