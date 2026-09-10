"""
Stage 8 - Cox diagnostics, robustness, and the weighted-vs-unweighted question.

WHY
---
The review flagged the multivariable model as fragile: three predictors on 21 events (~7 events per
parameter), reported with only HR and p, no confidence intervals, no assumption checks, and an
undeclared complete-case restriction. It also noted that the four Bayesian coefficients are not
actually used by the final predictor, which is an unweighted mean.

Everything here is a DIAGNOSTIC OF AN ALREADY-UNBLINDED COHORT, not a new validation. DK has been
seen; nothing below re-selects genes, re-fits the panel, or changes the locked score.

Reads only; writes to work/outputs_coxdx/.
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
from lifelines.statistics import proportional_hazard_test
from sksurv.metrics import concordance_index_censored
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_coxdx"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3)

RNG = np.random.default_rng(20260815)
PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


print("=" * 74)
print("Stage 8 - Cox diagnostics (sensitivity analysis on an unblinded cohort)")
print("=" * 74)

D = d3.load("all")
genes_all, T, E, clin = D["genes"], D["T"], D["E"], D["clin"]
dk_raw = D["dk_raw"]
Zd = d3.rint(dk_raw).values
gi = {g: i for i, g in enumerate(genes_all)}
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta = dict(zip(LOCK["genes"], LOCK["beta"]))
idx4 = np.array([gi[g] for g in PANEL4])
w4 = np.array([beta[g] for g in PANEL4])

score_sign = (Zd[:, idx4] * np.sign(w4)).mean(1)      # the deployed predictor
score_beta = Zd[:, idx4] @ w4                          # the Bayesian-weighted alternative

ref = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel4_score.npy")
assert np.allclose(score_sign, ref, atol=1e-9), "locked score did not reproduce"
print(f"locked 4-gene score reproduced exactly (max |diff| = {np.abs(score_sign-ref).max():.2e})")
print(f"DK n = {len(T)}, events at 5y = {int(E5.sum())}, events full = {int(E.sum())}")

res = {"note": "sensitivity analyses on an already-unblinded cohort; not new validation",
       "n_dk": int(len(T)), "events_5y": int(E5.sum()), "events_full": int(E.sum())}


# ================= 1. full multivariable table with CIs =================
print("\n" + "=" * 74)
print("1. multivariable Cox, 5-year OS - full table")
print("=" * 74)

df = pd.DataFrame({"time": T5, "event": E5, "score": stats.zscore(score_sign)}, index=clin.index)
df["stage_III_IV"] = (clin["STAGE"].astype(str) == "III-IV").astype(int)
df["solid"] = (clin["FORM"].astype(str) == "Solid").astype(int)
u = clin["FORM"].astype(str).isin(["Solid", "Tubulocribriform"]).values
cc = df.loc[u, ["time", "event", "score", "stage_III_IV", "solid"]].dropna()

n_cc, ev_cc = len(cc), int(cc["event"].sum())
epp = ev_cc / 3.0
print(f"complete-case n = {n_cc} of {len(df)}  (histology restricted to Solid/Tubulocribriform)")
print(f"events = {ev_cc}   parameters = 3   events per parameter = {epp:.1f}")
if epp < 10:
    print("  -> below the conventional 10 EPP guideline; coefficients are unstable by construction")

cph = CoxPHFitter().fit(cc, "time", "event")
tab = cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
tab.columns = ["HR", "CI_low", "CI_high", "p"]
print("\n" + tab.round(4).to_string())
tab.to_csv(f"{OUT}/multivariable_cox_full_table.csv")
res["multivariable"] = {"n_complete_case": n_cc, "n_total": int(len(df)), "events": ev_cc,
                        "n_parameters": 3, "events_per_parameter": round(epp, 2),
                        "table": {k: {"HR": float(tab.loc[k, "HR"]),
                                      "ci": [float(tab.loc[k, "CI_low"]), float(tab.loc[k, "CI_high"])],
                                      "p": float(tab.loc[k, "p"])} for k in tab.index}}

# univariable, for contrast
uni = {}
for v in ["score", "stage_III_IV", "solid"]:
    m = CoxPHFitter().fit(cc[["time", "event", v]], "time", "event")
    s = m.summary.iloc[0]
    uni[v] = {"HR": float(s["exp(coef)"]),
              "ci": [float(s["exp(coef) lower 95%"]), float(s["exp(coef) upper 95%"])],
              "p": float(s["p"])}
res["univariable"] = uni
print("\nunivariable, same complete-case set:")
for k, v in uni.items():
    print(f"  {k:<14s} HR {v['HR']:.3f} [{v['ci'][0]:.3f}-{v['ci'][1]:.3f}]  p={v['p']:.4g}")


# ================= 2. proportional hazards =================
print("\n" + "=" * 74)
print("2. proportional hazards (scaled Schoenfeld residuals)")
print("=" * 74)

ph = proportional_hazard_test(cph, cc, time_transform="rank")
ph_tab = ph.summary
print(ph_tab.round(4).to_string())
res["proportional_hazards"] = {str(i): {"test_stat": float(ph_tab.loc[i, "test_statistic"]),
                                        "p": float(ph_tab.loc[i, "p"])}
                               for i in ph_tab.index}
res["proportional_hazards"]["any_violation_p_lt_0.05"] = bool((ph_tab["p"] < 0.05).any())
print("  -> " + ("VIOLATION detected" if (ph_tab["p"] < 0.05).any()
                 else "no evidence of violation"))

sch = cph.compute_residuals(cc, "scaled_schoenfeld")
fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
for ax, v in zip(axes, ["score", "stage_III_IV", "solid"]):
    ax.scatter(np.arange(len(sch)), sch[v], s=18, color="#B4436C", alpha=.8, edgecolor="none")
    ax.axhline(0, ls=":", c="0.5", lw=1)
    z = np.polyfit(np.arange(len(sch)), sch[v], 1)
    ax.plot(np.arange(len(sch)), np.polyval(z, np.arange(len(sch))), c="0.35", lw=1.3, ls="--")
    ax.set_title(f"{v}  (p = {ph_tab.loc[v, 'p']:.3f})", fontsize=9)
    ax.set_xlabel("event order"); ax.set_ylabel("scaled Schoenfeld")
    ax.spines[["top", "right"]].set_visible(False)
save(fig, "COXDX_schoenfeld")


# ================= 3. functional form =================
print("\n" + "=" * 74)
print("3. functional form of the continuous score")
print("=" * 74)


def rcs_basis(x, knots):
    """Restricted cubic spline basis (Harrell), returns k-2 columns beyond the linear term."""
    k = len(knots)
    out = []
    for j in range(k - 2):
        def cub(t):
            return np.maximum(t, 0) ** 3
        num = (cub(x - knots[j])
               - cub(x - knots[k - 2]) * (knots[k - 1] - knots[j]) / (knots[k - 1] - knots[k - 2])
               + cub(x - knots[k - 1]) * (knots[k - 2] - knots[j]) / (knots[k - 1] - knots[k - 2]))
        out.append(num / (knots[k - 1] - knots[0]) ** 2)
    return np.column_stack(out)


s = cc["score"].values
kn = np.quantile(s, [0.10, 0.50, 0.90])
sp = rcs_basis(s, kn)
cc_sp = cc[["time", "event", "score", "stage_III_IV", "solid"]].copy()
for j in range(sp.shape[1]):
    cc_sp[f"spline{j}"] = sp[:, j]
m_lin = CoxPHFitter().fit(cc, "time", "event")
m_sp = CoxPHFitter().fit(cc_sp, "time", "event")
lr = 2 * (m_sp.log_likelihood_ - m_lin.log_likelihood_)
p_lr = stats.chi2.sf(lr, sp.shape[1])
print(f"linear vs restricted cubic spline (3 knots): LR = {lr:.3f}, df = {sp.shape[1]}, p = {p_lr:.4f}")
print("  -> " + ("non-linearity suggested" if p_lr < 0.05 else "linear form adequate"))
res["functional_form"] = {"lr_stat": float(lr), "df": int(sp.shape[1]), "p": float(p_lr),
                          "linear_adequate": bool(p_lr >= 0.05)}

cc_null = cc[["time", "event", "stage_III_IV", "solid"]]
null = CoxPHFitter().fit(cc_null, "time", "event")
mart = null.compute_residuals(cc_null, "martingale")["martingale"].reindex(cc.index).values
fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.4))
ax[0].scatter(s, mart, s=22, color="#3C6E8F", alpha=.8, edgecolor="none")
lo = np.polyfit(s, mart, 2)
xs = np.linspace(s.min(), s.max(), 100)
ax[0].plot(xs, np.polyval(lo, xs), c="0.35", lw=1.4)
ax[0].axhline(0, ls=":", c="0.5", lw=1)
ax[0].set_xlabel("4-gene score (z)"); ax[0].set_ylabel("martingale residual")
ax[0].set_title("functional form", fontsize=9)
ax[1].scatter(score_sign, score_beta, s=22, color="#B4436C", alpha=.8, edgecolor="none")
rr = stats.spearmanr(score_sign, score_beta)
ax[1].set_xlabel("sign-mean score (deployed)"); ax[1].set_ylabel("beta-weighted score")
ax[1].set_title(f"weighting schemes, rho = {rr.correlation:.3f}", fontsize=9)
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
save(fig, "COXDX_functional_form")


# ================= 4. ridge-penalised sensitivity =================
print("\n" + "=" * 74)
print("4. ridge-penalised Cox (coefficient stability)")
print("=" * 74)

pen_path = []
for pen in [0.0, 0.01, 0.05, 0.1, 0.5, 1.0]:
    m = CoxPHFitter(penalizer=pen, l1_ratio=0.0).fit(cc, "time", "event")
    row = {"penalizer": pen}
    for v in ["score", "stage_III_IV", "solid"]:
        row[f"HR_{v}"] = float(np.exp(m.params_[v]))
    row["C"] = float(m.concordance_index_)
    pen_path.append(row)
    print(f"  pen={pen:<5} " + "  ".join(f"HR_{v}={row[f'HR_{v}']:.3f}" for v in
                                         ["score", "stage_III_IV", "solid"])
          + f"  C={row['C']:.3f}")
pd.DataFrame(pen_path).to_csv(f"{OUT}/ridge_path.csv", index=False)
res["ridge_path"] = pen_path


# ================= 5. optimism-corrected C and calibration slope =================
print("\n" + "=" * 74)
print("5. optimism correction (Harrell bootstrap, B = 1000)")
print("=" * 74)

B = 1000
COVS = ["score", "stage_III_IV", "solid"]      # covariates only: lifelines misaligns if the
c_app = cph.concordance_index_                 # duration/event columns are passed to predict
opt_c, slopes = [], []
n = len(cc)
for _ in range(B):
    bi = RNG.choice(n, n, replace=True)
    boot = cc.iloc[bi].reset_index(drop=True)      # duplicated index breaks lifelines' internals
    if boot["event"].sum() < 5:
        continue
    try:
        mb = CoxPHFitter().fit(boot, "time", "event")
    except Exception:
        continue
    # NB: mb.concordance_index_ is unreliable on a resampled frame (lifelines silently reindexes
    # and misaligns predictions with outcomes, returning ~0.5). Compute it explicitly instead.
    lp_boot = np.log(mb.predict_partial_hazard(boot[COVS]).values)
    c_boot = concordance_index_censored(boot["event"].astype(bool).values,
                                        boot["time"].values, lp_boot)[0]
    lp_orig = np.log(mb.predict_partial_hazard(cc[COVS]).values)
    c_orig = concordance_index_censored(cc["event"].astype(bool).values,
                                        cc["time"].values, lp_orig)[0]
    opt_c.append(c_boot - c_orig)
    cal = pd.DataFrame({"time": cc["time"].values, "event": cc["event"].values, "lp": lp_orig})
    try:
        slopes.append(float(CoxPHFitter().fit(cal, "time", "event").params_["lp"]))
    except Exception:
        pass

optimism = float(np.mean(opt_c))
c_corr = float(c_app - optimism)
slope = float(np.mean(slopes))
print(f"apparent C          = {c_app:.3f}")
print(f"optimism            = {optimism:+.3f}   (B = {len(opt_c)} usable resamples)")
print(f"optimism-corrected C= {c_corr:.3f}")
print(f"calibration slope   = {slope:.3f}   (1.0 = perfectly calibrated; <1 = overfit)")
res["optimism"] = {"apparent_C": float(c_app), "optimism": optimism,
                   "corrected_C": c_corr, "calibration_slope": slope, "B": len(opt_c)}


# ================= 6. weighted vs unweighted =================
print("\n" + "=" * 74)
print("6. deployed sign-mean vs Bayesian beta-weighted score")
print("=" * 74)
print("Neither is chosen on DK: the scheme was fixed by the discovery parsimony rule.")


def summarise(sc, T_, E_):
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    m = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "s": stats.zscore(sc)}),
                          "time", "event")
    s_ = m.summary.iloc[0]
    return {"C": float(c), "hr_per_sd": float(s_["exp(coef)"]),
            "ci": [float(s_["exp(coef) lower 95%"]), float(s_["exp(coef) upper 95%"])],
            "p": float(s_["p"])}


cmp = {}
for nm, sc in [("sign_mean_deployed", score_sign), ("beta_weighted", score_beta)]:
    cmp[nm] = {"5yr": summarise(sc, T5, E5), "full": summarise(sc, T, E)}
    a = cmp[nm]["5yr"]; b = cmp[nm]["full"]
    print(f"  {nm:<20s} 5yr C={a['C']:.3f} HR/SD={a['hr_per_sd']:.3f} "
          f"[{a['ci'][0]:.3f}-{a['ci'][1]:.3f}] p={a['p']:.4g}  |  "
          f"full C={b['C']:.3f} p={b['p']:.4g}")
cmp["spearman_between_schemes"] = float(rr.correlation)
res["weighting_comparison"] = cmp
res["weighting_statement"] = (
    "All four locked coefficients are positive, so the deployed score is an unweighted mean of "
    "rank-normalised expression. The Bayesian multi-level model therefore functions as a "
    "gene-selection and direction-assignment device, not as the final predictive model.")

json.dump(res, open(f"{OUT}/COXDX_results.json", "w"), indent=1)
print(f"\nwrote {OUT}/COXDX_results.json")
