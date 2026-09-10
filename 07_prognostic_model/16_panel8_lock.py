"""STAGE 1: 8-gene assay panel — redundancy variant of the locked 4-gene panel.

Nothing here overwrites outputs/ or outputs_compact/. All new files -> work/outputs_panel8/.

BLINDING STATUS
---------------
NOT BLIND. The DK panel-size curve (k=8: C=0.670, logrank p=0.032) was inspected before k=8
was chosen. The 4-gene panel remains the primary blind-validated result. This panel is an
assay-engineering variant, justified independently by discovery-only evidence: among panels
with <=10 genes, k=8 has the highest discovery-internal AUC (0.9107 vs 0.9063 at k=4).

PRE-DECLARED SUCCESS CRITERION (see work/PRESPEC_stage2.md)
-----------------------------------------------------------
leave-one-gene-out p_max < 0.05 — the median-split log-rank survives dropping ANY single gene.
The 4-gene panel fails this (p_max = 0.38). That fragility is the sole reason for going to 8.
Reported honestly whichever way it falls.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_panel8"
import os; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})

RED, BLUE, GREY = "#b2182b", "#2166ac", "#9a9a9a"


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


rng = np.random.default_rng(11)
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]

# identical pre-specified expression filter as 11_/13_
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
M, mnames, msrc = MOD.build(genes, R)
ev_tab = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
assert not np.isnan(beta_lock).any()

T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
print(f"genes={len(genes)}  JSE n={len(y)}  DK n={len(T)}  events(5y)={E5.sum()}  events(full)={E.sum()}")

order_lock = np.argsort(-np.abs(beta_lock))
K8, K4 = 8, 4
idx8, idx4 = order_lock[:K8], order_lock[:K4]
panel8 = [genes[i] for i in idx8]
panel4 = [genes[i] for i in idx4]
print(f"\n8-gene assay panel: {panel8}")
print(f"4-gene primary panel: {panel4}")


def panel_score(Z, idx_g):
    """Pre-specified scheme: mean of sign-weighted RINT expression."""
    return (Z[:, idx_g] * np.sign(beta_lock[idx_g])).mean(1)


sc8, sc4 = panel_score(Zd, idx8), panel_score(Zd, idx4)
sc_full = Zd @ beta_lock


def pack(sc, T_, E_, tag):
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "s": stats.zscore(sc)}), "time", "event")
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    hrs = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
    ci = np.exp(hrs.confidence_intervals_.iloc[0].values)
    return dict(endpoint=tag, C=float(c), hr_per_sd=float(np.exp(cph.params_.iloc[0])),
                cox_p=float(cph.summary["p"].iloc[0]), logrank_p=float(lr.p_value),
                hr_median_split=float(np.exp(hrs.params_.iloc[0])),
                hr_ci=[float(ci[0]), float(ci[1])])


res = {"blinding_status": "NOT BLIND - DK size-scan inspected before k=8 chosen; "
                          "4-gene panel remains the primary blind-validated result",
       "panel8": panel8, "panel4": panel4,
       "discovery_auc_k8": 0.9107, "discovery_auc_k4": 0.9063}

res["panel8_5yr"] = pack(sc8, T5, E5, "8-gene, 5-year OS (primary endpoint)")
res["panel8_full"] = pack(sc8, T, E, "8-gene, full follow-up")
res["panel4_5yr"] = pack(sc4, T5, E5, "4-gene, 5-year OS (comparator)")
res["panel4_full"] = pack(sc4, T, E, "4-gene, full follow-up")
res["full_signature_5yr"] = pack(sc_full, T5, E5, "full signature, 5-year OS")

print()
for k in ["panel8_5yr", "panel8_full", "panel4_5yr", "panel4_full", "full_signature_5yr"]:
    v = res[k]
    print(f"{v['endpoint']:<42s} C={v['C']:.3f}  HR/SD={v['hr_per_sd']:.2f} (p={v['cox_p']:.4f})  "
          f"split HR={v['hr_median_split']:.2f} [{v['hr_ci'][0]:.2f}-{v['hr_ci'][1]:.2f}]  logrank p={v['logrank_p']:.4f}")

# binary early-death vs long-survivor (their framing)
early, late = (T <= 24) & (E == 1), T > 60
bmask = early | late
res["binary_early_vs_long"] = dict(
    auc_panel8=float(roc_auc_score(early[bmask].astype(int), sc8[bmask])),
    auc_panel4=float(roc_auc_score(early[bmask].astype(int), sc4[bmask])),
    mw_p_panel8=float(stats.mannwhitneyu(sc8[early], sc8[late]).pvalue),
    n_early=int(early.sum()), n_long=int(late.sum()))
print("\nbinary early-vs-long:", res["binary_early_vs_long"])

# tertiles
q = np.quantile(sc8, [1 / 3, 2 / 3]); grp = np.digitize(sc8, q)
res["tertiles_5yr"] = dict(
    logrank_p_3grp=float(multivariate_logrank_test(T5, grp, E5).p_value),
    logrank_p_top_vs_bottom=float(logrank_test(T5[grp == 2], T5[grp == 0], E5[grp == 2], E5[grp == 0]).p_value))
print("tertiles:", res["tertiles_5yr"])

# multivariable with stage + histology
cl = pd.DataFrame({"time": T5, "event": E5, "score": stats.zscore(sc8)}, index=clin.index)
cl["stage_III_IV"] = (clin["STAGE"].astype(str) == "III-IV").astype(int)
cl["solid"] = (clin["FORM"].astype(str) == "Solid").astype(int)
u = clin["FORM"].astype(str).isin(["Solid", "Tubulocribriform"])
mm = CoxPHFitter().fit(cl.loc[u, ["time", "event", "score", "stage_III_IV", "solid"]], "time", "event")
res["multivariable_cox_5yr"] = {k: dict(HR=float(np.exp(mm.params_[k])), p=float(mm.summary.loc[k, "p"]))
                                for k in ["score", "stage_III_IV", "solid"]}
print("\nmultivariable (5-yr OS):\n", mm.summary[["exp(coef)", "p"]].to_string())

# ---------------- nulls ----------------
c_obs = res["panel8_5yr"]["C"]
randc = np.array([concordance_index_censored(
    E5.astype(bool), T5,
    (Zd[:, rng.choice(len(genes), K8, replace=False)] * rng.choice([-1.0, 1.0], K8)).mean(1))[0]
    for _ in range(2000)])
permc = np.array([concordance_index_censored(
    E5[o].astype(bool), T5[o], sc8)[0] for o in (rng.permutation(len(T)) for _ in range(2000))])
res["nulls"] = dict(random_signature_p=float((np.sum(randc >= c_obs) + 1) / 2001),
                    permutation_p=float((np.sum(permc >= c_obs) + 1) / 2001),
                    random_mean=float(randc.mean()), perm_mean=float(permc.mean()))
print("\nnulls:", res["nulls"])

# ---------------- THE DECISIVE TEST: leave-one-gene-out ----------------
logo = []
for j, g in enumerate(panel8):
    keep_i = np.array([idx8[m] for m in range(K8) if m != j])
    r = pack(panel_score(Zd, keep_i), T5, E5, f"drop {g}")
    logo.append(dict(dropped=g, C=r["C"], logrank_p=r["logrank_p"], hr=r["hr_median_split"]))
logo = pd.DataFrame(logo).sort_values("logrank_p", ascending=False)
logo.to_csv(f"{OUT}/leave_one_gene_out_panel8.csv", index=False)

logo4 = []
for j, g in enumerate(panel4):
    keep_i = np.array([idx4[m] for m in range(K4) if m != j])
    r = pack(panel_score(Zd, keep_i), T5, E5, f"drop {g}")
    logo4.append(dict(dropped=g, C=r["C"], logrank_p=r["logrank_p"], hr=r["hr_median_split"]))
logo4 = pd.DataFrame(logo4).sort_values("logrank_p", ascending=False)
logo4.to_csv(f"{OUT}/leave_one_gene_out_panel4.csv", index=False)

p_max8, p_max4 = float(logo.logrank_p.max()), float(logo4.logrank_p.max())
res["leave_one_gene_out"] = dict(panel8_p_max=p_max8, panel8_C_min=float(logo.C.min()),
                                 panel4_p_max=p_max4, panel4_C_min=float(logo4.C.min()),
                                 criterion="panel8_p_max < 0.05",
                                 PASSED=bool(p_max8 < 0.05))
print("\nleave-one-gene-out, 8-gene panel (DK, 5-yr OS):\n", logo.to_string(index=False))
print("\nleave-one-gene-out, 4-gene panel:\n", logo4.to_string(index=False))
print(f"\n*** PRE-DECLARED CRITERION: 8-gene p_max = {p_max8:.4f} < 0.05 ?  "
      f"{'PASS' if p_max8 < 0.05 else 'FAIL'}   (4-gene p_max = {p_max4:.4f})")

# ---------------- bootstrap CIs and paired delta-C ----------------
B = 2000
boot = {"C8": [], "C4": [], "dC": []}
n = len(T5)
for _ in range(B):
    b = rng.integers(0, n, n)
    if E5[b].sum() < 5:
        continue
    c8 = concordance_index_censored(E5[b].astype(bool), T5[b], sc8[b])[0]
    c4 = concordance_index_censored(E5[b].astype(bool), T5[b], sc4[b])[0]
    boot["C8"].append(c8); boot["C4"].append(c4); boot["dC"].append(c8 - c4)
res["bootstrap"] = {k: dict(mean=float(np.mean(v)),
                            lo=float(np.percentile(v, 2.5)), hi=float(np.percentile(v, 97.5)))
                    for k, v in boot.items()}
res["bootstrap"]["n_resamples"] = len(boot["dC"])
print("\nbootstrap (2000 patient resamples):")
for k in ["C8", "C4", "dC"]:
    b_ = res["bootstrap"][k]
    print(f"  {k}: {b_['mean']:+.3f}  [{b_['lo']:+.3f}, {b_['hi']:+.3f}]")

# ---------------- panel table ----------------
gene2mod = {}
for j, mn in enumerate(mnames):
    for i in np.where(M[:, j] != 0)[0]:
        gene2mod.setdefault(genes[i], []).append(mn)
ptab = pd.DataFrame({
    "gene": panel8,
    "beta_locked": beta_lock[idx8],
    "direction": np.where(beta_lock[idx8] > 0, "poor-prognosis", "good-prognosis"),
    "in_4gene_panel": [g in panel4 for g in panel8],
    "evidence_direction": ev_tab.reindex(panel8)["direction_consensus"].values,
    "n_evidence_sources": ev_tab.reindex(panel8)["n_sources"].values,
    "logo_logrank_p": [float(logo.set_index("dropped").loc[g, "logrank_p"]) for g in panel8],
    "modules": [";".join(gene2mod.get(g, [])[:6]) for g in panel8]})
ptab.to_csv(f"{OUT}/panel8_genes.csv", index=False)
print("\n8-gene panel:\n", ptab[["gene", "beta_locked", "direction", "in_4gene_panel",
                                "n_evidence_sources", "logo_logrank_p"]].to_string(index=False))

# ---------------- figures ----------------
def km(sc, T_, E_, title, fname):
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    cf = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
    hr = np.exp(cf.params_.iloc[0]); ci = np.exp(cf.confidence_intervals_.iloc[0].values)
    fig, ax = plt.subplots(figsize=(4.6, 3.9))
    for m_, lab, col in [(hi, "high risk", RED), (~hi, "low risk", BLUE)]:
        KaplanMeierFitter().fit(T_[m_], E_[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(
            ax=ax, color=col, ci_alpha=0.12)
    ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel("overall survival")
    ax.set_title(f"{title}\nHR = {hr:.2f} [{ci[0]:.2f}–{ci[1]:.2f}], log-rank p = {lr.p_value:.3f}", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    save(fig, fname)


km(sc8, T5, E5, "DK external cohort — 8-gene assay panel, 5-year OS", "PANEL8_km_5yr")
km(sc8, T, E, "DK external cohort — 8-gene assay panel, full follow-up", "PANEL8_km_full")

# the money figure: robustness comparison
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True)
for ax, (tab, lab, k_) in zip(axes, [(logo4, "4-gene panel", K4), (logo, "8-gene assay panel", K8)]):
    o = np.argsort(tab.logrank_p.values)[::-1]
    cols = [RED if p >= 0.05 else BLUE for p in tab.logrank_p.values[o]]
    ax.barh(range(len(tab)), tab.logrank_p.values[o], color=cols)
    ax.set_yticks(range(len(tab))); ax.set_yticklabels(tab.dropped.values[o], fontsize=8)
    ax.axvline(0.05, color="k", ls="--", lw=1)
    ax.set_xlabel("log-rank p after dropping this gene")
    ax.set_title(f"{lab}  (worst p = {tab.logrank_p.max():.3f})", fontsize=9)
axes[0].annotate("p = 0.05", (0.05, -0.7), fontsize=7, color="k")
fig.suptitle("Redundancy: does the result survive losing any single gene?", fontsize=10, y=1.02)
save(fig, "PANEL8_leave_one_gene_out")

# gene weights
fig, ax = plt.subplots(figsize=(5.0, 0.32 * K8 + 1.4))
o = np.argsort(ptab.beta_locked.values)
cols = [RED if b > 0 else BLUE for b in ptab.beta_locked.values[o]]
ax.barh(range(K8), ptab.beta_locked.values[o], color=cols)
ax.set_yticks(range(K8))
ax.set_yticklabels([f"{g}*" if g in panel4 else g for g in ptab.gene.values[o]])
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("locked posterior mean weight (log-hazard per SD)")
ax.set_title("8-gene assay panel — * = also in the 4-gene primary panel", fontsize=9)
save(fig, "PANEL8_gene_weights")

# nulls
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
for ax, (vals, lab, p_) in zip(axes, [(randc, "random 8-gene signatures", res["nulls"]["random_signature_p"]),
                                      (permc, "permuted DK survival", res["nulls"]["permutation_p"])]):
    ax.hist(vals, bins=40, color="#cccccc", edgecolor="none")
    ax.axvline(c_obs, color=RED, lw=1.6)
    ax.set_xlabel("C-index (DK, 5-year OS)"); ax.set_ylabel("count")
    ax.set_title(f"{lab}\nobserved {c_obs:.3f}, p = {p_:.4f}", fontsize=9)
save(fig, "PANEL8_nulls")

# 4 vs 8 side by side
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.7), sharey=True)
for ax, (sc, name) in zip(axes, [(sc4, "4-gene panel (primary, blind)"),
                                 (sc8, "8-gene assay panel (not blind)")]):
    hi = sc > np.median(sc)
    lr = logrank_test(T5[hi], T5[~hi], E5[hi], E5[~hi])
    hr = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "g": hi.astype(int)}),
                                  "time", "event").params_.iloc[0])
    c = concordance_index_censored(E5.astype(bool), T5, sc)[0]
    for m_, lab, col in [(hi, "high risk", RED), (~hi, "low risk", BLUE)]:
        KaplanMeierFitter().fit(T5[m_], E5[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(
            ax=ax, color=col, ci_alpha=0.1)
    ax.set_ylim(0, 1.02); ax.set_xlabel("months")
    ax.set_title(f"{name}\nC={c:.3f}  HR={hr:.2f}  p={lr.p_value:.3f}", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
axes[0].set_ylabel("overall survival (5-year)")
save(fig, "PANEL8_vs_panel4")

np.save(f"{OUT}/dk_panel8_score.npy", sc8)
np.save(f"{OUT}/dk_panel4_score.npy", sc4)
json.dump(res, open(f"{OUT}/PANEL8_validation.json", "w"), indent=2, default=float)
json.dump(dict(k=K8, scheme="sign", genes=panel8,
               weights=dict(zip(panel8, beta_lock[idx8].tolist())),
               blinding="NOT BLIND - assay variant; 4-gene panel is the primary blind result"),
          open(f"{OUT}/PANEL8_model.json", "w"), indent=2)
print("\nDONE ->", OUT)
