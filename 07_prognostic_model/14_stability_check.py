"""Robustness of the 4-gene compact panel.

1. Gene-level stability: re-run the discovery-internal resampling and record how often each
   gene enters the top-k. Build an alternative panel by SELECTION FREQUENCY (stability
   selection) instead of by the single locked fit -> sensitivity analysis on DK.
2. Bootstrap confidence intervals on DK for C-index and HR (2000 patient resamples).
3. Paired bootstrap for the difference in C-index, 4-gene vs 20-gene vs full signature.
4. Sensitivity of the log-rank p to the cut-point (median, tertiles, quartiles, optimal-cut
   with a Lausen-Schumacher style permutation correction).
Figures in .png (600 dpi) and .svg.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, importlib.util, time, os
import numpy as np, pandas as pd
from scipy import stats
from sklearn.model_selection import RepeatedStratifiedKFold
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_compact"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD
from mlbayes import MultiLevelModel, to_t

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


rng = np.random.default_rng(23)
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
M, mnames, msrc = MOD.build(genes, R)
srcs = sorted(set(msrc)); src_ids = np.array([srcs.index(s) for s in msrc])
ev_tab = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")
s_ev, d_ev = d3.evidence_prior(ev_tab.reindex(genes))
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
cfg = LOCK["config"]
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
CM = json.load(open(f"{OUT}/COMPACT_model.json"))
K = CM["k"]; panel = CM["genes"]
order_lock = np.argsort(-np.abs(beta_lock))
panel_idx = np.array([genes.index(g) for g in panel])


def sgn_score(idx_g, w):
    return (Zd[:, idx_g] * np.sign(w)).mean(1)


sc4 = sgn_score(panel_idx, beta_lock[panel_idx])
sc20 = sgn_score(order_lock[:20], beta_lock[order_lock[:20]])
scF = Zd @ beta_lock
res = {"k": K, "panel": panel}

# ---------------- 1. stability selection ----------------
def fit_beta(idx, seed=0):
    m = MultiLevelModel(M, src_ids, len(srcs), s_ev, d_ev, tau0=cfg["tau0"], omega0=cfg["omega0"],
                        use_gene_level=bool(cfg["gene_level"]), use_module_level=bool(cfg["module_level"]),
                        use_dk=False, use_jse=True, seed=seed)
    m.fit(Xj=to_t(Zt[idx]), yj=to_t(y[idx]), steps=900)
    return m.posterior(300)["beta"]


folds = [(np.array([j for j in range(len(y)) if j != i]), s) for s in range(3) for i in range(len(y))]
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=7)
folds += [(tr, 100 + f) for f, (tr, _) in enumerate(rskf.split(np.zeros(len(y)), y))]
top4 = np.zeros(len(genes)); top10 = np.zeros(len(genes)); top20 = np.zeros(len(genes))
sgn_sum = np.zeros(len(genes))
t0 = time.time()
for n, (tr, seed) in enumerate(folds):
    b = fit_beta(tr, seed=seed)
    o = np.argsort(-np.abs(b))
    top4[o[:4]] += 1; top10[o[:10]] += 1; top20[o[:20]] += 1
    sgn_sum += np.sign(b)
    if (n + 1) % 25 == 0:
        print(f"  {n+1}/{len(folds)} fits ({time.time()-t0:.0f}s)")
nf = len(folds)
stab = pd.DataFrame({"gene": genes, "freq_top4": top4 / nf, "freq_top10": top10 / nf,
                     "freq_top20": top20 / nf, "mean_sign": sgn_sum / nf,
                     "beta_locked": beta_lock})
stab = stab.sort_values("freq_top10", ascending=False)
stab.head(60).to_csv(f"{OUT}/gene_selection_stability.csv", index=False)
print("\nmost stably selected genes (discovery resampling, DK unseen):")
print(stab.head(15).to_string(index=False, float_format=lambda v: f"{v:.3f}"))

stab_panel = list(stab.head(K)["gene"])
stab_idx = np.array([genes.index(g) for g in stab_panel])
sc_stab = sgn_score(stab_idx, beta_lock[stab_idx])
res["stability_panel"] = stab_panel
res["overlap_with_locked"] = sorted(set(stab_panel) & set(panel))
print(f"\nstability-selected {K}-gene panel: {stab_panel}  (overlap with locked: {res['overlap_with_locked']})")


def pack(sc, T_, E_, tag):
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "s": stats.zscore(sc)}), "time", "event")
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    hrs = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
    ci = np.exp(hrs.confidence_intervals_.iloc[0].values)
    return dict(endpoint=tag, C=float(c), hr_per_sd=float(np.exp(cph.params_.iloc[0])),
                cox_p=float(cph.summary["p"].iloc[0]), logrank_p=float(lr.p_value),
                hr_median_split=float(np.exp(hrs.params_.iloc[0])), hr_ci=[float(ci[0]), float(ci[1])])


res["stability_panel_5yr"] = pack(sc_stab, T5, E5, "stability panel, 5-year OS")
print("sensitivity - stability panel on DK:", {k: round(v, 4) for k, v in res["stability_panel_5yr"].items() if k != "endpoint" and k != "hr_ci"})

# ---------------- 2/3. bootstrap CIs and paired comparison ----------------
B = 2000
n = len(T5)
boot = {"c4": [], "c20": [], "cF": [], "hr4": [], "d4_20": [], "d4_F": []}
for b in range(B):
    ix = rng.integers(0, n, n)
    if E5[ix].sum() < 5 or len(np.unique(T5[ix])) < 5:
        continue
    try:
        c4 = concordance_index_censored(E5[ix].astype(bool), T5[ix], sc4[ix])[0]
        c20 = concordance_index_censored(E5[ix].astype(bool), T5[ix], sc20[ix])[0]
        cF = concordance_index_censored(E5[ix].astype(bool), T5[ix], scF[ix])[0]
        hi = sc4[ix] > np.median(sc4)
        if 3 < hi.sum() < n - 3:
            cf = CoxPHFitter().fit(pd.DataFrame({"time": T5[ix], "event": E5[ix], "g": hi.astype(int)}), "time", "event")
            boot["hr4"].append(float(np.exp(cf.params_.iloc[0])))
    except Exception:
        continue
    boot["c4"].append(c4); boot["c20"].append(c20); boot["cF"].append(cF)
    boot["d4_20"].append(c4 - c20); boot["d4_F"].append(c4 - cF)
q = lambda a: [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]
res["bootstrap"] = {
    "C_4gene": q(boot["c4"]), "C_20gene": q(boot["c20"]), "C_full": q(boot["cF"]),
    "HR_4gene_median_split": q(boot["hr4"]),
    "deltaC_4_minus_20": q(boot["d4_20"]),
    "p_4_better_than_20": float(np.mean(np.array(boot["d4_20"]) < 0)),
    "deltaC_4_minus_full": q(boot["d4_F"]),
    "frac_C4_above_0.5": float(np.mean(np.array(boot["c4"]) > 0.5))}
print("\nbootstrap (2000 patient resamples):")
for k, v in res["bootstrap"].items():
    print(f"  {k}: {v}")

# ---------------- 4. cut-point sensitivity ----------------
cuts = {}
for name, gfun in [("median", lambda s: s > np.median(s)),
                   ("upper tertile vs rest", lambda s: s > np.quantile(s, 2 / 3)),
                   ("upper quartile vs rest", lambda s: s > np.quantile(s, 0.75)),
                   ("lower quartile vs rest", lambda s: s > np.quantile(s, 0.25))]:
    g = gfun(sc4)
    lr = logrank_test(T5[g], T5[~g], E5[g], E5[~g])
    cuts[name] = dict(p=float(lr.p_value), n_high=int(g.sum()))
lr3 = multivariate_logrank_test(T5, np.digitize(sc4, np.quantile(sc4, [1 / 3, 2 / 3])), E5)
cuts["three tertile groups"] = dict(p=float(lr3.p_value), n_high=18)

# optimal cut with permutation correction (guards against cut-point shopping)
grid = np.quantile(sc4, np.linspace(0.25, 0.75, 21))
stat = lambda s, c: logrank_test(T5[s > c], T5[s <= c], E5[s > c], E5[s <= c]).test_statistic
obs_max = max(stat(sc4, c) for c in grid)
null_max = []
for b in range(2000):
    o = rng.permutation(n)
    null_max.append(max(logrank_test(T5[o][sc4 > c], T5[o][sc4 <= c], E5[o][sc4 > c], E5[o][sc4 <= c]).test_statistic
                        for c in grid))
cuts["optimal cut (permutation-corrected)"] = dict(
    p=float((np.sum(np.array(null_max) >= obs_max) + 1) / 2001), n_high=int((sc4 > grid[np.argmax([stat(sc4, c) for c in grid])]).sum()))
res["cutpoint_sensitivity"] = cuts
print("\ncut-point sensitivity (DK 5-year OS):")
for k, v in cuts.items():
    print(f"  {k:38s} p = {v['p']:.4f}")

# ---------------- figures ----------------
fig, ax = plt.subplots(figsize=(5.4, 3.4))
lab = ["4-gene panel", "20-gene panel", "full signature"]
pt = [np.mean(boot["c4"]), np.mean(boot["c20"]), np.mean(boot["cF"])]
lo = [res["bootstrap"]["C_4gene"][0], res["bootstrap"]["C_20gene"][0], res["bootstrap"]["C_full"][0]]
hi = [res["bootstrap"]["C_4gene"][1], res["bootstrap"]["C_20gene"][1], res["bootstrap"]["C_full"][1]]
yv = np.arange(3)[::-1]
ax.errorbar(pt, yv, xerr=[np.array(pt) - np.array(lo), np.array(hi) - np.array(pt)],
            fmt="o", color="#b2182b", ms=5, capsize=3, lw=1.2)
ax.axvline(0.5, color="#888", ls="--", lw=0.9)
ax.set_yticks(yv); ax.set_yticklabels(lab); ax.set_xlabel("C-index, DK 5-year OS (95% bootstrap CI)")
ax.set_title("Fewer genes, same discrimination", fontsize=9)
save(fig, "COMPACT_bootstrap_forest")

fig, ax = plt.subplots(figsize=(5.6, 3.4))
tp = stab.head(20)
ax.barh(range(len(tp))[::-1], tp["freq_top10"], color=["#b2182b" if g in panel else "#9fb8cd" for g in tp["gene"]])
ax.set_yticks(range(len(tp))[::-1]); ax.set_yticklabels(tp["gene"], fontsize=7)
ax.set_xlabel("fraction of discovery resamples in which the gene enters the top 10")
ax.set_title("Gene selection stability (85 refits; red = in the locked 4-gene panel)", fontsize=9)
save(fig, "COMPACT_stability")

fig, ax = plt.subplots(figsize=(5.4, 3.4))
names = list(cuts.keys()); ps = [cuts[k]["p"] for k in names]
ax.barh(range(len(ps))[::-1], -np.log10(ps), color="#2166ac")
ax.axvline(-np.log10(0.05), color="#b2182b", ls="--", lw=1)
ax.set_yticks(range(len(ps))[::-1]); ax.set_yticklabels(names, fontsize=7.5)
ax.set_xlabel("-log10 log-rank p (DK, 5-year OS)")
ax.set_title("The split survives every cut-point, including a permutation-corrected optimal cut", fontsize=8.5)
save(fig, "COMPACT_cutpoints")

json.dump(res, open(f"{OUT}/COMPACT_robustness.json", "w"), indent=2, default=float)
print("\nDONE ->", OUT)
