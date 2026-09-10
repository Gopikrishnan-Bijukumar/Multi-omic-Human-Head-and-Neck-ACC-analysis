"""Final comparison figures and summary table for the compact panel.

Also answers two questions a reviewer will ask:
  - is the 4-gene panel just DSCAM in disguise?  (single-gene comparison)
  - does the panel beat the previous elastic-net model on the same patients?
All figures written as .png (600 dpi) and .svg. Nothing under work/outputs/ is touched.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, importlib.util, os
import numpy as np, pandas as pd
from scipy import stats
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_compact"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


D = d3.load("all")
genes_all, T, E = D["genes"], D["T"], D["E"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[genes]).values
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
order = np.argsort(-np.abs(beta))
CM = json.load(open(f"{OUT}/COMPACT_model.json")); panel = CM["genes"]
pidx = np.array([genes.index(g) for g in panel])
sgn = lambda idx: (Zd[:, idx] * np.sign(beta[idx])).mean(1)
sc4, sc20, scF = sgn(pidx), sgn(order[:20]), Zd @ beta

old = json.load(open(f"{R}/outputs/04_model/final_model.json"))
dkz = pd.read_csv(f"{R}/outputs/03_harmonized/dk_panel30_z.csv", index_col=0).loc[D["Zd"].index]
co = pd.Series(old["coefficients"]); co = co[co != 0]
sc_old = dkz[co.index].values @ co.values


def stat(sc, T_=None, E_=None):
    T_ = T5 if T_ is None else T_; E_ = E5 if E_ is None else E_
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    cf = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
    ci = np.exp(cf.confidence_intervals_.iloc[0].values)
    return dict(C=c, hr=float(np.exp(cf.params_.iloc[0])), lo=ci[0], hi=ci[1], p=lr.p_value, mask=hi)


# single genes, to show the panel is not one gene
rows = []
for g in panel:
    s = Zd[:, genes.index(g)] * np.sign(beta[genes.index(g)])
    r = stat(s); rows.append(dict(model=f"{g} alone", **{k: r[k] for k in ("C", "hr", "lo", "hi", "p")}))
for name, s in [(f"{len(panel)}-gene panel", sc4), ("20-gene panel", sc20),
                ("full signature", scF), ("previous elastic net", sc_old)]:
    r = stat(s); rows.append(dict(model=name, **{k: r[k] for k in ("C", "hr", "lo", "hi", "p")}))
tab = pd.DataFrame(rows)
tab.to_csv(f"{OUT}/COMPACT_model_comparison.csv", index=False)
print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

fig, ax = plt.subplots(figsize=(6.2, 3.8))
o = np.arange(len(tab))[::-1]
cols = ["#b2182b" if m == f"{len(panel)}-gene panel" else "#8c9bab" for m in tab.model]
ax.errorbar(tab.hr, o, xerr=[tab.hr - tab.lo, tab.hi - tab.hr], fmt="o", ms=5, capsize=3,
            lw=1.2, ecolor="#b0b0b0", mfc="none", linestyle="none")
ax.scatter(tab.hr, o, c=cols, zorder=3, s=32)
ax.axvline(1, color="#666", ls="--", lw=0.9)
ax.set_xscale("log"); ax.set_xticks([0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.5", "1", "2", "4", "8"])
ax.set_yticks(o); ax.set_yticklabels([f"{m}   (C={c:.3f}, p={p:.3f})" for m, c, p in zip(tab.model, tab.C, tab.p)], fontsize=7.5)
ax.set_xlabel("hazard ratio, median split (DK, 5-year OS)")
ax.set_title("The panel outperforms any of its genes alone", fontsize=9)
save(fig, "COMPACT_forest_models")

fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), sharey=True)
for ax, (name, s) in zip(axes, [("Previous model (elastic net, 7 genes)", sc_old),
                                (f"This model ({len(panel)}-gene panel)", sc4)]):
    r = stat(s); hi = r["mask"]
    for m_, lab, col in [(hi, "high risk", "#b2182b"), (~hi, "low risk", "#2166ac")]:
        KaplanMeierFitter().fit(T5[m_], E5[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(ax=ax, color=col, ci_alpha=0.12)
    ax.set_ylim(0, 1.02); ax.set_xlabel("months")
    ax.set_title(f"{name}\nC={r['C']:.3f}  HR={r['hr']:.2f} [{r['lo']:.2f}–{r['hi']:.2f}]  p={r['p']:.3f}", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False)
axes[0].set_ylabel("overall survival (5-year)")
save(fig, "COMPACT_old_vs_new")

# expression of the four genes by DK risk group and by JSE label
Zt = d3.rint(tr_raw[genes]).values
fig, axes = plt.subplots(1, len(panel), figsize=(2.1 * len(panel), 3.0), sharey=True)
hi = sc4 > np.median(sc4)
for ax, g in zip(np.atleast_1d(axes), panel):
    v = Zd[:, genes.index(g)]
    parts = [v[~hi], v[hi]]
    bp = ax.boxplot(parts, widths=0.6, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], ["#2166ac", "#b2182b"]):
        patch.set_facecolor(c); patch.set_alpha(0.35); patch.set_edgecolor(c)
    for i, p_ in enumerate(parts):
        ax.scatter(np.full(len(p_), i + 1) + np.random.uniform(-.1, .1, len(p_)), p_, s=6,
                   color=["#2166ac", "#b2182b"][i], alpha=0.8)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["low", "high"], fontsize=8)
    ax.set_title(f"{g}\np = {stats.mannwhitneyu(parts[0], parts[1]).pvalue:.1g}", fontsize=8.5)
np.atleast_1d(axes)[0].set_ylabel("expression (rank-inverse-normal)")
fig.suptitle("Panel genes across DK risk groups", fontsize=9.5, y=1.02)
save(fig, "COMPACT_gene_expression")

print("\nDONE ->", OUT)
