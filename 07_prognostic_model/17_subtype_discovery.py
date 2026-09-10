"""STAGE 2: subtype discovery in DK — a data-defined group instead of an imposed median split.

Pre-specified in work/PRESPEC_stage2.md BEFORE this script was run.

Every group assignment is computed from EXPRESSION ONLY. The script writes all assignments to
disk (step 1) before reading any DK survival (step 2). DK survival is never fitted on.

Feature sets : A = locked 8-gene panel ; B = top 100 genes by |beta|
Methods      : 1 nearest-centroid (PRIMARY)  2 GMM k=2 and k=3  3 JSE-locked Youden threshold
Primary test : Method 1 / set A / 5-year OS log-rank
Secondary    : all remaining combinations, Holm-corrected across the pre-declared family
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
from sklearn.mixture import GaussianMixture
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_subtype"
import os; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
RED, BLUE, GREY = "#b2182b", "#2166ac", "#9a9a9a"


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
dk_ids = list(dk_raw.index)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
order_lock = np.argsort(-np.abs(beta_lock))
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

FEATURES = {"A_panel8": order_lock[:8], "B_top100": order_lock[:100]}
print(f"DK n={len(T)} events(5y)={E5.sum()}   feature sets: "
      + ", ".join(f"{k}({len(v)} genes)" for k, v in FEATURES.items()))

# =====================================================================
# STEP 1 — group assignment from EXPRESSION ONLY. No survival is read here.
# =====================================================================
assign = pd.DataFrame(index=dk_ids)
meta = {}

for fname, idx in FEATURES.items():
    Xd, Xt = Zd[:, idx], Zt[:, idx]

    # ---- Method 1: nearest centroid to JSE Poor / Good ----
    c_poor, c_good = Xt[y == 1].mean(0), Xt[y == 0].mean(0)
    r_poor = np.array([stats.pearsonr(x, c_poor)[0] for x in Xd])
    r_good = np.array([stats.pearsonr(x, c_good)[0] for x in Xd])
    m1 = (r_poor > r_good).astype(int)
    assign[f"M1_centroid_{fname}"] = m1
    assign[f"M1_delta_r_{fname}"] = r_poor - r_good
    meta[f"M1_centroid_{fname}"] = dict(method="nearest centroid to JSE Poor/Good",
                                        features=fname, n_poor_like=int(m1.sum()))

    # ---- Method 2: unsupervised GMM on DK alone, k = 2 and 3 ----
    for k in (2, 3):
        gm = GaussianMixture(n_components=k, covariance_type="full", random_state=0,
                             n_init=10, reg_covar=1e-4).fit(Xd)
        lab = gm.predict(Xd)
        # anchor identity: which cluster centroid correlates most with the JSE Poor centroid
        rr = [stats.pearsonr(Xd[lab == c].mean(0), c_poor)[0] if (lab == c).sum() > 1 else -np.inf
              for c in range(k)]
        poor_c = int(np.argmax(rr))
        mk = (lab == poor_c).astype(int)
        assign[f"M2_gmm{k}_{fname}"] = mk
        assign[f"M2_gmm{k}_raw_{fname}"] = lab
        meta[f"M2_gmm{k}_{fname}"] = dict(method=f"GMM k={k} on DK, poor cluster anchored to JSE",
                                          features=fname, n_poor_like=int(mk.sum()),
                                          cluster_sizes=np.bincount(lab, minlength=k).tolist(),
                                          corr_to_jse_poor=[float(v) for v in rr])

    # ---- Method 3: JSE-locked Youden threshold on the panel score ----
    # score is the sign-weighted mean of RINT expression; with genes and signs fixed by the
    # locked model, no refitting is involved, so this cut is LOO-invariant inside JSE.
    sgn = np.sign(beta_lock[idx])
    s_jse = (Xt * sgn).mean(1)
    s_dk = (Xd * sgn).mean(1)
    cuts = np.unique(s_jse)
    J = [( (s_jse[y == 1] > c).mean() - (s_jse[y == 0] > c).mean() ) for c in cuts]
    thr = float(cuts[int(np.argmax(J))])
    m3 = (s_dk > thr).astype(int)
    assign[f"M3_youden_{fname}"] = m3
    meta[f"M3_youden_{fname}"] = dict(method="JSE-locked Youden cut-point", features=fname,
                                      threshold=thr, youden_J=float(np.max(J)),
                                      n_poor_like=int(m3.sum()))
    assign[f"score_{fname}"] = s_dk

assign.to_csv(f"{OUT}/dk_group_assignments.csv")
json.dump(meta, open(f"{OUT}/group_definitions.json", "w"), indent=2, default=float)
print("\nSTEP 1 complete — assignments written from expression alone:")
for k, v in meta.items():
    print(f"  {k:<28s} poor-like n = {v['n_poor_like']:>2d} / {len(dk_ids)}")

# =====================================================================
# STEP 2 — unblind DK survival ONCE and evaluate the pre-declared family.
# =====================================================================
GROUPS = [c for c in assign.columns if c.startswith(("M1_centroid", "M2_gmm", "M3_youden"))
          and "raw" not in c and "delta" not in c]
PRIMARY = "M1_centroid_A_panel8"
assert PRIMARY in GROUPS

rows = []
for g in GROUPS:
    hi = assign[g].values.astype(bool)
    if hi.sum() < 3 or (~hi).sum() < 3:
        rows.append(dict(group=g, n_poor=int(hi.sum()), C=np.nan, hr=np.nan,
                         lo=np.nan, up=np.nan, logrank_p=np.nan,
                         med_poor=np.nan, med_good=np.nan, note="group too small"))
        continue
    lr = logrank_test(T5[hi], T5[~hi], E5[hi], E5[~hi])
    cf = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "g": hi.astype(int)}), "time", "event")
    ci = np.exp(cf.confidence_intervals_.iloc[0].values)
    c = concordance_index_censored(E5.astype(bool), T5, hi.astype(float))[0]
    kp = KaplanMeierFitter().fit(T[hi], E[hi]); kg = KaplanMeierFitter().fit(T[~hi], E[~hi])
    rows.append(dict(group=g, n_poor=int(hi.sum()), C=float(c),
                     hr=float(np.exp(cf.params_.iloc[0])), lo=float(ci[0]), up=float(ci[1]),
                     logrank_p=float(lr.p_value),
                     med_poor=float(kp.median_survival_time_), med_good=float(kg.median_survival_time_),
                     note=""))
tab = pd.DataFrame(rows)

# Holm correction across the pre-declared secondary family (primary excluded)
sec = tab[tab.group != PRIMARY].copy().sort_values("logrank_p")
m = sec.logrank_p.notna().sum()
holm, run = [], 0.0
for i, p in enumerate(sec.logrank_p.values):
    adj = min(1.0, max(run, (m - i) * p)) if np.isfinite(p) else np.nan
    run = adj if np.isfinite(adj) else run
    holm.append(adj)
sec["holm_p"] = holm
tab = tab.merge(sec[["group", "holm_p"]], on="group", how="left")
tab.loc[tab.group == PRIMARY, "holm_p"] = np.nan
tab = tab.sort_values("logrank_p")
tab.to_csv(f"{OUT}/subtype_survival_results.csv", index=False)

print("\nSTEP 2 — DK survival, 5-year OS (primary marked *, others Holm-corrected):")
for _, r in tab.iterrows():
    star = " *PRIMARY*" if r.group == PRIMARY else ""
    hp = "" if not np.isfinite(r.holm_p) else f"  holm={r.holm_p:.3f}"
    if not np.isfinite(r.C):
        print(f"  {r.group:<28s} n={r.n_poor:>2d}  {r.note}{star}")
    else:
        print(f"  {r.group:<28s} n={r.n_poor:>2d}  HR={r.hr:5.2f} [{r.lo:.2f}-{r.up:.2f}]  "
              f"p={r.logrank_p:.4f}{hp}  medOS {r.med_poor:.0f} vs {r.med_good:.0f} mo{star}")

pr = tab[tab.group == PRIMARY].iloc[0]
res = dict(primary_test=PRIMARY,
           primary=dict(n_poor_like=int(pr.n_poor), n_total=len(dk_ids),
                        fraction=float(pr.n_poor / len(dk_ids)),
                        HR=float(pr.hr), ci=[float(pr.lo), float(pr.up)],
                        logrank_p=float(pr.logrank_p),
                        median_os_poor=float(pr.med_poor), median_os_good=float(pr.med_good),
                        significant=bool(pr.logrank_p < 0.05)),
           family_size=int(m), all_results=tab.to_dict("records"),
           median_split_comparator=None)

# comparator: the imposed 50/50 split on the same score
sc8 = assign["score_A_panel8"].values
hi = sc8 > np.median(sc8)
lr = logrank_test(T5[hi], T5[~hi], E5[hi], E5[~hi])
cf = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "g": hi.astype(int)}), "time", "event")
res["median_split_comparator"] = dict(n_poor_like=int(hi.sum()), HR=float(np.exp(cf.params_.iloc[0])),
                                      logrank_p=float(lr.p_value))
print(f"\ncomparator — imposed 50/50 median split: n={hi.sum()}  HR={np.exp(cf.params_.iloc[0]):.2f}  p={lr.p_value:.4f}")

json.dump(res, open(f"{OUT}/SUBTYPE_results.json", "w"), indent=2, default=float)

# =====================================================================
# figures
# =====================================================================
def km_group(mask, title, fname, T_=None, E_=None, ylab="overall survival (5-year)"):
    T_ = T5 if T_ is None else T_; E_ = E5 if E_ is None else E_
    lr = logrank_test(T_[mask], T_[~mask], E_[mask], E_[~mask])
    cf = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": mask.astype(int)}), "time", "event")
    hr = np.exp(cf.params_.iloc[0]); ci = np.exp(cf.confidence_intervals_.iloc[0].values)
    fig, ax = plt.subplots(figsize=(4.7, 3.9))
    for m_, lab, col in [(mask, "poor-like subtype", RED), (~mask, "rest", BLUE)]:
        KaplanMeierFitter().fit(T_[m_], E_[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(
            ax=ax, color=col, ci_alpha=0.12)
    ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nHR = {hr:.2f} [{ci[0]:.2f}–{ci[1]:.2f}], log-rank p = {lr.p_value:.4f}", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    save(fig, fname)


km_group(assign[PRIMARY].values.astype(bool),
         "DK — nearest-centroid poor-like subtype (PRIMARY)", "SUBTYPE_km_primary")
km_group(assign[PRIMARY].values.astype(bool),
         "DK — nearest-centroid poor-like subtype, full follow-up",
         "SUBTYPE_km_primary_full", T_=T, E_=E, ylab="overall survival")
km_group(hi, "DK — imposed 50/50 median split (comparator)", "SUBTYPE_km_median_comparator")

# forest of every pre-declared definition
f_ = tab[np.isfinite(tab.hr)].sort_values("hr")
fig, ax = plt.subplots(figsize=(7.0, 0.34 * len(f_) + 1.6))
yy = np.arange(len(f_))
cols = [RED if p < 0.05 else GREY for p in f_.logrank_p.values]
ax.errorbar(f_.hr.values, yy, xerr=[f_.hr.values - f_.lo.values, f_.up.values - f_.hr.values],
            fmt="o", ms=5, color="k", ecolor=cols, elinewidth=2.2, capsize=0, ls="none")
for i, (h, c) in enumerate(zip(f_.hr.values, cols)):
    ax.plot(h, i, "o", ms=5.5, color=c)
ax.axvline(1, color="k", lw=0.9, ls="--")
ax.set_yticks(yy)
ax.set_yticklabels([f"{g}  (n={n})" + ("  *PRIMARY*" if g == PRIMARY else "")
                    for g, n in zip(f_.group.values, f_.n_poor.values)], fontsize=8)
ax.set_xscale("log"); ax.set_xlabel("hazard ratio, poor-like vs rest (DK, 5-year OS)")
ax.set_title("Every pre-declared group definition\nred = nominal p < 0.05", fontsize=9)
save(fig, "SUBTYPE_forest")

# group size vs effect: the dilution question
fig, ax = plt.subplots(figsize=(5.2, 3.7))
ok = np.isfinite(f_.hr.values)
ax.scatter(f_.n_poor.values[ok] / len(dk_ids) * 100, f_.hr.values[ok],
           c=[RED if p < 0.05 else GREY for p in f_.logrank_p.values[ok]], s=45, zorder=3)
ax.axhline(res["median_split_comparator"]["HR"], color=BLUE, ls="--", lw=1.2,
           label=f"imposed 50/50 split (HR={res['median_split_comparator']['HR']:.2f})")
ax.axhline(1, color="k", lw=0.8)
ax.set_xlabel("size of the poor-like group (% of DK cohort)")
ax.set_ylabel("hazard ratio vs rest"); ax.set_yscale("log")
ax.set_title("Does a smaller, data-defined group separate better?", fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "SUBTYPE_size_vs_effect")

print("\nDONE ->", OUT)
