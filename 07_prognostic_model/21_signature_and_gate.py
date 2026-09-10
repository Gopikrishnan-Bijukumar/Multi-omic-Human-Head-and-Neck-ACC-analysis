"""STAGE 4b: build the signature matrix and run the SANITY GATE.

The gate, from work/deconvolv_plan.md, written before this was run:

    "Deconvolve JSE too. If the estimated proportions do not reproduce the scCODA-observed
     Poor vs Good shift in the cohort where we know the single-cell truth, the deconvolution
     is not trustworthy on DK and the analysis stops there."

We can do better than the plan asked, because 19 of the 20 JSE bulk patients have matched
single-cell data. So the gate is not just "does the group shift replicate" but the far stricter
"does the estimated composition match that same patient's TRUE composition".

Leave-one-patient-out: the signature matrix used to deconvolve patient i is built from the other
18 patients only. No patient contributes to their own reference. DK is not touched in this script.

PRE-SPECIFIED GATE CRITERIA (all three must hold to proceed to DK):
  G1  median per-cell-type Spearman r (estimated vs true, across patients) >= 0.40
  G2  each of the three dominant compartments (Myoepithelial_Tumor, Epithelial_Tumor, Dividing)
      individually >= 0.40
  G3  mean per-patient composition correlation >= 0.70
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import os, json, itertools
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import nnls
from sklearn.svm import NuSVR
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_deconv"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
RED, BLUE, GREY, GREEN = "#b2182b", "#2166ac", "#9a9a9a", "#1b7837"


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


META = json.load(open(f"{OUT}/reference_meta.json"))
COARSE = META["coarse_map"]
pb = pd.read_parquet(f"{OUT}/pseudobulk_ct_by_patient.parquet")
ncell = pd.read_parquet(f"{OUT}/ncells_ct_by_patient.parquet")["n_cells"]
true_c = pd.read_csv(f"{OUT}/sc_true_composition_coarse.csv", index_col=0)

# bulk (log2 -> linear); read-only
tr_log = pd.read_csv(f"{R}/outputs/03_harmonized/train_norm_log2.csv", index_col=0)
dk_log = pd.read_csv(f"{R}/outputs/03_harmonized/dk_norm_log2.csv", index_col=0)
shared = sorted(set(pb.columns) & set(tr_log.columns) & set(dk_log.columns))
print(f"genes shared by single-cell, JSE bulk and DK bulk: {len(shared)}")

# collapse pseudobulk to coarse cell types, keep the patient axis
pb = pb[shared]
pb.index = pd.MultiIndex.from_arrays(
    [[COARSE.get(c, "Other") for c in pb.index.get_level_values("cell_type")],
     pb.index.get_level_values("patient")], names=["cell_type", "patient"])
pb = pb.groupby(level=["cell_type", "patient"]).sum()
CT = sorted(pb.index.get_level_values("cell_type").unique())
PT_SC = sorted(pb.index.get_level_values("patient").unique())
print(f"coarse cell types ({len(CT)}): {CT}")

shared_pts = sorted(set(PT_SC) & set(tr_log.index))
print(f"patients with BOTH single-cell and bulk: {len(shared_pts)}")


def cpm(df, axis=1):
    return df.div(df.sum(axis=axis), axis=0) * 1e6


def build_signature(exclude_patient=None, n_markers=60):
    """Signature matrix (genes x cell types) in CPM, plus the marker gene list."""
    sub = pb
    if exclude_patient is not None:
        sub = pb[pb.index.get_level_values("patient") != exclude_patient]
    prof = sub.groupby(level="cell_type").sum()          # cell type x gene raw counts
    prof = cpm(prof)                                      # CPM within cell type
    # marker selection: specificity = log2( own / max(other) )
    markers = []
    lv = np.log2(prof + 1.0)
    for ct in prof.index:
        others = lv.drop(index=ct).max(axis=0)
        spec = lv.loc[ct] - others
        ok = (prof.loc[ct] > 20)                          # must be expressed in that cell type
        cand = spec[ok].sort_values(ascending=False)
        markers += list(cand.index[:n_markers])
    markers = sorted(set(markers))
    return prof.T.loc[markers], markers                   # genes x cell types


def deconvolve(bulk_vec, S, method="nusvr"):
    """bulk_vec: genes (linear, marker-subset). S: genes x cell types. Returns fractions."""
    b = bulk_vec.values.astype(float)
    A = S.values.astype(float)
    # scale both to unit norm per CIBERSORT convention
    bz = (b - b.mean()) / (b.std() + 1e-12)
    Az = (A - A.mean(0)) / (A.std(0) + 1e-12)
    if method == "nusvr":
        best, best_r = None, -np.inf
        for nu in (0.25, 0.5, 0.75):
            m = NuSVR(nu=nu, C=1.0, kernel="linear").fit(Az, bz)
            w = m.coef_.ravel()
            pred = Az @ w
            r = np.corrcoef(pred, bz)[0, 1]
            if np.isfinite(r) and r > best_r:
                best_r, best = r, w
        w = best
    elif method == "nnls":
        w = nnls(A, b)[0]
    elif method == "dwls":
        w = nnls(A, b)[0]
        for _ in range(12):
            fit = A @ w
            wt = 1.0 / np.maximum(fit, np.percentile(fit[fit > 0], 5)) ** 2
            sw = np.sqrt(wt)
            w_new = nnls(A * sw[:, None], b * sw)[0]
            if np.max(np.abs(w_new - w)) < 1e-8:
                w = w_new; break
            w = w_new
    w = np.clip(w, 0, None)
    return pd.Series(w / w.sum() if w.sum() > 0 else np.full(len(w), 1 / len(w)), index=S.columns)


# ---------------------------------------------------------------
# SANITY GATE: leave-one-patient-out deconvolution of JSE bulk
# ---------------------------------------------------------------
tr_lin = (2.0 ** tr_log[shared]) - 1.0
tr_lin[tr_lin < 0] = 0.0

results = {}
for method in ["nusvr", "dwls", "nnls"]:
    est = {}
    for p in shared_pts:
        S, mk = build_signature(exclude_patient=p)
        est[p] = deconvolve(tr_lin.loc[p, S.index], S, method=method)
    est = pd.DataFrame(est).T
    est = est.reindex(columns=true_c.columns.intersection(est.columns))
    tr_true = true_c.loc[est.index, est.columns]

    per_ct = {}
    for ct in est.columns:
        if tr_true[ct].std() < 1e-9:
            continue
        per_ct[ct] = dict(
            spearman=float(stats.spearmanr(est[ct], tr_true[ct]).statistic),
            pearson=float(stats.pearsonr(est[ct], tr_true[ct])[0]),
            mean_est=float(est[ct].mean()), mean_true=float(tr_true[ct].mean()))
    per_pt = [float(stats.pearsonr(est.loc[p], tr_true.loc[p])[0]) for p in est.index]
    rmse = float(np.sqrt(((est.values - tr_true.values) ** 2).mean()))

    med = float(np.median([v["spearman"] for v in per_ct.values()]))
    dom = {c: per_ct[c]["spearman"] for c in ["Myoepithelial_Tumor", "Epithelial_Tumor", "Dividing"]
           if c in per_ct}
    results[method] = dict(per_celltype=per_ct, median_spearman=med,
                           dominant=dom, mean_patient_r=float(np.mean(per_pt)), rmse=rmse,
                           G1=bool(med >= 0.40), G2=bool(all(v >= 0.40 for v in dom.values())),
                           G3=bool(np.mean(per_pt) >= 0.70))
    est.to_csv(f"{OUT}/jse_estimated_composition_{method}.csv")

    print(f"\n=== {method.upper()} — leave-one-patient-out on {len(est)} JSE patients ===")
    for ct, v in sorted(per_ct.items(), key=lambda kv: -kv[1]["spearman"]):
        print(f"  {ct:<22s} spearman={v['spearman']:+.3f}  pearson={v['pearson']:+.3f}  "
              f"est={v['mean_est']:.3f} true={v['mean_true']:.3f}")
    print(f"  median spearman = {med:.3f}   mean per-patient r = {np.mean(per_pt):.3f}   RMSE = {rmse:.4f}")
    print(f"  G1 {'PASS' if results[method]['G1'] else 'FAIL'}   "
          f"G2 {'PASS' if results[method]['G2'] else 'FAIL'}   "
          f"G3 {'PASS' if results[method]['G3'] else 'FAIL'}")

# choose the method that passes the most gates, tie-broken by median spearman
rank = sorted(results.items(),
              key=lambda kv: (kv[1]["G1"] + kv[1]["G2"] + kv[1]["G3"], kv[1]["median_spearman"]),
              reverse=True)
BEST = rank[0][0]
GATE_PASSED = all([results[BEST]["G1"], results[BEST]["G2"], results[BEST]["G3"]])
print(f"\n*** best method: {BEST}   GATE {'PASSED' if GATE_PASSED else 'FAILED'} ***")

# ---------------------------------------------------------------
# Poor vs Good compositional shift, from TRUE single-cell composition
# (this is the direction the DK analysis will be pre-specified against)
# ---------------------------------------------------------------
outc = pd.Series(META["patient_outcome"])
tc = true_c.loc[[p for p in true_c.index if p in outc.index]]
lab = outc.loc[tc.index]
shift = []
for ct in tc.columns:
    a, b = tc.loc[lab == "Poor", ct], tc.loc[lab == "Good", ct]
    shift.append(dict(cell_type=ct, mean_poor=float(a.mean()), mean_good=float(b.mean()),
                      log2fc=float(np.log2((a.mean() + 1e-4) / (b.mean() + 1e-4))),
                      mw_p=float(stats.mannwhitneyu(a, b).pvalue)))
shift = pd.DataFrame(shift).sort_values("log2fc", ascending=False)
shift.to_csv(f"{OUT}/sc_true_poor_vs_good_shift.csv", index=False)
print("\nTRUE single-cell composition shift, Poor vs Good (24 sc patients):")
print(shift.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

json.dump({"gate": results, "best_method": BEST, "gate_passed": bool(GATE_PASSED),
           "n_patients_gate": len(shared_pts), "n_shared_genes": len(shared),
           "poor_vs_good_shift": shift.to_dict("records")},
          open(f"{OUT}/GATE_results.json", "w"), indent=2, default=float)

# ---------------------------------------------------------------
# figures
# ---------------------------------------------------------------
est = pd.read_csv(f"{OUT}/jse_estimated_composition_{BEST}.csv", index_col=0)
tr_true = true_c.loc[est.index, est.columns]
ncols = 3; nrows = int(np.ceil(len(est.columns) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.7 * nrows))
for ax, ct in zip(axes.ravel(), est.columns):
    r = stats.spearmanr(est[ct], tr_true[ct]).statistic
    ax.scatter(tr_true[ct], est[ct], s=26, color=RED if r >= 0.4 else GREY, zorder=3)
    lim = max(tr_true[ct].max(), est[ct].max()) * 1.1 + 1e-3
    ax.plot([0, lim], [0, lim], "k:", lw=.8)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_title(f"{ct}\nρ = {r:+.2f}", fontsize=8)
    ax.set_xlabel("true (single-cell)", fontsize=7.5); ax.set_ylabel("estimated (bulk)", fontsize=7.5)
for ax in axes.ravel()[len(est.columns):]:
    ax.axis("off")
fig.suptitle(f"Sanity gate — {BEST.upper()}, leave-one-patient-out on {len(est)} JSE patients",
             fontsize=10, y=1.005)
fig.tight_layout()
save(fig, "GATE_estimated_vs_true")

fig, ax = plt.subplots(figsize=(6.2, 3.6))
w = 0.38; xx = np.arange(len(est.columns))
ax.bar(xx - w / 2, tr_true.mean(0).values, w, color=BLUE, label="true (single-cell)")
ax.bar(xx + w / 2, est.mean(0).values, w, color=RED, label="estimated (bulk deconvolution)")
ax.set_xticks(xx); ax.set_xticklabels(est.columns, rotation=40, ha="right", fontsize=8)
ax.set_ylabel("mean fraction across patients")
ax.set_title("Does bulk deconvolution recover the right average composition?", fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "GATE_mean_composition")

print(f"\nDONE -> {OUT}")
