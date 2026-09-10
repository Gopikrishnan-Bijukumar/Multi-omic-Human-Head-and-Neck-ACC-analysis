"""Diagnostic: how much prognostic signal exists in DK, and does JSE evidence enrich for it?

READ-ONLY on inputs. Writes only work/outputs.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from scipy import stats

ROOT = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"

dk = pd.read_csv(f"{ROOT}/outputs/03_harmonized/dk_norm_log2.csv", index_col=0)
tr = pd.read_csv(f"{ROOT}/outputs/03_harmonized/train_norm_log2.csv", index_col=0)
lab = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/sample_labels.csv").set_index("sample")["outcome"]
surv = pd.read_csv(f"{ROOT}/external_dk/DK_sample_data_table2.csv").set_index("UID")

s = surv.loc[dk.index]
ok = pd.to_numeric(s["SURV_MO"], errors="coerce").notna() & pd.to_numeric(s["SURV_CENS"], errors="coerce").notna()
s = s[ok]
X = dk.loc[s.index]
T = pd.to_numeric(s["SURV_MO"]).values
E = pd.to_numeric(s["SURV_CENS"]).values.astype(int)
print(f"DK usable n={len(T)} events={E.sum()}")
print("\nGROUP_Letter x status:\n", pd.crosstab(s["GROUP_Letter"], s["STATUS"]))
print("\nSTAGE:\n", s["STAGE"].value_counts(), "\nFORM:\n", s["FORM"].value_counts())

# clinical baselines
for col in ["STAGE", "FORM", "GROUP_Letter", "GENDER"]:
    v = s[col].astype(str)
    keep = ~v.isin(["unk", "dataMissing", "nan", "NA"])
    if keep.sum() < 10: continue
    g = v[keep]
    cats = g.unique()
    if len(cats) == 2:
        m = (g == cats[0]).values
        lr = logrank_test(T[keep.values][m], T[keep.values][~m], E[keep.values][m], E[keep.values][~m])
        print(f"{col}: {cats} n={keep.sum()} logrank p={lr.p_value:.4g}")

# ---- univariate Cox screen over all measurable genes (fast, vectorized score test) ----
# Use lifelines only for a subset; use the Cox score test approximation for speed.
def cox_score_test(Xm, T, E):
    """Vectorized log-rank/score test for each column (standardized). Returns z."""
    order = np.argsort(T)
    Xs = (Xm - Xm.mean(0)) / (Xm.std(0) + 1e-9)
    Xs = Xs[order]; Ts = T[order]; Es = E[order]
    n, p = Xs.shape
    # risk set = suffix sums
    csum = np.cumsum(Xs[::-1], axis=0)[::-1]
    csum2 = np.cumsum((Xs**2)[::-1], axis=0)[::-1]
    nrisk = np.arange(n, 0, -1)[:, None]
    mean_r = csum / nrisk
    var_r = csum2 / nrisk - mean_r**2
    ev = Es == 1
    U = (Xs[ev] - mean_r[ev]).sum(0)
    V = var_r[ev].sum(0)
    return U / np.sqrt(V + 1e-12)

Xm = X.values.astype(np.float64)
z = cox_score_test(Xm, T.astype(float), E)
p = 2 * stats.norm.sf(np.abs(z))
res = pd.DataFrame({"gene": X.columns, "z": z, "p": p}).sort_values("p")
from statsmodels.stats.multitest import multipletests
res["fdr"] = multipletests(res["p"], method="fdr_bh")[1]
res.to_csv(f"{OUT}/dk_univariate_cox.csv", index=False)
print(f"\nDK univariate: genes p<0.01 = {(res.p<0.01).sum()} / {len(res)} (expect {0.01*len(res):.0f}); FDR<0.10 = {(res.fdr<0.10).sum()}")
print(res.head(25).to_string(index=False))

# ---- JSE binary signal: t-stat per gene ----
y = (lab.loc[tr.index] == "Poor").values.astype(int)
A = tr.values[y == 1]; B = tr.values[y == 0]
t, pt = stats.ttest_ind(A, B, axis=0)
jse = pd.DataFrame({"gene": tr.columns, "t": t, "p": pt})
jse.to_csv(f"{OUT}/jse_univariate_t.csv", index=False)

# ---- Do JSE-significant genes enrich for DK prognostic signal? (direction-aware) ----
m = res.set_index("gene").loc[jse.gene]
agree = np.sign(m["z"].values) == np.sign(jse["t"].values)  # higher in Poor -> higher hazard
for thr in [0.05, 0.01, 0.001]:
    sel = jse.p.values < thr
    print(f"JSE p<{thr}: n={sel.sum()}  direction agreement with DK hazard = {agree[sel].mean():.3f} "
          f"(all genes {agree.mean():.3f})  binom p={stats.binomtest(int(agree[sel].sum()), int(sel.sum()), 0.5).pvalue:.3g}")
print("corr(JSE t, DK z) =", stats.spearmanr(jse["t"], m["z"]).statistic)

# candidate panel enrichment
cand = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/candidate_genes.csv")
cg = [g for g in cand.gene if g in res.gene.values]
print(f"\n436-candidate set: n present={len(cg)}; median |z| = {res.set_index('gene').loc[cg,'z'].abs().median():.3f} vs all {res.z.abs().median():.3f}")
print("frac p<0.05:", (res.set_index('gene').loc[cg,'p'] < .05).mean(), "vs all", (res.p<.05).mean())
