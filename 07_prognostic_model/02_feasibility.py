"""Feasibility scan: how well can a JSE-derived, direction-consensus score split DK?

No DK survival is used to build any score here; DK survival is only read to score
the resulting fixed signatures. Writes only work/outputs.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sksurv.metrics import concordance_index_censored

ROOT = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"
rng = np.random.default_rng(6)

dk = pd.read_csv(f"{ROOT}/outputs/03_harmonized/dk_norm_log2.csv", index_col=0)
tr = pd.read_csv(f"{ROOT}/outputs/03_harmonized/train_norm_log2.csv", index_col=0)
lab = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/sample_labels.csv").set_index("sample")["outcome"]
surv = pd.read_csv(f"{ROOT}/external_dk/DK_sample_data_table2.csv").set_index("UID")
ev = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")

s = surv.loc[dk.index]
ok = pd.to_numeric(s["SURV_MO"], errors="coerce").notna() & pd.to_numeric(s["SURV_CENS"], errors="coerce").notna()
s = s[ok]
T = pd.to_numeric(s["SURV_MO"]).values.astype(float)
E = pd.to_numeric(s["SURV_CENS"]).values.astype(int)
DK = dk.loc[s.index]
y = (lab.loc[tr.index] == "Poor").values.astype(int)

# --- robust per-cohort standardisation: rank-inverse-normal within cohort (platform-free) ---
def rint(df):
    r = df.rank(axis=0)
    n = df.shape[0]
    return pd.DataFrame(stats.norm.ppf((r.values - 0.375) / (n + 0.25)), index=df.index, columns=df.columns)

TRz, DKz = rint(tr), rint(DK)

# JSE per-gene statistics
t, pt = stats.ttest_ind(TRz.values[y == 1], TRz.values[y == 0], axis=0)
jse = pd.DataFrame({"t": t, "p": pt}, index=tr.columns)
# also raw-scale effect (auc-like) for robustness
auc_j = np.array([stats.mannwhitneyu(TRz.values[y == 1, i], TRz.values[y == 0, i]).statistic for i in range(TRz.shape[1])]) / 100.0

def evaluate(score, name, extra=""):
    score = np.asarray(score, float)
    c = concordance_index_censored(E.astype(bool), T, score)[0]
    df = pd.DataFrame({"time": T, "event": E, "score": stats.zscore(score)})
    cph = CoxPHFitter().fit(df, "time", "event")
    hi = score > np.median(score)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    hr_split = np.nan
    d2 = df.copy(); d2["grp"] = hi.astype(int)
    try:
        hr_split = np.exp(CoxPHFitter().fit(d2[["time", "event", "grp"]], "time", "event").params_.iloc[0])
    except Exception: pass
    print(f"{name:52s} C={c:.3f} coxp={cph.summary['p'].iloc[0]:.3g} HRperSD={np.exp(cph.params_.iloc[0]):.2f} "
          f"logrank={lr.p_value:.3g} HRsplit={hr_split:.2f} {extra}")
    return dict(name=name, C=c, cox_p=cph.summary["p"].iloc[0], hr_sd=np.exp(cph.params_.iloc[0]),
                logrank_p=lr.p_value, hr_split=hr_split)

rows = []
# 1. sign-consensus score over JSE-significant genes, various thresholds / sizes
for thr in [0.05, 0.01, 0.005, 0.001]:
    g = jse.index[(jse.p < thr)]
    g = [x for x in g if x in DKz.columns]
    w = np.sign(jse.loc[g, "t"].values)
    rows.append(evaluate(DKz[g].values @ w / len(g), f"sign-consensus JSE p<{thr} (n={len(g)})"))

for k in [25, 50, 100, 200, 400, 800]:
    g = jse.reindex(DKz.columns).dropna().sort_values("p").index[:k].tolist()
    w = np.sign(jse.loc[g, "t"].values)
    rows.append(evaluate(DKz[g].values @ w / k, f"sign-consensus top-{k} JSE"))
    wt = jse.loc[g, "t"].values
    rows.append(evaluate(DKz[g].values @ wt / np.abs(wt).sum(), f"t-weighted top-{k} JSE"))

# 2. restricted to multi-omic candidate genes (evidence-filtered)
cand = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/candidate_genes.csv")
for pool_name, pool in [("436 candidates", cand.gene.tolist()),
                        ("cross-modal tiers1-3", cand.loc[cand.cross_modal == True, "gene"].tolist()),
                        ("evidence 1803", ev.index.tolist())]:
    g = [x for x in pool if x in DKz.columns and x in jse.index]
    d = jse.loc[g].sort_values("p")
    for k in [30, 60, 120, len(g)]:
        k = min(k, len(g)); gg = d.index[:k].tolist()
        w = np.sign(jse.loc[gg, "t"].values)
        rows.append(evaluate(DKz[gg].values @ w / k, f"{pool_name}: sign top-{k}"))

# 3. direction from evidence consensus only (no JSE t-test at all)
g = [x for x in ev.index if x in DKz.columns and ev.loc[x, "direction_consensus"] in ("Poor", "Good")]
w = np.array([1.0 if ev.loc[x, "direction_consensus"] == "Poor" else -1.0 for x in g])
rows.append(evaluate(DKz[g].values @ w / len(g), f"evidence-direction only (n={len(g)})"))

pd.DataFrame(rows).to_csv(f"{OUT}/feasibility_scan.csv", index=False)
