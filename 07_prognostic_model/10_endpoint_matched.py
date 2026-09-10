"""Endpoint-matched validation.

The discovery labels are 'dead or alive <2 years' (Poor) vs 'alive >5 years' (Good), so a
signature trained on them predicts EARLY death, not lifetime hazard. DK follow-up runs to
322 months. This script therefore evaluates the pre-specified, DK-independent signatures
against the matched endpoint - 5-year overall survival, administratively censored at 60
months - and against the discovery-identical binary contrast (dead <=24 mo vs alive >60 mo).

It also tests an unsupervised expression filter: DK libraries are shallow (~2.4M counts),
so genes measured with few reads contribute noise to every score.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, importlib.util, json
import numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sksurv.metrics import concordance_index_censored
from sklearn.metrics import roc_auc_score
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD
from mlbayes import MultiLevelModel, prep_cox, to_t

D = d3.load("all")
genes, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]

# ---- unsupervised expression filter (DK depth + training detection; no survival) ----
dk_mean = dk_raw.mean(0).values
dk_det = (dk_raw > 0).mean(0).values
tr_mean = tr_raw.mean(0).values
FILTERS = {
    "none":     np.ones(len(genes), bool),
    "moderate": (dk_mean > 1.0) & (dk_det > 0.6) & (tr_mean > 1.0),
    "strict":   (dk_mean > 2.0) & (dk_det > 0.9) & (tr_mean > 2.0),
}
for k, v in FILTERS.items(): print(f"filter {k}: {v.sum()} genes")

ev_tab = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")

ENDPOINTS = {
    "OS_full":  (T, E),
    "OS_5yr":   (np.minimum(T, 60.0), np.where(T > 60, 0, E)),
    "OS_10yr":  (np.minimum(T, 120.0), np.where(T > 120, 0, E)),
}
# discovery-identical binary contrast
early = (T <= 24) & (E == 1)
late = T > 60
bin_mask = early | late
bin_y = early[bin_mask].astype(int)
print(f"\nbinary contrast: {early.sum()} early deaths vs {late.sum()} long survivors")


def ev_score(sc, T_, E_, tag):
    sc = np.asarray(sc, float)
    if np.std(sc) < 1e-12: return None
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "s": stats.zscore(sc)}), "time", "event")
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    hrs = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event").params_.iloc[0])
    auc = roc_auc_score(bin_y, sc[bin_mask]) if len(np.unique(bin_y)) == 2 else np.nan
    return dict(score=tag, C=c, hr_per_sd=np.exp(cph.params_.iloc[0]), cox_p=cph.summary["p"].iloc[0],
                logrank_p=lr.p_value, hr_split=hrs, binary_auc=auc)


rows = []
for fname, fmask in FILTERS.items():
    gsub = [g for g, m in zip(genes, fmask) if m]
    Zd = d3.rint(dk_raw[gsub]).values
    Zt = d3.rint(tr_raw[gsub]).values
    M, mnames, msrc = MOD.build(gsub, R)
    evg = ev_tab.reindex(gsub)
    s_ev, d_ev = d3.evidence_prior(evg)
    tst = stats.ttest_ind(Zt[y == 1], Zt[y == 0], axis=0)

    SCORES = {
        "evidence_direction": Zd @ (d_ev * s_ev),
        "jse_sign_top50":  Zd @ _b(np.argsort(tst.pvalue)[:50], np.sign(tst.statistic), len(gsub)) if False else None,
    }
    def sign_top(k):
        o = np.argsort(tst.pvalue)[:k]
        b = np.zeros(len(gsub)); b[o] = np.sign(tst.statistic[o])
        return Zd @ b
    SCORES["jse_sign_top50"] = sign_top(50)
    SCORES["jse_sign_top100"] = sign_top(100)
    SCORES["jse_sign_top200"] = sign_top(200)
    for mm in ["hypoxia", "MYC_translation", "cell_cycle_G2M", "angiogenesis",
               "sc_Actively_Dividing_cells_up_poor"]:
        if mm in mnames: SCORES[mm] = Zd @ M[:, mnames.index(mm)]
    if "hypoxia" in mnames and "MYC_translation" in mnames:
        SCORES["hyp+myc"] = stats.zscore(Zd @ M[:, mnames.index("hypoxia")]) + \
                            stats.zscore(Zd @ M[:, mnames.index("MYC_translation")])

    # evidence-informed Bayesian shrinkage fitted on the JSE cohort ONLY (DK never seen)
    m = MultiLevelModel(M, np.array([sorted(set(msrc)).index(s) for s in msrc]), len(set(msrc)),
                        s_ev, d_ev, tau0=0.05, omega0=0.3, use_dk=False, use_jse=True)
    m.fit(Xj=to_t(Zt), yj=to_t(y), steps=1500)
    post = m.posterior(400)
    SCORES["bayes_JSEonly_multilevel"] = Zd @ post["beta"]
    np.save(f"{OUT}/bayes_beta_{fname}.npy", post["beta"])
    pd.DataFrame({"gene": gsub, "beta": post["beta"], "beta_sd": post["beta_sd"]}).to_csv(
        f"{OUT}/bayes_beta_{fname}.csv", index=False)

    for ename, (T_, E_) in ENDPOINTS.items():
        for sname, sc in SCORES.items():
            if sc is None: continue
            r = ev_score(sc, T_, E_, f"{sname} | filter={fname} | {ename}")
            if r: rows.append({**r, "filter": fname, "endpoint": ename, "signature": sname})

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/endpoint_matched.csv", index=False)
pd.set_option("display.width", 240)
for e in ENDPOINTS:
    print(f"\n===== endpoint {e} =====")
    print(res[res.endpoint == e].sort_values("cox_p")[
        ["signature", "filter", "C", "hr_per_sd", "cox_p", "logrank_p", "hr_split", "binary_auc"]
    ].head(12).to_string(index=False))
