"""STAGE 3 addendum: give the Brayer classifier its best reasonable shot.

POST-HOC by construction, and labelled as such throughout. The pre-specified test (median split,
5-year OS) is in 18_brayer_benchmark.py. This script asks whether the negative result there is an
artefact of evaluating their classifier in a framing they never intended:

  * they describe a MINORITY poor-prognosis group (median OS ~20 mo vs >120 mo), not a 50/50 split
  * their published survival analysis used FULL overall survival, not 5-year censored OS
  * their classifier is a binary class assignment, so a top-quantile cut is the closer analogue

Every cut-point below is scanned, and the best one is reported WITH a permutation correction for
having scanned it. Our own panel is put through the identical scan so the comparison is fair.
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
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_brayer"
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
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[genes]).values
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
order_lock = np.argsort(-np.abs(beta_lock))
idx8 = order_lock[:8]
ours = (Zd[:, idx8] * np.sign(beta_lock[idx8])).mean(1)
theirs = np.load(f"{OUT}/dk_brayer14_score.npy")

rng = np.random.default_rng(17)
QS = np.round(np.arange(0.50, 0.96, 0.05), 2)      # top 50% down to top 5%
res = {"note": "POST-HOC cut-point scan. p_perm corrects for scanning all cut-points."}


def scan_p(sc, T_, E_):
    """Log-rank p at every pre-listed cut-point. No Cox fit — safe under permutation."""
    ps = []
    for q in QS:
        hi = sc > np.quantile(sc, q)
        if hi.sum() < 4 or (~hi).sum() < 4 or E_[hi].sum() + E_[~hi].sum() == 0:
            continue
        ps.append(logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi]).p_value)
    return np.array(ps, dtype=float)


def scan(sc, T_, E_):
    out = []
    for q in QS:
        hi = sc > np.quantile(sc, q)
        if hi.sum() < 4 or (~hi).sum() < 4:
            continue
        lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
        try:
            cf = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}),
                                   "time", "event")
            hr = float(np.exp(cf.params_.iloc[0]))
        except Exception:
            hr = float("nan")
        out.append(dict(top_frac=float(1 - q), n=int(hi.sum()), hr=hr, p=float(lr.p_value)))
    return pd.DataFrame(out)


def perm_corrected(sc, T_, E_, best_p, B=2000):
    """Probability that ANY cut-point in the scan beats best_p under permuted survival."""
    cnt = 0
    for _ in range(B):
        o = rng.permutation(len(T_))
        ps = scan_p(sc, T_[o], E_[o])
        if ps.size and np.nanmin(ps) <= best_p:
            cnt += 1
    return (cnt + 1) / (B + 1)


for endpoint, (T_, E_) in [("5yr", (T5, E5)), ("full", (T, E))]:
    for who, sc in [("brayer14", theirs), ("ours_8gene", ours)]:
        s = scan(sc, T_, E_)
        s.to_csv(f"{OUT}/cutscan_{who}_{endpoint}.csv", index=False)
        b = s.loc[s.p.idxmin()]
        pp = perm_corrected(sc, T_, E_, float(b.p), B=1000)
        res[f"{who}_{endpoint}"] = dict(best_top_frac=float(b.top_frac), best_n=int(b.n),
                                        best_hr=float(b.hr), best_p_raw=float(b.p),
                                        p_permutation_corrected=float(pp),
                                        scan=s.to_dict("records"))
        print(f"{who:<12s} {endpoint:>4s}  best cut = top {b.top_frac*100:.0f}% (n={int(b.n):>2d})  "
              f"HR={b.hr:5.2f}  raw p={b.p:.4f}  perm-corrected p={pp:.4f}")

# continuous C-index on full follow-up too
for who, sc in [("brayer14", theirs), ("ours_8gene", ours)]:
    c5 = concordance_index_censored(E5.astype(bool), T5, sc)[0]
    cf_ = concordance_index_censored(E.astype(bool), T, sc)[0]
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T, "event": E, "s": stats.zscore(sc)}), "time", "event")
    res[f"{who}_continuous"] = dict(C_5yr=float(c5), C_full=float(cf_),
                                    hr_per_sd_full=float(np.exp(cph.params_.iloc[0])),
                                    cox_p_full=float(cph.summary["p"].iloc[0]))
    print(f"{who:<12s} continuous: C(5yr)={c5:.3f}  C(full)={cf_:.3f}  "
          f"HR/SD(full)={np.exp(cph.params_.iloc[0]):.2f} (p={cph.summary['p'].iloc[0]:.4f})")

json.dump(res, open(f"{OUT}/BRAYER_fairness_posthoc.json", "w"), indent=2, default=float)

# figure: cut-point scan, both signatures, both endpoints
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True)
for ax, endpoint in zip(axes, ["5yr", "full"]):
    for who, col in [("ours_8gene", RED), ("brayer14", BLUE)]:
        s = pd.DataFrame(res[f"{who}_{endpoint}"]["scan"])
        ax.plot(s.top_frac * 100, s.p, "o-", ms=4, color=col, label=who)
    ax.axhline(0.05, color="k", ls="--", lw=1)
    ax.set_xlabel("size of high-risk group (top % of cohort)")
    ax.set_yscale("log")
    ax.set_title(f"{'5-year OS' if endpoint == '5yr' else 'full follow-up'}", fontsize=9)
axes[0].set_ylabel("log-rank p (uncorrected)")
axes[0].legend(fontsize=8, frameon=False)
fig.suptitle("Post-hoc cut-point scan — giving every cut its best chance", fontsize=10, y=1.02)
save(fig, "BRAYER_cutpoint_scan")
print("\nDONE ->", OUT)
