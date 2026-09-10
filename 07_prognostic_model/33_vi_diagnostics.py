"""Stage 33 - inference diagnostics for the variational fit (review item #1).

The first methods review asked for optimization traces, multiple initializations,
coefficient/rank stability across fits, sensitivity to tau0/omega0/module structure,
and a justification for variational inference over MCMC or a penalized model.

Nothing in stages 00-32 persisted any of it: mlbayes.py printed the ELBO only when
verbose=True and no caller ever passed it.

One structural fact drives this whole stage and was found while building it:

    mlbayes.MultiLevelModel(seed=...) seeds ONLY the variational parameter
    initialisation (mlbayes.py:47, a CPU torch.Generator). The ADVI gradient noise
    at mlbayes.py:74 is drawn from the GLOBAL torch RNG, which no production stage
    ever seeds. torch.manual_seed appears only in 31d/31e.

So every locked artefact was produced under unseeded stochastic gradient noise.
Part A measures that directly; Parts B-D measure what it costs.

Writes to work/outputs_vi/. Reads only; overwrites nothing.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, time, importlib.util, itertools
import numpy as np, pandas as pd
import torch
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sksurv.metrics import concordance_index_censored
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_vi"
import os; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD
from mlbayes import MultiLevelModel, to_t

PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]
STEPS = 900
LOCKED = dict(tau0=0.02, omega0=0.1, gene_level=True, module_level=True)

# ---------------------------------------------------------------- data (mirrors stage 11)
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]
keep = (dk_raw.mean(axis=0).values > 1.0) & ((dk_raw > 0).mean(axis=0).values > 0.6) & (tr_raw.mean(axis=0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
gidx = {g: i for i, g in enumerate(genes)}
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
M, mnames, msrc = MOD.build(genes, R)
srcs = sorted(set(msrc)); src_ids = np.array([srcs.index(s) for s in msrc])
ev_tab = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")
s_ev, d_ev = d3.evidence_prior(ev_tab.reindex(genes))
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
print(f"genes={len(genes)} modules={M.shape[1]} JSE n={len(y)} DK n={len(T)} events5y={E5.sum()}", flush=True)


def fit_beta(idx, tau0, omega0, gene_level, module_level, seed=0, steps=STEPS,
             torch_seed=None, trace=False):
    """Stage-11's fit_beta plus explicit control of the global RNG that drives ADVI noise."""
    if torch_seed is not None:
        torch.manual_seed(torch_seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(torch_seed)
    m = MultiLevelModel(M, src_ids, len(srcs), s_ev, d_ev, tau0=tau0, omega0=omega0,
                        use_gene_level=gene_level, use_module_level=module_level,
                        use_dk=False, use_jse=True, seed=seed)
    m.fit(Xj=to_t(Zt[idx]), yj=to_t(y[idx]), steps=steps, trace=trace)
    post = m.posterior(300)
    post["elbo_trace"] = m.elbo_trace
    return post


def top_panel(beta, k=4):
    return [genes[i] for i in np.argsort(-np.abs(beta))[:k]]


def dk_c(panel, beta):
    """DK 5-year C for the locked scoring rule: mean of sign-weighted RINT expression."""
    ii = [gidx[g] for g in panel]
    sc = (Zd[:, ii] * np.sign(beta[ii])).mean(1)
    return float(concordance_index_censored(E5.astype(bool), T5, sc)[0])


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b)


res = {"_note": "Stage 33 inference diagnostics. Reads only; no locked artefact modified.",
       "n_genes": len(genes), "n_modules": int(M.shape[1]), "steps": STEPS, "locked_config": LOCKED}
t00 = time.time()

# =============================================================== A. determinism probe
print("\n[A] determinism probe", flush=True)
bA1 = fit_beta(np.arange(len(y)), **LOCKED, torch_seed=123)["beta"]
bA2 = fit_beta(np.arange(len(y)), **LOCKED, torch_seed=123)["beta"]
bB1 = fit_beta(np.arange(len(y)), **LOCKED, torch_seed=None)["beta"]
bB2 = fit_beta(np.arange(len(y)), **LOCKED, torch_seed=None)["beta"]
res["determinism"] = {
    "seeded_global_rng": {
        "max_abs_diff": float(np.max(np.abs(bA1 - bA2))),
        "pearson": float(np.corrcoef(bA1, bA2)[0, 1]),
        "spearman_abs": float(stats.spearmanr(np.abs(bA1), np.abs(bA2)).statistic),
        "top4_a": top_panel(bA1), "top4_b": top_panel(bA2)},
    "unseeded_global_rng": {
        "max_abs_diff": float(np.max(np.abs(bB1 - bB2))),
        "pearson": float(np.corrcoef(bB1, bB2)[0, 1]),
        "spearman_abs": float(stats.spearmanr(np.abs(bB1), np.abs(bB2)).statistic),
        "top4_a": top_panel(bB1), "top4_b": top_panel(bB2)},
    "interpretation": ("mlbayes seed= controls initialisation only; ADVI gradient noise comes from "
                       "the global torch RNG, unseeded in every production stage (11/13/16/25/25b/25c).")}
print("  seeded   max|d|", round(res["determinism"]["seeded_global_rng"]["max_abs_diff"], 6),
      "| unseeded max|d|", round(res["determinism"]["unseeded_global_rng"]["max_abs_diff"], 6), flush=True)

# =============================================================== B. multi-seed ensemble
NSEED = 20
print(f"\n[B] {NSEED}-seed ensemble at the locked configuration", flush=True)
betas, traces = [], []
for i in range(NSEED):
    p = fit_beta(np.arange(len(y)), **LOCKED, seed=i, torch_seed=1000 + i, trace=True)
    betas.append(p["beta"]); traces.append(p["elbo_trace"])
    if (i + 1) % 5 == 0:
        print(f"  {i+1}/{NSEED} fits ({time.time()-t00:.0f}s)", flush=True)
B = np.vstack(betas)
np.savez_compressed(f"{OUT}/multiseed_beta.npz", beta=B, genes=np.array(genes),
                    elbo=np.array(traces, dtype=np.float32))

pear = np.corrcoef(B)
sp = np.array([[stats.spearmanr(np.abs(B[i]), np.abs(B[j])).statistic for j in range(NSEED)] for i in range(NSEED)])
off = ~np.eye(NSEED, dtype=bool)
panels = {k: [top_panel(b, k) for b in B] for k in (4, 8, 20)}
jac = {k: [jaccard(panels[k][i], panels[k][j]) for i in range(NSEED) for j in range(i + 1, NSEED)] for k in (4, 8, 20)}
freq4 = pd.Series([g for p in panels[4] for g in p]).value_counts()
freq20 = pd.Series([g for p in panels[20] for g in p]).value_counts()
c_seed = [dk_c(panels[4][i], B[i]) for i in range(NSEED)]
c_locked_panel = [dk_c(PANEL4, B[i]) for i in range(NSEED)]

res["multiseed"] = {
    "n_seeds": NSEED,
    "beta_pearson_offdiag": {"mean": float(pear[off].mean()), "min": float(pear[off].min()), "max": float(pear[off].max())},
    "absbeta_spearman_offdiag": {"mean": float(sp[off].mean()), "min": float(sp[off].min()), "max": float(sp[off].max())},
    "topk_jaccard": {f"top{k}": {"mean": float(np.mean(v)), "min": float(np.min(v)), "max": float(np.max(v))}
                     for k, v in jac.items()},
    "locked_panel_recovery_top4": {g: int(freq4.get(g, 0)) for g in PANEL4},
    "locked_panel_recovery_top20": {g: int(freq20.get(g, 0)) for g in PANEL4},
    "most_frequent_top4_genes": freq4.head(15).to_dict(),
    "most_frequent_top20_genes": freq20.head(25).to_dict(),
    "dk_c_seed_own_panel": {"mean": float(np.mean(c_seed)), "sd": float(np.std(c_seed, ddof=1)),
                            "min": float(np.min(c_seed)), "max": float(np.max(c_seed)), "all": c_seed},
    "dk_c_locked_panel_seed_signs": {"mean": float(np.mean(c_locked_panel)),
                                     "sd": float(np.std(c_locked_panel, ddof=1))},
    "dk_c_locked_reference": 0.67149}
pd.DataFrame({"seed": range(NSEED), "top4": [";".join(p) for p in panels[4]],
              "dk_c_5yr": c_seed}).to_csv(f"{OUT}/multiseed_panels.csv", index=False)
print("  mean |beta| rank spearman across seeds:", round(res["multiseed"]["absbeta_spearman_offdiag"]["mean"], 3), flush=True)
print("  top4 Jaccard mean:", round(res["multiseed"]["topk_jaccard"]["top4"]["mean"], 3), flush=True)
print("  DK C of each seed's own top4: %.3f +- %.3f" % (np.mean(c_seed), np.std(c_seed, ddof=1)), flush=True)

# =============================================================== C. convergence
print("\n[C] convergence: 900 vs 3000 steps", flush=True)
p900 = fit_beta(np.arange(len(y)), **LOCKED, seed=0, torch_seed=77, trace=True)
p3000 = fit_beta(np.arange(len(y)), **LOCKED, seed=0, torch_seed=77, steps=3000, trace=True)
e9, e3 = np.array(p900["elbo_trace"]), np.array(p3000["elbo_trace"])
res["convergence"] = {
    "elbo_900": {"first": float(e9[0]), "last": float(e9[-1]),
                 "mean_last100": float(e9[-100:].mean()), "sd_last100": float(e9[-100:].std())},
    "elbo_3000": {"last": float(e3[-1]), "mean_last100": float(e3[-100:].mean()),
                  "sd_last100": float(e3[-100:].std()),
                  "mean_steps_800_900": float(e3[800:900].mean())},
    "beta_900_vs_3000": {"pearson": float(np.corrcoef(p900["beta"], p3000["beta"])[0, 1]),
                         "spearman_abs": float(stats.spearmanr(np.abs(p900["beta"]), np.abs(p3000["beta"])).statistic),
                         "top4_900": top_panel(p900["beta"]), "top4_3000": top_panel(p3000["beta"]),
                         "dk_c_900": dk_c(top_panel(p900["beta"]), p900["beta"]),
                         "dk_c_3000": dk_c(top_panel(p3000["beta"]), p3000["beta"])}}
np.savez_compressed(f"{OUT}/convergence_traces.npz", elbo900=e9, elbo3000=e3)
print("  ELBO last100: 900-step %.1f | 3000-step %.1f" % (e9[-100:].mean(), e3[-100:].mean()), flush=True)

# =============================================================== D. hyperparameter sensitivity
print("\n[D] hyperparameter sensitivity, 12 configs x 5 seeds", flush=True)
GRID = [dict(tau0=t, omega0=o, gene_level=g, module_level=mo)
        for t in [0.02, 0.1] for o in [0.1, 0.5] for g, mo in [(True, True), (False, True), (True, False)]]
rows = []
for ci, cfg in enumerate(GRID):
    for s in range(5):
        p = fit_beta(np.arange(len(y)), **cfg, seed=s, torch_seed=2000 + 100 * ci + s)
        b = p["beta"]; pan = top_panel(b)
        loo_like = roc_auc_score(y, Zt @ b)   # in-sample separation, not a performance estimate
        rows.append({**cfg, "seed": s, "top4": ";".join(pan), "dk_c_5yr": dk_c(pan, b),
                     "dk_c_locked_panel": dk_c(PANEL4, b), "insample_auc": float(loo_like),
                     "n_locked_in_top4": len(set(pan) & set(PANEL4))})
    print(f"  config {ci+1}/12 done ({time.time()-t00:.0f}s)", flush=True)
sens = pd.DataFrame(rows)
sens.to_csv(f"{OUT}/hyperparameter_sensitivity.csv", index=False)
res["sensitivity"] = {
    "note": ("All 12 configs tied at LOO AUC 0.96 in work/outputs/jse_loo_tuning.csv, so the locked "
             "config was an argmax over a 12-way tie. Sensitivity is therefore assessed on the "
             "selected panel and its external performance, not on the tuning metric."),
    "dk_c_own_panel": {"mean": float(sens.dk_c_5yr.mean()), "sd": float(sens.dk_c_5yr.std()),
                       "min": float(sens.dk_c_5yr.min()), "max": float(sens.dk_c_5yr.max())},
    "dk_c_locked_panel": {"mean": float(sens.dk_c_locked_panel.mean()), "sd": float(sens.dk_c_locked_panel.std())},
    "n_distinct_panels": int(sens.top4.nunique()),
    "locked_genes_in_top4_mean": float(sens.n_locked_in_top4.mean()),
    "by_config": sens.groupby(["tau0", "omega0", "gene_level", "module_level"]).agg(
        dk_c_mean=("dk_c_5yr", "mean"), dk_c_sd=("dk_c_5yr", "std"),
        locked_hits=("n_locked_in_top4", "mean")).reset_index().to_dict("records")}
print("  DK C across all 60 fits: %.3f +- %.3f" % (sens.dk_c_5yr.mean(), sens.dk_c_5yr.std()), flush=True)

# =============================================================== E. penalized comparator
print("\n[E] penalized-regression comparator (LOO inside discovery)", flush=True)
def loo_penalized(kind):
    pred = np.zeros(len(y))
    for i in range(len(y)):
        tr = np.array([j for j in range(len(y)) if j != i])
        if kind == "l1":
            clf = LogisticRegression(penalty="l1", C=0.1, solver="liblinear", max_iter=5000)
        elif kind == "l2":
            clf = LogisticRegression(penalty="l2", C=0.1, solver="lbfgs", max_iter=5000)
        else:
            clf = LogisticRegression(penalty="elasticnet", C=0.1, l1_ratio=0.5,
                                     solver="saga", max_iter=5000)
        clf.fit(Zt[tr], y[tr])
        pred[i] = clf.decision_function(Zt[i:i + 1])[0]
    return float(roc_auc_score(y, pred))

comp = {}
for kind in ("l1", "l2", "enet"):
    auc = loo_penalized(kind)
    if kind == "l1":
        clf = LogisticRegression(penalty="l1", C=0.1, solver="liblinear", max_iter=5000)
    elif kind == "l2":
        clf = LogisticRegression(penalty="l2", C=0.1, solver="lbfgs", max_iter=5000)
    else:
        clf = LogisticRegression(penalty="elasticnet", C=0.1, l1_ratio=0.5, solver="saga", max_iter=5000)
    clf.fit(Zt, y)
    co = clf.coef_.ravel()
    pan = top_panel(co)
    comp[kind] = {"loo_auc": auc, "top4": pan, "dk_c_5yr": dk_c(pan, co),
                  "n_nonzero": int((co != 0).sum())}
    print(f"  {kind:5s} LOO AUC {auc:.3f}  DK C {comp[kind]['dk_c_5yr']:.3f}  top4 {pan}", flush=True)
res["comparator"] = {
    "penalized": comp,
    "bayesian_loo_auc_reported": 0.96,
    "why_not_mcmc": ("The model has 12,700 gene-level and 52 module-level latent parameters plus "
                     "hierarchical scales, i.e. >25,000 latents fitted to 20 binary outcomes. "
                     "Full MCMC was not attempted; the 3000-step run in part C is used as the "
                     "convergence reference instead. This is a stated limitation, not a claim of "
                     "posterior accuracy.")}

res["runtime_s"] = round(time.time() - t00, 1)
json.dump(res, open(f"{OUT}/VI_DIAGNOSTICS.json", "w"), indent=2)
print(f"\nwrote {OUT}/VI_DIAGNOSTICS.json  ({res['runtime_s']}s)", flush=True)
