"""Stage 34 - is the biology stable even though the gene membership is not? (review item #2)

The review's second objection: under leakage-free refitting no gene was selected in more
than 5 of 25 folds (DSCAM 5, ODC1 5, NCAPG 3, CCNB2 0), so "the evidence supports a
proliferation/MYC-related expression axis more strongly than these exact four transcripts".

That is a claim about programmes, and nothing in stages 00-32 tested it. This stage does,
using the per-fold gene lists already stored in outputs_nested/fold_selections.csv (the
nested-k leakage-free arm), and asks three questions:

  1. Do the fold panels concentrate on the same curated programmes, against a
     size-matched random-panel null?
  2. Do the fold panels produce DK scores that agree with the locked 4-gene score?
  3. Do they reach comparable external discrimination?

If (1)-(3) hold, membership instability is cosmetic and the axis is the real finding,
which is what the manuscript should claim. Signs are taken from the locked discovery
beta (work/outputs/FINAL_model.json) - DK survival is used nowhere in constructing any
score here.

Writes to work/outputs_axis/. Reads only; overwrites nothing.
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
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_axis"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD

PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]
PROLIF = ["cell_cycle_G2M", "cell_cycle_S", "E2F_targets", "MYC_translation", "DNA_repair"]
rng = np.random.default_rng(34)

# ---------------------------------------------------------------- data (mirrors stage 11)
D = d3.load("all")
genes_all, T, E = D["genes"], D["T"], D["E"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]
keep = (dk_raw.mean(axis=0).values > 1.0) & ((dk_raw > 0).mean(axis=0).values > 0.6) & (tr_raw.mean(axis=0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
gidx = {g: i for i, g in enumerate(genes)}
Zd = d3.rint(dk_raw[genes]).values
M, mnames, msrc = MOD.build(genes, R)
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

beta_locked = np.asarray(json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))["beta"], float)
assert len(beta_locked) == len(genes), "gene universe does not match the locked beta"
sgn = np.sign(beta_locked)

locked_score = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel4_score.npy")
print(f"genes={len(genes)} modules={len(mnames)} DK n={len(T)} events5y={E5.sum()}", flush=True)


def dk_score(panel):
    ii = [gidx[g] for g in panel if g in gidx]
    return (Zd[:, ii] * sgn[ii]).mean(1)


def dk_c(sc):
    return float(concordance_index_censored(E5.astype(bool), T5, sc)[0])


# sanity: our reconstruction of the locked score must match the stored one
ours = dk_score(PANEL4)
assert np.corrcoef(ours, locked_score)[0, 1] > 0.999, "locked score reconstruction failed"
C_LOCKED = dk_c(locked_score)
print(f"locked 4-gene DK 5-yr C = {C_LOCKED:.5f} (stored score reproduced)", flush=True)

# ---------------------------------------------------------------- fold panels
fs = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_nested/fold_selections.csv")
lf = fs[fs.arm == "leakage_free"].reset_index(drop=True)
lf["panel"] = lf.genes.str.split(";")
print(f"leakage-free folds: {len(lf)}  panel sizes {lf.panel.map(len).min()}-{lf.panel.map(len).max()}", flush=True)

res = {"_note": "Stage 34 axis-level stability. Reads only; no locked artefact modified.",
       "n_folds": int(len(lf)), "locked_dk_c_5yr": C_LOCKED,
       "sign_source": "work/outputs/FINAL_model.json (discovery-only fit); DK survival not used"}

# =============================================================== 1. programme concentration
mod_of = {}
for k, nm in enumerate(mnames):
    if msrc[k] != "curated":
        continue
    mod_of[nm] = set(np.nonzero(M[:, k])[0])

def prog_hits(panel):
    ii = {gidx[g] for g in panel if g in gidx}
    return {nm: len(ii & mem) for nm, mem in mod_of.items()}

obs = pd.DataFrame([prog_hits(p) for p in lf.panel])
fold_has = (obs > 0).mean(0)                      # fraction of folds hitting each programme
fold_genes = obs.sum(0)                            # total genes drawn from each programme

# size-matched random-panel null over the same gene universe
B = 10000
sizes = lf.panel.map(len).tolist()
null_has = np.zeros((B, len(mod_of)))
null_genes = np.zeros((B, len(mod_of)))
mem_arr = [np.zeros(len(genes), bool) for _ in mod_of]
for j, nm in enumerate(mod_of):
    mem_arr[j][list(mod_of[nm])] = True
for b in range(B):
    hits = np.zeros(len(mod_of)); tot = np.zeros(len(mod_of))
    for s in sizes:
        pick = rng.choice(len(genes), size=s, replace=False)
        for j in range(len(mod_of)):
            n = int(mem_arr[j][pick].sum())
            tot[j] += n; hits[j] += (n > 0)
    null_has[b] = hits / len(sizes); null_genes[b] = tot

prog = pd.DataFrame({
    "programme": list(mod_of),
    "n_genes_in_universe": [len(v) for v in mod_of.values()],
    "folds_hit": [int((obs[nm] > 0).sum()) for nm in mod_of],
    "frac_folds_hit": [float(fold_has[nm]) for nm in mod_of],
    "null_frac_folds_hit": null_has.mean(0),
    "genes_drawn": [int(fold_genes[nm]) for nm in mod_of],
    "null_genes_drawn": null_genes.mean(0),
    "p_empirical": [float((null_genes[:, j] >= fold_genes[list(mod_of)[j]]).mean()) for j in range(len(mod_of))],
})
prog["enrichment"] = prog.genes_drawn / prog.null_genes_drawn.clip(lower=1e-9)
prog = prog.sort_values("genes_drawn", ascending=False)
prog.to_csv(f"{OUT}/programme_frequency.csv", index=False)

prolif_obs = int(fold_genes[PROLIF].sum())
prolif_null = null_genes[:, [list(mod_of).index(p) for p in PROLIF if p in mod_of]].sum(1)
prolif_folds = int((obs[PROLIF].sum(1) > 0).sum())
res["programme_concentration"] = {
    "top_programmes": prog.head(8)[["programme", "folds_hit", "genes_drawn",
                                    "null_genes_drawn", "enrichment", "p_empirical"]].to_dict("records"),
    "proliferation_union": {
        "programmes": PROLIF,
        "genes_drawn": prolif_obs,
        "null_mean": float(prolif_null.mean()),
        "enrichment": float(prolif_obs / prolif_null.mean()),
        "p_empirical": float((prolif_null >= prolif_obs).mean()),
        "folds_hitting_at_least_one": prolif_folds,
        "frac_folds": float(prolif_folds / len(lf))},
    "contrast_with_gene_level": {
        "max_single_gene_folds_out_of_25": 5,
        "note": "compare 5/25 for the best single gene with the proliferation-union fold coverage"}}
print("  proliferation union: %d genes drawn vs null %.1f (enrich %.2fx, p=%.4f), %d/%d folds" % (
    prolif_obs, prolif_null.mean(), prolif_obs / prolif_null.mean(),
    (prolif_null >= prolif_obs).mean(), prolif_folds, len(lf)), flush=True)

# =============================================================== 2/3. score agreement + external C
rows = []
for i, r in lf.iterrows():
    sc = dk_score(r.panel)
    rows.append({"fold": int(r.fold), "repeat": int(r["repeat"]), "k": int(r.k), "scheme": r.scheme,
                 "n_genes": len(r.panel),
                 "overlap_locked": len(set(r.panel) & set(PANEL4)),
                 "pearson_vs_locked": float(np.corrcoef(sc, locked_score)[0, 1]),
                 "spearman_vs_locked": float(stats.spearmanr(sc, locked_score).statistic),
                 "dk_c_5yr": dk_c(sc)})
fold_df = pd.DataFrame(rows)
fold_df.to_csv(f"{OUT}/fold_panel_agreement.csv", index=False)

# size-matched random panels as the reference for "any 4-25 genes would do this"
null_c, null_r = [], []
for b in range(2000):
    s = sizes[b % len(sizes)]
    pick = rng.choice(len(genes), size=s, replace=False)
    sc = (Zd[:, pick] * sgn[pick]).mean(1)
    null_c.append(dk_c(sc)); null_r.append(float(stats.spearmanr(sc, locked_score).statistic))

res["score_agreement"] = {
    "spearman_vs_locked": {"mean": float(fold_df.spearman_vs_locked.mean()),
                           "sd": float(fold_df.spearman_vs_locked.std()),
                           "min": float(fold_df.spearman_vs_locked.min()),
                           "max": float(fold_df.spearman_vs_locked.max())},
    "null_spearman_vs_locked": {"mean": float(np.mean(null_r)), "sd": float(np.std(null_r))},
    "folds_with_zero_locked_genes": int((fold_df.overlap_locked == 0).sum())}
res["external_discrimination"] = {
    "fold_panel_dk_c": {"mean": float(fold_df.dk_c_5yr.mean()), "sd": float(fold_df.dk_c_5yr.std()),
                        "min": float(fold_df.dk_c_5yr.min()), "max": float(fold_df.dk_c_5yr.max()),
                        "median": float(fold_df.dk_c_5yr.median())},
    "locked_dk_c": C_LOCKED,
    "folds_beating_locked": int((fold_df.dk_c_5yr >= C_LOCKED).sum()),
    "random_panel_dk_c": {"mean": float(np.mean(null_c)), "sd": float(np.std(null_c)),
                          "q95": float(np.quantile(null_c, 0.95))},
    "p_fold_mean_vs_random": float(np.mean(np.array(null_c) >= fold_df.dk_c_5yr.mean())),
    "zero_overlap_subset": {
        "n": int((fold_df.overlap_locked == 0).sum()),
        "dk_c_mean": float(fold_df.loc[fold_df.overlap_locked == 0, "dk_c_5yr"].mean()),
        "spearman_mean": float(fold_df.loc[fold_df.overlap_locked == 0, "spearman_vs_locked"].mean())}}
print("  fold-panel DK C: %.3f +- %.3f  (locked %.3f, random %.3f)" % (
    fold_df.dk_c_5yr.mean(), fold_df.dk_c_5yr.std(), C_LOCKED, np.mean(null_c)), flush=True)
print("  spearman vs locked score: %.3f (random %.3f)" % (
    fold_df.spearman_vs_locked.mean(), np.mean(null_r)), flush=True)

json.dump(res, open(f"{OUT}/AXIS_stability.json", "w"), indent=2)
print(f"\nwrote {OUT}/AXIS_stability.json", flush=True)
