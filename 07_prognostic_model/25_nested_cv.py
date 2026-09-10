"""
Stage 7 - leakage-free nested cross-validation of the discovery procedure.

WHY
---
The reported discovery AUC of 0.906 is optimistic. `13_compact_panel.py` refits the Bayesian
model inside each fold, but the things the model is *conditioned on* were built once from all 20
patients using the outcome:

  * `gene_evidence_wide.csv`  - a 1,803-gene candidate list whose membership is set by
    outcome-derived bulk DE, single-cell DE and WGCNA gene-significance;
  * `evidence_prior()`        - strength `s_g` and direction `d_g` (direction_consensus = Poor/Good);
  * `modules.build()`         - the `sc_*_up_poor / up_good` and `evidence_poor / good` module
    columns, i.e. outcome-derived *features*;
  * `rint()`                  - fitted across all 20 patients at once;
  * panel size k and scoring scheme - chosen by inspecting the very CV curve being reported.

This script re-runs the whole procedure with an outer loop that never sees held-out information,
and quantifies the leak by running the same harness in two arms.

  ARM A  "as_built"      - reproduces the current pipeline. Sanity check: should land near 0.906.
  ARM B  "leakage_free"  - evidence, priors, outcome-derived modules, RINT and (k, scheme) all
                           rebuilt from the training patients of that fold alone.

Held fixed in both arms, and declared as such in the paper:
  * curated module membership          - literature-derived, no data
  * WGCNA *module assignment* and MM   - co-expression topology, outcome-free
    (its supervised gene-significance GS *is* refit in Arm B)
  * the gene-universe expression filter

Disclosed approximation: in-fold DE uses a moderated t-test (limma-style eBayes shrinkage) rather
than PyDESeq2, which is impractical at 150+ folds. The agreement between the two on the full data
is measured and reported so the substitution can be judged.

DK is not touched anywhere in this file.
Reads only; writes to work/outputs_nested/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import sys
import time
import importlib.util
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_nested"
CT = f"{R}/outputs/01_candidate_table"
WGCNA = f"{ACC_DATA_ROOT}/bulk_rna_data/outputs/images/clinical_outcome/wgcna/06_gene_info_complete.csv"

sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3)
import modules as MOD
from mlbayes import MultiLevelModel, to_t

STEPS = 900
KS = [3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30]
SCHEMES = ("sign", "beta")
PARSIMONY_TOL = 0.02


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


# ===================== transforms =====================
def rint_train(X):
    """Standard RINT on the training matrix (samples x genes)."""
    r = stats.rankdata(X, axis=0)
    n = X.shape[0]
    return stats.norm.ppf((r - 0.375) / (n + 0.25))


def rint_apply(X_new, X_ref):
    """Map new samples into the training cohort's rank space (frozen reference quantiles).

    For each gene, a new value is placed among the sorted training values and converted to a
    normal deviate. This is what a deployed assay would have to do: a held-out sample never
    influences its own transform, and never sees the other held-out samples.
    """
    ref = np.sort(X_ref, axis=0)
    n = ref.shape[0]
    out = np.empty_like(X_new, dtype=float)
    for j in range(X_new.shape[1]):
        pos = np.searchsorted(ref[:, j], X_new[:, j], side="left")
        p = (pos + 0.5) / (n + 1.0)
        out[:, j] = stats.norm.ppf(p)
    return out


def log2cpm(counts):
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log2(counts / lib * 1e6 + 1.0)


def moderated_t(X, grp, df_prior=4.0):
    """limma-style moderated t-test. X samples x genes (log2 scale), grp boolean = 'Poor'.

    Returns (log2fc, pval, padj) with log2fc oriented Poor - Good.
    """
    a, b = X[grp], X[~grp]
    n1, n2 = a.shape[0], b.shape[0]
    if n1 < 2 or n2 < 2:
        z = np.zeros(X.shape[1])
        return z, np.ones(X.shape[1]), np.ones(X.shape[1])
    diff = a.mean(0) - b.mean(0)
    df = n1 + n2 - 2
    s2 = ((n1 - 1) * a.var(0, ddof=1) + (n2 - 1) * b.var(0, ddof=1)) / df
    s2_prior = np.median(s2[s2 > 0]) if np.any(s2 > 0) else 1.0
    s2_post = (df_prior * s2_prior + df * s2) / (df_prior + df)
    se = np.sqrt(s2_post * (1.0 / n1 + 1.0 / n2))
    se[se <= 0] = np.inf
    t = diff / se
    p = 2 * stats.t.sf(np.abs(t), df + df_prior)
    order = np.argsort(p)
    m = len(p)
    padj = np.empty(m)
    padj[order] = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    return diff, p, np.clip(padj, 0, 1)


# ===================== data =====================
print("=" * 74)
print("Stage 7 - leakage-free nested CV")
print("=" * 74)

D = d3.load("all")
genes_all = D["genes"]
tr_raw, dk_raw = D["tr_raw"], D["dk_raw"]
y = D["y"].astype(int)

# gene universe: identical filter to 13_compact_panel.py:57 (uses DK *expression*, no survival)
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
gi = {g: i for i, g in enumerate(genes)}
Xtr_raw = tr_raw[genes].values            # log2 normalised expression, 20 x P
samples = list(tr_raw.index)
print(f"patients={len(y)}  genes={len(genes)}  poor={int(y.sum())}")

# ---- raw sources for in-fold recomputation ----
bulk_counts = pd.read_csv(f"{CT}/bulk_expression/counts.csv", index_col=0)
bulk_counts = bulk_counts.reindex(samples)
bulk_shared = [g for g in genes if g in bulk_counts.columns]
bulk_l2 = log2cpm(bulk_counts[bulk_shared].values.astype(float))
bulk_gi = {g: i for i, g in enumerate(bulk_shared)}
print(f"bulk counts for in-fold DE: {bulk_counts.shape} -> {len(bulk_shared)} genes in universe")

sc_meta = pd.read_csv(f"{CT}/sc_pseudobulk/pseudobulk_meta.csv")
sc_counts = pd.read_csv(f"{CT}/sc_pseudobulk/pseudobulk_counts_donor_x_celltype.csv", index_col=0)
SC_TYPES = ["Epithelial cells - Tumor", "Myoepithelial cells - Tumor", "Actively Dividing cells"]
sc_meta = sc_meta[sc_meta.cell_type.isin(SC_TYPES) & sc_meta.kept]
sc_shared = [g for g in genes if g in sc_counts.columns]
sc_gi = {g: i for i, g in enumerate(sc_shared)}
print(f"sc pseudobulk: {sc_counts.shape}, donors={sc_meta.donor.nunique()}, "
      f"{len(sc_shared)} genes in universe")
print(f"bulk patients with matched single cell: {len(set(samples) & set(sc_meta.donor))}/20")

wg = pd.read_csv(WGCNA)
wg.columns = [c.strip('"') for c in wg.columns]
wg_gene_col = wg.columns[0]
wg = wg.set_index(wg_gene_col)
wg_module = wg["Module"].reindex(genes)
wg_mm = wg.filter(like="MM.").abs().max(axis=1).reindex(genes)     # unsupervised
print(f"WGCNA: {wg_module.notna().sum()} genes mapped, "
      f"{wg_module.dropna().nunique()} modules")

# ---- fixed (as-built) artefacts ----
EV_FULL = pd.read_csv(f"{CT}/evidence/gene_evidence_wide.csv").set_index("gene")

# Evidence-table sizes, calibrated once from the stored full-data table restricted to this gene
# universe. In-fold tables are built by taking the top-N genes by significance rather than by
# applying a p-value threshold: PyDESeq2 ran on the unfiltered ~59k matrix with its own independent
# filtering, so a padj cut-off is not transferable to a moderated t-test on 12,700 genes. Fixing N
# keeps the *strength* of the prior identical across folds and arms, so the Arm A / Arm B contrast
# isolates the information the prior is built from rather than how much prior there is.
_inu = EV_FULL.index.isin(set(genes))
N_BULK = int((EV_FULL["in_bulk_de"][_inu] == True).sum())
N_SC = int((EV_FULL["in_sc_de"][_inu] == True).sum())
N_HUB = int((EV_FULL["in_wgcna_hub"][_inu] == True).sum())
print(f"in-fold evidence sizes calibrated to stored table: "
      f"bulk={N_BULK}  sc={N_SC}  wgcna_hub={N_HUB}")
M_FULL, MN_FULL, MS_FULL = MOD.build(genes, R)
SRC_FULL = sorted(set(MS_FULL))
SRCID_FULL = np.array([SRC_FULL.index(s) for s in MS_FULL])
S_FULL, DIR_FULL = d3.evidence_prior(EV_FULL.reindex(genes))
Z_FULL = d3.rint(tr_raw[genes]).values
LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
CFG = LOCK["config"]
print(f"as-built artefacts: modules={M_FULL.shape[1]}, evidence rows={len(EV_FULL)}, "
      f"locked cfg={CFG}")


# ===================== in-fold evidence + modules =====================
def evidence_infold(tr_idx):
    """Rebuild the candidate/evidence table using training patients only.

    Selection is top-N by significance, with N fixed to the stored table's size in this universe
    (see the N_BULK / N_SC / N_HUB calibration above).
    """
    tr_samples = [samples[i] for i in tr_idx]
    ytr = y[tr_idx].astype(bool)

    # --- bulk DE, training patients only ---
    lfc_b, p_b, _ = moderated_t(bulk_l2[tr_idx], ytr)
    bulk_lfc = np.zeros(len(genes))
    bulk_p = np.ones(len(genes))
    idx = np.array([gi[g] for g in bulk_shared])
    bulk_lfc[idx] = lfc_b
    bulk_p[idx] = p_b
    in_bulk = np.zeros(len(genes), bool)
    in_bulk[np.argsort(bulk_p)[:N_BULK]] = True

    # --- single-cell DE, dropping donors of held-out bulk patients ---
    drop = set(samples) - set(tr_samples)
    sc_best_padj = np.ones(len(genes))
    sc_dir_score = np.zeros(len(genes))
    sc_sets = {}
    for ct in SC_TYPES:
        mm = sc_meta[(sc_meta.cell_type == ct) & (~sc_meta.donor.isin(drop))]
        if mm.donor.nunique() < 6:
            continue
        rows = [f"{d}||{ct}" for d in mm.donor]
        rows = [r for r in rows if r in sc_counts.index]
        if len(rows) < 6:
            continue
        sub = sc_counts.loc[rows, sc_shared].values.astype(float)
        grp = mm.set_index("donor").loc[[r.split("||")[0] for r in rows], "outcome"].values == "Poor"
        if grp.sum() < 2 or (~grp).sum() < 2:
            continue
        lfc_s, p_s, _ = moderated_t(log2cpm(sub), grp)
        ii = np.array([gi[g] for g in sc_shared])
        sc_best_padj[ii] = np.minimum(sc_best_padj[ii], p_s)
        # per-cell-type gene sets keep the same size as the stored single-cell modules
        n_set = max(int(N_SC / len(SC_TYPES)), 20)
        sig = np.zeros(len(p_s), bool)
        sig[np.argsort(p_s)[:n_set * 2]] = True
        sc_dir_score[ii] += np.sign(lfc_s) * sig
        tag = ct.replace(" ", "_").replace("-", "")
        sc_sets[f"sc_{tag}_up_poor"] = [sc_shared[j] for j in np.where(sig & (lfc_s > 0))[0]]
        sc_sets[f"sc_{tag}_up_good"] = [sc_shared[j] for j in np.where(sig & (lfc_s < 0))[0]]
    in_sc = np.zeros(len(genes), bool)
    in_sc[np.argsort(sc_best_padj)[:N_SC]] = True

    # --- WGCNA gene significance, training patients only (module/MM stay fixed) ---
    Xc = Xtr_raw[tr_idx]
    Xc = Xc - Xc.mean(0)
    yc = ytr.astype(float) - ytr.mean()
    denom = np.sqrt((Xc ** 2).sum(0) * (yc ** 2).sum())
    denom[denom == 0] = np.inf
    gs = (Xc * yc[:, None]).sum(0) / denom
    mm_abs = np.nan_to_num(wg_mm.values.astype(float))
    elig = mm_abs > 0.8                       # module membership: unsupervised, fixed
    in_hub = np.zeros(len(genes), bool)
    cand = np.where(elig)[0]
    in_hub[cand[np.argsort(-np.abs(gs[cand]))[:N_HUB]]] = True

    # --- assemble, same schema evidence_prior() consumes ---
    dir_bulk = np.where(in_bulk, np.where(bulk_lfc > 0, 1, -1), 0)
    dir_sc = np.where(in_sc, np.sign(sc_dir_score), 0)
    dir_wg = np.where(in_hub, np.where(gs > 0, 1, -1), 0)
    n_sources = in_bulk.astype(int) + in_sc.astype(int) + in_hub.astype(int)
    consensus = np.sign(dir_bulk + dir_sc + dir_wg)

    ev = pd.DataFrame({
        "in_wgcna_hub": in_hub, "wgcna_gs": gs,
        "in_sc_de": in_sc, "sc_best_padj": sc_best_padj,
        "in_bulk_de": in_bulk, "bulk_de_log2fc": bulk_lfc,
        "n_sources": np.maximum(n_sources, 1),
        "direction_consensus": np.where(consensus > 0, "Poor",
                                        np.where(consensus < 0, "Good", "None")),
    }, index=genes)
    # genes with no evidence at all are absent from the stored table -> neutral prior
    ev.loc[n_sources == 0, ["wgcna_gs", "bulk_de_log2fc"]] = 0.0
    ev.loc[n_sources == 0, "sc_best_padj"] = 1.0
    return ev, sc_sets


def modules_infold(sc_sets, ev):
    """Module matrix with the outcome-derived columns rebuilt from this fold's data."""
    cols, names, src = [], [], []

    def add(name, gset, source):
        hit = [gi[g] for g in gset if g in gi]
        if len(hit) < 8:
            return
        v = np.zeros(len(genes), np.float32)
        v[hit] = 1.0 / np.sqrt(len(hit))
        cols.append(v); names.append(name); src.append(source)

    for k, v in MOD.CURATED.items():                       # literature, no data
        add(k, v.split(), "curated")
    for mod, sub in wg_module.dropna().groupby(wg_module.dropna()):   # co-expression topology
        if mod == "grey":
            continue
        gl = list(sub.index)
        if 20 <= len(gl) <= 2500:
            add(f"wgcna_{mod}", gl, "wgcna")
    for nm, gl in sc_sets.items():                          # in-fold single-cell DE
        add(nm, gl, "single_cell")
    for tag, d in [("poor", "Poor"), ("good", "Good")]:     # in-fold evidence consensus
        add(f"evidence_{tag}", ev.index[ev.direction_consensus == d].tolist(), "evidence")
        sub = ev[(ev.direction_consensus == d) & (ev.n_sources >= 3)]
        add(f"evidence3_{tag}", sub.index.tolist(), "evidence")

    M = np.stack(cols, 1)
    srcs = sorted(set(src))
    return M, np.array([srcs.index(s) for s in src]), len(srcs)


# ===================== model =====================
def fit_beta(Z, yy, M, src_ids, n_src, s_ev, d_ev, seed):
    m = MultiLevelModel(M, src_ids, n_src, s_ev, d_ev,
                        tau0=CFG["tau0"], omega0=CFG["omega0"],
                        use_gene_level=bool(CFG["gene_level"]),
                        use_module_level=bool(CFG["module_level"]),
                        use_dk=False, use_jse=True, seed=seed)
    m.fit(Xj=to_t(Z), yj=to_t(yy), steps=STEPS)
    return m.posterior(300)["beta"]


def panel_score(Z, idx_g, w, scheme):
    if scheme == "sign":
        return (Z[:, idx_g] * np.sign(w)).mean(1)
    return Z[:, idx_g] @ w


def fold_artifacts(tr_idx, arm):
    """Everything the model conditions on, built appropriately for the arm."""
    if arm == "as_built":
        Ztr, Zte_fn = Z_FULL[tr_idx], (lambda te: Z_FULL[te])
        return Ztr, Zte_fn, M_FULL, SRCID_FULL, len(SRC_FULL), S_FULL, DIR_FULL
    ev, sc_sets = evidence_infold(tr_idx)
    M, src_ids, n_src = modules_infold(sc_sets, ev)
    s_ev, d_ev = d3.evidence_prior(ev)
    Ztr = rint_train(Xtr_raw[tr_idx])
    Zte_fn = lambda te: rint_apply(Xtr_raw[te], Xtr_raw[tr_idx])
    return Ztr, Zte_fn, M, src_ids, n_src, s_ev, d_ev


def select_k_scheme(tr_idx, arm, seed):
    """Inner loop: parsimony rule applied on training patients only."""
    ytr = y[tr_idx]
    inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    pred = {(k, s): np.full(len(tr_idx), np.nan) for k in KS for s in SCHEMES}
    for itr, ite in inner.split(np.zeros(len(ytr)), ytr):
        gtr = tr_idx[itr]
        Ztr, Zte_fn, M, sid, nsrc, s_ev, d_ev = fold_artifacts(gtr, arm)
        b = fit_beta(Ztr, y[gtr], M, sid, nsrc, s_ev, d_ev, seed)
        Zte = Zte_fn(tr_idx[ite])
        order = np.argsort(-np.abs(b))
        for k in KS:
            idx_g = order[:k]
            for sch in SCHEMES:
                pred[(k, sch)][ite] = panel_score(Zte, idx_g, b[idx_g], sch)
    rows = []
    for k in KS:
        for sch in SCHEMES:
            p = pred[(k, sch)]
            ok = ~np.isnan(p)
            auc = roc_auc_score(ytr[ok], p[ok]) if len(np.unique(ytr[ok])) > 1 else 0.5
            rows.append((k, sch, auc))
    tab = pd.DataFrame(rows, columns=["k", "scheme", "auc"])
    best = tab.auc.max()
    elig = tab[tab.auc >= best - PARSIMONY_TOL].sort_values(["k", "auc"], ascending=[True, False])
    return int(elig.iloc[0]["k"]), str(elig.iloc[0]["scheme"]), float(best)


# ===================== outer loop =====================
def run_arm(arm, n_repeats=5, tag=""):
    print("\n" + "=" * 74)
    print(f"ARM: {arm}{tag}")
    print("=" * 74)
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=7)
    oof = {r: np.full(len(y), np.nan) for r in range(n_repeats)}
    picks = []
    t0 = time.time()
    for f, (tr_idx, te_idx) in enumerate(rskf.split(np.zeros(len(y)), y)):
        rep = f // 5
        k, sch, inner_auc = select_k_scheme(tr_idx, arm, seed=100 + f)
        Ztr, Zte_fn, M, sid, nsrc, s_ev, d_ev = fold_artifacts(tr_idx, arm)
        b = fit_beta(Ztr, y[tr_idx], M, sid, nsrc, s_ev, d_ev, seed=100 + f)
        order = np.argsort(-np.abs(b))
        idx_g = order[:k]
        oof[rep][te_idx] = panel_score(Zte_fn(te_idx), idx_g, b[idx_g], sch)
        picks.append({"fold": f, "repeat": rep, "k": k, "scheme": sch,
                      "inner_auc": inner_auc,
                      "genes": ";".join(genes[i] for i in idx_g)})
        print(f"  fold {f+1:2d}/{n_repeats*5}  k={k:2d} scheme={sch:4s} "
              f"inner_auc={inner_auc:.3f}  ({time.time()-t0:.0f}s)")
    aucs = [roc_auc_score(y, oof[r]) for r in range(n_repeats)]
    return {"arm": arm + tag, "auc_mean": float(np.mean(aucs)), "auc_sd": float(np.std(aucs)),
            "auc_per_repeat": [float(a) for a in aucs],
            "k_selected": pd.Series([p["k"] for p in picks]).value_counts().to_dict(),
            "scheme_selected": pd.Series([p["scheme"] for p in picks]).value_counts().to_dict()}, picks


# ---- validate the DE surrogate before using it ----
print("\n" + "-" * 74)
print("validating the moderated-t surrogate against the stored PyDESeq2 ranking")
lfc_all, p_all, _ = moderated_t(bulk_l2, y.astype(bool))
ref = EV_FULL.reindex(bulk_shared)
m_ok = ref["bulk_de_log2fc"].notna().values
# the stored log2FC is oriented Good-Poor; evidence_prior() uses only |log2FC|, so compare magnitudes
rho_lfc = stats.spearmanr(np.abs(lfc_all[m_ok]), np.abs(ref["bulk_de_log2fc"].values[m_ok]))
rho_signed = stats.spearmanr(lfc_all[m_ok], -ref["bulk_de_log2fc"].values[m_ok])
sig_ref = set(ref.index[(ref.in_bulk_de == True)])
sig_new = set(np.array(bulk_shared)[np.argsort(p_all)[:len(sig_ref)]])
jac = len(sig_ref & sig_new) / max(len(sig_ref | sig_new), 1)
print(f"  |log2FC| Spearman rho vs PyDESeq2 = {rho_lfc.correlation:.3f} (n={int(m_ok.sum())})")
print(f"  signed log2FC (orientation-corrected) rho = {rho_signed.correlation:.3f}")
print(f"  top-{len(sig_ref)} DE-gene overlap with stored set: "
      f"{len(sig_ref & sig_new)}/{len(sig_ref)}  (Jaccard {jac:.3f})")
surrogate = {"abs_lfc_spearman_vs_pydeseq2": float(rho_lfc.correlation),
             "signed_lfc_spearman_orientation_corrected": float(rho_signed.correlation),
             "topN_de_overlap": len(sig_ref & sig_new), "topN": len(sig_ref),
             "de_set_jaccard": float(jac),
             "note": "stored log2FC is oriented Good-Poor; evidence_prior uses |log2FC| only"}

res = {"n_patients": int(len(y)), "n_genes": len(genes),
       "reported_disc_auc_13_compact": 0.906,
       "de_surrogate_validation": surrogate,
       "design": {"outer": "RepeatedStratifiedKFold 5x5",
                  "inner": "StratifiedKFold 5 on training patients, parsimony rule",
                  "parsimony_tol": PARSIMONY_TOL,
                  "fixed_across_folds": ["curated modules", "WGCNA module assignment + MM",
                                         "gene-universe expression filter"],
                  "refit_in_fold_armB": ["bulk DE", "single-cell DE (held-out donors dropped)",
                                         "evidence table + n_sources + direction_consensus",
                                         "WGCNA gene-significance", "evidence_prior s and d",
                                         "outcome-derived modules", "RINT", "k and scheme"]}}

armA, picksA = run_arm("as_built")
armB, picksB = run_arm("leakage_free")
res["arm_as_built"] = armA
res["arm_leakage_free"] = armB
res["measured_leak_auc"] = float(armA["auc_mean"] - armB["auc_mean"])

print("\n" + "=" * 74)
print(f"as_built      AUC = {armA['auc_mean']:.3f} +/- {armA['auc_sd']:.3f}")
print(f"leakage_free  AUC = {armB['auc_mean']:.3f} +/- {armB['auc_sd']:.3f}")
print(f"measured leak     = {res['measured_leak_auc']:+.3f} AUC")
print("=" * 74)

pd.DataFrame(picksA + picksB).assign(
    arm=["as_built"] * len(picksA) + ["leakage_free"] * len(picksB)
).to_csv(f"{OUT}/fold_selections.csv", index=False)

# ---- figure ----
fig, ax = plt.subplots(1, 2, figsize=(9, 3.8))
C1, C2 = "#B4436C", "#3C6E8F"
ax[0].bar([0, 1], [armA["auc_mean"], armB["auc_mean"]],
          yerr=[armA["auc_sd"], armB["auc_sd"]], color=[C2, C1], width=.55, capsize=4)
ax[0].axhline(0.906, ls="--", c="0.4", lw=1)
ax[0].text(1.45, 0.908, "reported 0.906", fontsize=8, ha="right", color="0.35")
ax[0].axhline(0.5, ls=":", c="0.6", lw=1)
ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["as-built", "leakage-free"])
ax[0].set_ylabel("nested CV AUC"); ax[0].set_ylim(0.4, 1.0)
ax[0].set_title("discovery AUC, honest vs as-built", fontsize=10)

kk = pd.DataFrame(picksB)["k"].value_counts().sort_index()
ax[1].bar(kk.index.astype(str), kk.values, color=C1, width=.7)
ax[1].set_xlabel("panel size k selected in the inner loop")
ax[1].set_ylabel("outer folds")
ax[1].set_title("leakage-free k selection stability", fontsize=10)
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
save(fig, "NESTED_cv_summary")

json.dump(res, open(f"{OUT}/NESTED_results.json", "w"), indent=1)
print(f"\nwrote {OUT}/NESTED_results.json")
