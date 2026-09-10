"""Shared data assembly for the t3 rebuild. Read-only on all project inputs."""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import numpy as np, pandas as pd
from scipy import stats

ROOT = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"


def rint(df):
    """Rank-inverse-normal transform within cohort: removes platform scale/shape entirely."""
    r = df.rank(axis=0)
    n = df.shape[0]
    return pd.DataFrame(stats.norm.ppf((r.values - 0.375) / (n + 0.25)),
                        index=df.index, columns=df.columns)


def load(pool="evidence"):
    dk = pd.read_csv(f"{ROOT}/outputs/03_harmonized/dk_norm_log2.csv", index_col=0)
    tr = pd.read_csv(f"{ROOT}/outputs/03_harmonized/train_norm_log2.csv", index_col=0)
    lab = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/sample_labels.csv").set_index("sample")["outcome"]
    surv = pd.read_csv(f"{ROOT}/external_dk/DK_sample_data_table2.csv").set_index("UID")
    ev = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")

    s = surv.loc[dk.index]
    ok = (pd.to_numeric(s["SURV_MO"], errors="coerce").notna()
          & pd.to_numeric(s["SURV_CENS"], errors="coerce").notna())
    s = s[ok]
    T = pd.to_numeric(s["SURV_MO"]).values.astype(float)
    E = pd.to_numeric(s["SURV_CENS"]).values.astype(int)
    dk = dk.loc[s.index]

    # gene pool (DK-survival independent by construction)
    if pool == "evidence":
        genes = [g for g in ev.index if g in dk.columns and g in tr.columns]
    elif pool == "candidates":
        cand = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/candidate_genes.csv")
        genes = [g for g in cand.gene if g in dk.columns and g in tr.columns]
    else:
        genes = [g for g in dk.columns if g in tr.columns]

    # drop degenerate genes (constant / near-constant in either cohort)
    d, t = dk[genes], tr[genes]
    keep = (d.std(axis=0) > 1e-6) & (t.std(axis=0) > 1e-6) & (d.gt(0).mean(axis=0) >= 0.3)
    genes = [g for g in genes if keep[g]]

    Zd = rint(dk[genes])
    Zt = rint(tr[genes])
    y = (lab.loc[tr.index] == "Poor").values.astype(int)

    evg = ev.reindex(genes)
    return dict(genes=genes, Zd=Zd, Zt=Zt, y=y, T=T, E=E, ev=evg, clin=s, dk_raw=dk[genes], tr_raw=tr[genes])


def evidence_prior(evg, floor=0.35):
    """Evidence multiplier s_g in (floor, ~1.6] and direction d_g in {-1,+1}.

    Built ONLY from discovery-cohort evidence (WGCNA / single-cell DE / bulk DE).
    DK survival contributes nothing.
    """
    def pct(x):
        x = pd.Series(x).fillna(0.0).values
        return stats.rankdata(x) / len(x)

    sc = pct(-np.log10(evg["sc_best_padj"].fillna(1.0).clip(1e-300)) * evg["in_sc_de"].fillna(False).astype(float))
    wg = pct(evg["wgcna_gs"].abs().fillna(0.0) * evg["in_wgcna_hub"].fillna(False).astype(float))
    bd = pct(evg["bulk_de_log2fc"].abs().fillna(0.0) * evg["in_bulk_de"].fillna(False).astype(float))
    nsrc = evg["n_sources"].fillna(1).values.astype(float)

    strength = 0.4 * sc + 0.3 * wg + 0.3 * bd
    s = floor + 1.25 * (0.5 * strength + 0.5 * (nsrc - 1) / 2.0)
    d = np.where(evg["direction_consensus"].values == "Poor", 1.0,
                 np.where(evg["direction_consensus"].values == "Good", -1.0, 0.0))
    return s.astype(np.float32), d.astype(np.float32)
