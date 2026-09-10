import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

f"""STAGE 4a: build the single-cell reference for deconvolution.

Reads {ACC_DATA_ROOT}/outputs/h5ad_files/inprogress_3.h5ad READ-ONLY
(h5py, mode='r'), streams the sparse `counts` layer in cell blocks, and aggregates raw counts
to a (cell type x patient x gene) tensor.

Everything downstream is derived from this cache, so the 8 GB object is touched exactly once.
Cached to work/outputs_deconv/ - nothing is written back to the h5ad.

Why per-PATIENT as well as per-cell-type: 19 of the 20 JSE bulk patients have matched single-cell
data, so the signature matrix can be rebuilt leaving one patient out. That makes the sanity gate
(does bulk deconvolution recover the true single-cell composition?) non-circular.
"""
import os, json
import numpy as np, pandas as pd, h5py

R = f"{ACC_DATA_ROOT}"
SC = f"{ACC_DATA_ROOT}/outputs/h5ad_files/inprogress_3.h5ad"
OUT = f"{ACC_DATA_ROOT}/model_outputs_deconv"
os.makedirs(OUT, exist_ok=True)

BLOCK = 10000

# coarse grouping of annot_5 - biologically motivated, fixed here before any survival analysis
COARSE = {
    "Myoepithelial cells - Tumor": "Myoepithelial_Tumor",
    "Actively Dividing Myoepithelial cells - Tumor": "Dividing",
    "Actively Dividing cells": "Dividing",
    "Epithelial cells - Tumor": "Epithelial_Tumor",
    "Epithelial cells - Basal - Tumor": "Epithelial_Tumor",
    "Basal cells - Tumor": "Epithelial_Tumor",
    "Epithelial cells - Basal": "Epithelial_Normal",
    "Epithelial cells - EA": "Epithelial_Normal",
    "Epithelial cells - Ductal": "Epithelial_Normal",
    "Epithelial cells - Glandular": "Epithelial_Normal",
    "Epithelial cells - Secretory": "Epithelial_Normal",
    "Epithelial cells - Ciliated": "Epithelial_Normal",
    "Fibroblast cells": "Fibroblast",
    "Fibroblast cells - CAF": "Fibroblast",
    "Fibroblast cells - Inflammatory": "Fibroblast",
    "Endothelial cells": "Vascular",
    "Mural cells": "Vascular",
    "Macrophage - M2": "Myeloid",
    "Macrophage - AP/TAM": "Myeloid",
    "Mast cells": "Myeloid",
    "CD4+ T cells": "Lymphoid",
    "CD8+ T cells": "Lymphoid",
    "B cells": "Lymphoid",
    "Memory B cells": "Lymphoid",
    "Plasma cells": "Lymphoid",
    "Muscle cells - Tongue": "Muscle",
}


def read_cat(f, name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group):
        cats = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
        return np.asarray(cats), g["codes"][:]
    v = g[:]
    v = np.array([x.decode() if isinstance(x, bytes) else x for x in v])
    cats, codes = np.unique(v, return_inverse=True)
    return cats, codes


with h5py.File(SC, "r") as f:
    genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var/_index"][:]])
    ct_cats, ct_codes = read_cat(f, "annot_5")
    pt_cats, pt_codes = read_cat(f, "batch")
    oc_cats, oc_codes = read_cat(f, "clinical_outcome")

    n_cells, n_genes = f["layers/counts"].attrs["shape"]
    n_ct, n_pt = len(ct_cats), len(pt_cats)
    print(f"cells={n_cells}  genes={n_genes}  cell types={n_ct}  patients={n_pt}")

    # patient -> outcome (constant within patient)
    pt_outcome = {}
    for pi in range(n_pt):
        m = pt_codes == pi
        pt_outcome[pt_cats[pi]] = str(oc_cats[np.bincount(oc_codes[m], minlength=len(oc_cats)).argmax()])

    # group key = cell type * n_pt + patient
    key = ct_codes.astype(np.int64) * n_pt + pt_codes.astype(np.int64)
    n_key = n_ct * n_pt
    agg = np.zeros((n_key, n_genes), dtype=np.float64)
    ncell = np.bincount(key, minlength=n_key)

    dat = f["layers/counts/data"]
    ind = f["layers/counts/indices"]
    iptr = f["layers/counts/indptr"][:]

    for s in range(0, n_cells, BLOCK):
        e = min(s + BLOCK, n_cells)
        lo, hi = int(iptr[s]), int(iptr[e])
        d = dat[lo:hi].astype(np.float64)
        c = ind[lo:hi].astype(np.int64)
        rowlen = np.diff(iptr[s:e + 1]).astype(np.int64)
        rows = np.repeat(np.arange(s, e, dtype=np.int64), rowlen)
        np.add.at(agg, (key[rows], c), d)
        print(f"  {e}/{n_cells} cells", end="\r")
print()

ct_names = np.array([f"{ct_cats[i]}" for i in range(n_ct)])
idx = pd.MultiIndex.from_product([ct_cats, pt_cats], names=["cell_type", "patient"])
pb = pd.DataFrame(agg, index=idx, columns=genes)
pb = pb.loc[ncell > 0]                       # drop empty (cell type, patient) combos
ncell_s = pd.Series(ncell, index=idx).loc[pb.index]

pb.to_parquet(f"{OUT}/pseudobulk_ct_by_patient.parquet")
ncell_s.rename("n_cells").to_frame().to_parquet(f"{OUT}/ncells_ct_by_patient.parquet")
json.dump({"cell_types": ct_cats.tolist(), "patients": pt_cats.tolist(),
           "patient_outcome": pt_outcome, "coarse_map": COARSE,
           "n_cells_total": int(n_cells), "n_genes_sc": int(n_genes)},
          open(f"{OUT}/reference_meta.json", "w"), indent=2)

# true single-cell composition per patient (fine and coarse) - the answer key for the sanity gate
comp = ncell_s.unstack("cell_type").fillna(0.0)
comp_frac = comp.div(comp.sum(1), axis=0)
comp_frac.to_csv(f"{OUT}/sc_true_composition_fine.csv")

coarse_of = pd.Series({c: COARSE.get(c, "Other") for c in comp.columns})
comp_c = comp.T.groupby(coarse_of).sum().T
comp_c_frac = comp_c.div(comp_c.sum(1), axis=0)
comp_c_frac.to_csv(f"{OUT}/sc_true_composition_coarse.csv")

print(f"\npseudobulk: {pb.shape[0]} (cell type x patient) groups x {pb.shape[1]} genes")
print(f"cells per patient: median {comp.sum(1).median():.0f}, min {comp.sum(1).min():.0f}")
print("\ntrue coarse composition (mean fraction across patients):")
print(comp_c_frac.mean(0).sort_values(ascending=False).to_string(float_format=lambda v: f"{v:.4f}"))
print(f"\ncached -> {OUT}")
