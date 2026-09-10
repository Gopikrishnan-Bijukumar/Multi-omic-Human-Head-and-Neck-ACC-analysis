<div align="center">

# Head & Neck Adenoid Cystic Carcinoma — Multi-omic Analysis

**Bulk RNA-seq · single-cell RNA-seq · spatial transcriptomics · a locked four-gene prognostic panel**

[![License: MIT](https://img.shields.io/badge/License-MIT-0b7285.svg)](LICENSE)
[![Python 3.11.15](https://img.shields.io/badge/Python-3.11.15-3776ab.svg)](requirements.txt)
[![R 4.4.3](https://img.shields.io/badge/R-4.4.3-276dc3.svg)](r_dependencies.txt)

</div>

---

## The study in one paragraph

Adenoid cystic carcinoma of the head and neck is indolent for years and then, in some patients,
is not. This repository holds the analysis code for a multi-omic comparison of tumours from
patients with **poor (< 2 years)** and **favourable (> 5 years)** survival, across bulk RNA-seq,
single-cell RNA-seq and spatial transcriptomics — and the derivation of a four-gene prognostic
panel validated in three independent external cohorts.

Every top-level folder corresponds to a numbered section of the paper's Supplementary Methods.

> ### The panel
>
> **`DSCAM` · `ODC1` · `NCAPG` · `CCNB2`** — scored as the unweighted mean of sign-weighted
> rank-inverse-normal expression, higher meaning worse prognosis.
>
> | | Danish cohort (n = 54) |
> |---|---|
> | 5-year overall survival, C-index | **0.671** (0.541–0.801) |
> | Hazard ratio per score SD | **1.74** (1.10–2.75), p = 0.017 |
> | Time-dependent AUC at 60 months | **0.688** (0.526–0.832) |
>
> Separately, the same score separates the MD Anderson ACC-I / ACC-II molecular subtypes at
> **AUC 0.884** (0.773–0.968), replicated in Frerich 2018.
>
> **No expression or outcome data from any external cohort contributed to gene selection or
> model fitting at any stage.** A three-level hierarchical Bayesian shrinkage procedure
> (source → module → gene, fitted by ADVI) chose the genes and their directions; it does not
> weight the deployed score, and no posterior quantity enters any reported C-index or hazard
> ratio.

---

## How the pieces fit together

```mermaid
flowchart TD
    subgraph DISC["Internal discovery cohort — 25 patients, pseudonymised P01–P25"]
        B["<b>01</b> Bulk RNA-seq<br/>20 tumours · 10 sites<br/><i>PyDESeq2 · GSEA · WGCNA · ORA</i>"]
        S["<b>02</b> Single-cell RNA-seq<br/>24 samples<br/><i>scVI · 3-way doublet consensus</i>"]
        C["<b>03</b> scCODA<br/><i>cell-type abundance</i>"]
        P["<b>04</b> CellPhoneDB<br/><i>cell–cell communication</i>"]
        T["<b>05</b> Trajectory<br/><i>CytoTRACE2 · Monocle3</i>"]
        X["<b>06</b> Spatial<br/>4 Visium 55µm · 2 HD 16µm<br/><i>Giotto · SpatialDWLS</i>"]
    end

    S --> C
    S --> P
    S --> T
    S --> X

    B ==> CAND["<b>Candidate gene nomination</b><br/>WGCNA modules · pseudobulk DE · spatial evidence"]
    S ==> CAND
    X ==> CAND

    CAND ==> BAYES["<b>07</b> Three-level hierarchical<br/>Bayesian shrinkage · ADVI<br/><i>source → module → gene</i>"]
    BAYES ==> PANEL["<b>LOCKED FOUR-GENE PANEL</b><br/>DSCAM · ODC1 · NCAPG · CCNB2"]

    PANEL ==> DK["<b>Danish cohort</b> · n = 54<br/>overall survival<br/><b>C = 0.671</b> · HR/SD 1.74<br/><i>the single primary claim</i>"]
    PANEL --> CCR["<b>CCR2020</b> MD Anderson · n = 54<br/>ACC-I vs ACC-II subtype<br/><b>AUC = 0.884</b>"]
    PANEL --> FRE["<b>Frerich 2018</b> · n = 66<br/>cross-platform replication<br/><b>OR = 7.5</b>"]

    %% ---- UNM palette: Cherry #BA0C2F, Turquoise #00A9B7, Zia Yellow #FFC72C, Silver #63666A
    style DISC fill:#f7f8f8,stroke:#63666A,stroke-width:2px,color:#2b2d2f
    style B    fill:#d4f0f3,stroke:#00A9B7,stroke-width:2px,color:#05343a
    style S    fill:#d4f0f3,stroke:#00A9B7,stroke-width:2px,color:#05343a
    style C    fill:#fff3d1,stroke:#D9A400,stroke-width:2px,color:#4a3600
    style P    fill:#fff3d1,stroke:#D9A400,stroke-width:2px,color:#4a3600
    style T    fill:#fff3d1,stroke:#D9A400,stroke-width:2px,color:#4a3600
    style X    fill:#fff3d1,stroke:#D9A400,stroke-width:2px,color:#4a3600
    style CAND  fill:#e8eaeb,stroke:#63666A,stroke-width:2px,color:#2b2d2f
    style BAYES fill:#63666A,stroke:#3f4143,stroke-width:2px,color:#ffffff
    style PANEL fill:#BA0C2F,stroke:#7A0820,stroke-width:4px,color:#ffffff
    style DK  fill:#fadde3,stroke:#BA0C2F,stroke-width:3px,color:#5c0618
    style CCR fill:#fdeef1,stroke:#BA0C2F,stroke-width:2px,color:#5c0618
    style FRE fill:#fdeef1,stroke:#BA0C2F,stroke-width:2px,color:#5c0618

    linkStyle 4,5,6,7,8,9,10,11 stroke:#BA0C2F,stroke-width:2.5px
```

---

## Repository map

| Folder | Supplementary Methods | Contents |
|---|---|---|
| [`01_bulk_rnaseq/`](01_bulk_rnaseq/) | Bulk RNA-seq 1–5 | PyDESeq2 differential expression, GSEA, WGCNA, over-representation analysis |
| [`02_single_cell/`](02_single_cell/) | Single cell 1–3 | QC, three-way doublet consensus, scVI integration, annotation, pseudobulk DE |
| [`03_sccoda/`](03_sccoda/) | Single cell 4 | scCODA compositional analysis of cell-type abundance |
| [`04_cellphonedb/`](04_cellphonedb/) | Single cell 5a–5d | CellPhoneDB cell–cell communication, patient-level and pooled analyses |
| [`05_trajectory/`](05_trajectory/) | Single cell 6 | CytoTRACE2 potency and Monocle3 trajectory inference |
| [`06_spatial/`](06_spatial/) | Single cell 7a–7e | Giotto pipelines for Visium (55 µm) and Visium HD (16 µm) |
| [`07_prognostic_model/`](07_prognostic_model/) | Section 8 | Hierarchical Bayesian shrinkage panel, external validation, frozen release |
| [`config/`](config/) | — | Data root definition and per-sample spatial parameters |

Each folder has its own README describing every script and the methods element it implements.
`figures/` subfolders contain the notebooks and R scripts that **generate** the manuscript
panels — no image files are distributed here.

---

## Cohorts

| Cohort | n | Role | Used for fitting? |
|---|---|---|---|
| Internal bulk RNA-seq | 20 tumours, 10 sites | Discovery; model fitting | ✅ yes |
| Internal single-cell RNA-seq | 24 samples (13 good, 11 poor) | Cell-type discovery, candidate nomination, deconvolution reference | ✅ yes |
| Internal spatial | 4 Visium (55 µm), 2 Visium HD (16 µm) | Spatial architecture | — |
| **Danish bulk RNA-seq** | 54 | Confirmatory external survival validation — the single primary claim | ❌ never entered any fit |
| **MD Anderson CCR2020** | 54 | Pre-specified external validation of molecular *subtype*, not survival | ❌ |
| **Frerich 2018** | 66 | Cross-platform replication of the subtype result (~18 patients shared with CCR2020, so **not** independent of it) | ❌ |

Full detail, including endpoints, event counts and survival-time availability, is in
[`07_prognostic_model/release/cohorts.csv`](07_prognostic_model/release/cohorts.csv).

---

## Getting started

The scripts read and write large data files (h5ad objects, Giotto objects, CellPhoneDB outputs)
that are **not** distributed with this code. Point `ACC_DATA_ROOT` at a local copy of the data
tree:

```bash
export ACC_DATA_ROOT=/path/to/acc_data
```

[`config/paths.py`](config/paths.py) and [`config/paths.R`](config/paths.R) are the single
definition of that root for Python and R. Scripts resolve sibling files relative to their own
location, so the repository can live anywhere.

Sample identifiers throughout are **pseudonyms** (`P01`–`P25`) and do not correspond to any
clinical, pathology or accession number.

> [!NOTE]
> **Two column names differ from the originals.** `ModelOutcome1` and `ModelOutcome2` were
> renamed when this repository was prepared. If you are running these scripts against the
> original `.h5ad` or trait table, rename the corresponding columns to match, or the scripts
> will raise a `KeyError`.

> [!IMPORTANT]
> **Do not install everything into one environment.** The pipelines ran in several separate
> conda environments — scCODA and CellPhoneDB pin conflicting dependencies. Dependencies are
> listed in [`requirements.txt`](requirements.txt) (Python 3.11.15) and
> [`r_dependencies.txt`](r_dependencies.txt) (R 4.4.3);
> [`07_prognostic_model/release/environment.txt`](07_prognostic_model/release/environment.txt)
> records the exact environment behind every published model number.

### Where to start reading

| If you want… | Go to |
|---|---|
| The methods for one analysis | that section's `README.md` |
| The model itself | [`07_prognostic_model/release/model_spec.json`](07_prognostic_model/release/model_spec.json) — genes, signs, frozen reference, locked threshold |
| To score new samples | [`07_prognostic_model/release/score_acc4.py`](07_prognostic_model/release/score_acc4.py) — the one scoring function, cohort-batch and frozen single-sample |
| A specific number from the paper | [`07_prognostic_model/release/results_dictionary.csv`](07_prognostic_model/release/results_dictionary.csv) — every manuscript-facing value with its source file, generating script and evidence tier |
| To know what **not** to quote | [`07_prognostic_model/release/DEPRECATED.md`](07_prognostic_model/release/DEPRECATED.md) — superseded artefacts and the wrong numbers they would produce |

---

## Reproducibility

The prognostic model ships as a frozen release rather than a re-runnable pipeline, because the
published numbers came from one specific environment against one specific set of inputs:

- [`release/MANIFEST.json`](07_prognostic_model/release/MANIFEST.json) — sha256 fingerprints of
  every raw and harmonised input, including the external cohorts
- [`release/environment.txt`](07_prognostic_model/release/environment.txt) and
  [`release/requirements-frozen.txt`](07_prognostic_model/release/requirements-frozen.txt) — the
  interpreter, GPU and package versions actually used
- [`release/results_dictionary.csv`](07_prognostic_model/release/results_dictionary.csv) — each
  number tagged `confirmatory`, `prespecified`, `secondary`, `exploratory`, `audit` or
  `withdrawn`, so post-hoc analyses are never mistaken for pre-specified ones

`release/build_release.py` verifies as it builds: `score_acc4.py` must reproduce the stored
Danish score and the reported frozen-reference C-index, or the build fails.

---

## Citation

If you use this code, please cite the associated manuscript.

## Licence

[MIT](LICENSE).
