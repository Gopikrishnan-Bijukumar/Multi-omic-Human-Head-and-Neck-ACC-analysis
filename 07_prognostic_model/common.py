"""Shared helpers for the prognostic-model build scripts.

Config loading, provenance accumulation, and the normalisation primitives that
stages 05 and 06 must agree on exactly. Nothing here writes outside
`build_root`.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

CONFIG_PATH = Path(__file__).resolve().parent.parent / "config.yaml"


def load_config(path: Path | str = CONFIG_PATH) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["out_dir"] = Path(cfg["out_dir"])
    return cfg


def sha256(path: str | Path, max_bytes: int | None = None) -> str:
    """Hash a file. The h5ads are 8 GB, so they get a head-hash plus size."""
    h = hashlib.sha256()
    read = 0
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 22)
            if not chunk:
                break
            h.update(chunk)
            read += len(chunk)
            if max_bytes is not None and read >= max_bytes:
                break
    tag = "sha256" if max_bytes is None else f"sha256-first{max_bytes}B"
    return f"{tag}:{h.hexdigest()}:size={os.path.getsize(path)}"


def input_fingerprint(path: str | Path) -> dict:
    """Full hash for small files, head-hash + size for anything over 1 GB."""
    size = os.path.getsize(path)
    limit = None if size < (1 << 30) else (1 << 26)
    return {"path": str(path), "size_bytes": size, "hash": sha256(path, limit)}


def package_versions() -> dict:
    import anndata
    import numpy
    import pandas
    import scanpy
    import scipy

    vers = {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
        "scanpy": scanpy.__version__,
        "anndata": anndata.__version__,
    }
    try:
        import pydeseq2

        vers["pydeseq2"] = pydeseq2.__version__
    except Exception:
        vers["pydeseq2"] = None
    return vers


def write_provenance(
    cfg: dict, stage: str, payload: dict, out_dir: Path | str | None = None
) -> Path:
    """Merge this stage's record into <out_dir>/provenance.json.

    Defaults to the step-1 output directory; stage 05 passes `panel_dir` so the
    step-2 provenance lives beside the panel it describes.
    """
    prov_path = Path(out_dir or cfg["out_dir"]) / "provenance.json"
    prov_path.parent.mkdir(parents=True, exist_ok=True)
    prov = {}
    if prov_path.exists():
        with open(prov_path) as fh:
            prov = json.load(fh)
    prov.setdefault("config_path", str(CONFIG_PATH))
    prov.setdefault("stages", {})
    prov["stages"][stage] = {
        "finished_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "package_versions": package_versions(),
        **payload,
    }
    with open(prov_path, "w") as fh:
        json.dump(prov, fh, indent=2, default=str)
    return prov_path


# ---------------------------------------------------------------------------
# Normalisation primitives (step 3)
#
# Every matrix below is SAMPLES x GENES. These live here rather than in a stage
# script because stage 05's training-expression gate and stage 06's harmonised
# universe must apply the identical rule — if they could drift, the panel could
# contain a gene the universe excludes.
# ---------------------------------------------------------------------------


def size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq2 median-of-ratios size factors. Samples x genes -> one per sample.

    The reference is the geometric mean across samples, taken over genes with a
    non-zero count in EVERY sample (DESeq2's behaviour — a single zero makes the
    geometric mean zero and the ratio undefined). Fitted within one cohort; it
    is never estimated across a train/validation boundary.
    """
    usable = counts.loc[:, (counts > 0).all(axis=0)]
    assert usable.shape[1] >= 100, (
        f"only {usable.shape[1]} genes are non-zero in every sample — "
        "median-of-ratios has no stable reference"
    )
    log_counts = np.log(usable)
    log_ref = log_counts.mean(axis=0)
    sf = np.exp((log_counts - log_ref).median(axis=1))
    assert (sf > 0).all() and np.isfinite(sf).all(), "non-positive size factor"
    return sf.rename("size_factor")


def effective_size_factors(counts: pd.DataFrame,
                           sf: pd.Series | None = None) -> pd.Series:
    """Median-of-ratios factors rescaled to a per-million basis.

    Median-of-ratios corrects library COMPOSITION but leaves the result on the
    cohort's own count scale, and the two cohorts are sequenced to very
    different depths (DK arrived pre-equalised at ~2.36M per library). A
    threshold or a mean shift expressed on that scale would measure depth, not
    biology. Dividing by the cohort's mean library size in millions puts both
    cohorts on a CPM-like axis while keeping the composition correction:

        effective_sf_i = sf_i * mean_j(library_size_j) / 1e6

    When every size factor is 1 this reduces exactly to CPM.
    """
    if sf is None:
        sf = size_factors(counts)
    return (sf * counts.sum(axis=1).mean() / 1e6).rename("effective_size_factor")


def normalize_sf_log2(counts: pd.DataFrame, sf: pd.Series | None = None) -> pd.DataFrame:
    """Samples x genes raw counts -> log2(count / effective_size_factor + 1).

    This is the primary transform for both cohorts. It has no fitted parameters
    beyond the size factors, so the same function applied to each cohort gives a
    genuinely identical procedure rather than merely a similar one.
    """
    return np.log2(counts.div(effective_size_factors(counts, sf), axis=0) + 1.0)


def training_expressed(norm: pd.DataFrame, min_norm_log2: float,
                       min_frac_samples: float) -> pd.Series:
    """Low-expression filter, DEFINED FROM THE TRAINING DATA ONLY.

    `norm` is the normalised training matrix (samples x genes). A gene passes
    when it clears `min_norm_log2` in at least `min_frac_samples` of the training
    samples. No validation expression and no outcome label enters this decision.
    """
    frac = (norm >= min_norm_log2).mean(axis=0)
    return (frac >= min_frac_samples).rename("train_expressed")


def banner(text: str) -> None:
    print("\n" + "=" * 78)
    print(text)
    print("=" * 78, flush=True)
