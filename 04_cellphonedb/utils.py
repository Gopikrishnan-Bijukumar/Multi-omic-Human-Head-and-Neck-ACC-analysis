"""
Shared utilities for the CellPhoneDB v5 clinical-outcome re-analysis
(clinical_outcome_v2).

This module is the only place where provenance, data derivation, endpoint
extraction and inference logic lives. Notebooks orchestrate; they must not
re-implement anything found here.

Design rules honoured throughout
--------------------------------
* pathlib everywhere; no string path arithmetic.
* Every artefact-producing write is atomic (write to a temp file in the same
  directory, then os.replace).
* Provenance is JSON / JSONL and hashes are SHA-256.
* The full count matrix is never densified. Sparse validation touches only
  the ``.data`` buffer of a CSR/CSC matrix.
* Patient identifiers are never hardcoded. They are always derived from
  ``obs[sample_column]``.
* Structural missingness (a cell type that a patient simply does not have
  enough cells of) is *never* silently converted to a zero score.

Orientation conventions (read this before using the abundance matrix)
---------------------------------------------------------------------
``celltype_abundance_matrix`` returns a DataFrame whose **rows are cell types**
and whose **columns are samples/patients**::

    abundance.index.name   == cell_type_column   (e.g. 'annot_6')
    abundance.columns.name == sample_column      (e.g. 'batch')

``core_celltypes`` consumes that orientation. If you hand it a transposed
matrix it raises rather than silently returning nonsense.

Endpoint key convention
-----------------------
An *endpoint* is one (interaction, sender, receiver) triple. Its key is::

    f"{id_cp_interaction}{KEY_SEP}{sender}{DIR_SEP}{receiver}"
    # e.g. "CPI-SS0123456789::Fibroblast cells->Macrophages - M2"

``KEY_SEP`` is ``'::'`` and ``DIR_SEP`` is ``'->'``. Neither substring occurs in
any annot_6 label (the labels contain ``' - '`` and ``'/'`` only), and no label
contains ``'|'``, which is what makes the CellPhoneDB ``"<sender>|<receiver>"``
column convention unambiguous. Direction is always preserved: A->B and B->A are
different endpoints.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
import time
import warnings
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = [
    # provenance
    "load_config",
    "config_hash",
    "sha256_file",
    "write_manifest",
    "write_success_json",
    "run_is_complete",
    "atomic_publish",
    "log_jsonl",
    # data derivation
    "derive_patient_outcomes",
    "celltype_abundance_matrix",
    "pooled_retained_celltypes",
    "balance_downsample",
    "patient_retained_celltypes",
    "core_celltypes",
    "write_cpdb_input",
    "validate_counts_matrix",
    # endpoint extraction
    "parse_cellpair_columns",
    "extract_patient_endpoints",
    "classify_missingness",
    "make_endpoint_key",
    # inference
    "mannwhitney_hybrid",
    "rank_biserial",
    "bootstrap_ci",
    "permutation_p",
    "bh",
    "fisher_availability",
    # safety
    "safe_unlink",
    # constants
    "KEY_SEP",
    "DIR_SEP",
    "PAIR_SEP",
    "SUCCESS_FILENAME",
    "MISSINGNESS_STATUSES",
]

KEY_SEP = "::"
DIR_SEP = "->"
PAIR_SEP = "|"
SUCCESS_FILENAME = "_SUCCESS.json"
MANIFEST_FILENAME = "manifest.json"

MISSINGNESS_STATUSES = (
    "evaluated",
    "structural_missing",
    "unexpected_missing",
    "invalid",
)

_PATH_KEYS = ("project_root", "input_h5ad", "cpdb_database", "code_root", "output_root")
# output_root is created if absent; everything else must already exist.
_MUST_EXIST_KEYS = ("project_root", "input_h5ad", "cpdb_database", "code_root")

_SCHEMA_PATH = Path(__file__).resolve().parent / "config_schema.json"


# =============================================================================
# internal helpers
# =============================================================================

def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _json_default(obj: Any) -> Any:
    """Make numpy / pathlib / pandas scalars JSON-serialisable."""
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        value = float(obj)
        return value if np.isfinite(value) else str(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (pd.Timestamp, datetime)):
        return obj.isoformat()
    if isinstance(obj, (set, frozenset)):
        return sorted(obj, key=str)
    if isinstance(obj, pd.Index):
        return list(obj)
    if isinstance(obj, pd.Series):
        return obj.to_dict()
    return str(obj)


def _atomic_write_bytes(path: str | os.PathLike, payload: bytes) -> Path:
    """Write ``payload`` to ``path`` atomically (same-directory temp + replace)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise
    return path


def _atomic_write_json(path: str | os.PathLike, obj: Any) -> Path:
    text = json.dumps(obj, indent=2, sort_keys=True, default=_json_default)
    return _atomic_write_bytes(path, (text + "\n").encode("utf-8"))


def _canonical_json(obj: Any) -> str:
    """Deterministic JSON dump used for hashing."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=_json_default)


def _as_1d_float(values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=float).ravel()
    return arr


@contextmanager
def _quiet_all_nan():
    """Silence the expected 'All-NaN slice' noise from nan-aware reductions.

    An endpoint that is structurally missing for every patient in a group is a
    legitimate, informative outcome (it comes back as NaN and is reported as
    such); it is not a numerical accident worth warning about once per chunk.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN", category=RuntimeWarning)
        warnings.filterwarnings("ignore", message="Mean of empty slice", category=RuntimeWarning)
        with np.errstate(invalid="ignore", divide="ignore"):
            yield


# =============================================================================
# provenance
# =============================================================================

def load_config(path: str | os.PathLike) -> dict:
    """Load, validate and path-resolve the analysis configuration.

    Steps
    -----
    1. Parse the YAML file.
    2. Validate it against ``config_schema.json`` (jsonschema, draft 2020-12).
    3. Resolve every path key to an absolute, symlink-resolved string.
    4. Assert the input paths exist (h5ad, cpdb zip, project root, code root);
       create ``output_root`` if it does not exist yet.

    Returns
    -------
    dict
        Plain, JSON-serialisable dict (path values are ``str``, not ``Path``,
        so that :func:`config_hash` is stable). Wrap in ``Path(...)`` at use
        sites.
    """
    import yaml
    from jsonschema import Draft202012Validator

    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"config not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise TypeError(f"config must parse to a mapping, got {type(cfg).__name__}")

    with _SCHEMA_PATH.open("r", encoding="utf-8") as handle:
        schema = json.load(handle)

    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(cfg), key=lambda e: list(e.path))
    if errors:
        detail = "\n".join(
            f"  - {'/'.join(str(p) for p in err.path) or '<root>'}: {err.message}"
            for err in errors
        )
        raise ValueError(f"config failed schema validation ({path}):\n{detail}")

    for key in _PATH_KEYS:
        cfg[key] = str(Path(cfg[key]).expanduser().resolve())

    for key in _MUST_EXIST_KEYS:
        target = Path(cfg[key])
        if not target.exists():
            raise FileNotFoundError(f"config key '{key}' points at a missing path: {target}")

    Path(cfg["output_root"]).mkdir(parents=True, exist_ok=True)

    if cfg["sample_column"] == "batch_1" or "batch_1" in (cfg.get("covariates_available") or []):
        raise ValueError(
            "obs['batch_1'] is known to mislabel a patient and must never be used."
        )

    return cfg


def config_hash(cfg: Mapping[str, Any]) -> str:
    """SHA-256 of the canonical JSON dump of ``cfg`` (key order irrelevant)."""
    return hashlib.sha256(_canonical_json(dict(cfg)).encode("utf-8")).hexdigest()


def sha256_file(path: str | os.PathLike, chunk_size: int = 1 << 20) -> str:
    """Chunked SHA-256 of a file. Never loads the whole file into memory."""
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(paths: Iterable[str | os.PathLike], out_path: str | os.PathLike) -> None:
    """Write a JSON manifest describing ``paths`` (size, mtime, sha256).

    Directories are walked recursively. Missing paths are recorded with
    ``"exists": false`` rather than raising, so a manifest can always be
    written for diagnosis.
    """
    entries: list[dict] = []
    for raw in paths:
        target = Path(raw)
        if target.is_dir():
            members = sorted(p for p in target.rglob("*") if p.is_file())
        else:
            members = [target]
        for member in members:
            if not member.exists():
                entries.append({"path": str(member), "exists": False})
                continue
            stat = member.stat()
            entries.append(
                {
                    "path": str(member),
                    "exists": True,
                    "bytes": stat.st_size,
                    "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                    "sha256": sha256_file(member),
                }
            )
    _atomic_write_json(
        out_path,
        {"created_utc": _utcnow(), "n_entries": len(entries), "entries": entries},
    )


def write_success_json(run_dir: str | os.PathLike, **fields: Any) -> None:
    """Write ``_SUCCESS.json`` into ``run_dir``. Must be the LAST write of a run.

    The presence of this file is the only signal that a run finished. Callers
    should pass at least ``config_hash=...`` so that :func:`run_is_complete`
    can reject a directory produced under a different configuration.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    record = {"completed_utc": _utcnow(), "pid": os.getpid()}
    record.update(fields)
    _atomic_write_json(run_dir / SUCCESS_FILENAME, record)


def run_is_complete(
    run_dir: str | os.PathLike,
    required_files: Sequence[str | os.PathLike],
    cfg_hash: str,
) -> bool:
    """True iff ``run_dir`` holds a finished run made under ``cfg_hash``.

    All three must hold:
      * ``_SUCCESS.json`` exists and parses;
      * its ``config_hash`` equals ``cfg_hash``;
      * every entry of ``required_files`` exists (relative paths are resolved
        against ``run_dir``).
    """
    run_dir = Path(run_dir)
    success = run_dir / SUCCESS_FILENAME
    if not success.is_file():
        return False
    try:
        with success.open("r", encoding="utf-8") as handle:
            record = json.load(handle)
    except (json.JSONDecodeError, OSError):
        return False
    if record.get("config_hash") != cfg_hash:
        return False
    for required in required_files:
        candidate = Path(required)
        if not candidate.is_absolute():
            candidate = run_dir / candidate
        if not candidate.exists():
            return False
    return True


def atomic_publish(tmp_dir: str | os.PathLike, final_dir: str | os.PathLike) -> None:
    """Atomically move ``tmp_dir`` into place as ``final_dir``.

    A pre-existing ``final_dir`` is first renamed aside and only deleted once
    the new directory is in place, so a crash can never leave the publish point
    empty. ``tmp_dir`` must live on the same filesystem as ``final_dir``'s
    parent (put it next to the destination).
    """
    tmp_dir = Path(tmp_dir)
    final_dir = Path(final_dir)
    if not tmp_dir.is_dir():
        raise NotADirectoryError(f"tmp_dir is not a directory: {tmp_dir}")
    if tmp_dir.resolve() == final_dir.resolve():
        raise ValueError("tmp_dir and final_dir must differ")

    final_dir.parent.mkdir(parents=True, exist_ok=True)
    displaced: Path | None = None
    if final_dir.exists():
        displaced = final_dir.with_name(f"{final_dir.name}.superseded.{int(time.time() * 1000)}")
        os.replace(final_dir, displaced)
    try:
        os.replace(tmp_dir, final_dir)
    except BaseException:
        if displaced is not None:
            os.replace(displaced, final_dir)
        raise
    if displaced is not None:
        shutil.rmtree(displaced, ignore_errors=True)


def log_jsonl(path: str | os.PathLike, record: Mapping[str, Any]) -> None:
    """Append one JSON record as a line, flushing and fsyncing immediately.

    Flushing on every write is deliberate: these logs are the forensic record
    of long CellPhoneDB runs and must survive a kill -9.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(record)
    payload.setdefault("ts_utc", _utcnow())
    line = json.dumps(payload, sort_keys=True, default=_json_default)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
        handle.flush()
        os.fsync(handle.fileno())


# =============================================================================
# data derivation
# =============================================================================

def derive_patient_outcomes(
    obs: pd.DataFrame,
    sample_col: str,
    outcome_col: str,
) -> pd.Series:
    """Map each patient to its single clinical outcome.

    Raises loudly if any patient carries more than one distinct outcome label,
    which would mean the sample/outcome join is corrupt (this is exactly the
    failure mode of obs['batch_1']).

    Returns
    -------
    pd.Series
        Index = patient id (sorted), values = outcome label, name = outcome_col.
    """
    for col in (sample_col, outcome_col):
        if col not in obs.columns:
            raise KeyError(f"column '{col}' not present in obs")

    frame = pd.DataFrame(
        {
            "sample": obs[sample_col].astype(str).to_numpy(),
            "outcome": obs[outcome_col].astype(str).to_numpy(),
        }
    )
    if frame["outcome"].isin({"nan", "None", ""}).any():
        bad = sorted(frame.loc[frame["outcome"].isin({"nan", "None", ""}), "sample"].unique())
        raise ValueError(f"missing {outcome_col} for samples: {bad}")

    counts = frame.drop_duplicates().groupby("sample")["outcome"].nunique()
    conflicted = counts[counts > 1]
    if len(conflicted):
        detail = {
            str(sample): sorted(frame.loc[frame["sample"] == sample, "outcome"].unique())
            for sample in conflicted.index
        }
        raise ValueError(
            f"{len(conflicted)} sample(s) carry more than one {outcome_col} value: {detail}"
        )

    series = frame.drop_duplicates().set_index("sample")["outcome"].sort_index()
    series.index.name = sample_col
    series.name = outcome_col
    return series


def celltype_abundance_matrix(
    obs: pd.DataFrame,
    sample_col: str,
    celltype_col: str,
) -> pd.DataFrame:
    """Cell counts per (cell type, sample).

    Orientation: **rows = cell types, columns = samples** (see module
    docstring). Unobserved combinations are 0, and every observed cell type /
    sample appears even if its count is 0 in some cells of the grid.
    """
    for col in (sample_col, celltype_col):
        if col not in obs.columns:
            raise KeyError(f"column '{col}' not present in obs")

    matrix = pd.crosstab(
        obs[celltype_col].astype(str),
        obs[sample_col].astype(str),
    )
    matrix = matrix.sort_index().sort_index(axis=1).astype(int)
    matrix.index.name = celltype_col
    matrix.columns.name = sample_col
    return matrix


def pooled_retained_celltypes(
    obs: pd.DataFrame,
    celltype_col: str,
    outcome_col: str,
    min_per_group: int,
) -> list[str]:
    """Cell types with at least ``min_per_group`` cells in BOTH outcome groups.

    Used for the pooled Good-vs-Poor arm, where a cell type must be present on
    both sides for a comparison to mean anything.
    """
    for col in (celltype_col, outcome_col):
        if col not in obs.columns:
            raise KeyError(f"column '{col}' not present in obs")

    table = pd.crosstab(obs[celltype_col].astype(str), obs[outcome_col].astype(str))
    if table.shape[1] < 2:
        raise ValueError(
            f"expected two outcome groups in '{outcome_col}', found {list(table.columns)}"
        )
    keep = table.min(axis=1) >= int(min_per_group)
    return sorted(table.index[keep].astype(str).tolist())


def balance_downsample(
    obs: pd.DataFrame,
    celltype_col: str,
    outcome_col: str,
    keep_types: Sequence[str],
    seed: int,
) -> pd.Index:
    """Per cell type, downsample the larger outcome group to the smaller one.

    Returns
    -------
    pd.Index
        Barcodes to keep, ordered as in ``obs``. Cell types outside
        ``keep_types`` are dropped entirely.

    Determinism: candidates are sorted before sampling and the generator is
    ``np.random.default_rng(seed)``, so the result depends only on ``seed`` and
    the input barcodes, not on row order or dict iteration.
    """
    for col in (celltype_col, outcome_col):
        if col not in obs.columns:
            raise KeyError(f"column '{col}' not present in obs")

    rng = np.random.default_rng(int(seed))
    keep_types = list(dict.fromkeys(str(t) for t in keep_types))
    celltypes = obs[celltype_col].astype(str)
    outcomes = obs[outcome_col].astype(str)

    selected: list[np.ndarray] = []
    for cell_type in keep_types:
        type_mask = (celltypes == cell_type).to_numpy()
        if not type_mask.any():
            continue
        groups = sorted(outcomes[type_mask].unique().tolist())
        per_group = {
            group: np.sort(obs.index[type_mask & (outcomes == group).to_numpy()].to_numpy().astype(str))
            for group in groups
        }
        target = min(len(v) for v in per_group.values())
        if target == 0:
            continue
        for group in groups:
            pool = per_group[group]
            if len(pool) == target:
                selected.append(pool)
            else:
                picked = rng.choice(len(pool), size=target, replace=False)
                selected.append(pool[np.sort(picked)])

    if not selected:
        return pd.Index([], name=obs.index.name)

    chosen = set(np.concatenate(selected).tolist())
    mask = obs.index.astype(str).isin(chosen)
    return obs.index[mask]


def patient_retained_celltypes(
    obs: pd.DataFrame,
    patient: str | None,
    min_cells: int,
    *,
    sample_col: str = "batch",
    celltype_col: str = "annot_6",
) -> list[str]:
    """Cell types this patient has at least ``min_cells`` of.

    ``sample_col`` / ``celltype_col`` are keyword-only with the project
    defaults so that the pinned positional signature
    ``(obs, patient, min_cells)`` still holds; pass them explicitly from
    ``cfg`` in notebooks.

    ``patient=None`` means ``obs`` is already subset to a single patient.
    """
    if celltype_col not in obs.columns:
        raise KeyError(f"column '{celltype_col}' not present in obs")

    if patient is None:
        subset = obs
    else:
        if sample_col not in obs.columns:
            raise KeyError(f"column '{sample_col}' not present in obs")
        subset = obs.loc[obs[sample_col].astype(str) == str(patient)]
        if subset.empty:
            raise ValueError(f"patient '{patient}' has no cells in obs[{sample_col!r}]")

    counts = subset[celltype_col].astype(str).value_counts()
    return sorted(counts.index[counts >= int(min_cells)].astype(str).tolist())


def core_celltypes(
    abundance: pd.DataFrame,
    min_cells: int,
    min_patients: int,
) -> list[str]:
    """Cell types present at ``>= min_cells`` in at least ``min_patients`` patients.

    ``abundance`` must be the matrix produced by
    :func:`celltype_abundance_matrix`: **rows = cell types, columns = patients**.
    Verified: at ``min_cells=10, min_patients=22`` over the 24-patient cohort
    this returns exactly six types (Fibroblast cells, Epithelial cells - Tumor,
    ADC - Tumor, Macrophages - M2, ADMEC - Tumor, Endothelial cells).
    """
    if not isinstance(abundance, pd.DataFrame):
        raise TypeError("abundance must be a DataFrame of cell types x patients")
    if min_patients > abundance.shape[1]:
        raise ValueError(
            f"min_patients={min_patients} exceeds the {abundance.shape[1]} columns of the "
            "abundance matrix -- did you pass it transposed? Rows must be cell types."
        )
    n_qualifying = (abundance.to_numpy() >= int(min_cells)).sum(axis=1)
    keep = n_qualifying >= int(min_patients)
    return sorted(abundance.index[keep].astype(str).tolist())


def validate_counts_matrix(matrix: Any, *, name: str = "counts") -> dict:
    """Assert a matrix holds raw integral counts WITHOUT densifying it.

    For a sparse matrix only the ``.data`` buffer (the stored non-zeros) is
    inspected -- never ``.toarray()``. Returns a small summary dict suitable
    for a provenance record.
    """
    import scipy.sparse as sp

    if sp.issparse(matrix):
        data = matrix.data
        n_stored = int(matrix.nnz)
        sparse_format = matrix.format
    else:
        data = np.asarray(matrix).ravel()
        n_stored = int(data.size)
        sparse_format = "dense"

    if data.size:
        finite = np.isfinite(data)
        if not finite.all():
            raise ValueError(f"{name}: matrix contains non-finite values")
        if float(data.min()) < 0:
            raise ValueError(f"{name}: matrix contains negative values -- not raw counts")
        if not np.array_equal(data, np.round(data)):
            raise ValueError(
                f"{name}: matrix values are not integral -- this looks normalised/z-scored, "
                "not raw counts"
            )
        data_max = float(data.max())
        data_sum = float(data.sum())
    else:
        data_max = 0.0
        data_sum = 0.0

    return {
        "name": name,
        "shape": tuple(int(x) for x in matrix.shape),
        "format": sparse_format,
        # NB: read dtype off the object; a getattr default would be evaluated
        # eagerly and densify the matrix, which is exactly what we forbid.
        "dtype": str(matrix.dtype) if hasattr(matrix, "dtype") else str(data.dtype),
        "n_stored": n_stored,
        "max": data_max,
        "sum": data_sum,
        "densified": False,
    }


def write_cpdb_input(
    adata: Any,
    cells: Sequence[str] | pd.Index,
    out_dir: str | os.PathLike,
    celltype_col: str,
    *,
    counts_layer: str = "counts",
    validate: bool = True,
) -> dict:
    """Write the two files CellPhoneDB v5 needs: ``counts.h5ad`` + ``meta.tsv``.

    * ``counts.h5ad`` carries the RAW counts layer as ``X`` (sparse preserved;
      ``adata.X`` itself is z-scored and unusable, so it is never written).
    * ``meta.tsv`` is tab separated with an index named ``Cell`` and exactly one
      metadata column named ``cell_type`` -- CellPhoneDB requires that name.
    * Row order of the two files matches barcode-for-barcode; this is asserted
      after writing by re-reading the h5ad obs names.

    Returns
    -------
    dict
        Paths, sha256 of both files, cell/gene counts, the cell-type tally and
        the sparse-safe counts validation summary.
    """
    import anndata as ad
    import scipy.sparse as sp

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cells = pd.Index([str(c) for c in pd.Index(cells)])
    if cells.has_duplicates:
        dupes = cells[cells.duplicated()].unique().tolist()[:5]
        raise ValueError(f"duplicate barcodes requested for CPDB input: {dupes} ...")
    if len(cells) == 0:
        raise ValueError("refusing to write an empty CellPhoneDB input")

    obs_names = pd.Index(adata.obs_names.astype(str))
    missing = cells.difference(obs_names)
    if len(missing):
        raise KeyError(f"{len(missing)} requested barcodes are absent from adata, e.g. {list(missing[:5])}")
    if celltype_col not in adata.obs.columns:
        raise KeyError(f"column '{celltype_col}' not present in adata.obs")

    positions = obs_names.get_indexer(cells)
    subset = adata[positions]

    if counts_layer not in subset.layers:
        raise KeyError(
            f"layer '{counts_layer}' not found (available: {list(subset.layers.keys())}). "
            "adata.X is z-scored and must never be used as counts."
        )
    counts = subset.layers[counts_layer]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    else:
        counts = counts.tocsr()
    counts_summary = validate_counts_matrix(counts, name=counts_layer) if validate else {}

    labels = pd.Series(
        np.asarray(subset.obs[celltype_col].astype(str)),
        index=cells,
        name="cell_type",
    )
    if labels.isna().any() or labels.isin({"nan", "None", ""}).any():
        raise ValueError(f"missing '{celltype_col}' labels among the selected cells")

    counts_adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame(index=cells.copy()),
        var=pd.DataFrame(index=pd.Index([str(g) for g in subset.var_names])),
    )
    counts_adata.obs_names_make_unique()
    counts_adata.obs[celltype_col] = labels.to_numpy()

    counts_path = out_dir / "counts.h5ad"
    meta_path = out_dir / "meta.tsv"

    tmp_counts = counts_path.with_suffix(".h5ad.tmp")
    counts_adata.write_h5ad(tmp_counts)
    os.replace(tmp_counts, counts_path)

    meta = pd.DataFrame({"cell_type": labels.to_numpy()}, index=cells.copy())
    meta.index.name = "Cell"
    _atomic_write_bytes(meta_path, meta.to_csv(sep="\t").encode("utf-8"))

    written = ad.read_h5ad(counts_path, backed="r")
    try:
        written_names = pd.Index(written.obs_names.astype(str))
    finally:
        if getattr(written, "file", None) is not None:
            written.file.close()
    if not written_names.equals(cells):
        raise AssertionError("counts.h5ad barcode order does not match the requested cell order")

    return {
        "counts_h5ad": str(counts_path),
        "meta_tsv": str(meta_path),
        "counts_h5ad_sha256": sha256_file(counts_path),
        "meta_tsv_sha256": sha256_file(meta_path),
        "n_cells": int(len(cells)),
        "n_genes": int(counts_adata.n_vars),
        "celltype_column": celltype_col,
        "cell_type_counts": labels.value_counts().sort_index().to_dict(),
        "counts_validation": counts_summary,
    }


# =============================================================================
# endpoint extraction
# =============================================================================

def make_endpoint_key(interaction: str, sender: str, receiver: str) -> str:
    """Directional endpoint key: ``"<interaction>::<sender>-><receiver>"``."""
    return f"{interaction}{KEY_SEP}{sender}{DIR_SEP}{receiver}"


def parse_cellpair_columns(
    cols: Iterable[str],
    known_celltypes: Iterable[str],
) -> list[tuple[str, str]]:
    """Parse CellPhoneDB ``"<sender>|<receiver>"`` columns into ordered pairs.

    Columns without ``'|'`` are metadata (id_cp_interaction, interacting_pair,
    classification, ...) and are skipped. Any column that *does* contain ``'|'``
    must split into exactly two halves that are both in ``known_celltypes``;
    otherwise a ValueError is raised. Silence here would silently drop or
    mis-assign endpoints, so failure is loud by design.

    Direction is preserved: ``(sender, receiver)`` is returned in column order
    and ``A|B`` is never conflated with ``B|A``.
    """
    known = set(str(t) for t in known_celltypes)
    if not known:
        raise ValueError("known_celltypes is empty; cannot validate cell-pair columns")

    pairs: list[tuple[str, str]] = []
    problems: list[str] = []
    for col in cols:
        col = str(col)
        if PAIR_SEP not in col:
            continue
        parts = col.split(PAIR_SEP)
        if len(parts) != 2:
            problems.append(f"{col!r}: split into {len(parts)} parts, expected 2")
            continue
        sender, receiver = parts[0], parts[1]
        unknown = [p for p in (sender, receiver) if p not in known]
        if unknown:
            problems.append(f"{col!r}: unknown cell type(s) {unknown}")
            continue
        pairs.append((sender, receiver))

    if problems:
        raise ValueError(
            "unparseable CellPhoneDB cell-pair column(s):\n  " + "\n  ".join(problems)
        )
    return pairs


def extract_patient_endpoints(
    scores_csv: str | os.PathLike,
    retained: Sequence[str],
    known_celltypes: Sequence[str],
    *,
    patient: str | None = None,
    restrict_to: Sequence[str] | None = None,
    id_col: str = "id_cp_interaction",
) -> pd.DataFrame:
    """Melt one patient's ``interaction_scores.csv`` into long endpoint rows.

    The returned frame is the union of
      * every cell-pair column actually present in the CSV, and
      * the full expected grid ``retained x retained`` (optionally intersected
        with ``restrict_to``, e.g. the six core cell types),
    so that endpoints the CSV *should* have had but does not are visible as
    rows rather than vanishing.

    Scores are NOT filtered or zero-filled here. Classification is the job of
    :func:`classify_missingness`; this function only records observable facts.

    Columns
    -------
    patient, id_cp_interaction, interacting_pair, classification,
    sender, receiver, cell_pair, endpoint_key, score_raw, score,
    score_present (was there a column for this pair?),
    sender_retained, receiver_retained, non_numeric, duplicate_id.
    """
    scores_csv = Path(scores_csv)
    if patient is None:
        patient = scores_csv.parent.name

    df = pd.read_csv(scores_csv)
    if id_col not in df.columns:
        raise KeyError(f"{scores_csv}: expected an '{id_col}' column, found {list(df.columns)[:8]}")

    known = [str(t) for t in known_celltypes]
    retained_set = {str(t) for t in retained}
    unknown_retained = retained_set - set(known)
    if unknown_retained:
        raise ValueError(f"retained contains cell types absent from known_celltypes: {sorted(unknown_retained)}")

    present_pairs = parse_cellpair_columns(df.columns, known)
    pair_to_col = {pair: f"{pair[0]}{PAIR_SEP}{pair[1]}" for pair in present_pairs}

    meta_cols = [c for c in df.columns if PAIR_SEP not in str(c)]
    keep_meta = [c for c in (id_col, "interacting_pair", "classification") if c in meta_cols]

    df = df.copy()
    df[id_col] = df[id_col].astype(str)
    dup_ids = set(df.loc[df[id_col].duplicated(keep=False), id_col].tolist())

    grid_types = sorted(retained_set)
    if restrict_to is not None:
        allowed = {str(t) for t in restrict_to}
        grid_types = sorted(retained_set & allowed)
    grid_pairs = [(s, r) for s in grid_types for r in grid_types]

    all_pairs: list[tuple[str, str]] = list(dict.fromkeys(list(present_pairs) + grid_pairs))
    if not all_pairs:
        raise ValueError(f"{scores_csv}: no cell-pair columns and an empty expected grid")

    blocks: list[pd.DataFrame] = []
    base = df[keep_meta].reset_index(drop=True)
    for sender, receiver in all_pairs:
        block = base.copy()
        block["sender"] = sender
        block["receiver"] = receiver
        col = pair_to_col.get((sender, receiver))
        if col is None:
            block["score_raw"] = np.nan
            block["score_present"] = False
        else:
            block["score_raw"] = df[col].to_numpy()
            block["score_present"] = True
        blocks.append(block)

    long_df = pd.concat(blocks, ignore_index=True)
    long_df.insert(0, "patient", str(patient))
    long_df["cell_pair"] = long_df["sender"] + PAIR_SEP + long_df["receiver"]
    long_df["endpoint_key"] = [
        make_endpoint_key(i, s, r)
        for i, s, r in zip(long_df[id_col], long_df["sender"], long_df["receiver"])
    ]
    long_df["score"] = pd.to_numeric(long_df["score_raw"], errors="coerce")
    long_df["non_numeric"] = (
        long_df["score_present"]
        & long_df["score_raw"].notna()
        & long_df["score"].isna()
    )
    long_df["sender_retained"] = long_df["sender"].isin(retained_set)
    long_df["receiver_retained"] = long_df["receiver"].isin(retained_set)
    long_df["duplicate_id"] = long_df[id_col].isin(dup_ids)

    for optional in ("interacting_pair", "classification"):
        if optional not in long_df.columns:
            long_df[optional] = pd.NA

    ordered = [
        "patient", id_col, "interacting_pair", "classification",
        "sender", "receiver", "cell_pair", "endpoint_key",
        "score_raw", "score", "score_present",
        "sender_retained", "receiver_retained", "non_numeric", "duplicate_id",
    ]
    return long_df[ordered]


def classify_missingness(
    long_df: pd.DataFrame,
    retained_by_patient: Mapping[str, Sequence[str]] | Sequence[str],
) -> pd.DataFrame:
    """Assign every endpoint row exactly one of the four missingness statuses.

    Statuses (mutually exclusive, exhaustive)
    -----------------------------------------
    ``evaluated``           both cell types passed the patient's abundance
                            threshold AND a finite numeric score exists.
    ``structural_missing``  sender or receiver failed the abundance threshold
                            for this patient. The endpoint is *unobservable*,
                            not zero.
    ``unexpected_missing``  both types retained but the row/score is absent
                            (missing column, or NaN in a present column).
    ``invalid``             duplicate endpoint key, non-numeric score, or a
                            present-but-non-finite (+/-inf) score.

    A structurally missing endpoint is NEVER converted to 0. A genuinely
    computed score of 0.0 stays 0.0 -- the two are distinguished by
    ``score_evaluated`` (NaN vs 0.0) plus ``status``.

    ``retained_by_patient`` may be a dict patient -> retained cell types (the
    frame must then carry a ``patient`` column) or a flat sequence applied to
    every row.
    """
    required = {"sender", "receiver", "endpoint_key"}
    missing_cols = required - set(long_df.columns)
    if missing_cols:
        raise KeyError(f"long_df is missing required column(s): {sorted(missing_cols)}")

    out = long_df.copy()
    if "score" not in out.columns:
        out["score"] = pd.to_numeric(out.get("score_raw"), errors="coerce")
    if "score_present" not in out.columns:
        out["score_present"] = out["score"].notna()
    if "non_numeric" not in out.columns:
        out["non_numeric"] = False

    if isinstance(retained_by_patient, Mapping):
        if "patient" not in out.columns:
            raise KeyError("retained_by_patient is a mapping, so long_df needs a 'patient' column")
        retained_sets = {str(k): {str(t) for t in v} for k, v in retained_by_patient.items()}
        unknown_patients = set(out["patient"].astype(str)) - set(retained_sets)
        if unknown_patients:
            raise KeyError(
                f"no retained cell-type list for patient(s): {sorted(unknown_patients)}"
            )
        patients = out["patient"].astype(str)
        sender_ok = np.fromiter(
            (s in retained_sets[p] for p, s in zip(patients, out["sender"].astype(str))),
            dtype=bool,
            count=len(out),
        )
        receiver_ok = np.fromiter(
            (r in retained_sets[p] for p, r in zip(patients, out["receiver"].astype(str))),
            dtype=bool,
            count=len(out),
        )
        key_frame = pd.Series(patients.to_numpy() + KEY_SEP + out["endpoint_key"].astype(str).to_numpy())
    else:
        flat = {str(t) for t in retained_by_patient}
        sender_ok = out["sender"].astype(str).isin(flat).to_numpy()
        receiver_ok = out["receiver"].astype(str).isin(flat).to_numpy()
        key_frame = out["endpoint_key"].astype(str)

    out["sender_retained"] = sender_ok
    out["receiver_retained"] = receiver_ok

    score = pd.to_numeric(out["score"], errors="coerce").to_numpy(dtype=float)
    present = out["score_present"].to_numpy(dtype=bool)
    non_numeric = out["non_numeric"].to_numpy(dtype=bool)
    duplicate_key = key_frame.duplicated(keep=False).to_numpy()
    if "duplicate_id" in out.columns:
        duplicate_key = duplicate_key | out["duplicate_id"].to_numpy(dtype=bool)

    raw = out["score_raw"] if "score_raw" in out.columns else out["score"]
    raw_numeric = pd.to_numeric(raw, errors="coerce").to_numpy(dtype=float)
    non_finite_present = present & ~np.isnan(raw_numeric) & ~np.isfinite(raw_numeric)

    both_retained = sender_ok & receiver_ok
    has_score = present & np.isfinite(score)

    status = np.full(len(out), "unexpected_missing", dtype=object)
    status[both_retained & has_score] = "evaluated"
    status[~both_retained] = "structural_missing"
    invalid = duplicate_key | non_numeric | non_finite_present
    status[invalid] = "invalid"

    out["status"] = pd.Categorical(status, categories=list(MISSINGNESS_STATUSES))
    if out["status"].isna().any():
        raise AssertionError("classify_missingness produced an unknown status")

    # score_evaluated is the ONLY column downstream inference may consume.
    # Structural rows are NaN here forever; a real 0.0 survives as 0.0.
    evaluated_mask = (out["status"] == "evaluated").to_numpy()
    out["score_evaluated"] = np.where(evaluated_mask, score, np.nan)
    out["tested"] = evaluated_mask
    return out


# =============================================================================
# inference
# =============================================================================

def mannwhitney_hybrid(good: Sequence[float], poor: Sequence[float]) -> tuple[float, float, str]:
    """Two-sided Mann-Whitney U with an exact/asymptotic hybrid null.

    Rule
    ----
    * No ties in ``concat(good, poor)``  -> ``method='exact'``
    * Any tie present                    -> ``method='asymptotic'``

    ``good`` is ALWAYS the first argument to ``scipy.stats.mannwhitneyu``, so
    the returned U is the Good-group statistic and :func:`rank_biserial` gets
    its sign right (positive = larger in Good).

    Why the hybrid, and not asymptotic everywhere
    ---------------------------------------------
    At 13 Good vs 11 Poor the smallest attainable p-value is
    ``3.90e-05`` asymptotically but ``8.01e-07`` exactly. Under
    Benjamini-Hochberg across a family of ~20,000 endpoints the smallest
    achievable q is ``p * m / rank``; with ``m ~ 2e4`` the best asymptotic
    p-value cannot reach q < 0.05 for any endpoint, so an asymptotic-always
    policy would make the primary analysis return zero discoveries *by
    arithmetic*, regardless of the biology. The exact null is the only one with
    the resolution this design can spend.

    Asymptotic remains the correct fallback under ties because scipy's
    ``method='exact'`` is invalid (and refuses to correct) when ties exist.

    Returns
    -------
    (U, p, method)
        ``method`` is one of ``'exact'``, ``'asymptotic'``, ``'insufficient'``.
        For ``'insufficient'`` (an empty group after dropping non-finite
        values) U and p are NaN.
    """
    from scipy.stats import mannwhitneyu

    good_arr = _as_1d_float(good)
    poor_arr = _as_1d_float(poor)
    good_arr = good_arr[np.isfinite(good_arr)]
    poor_arr = poor_arr[np.isfinite(poor_arr)]

    if good_arr.size == 0 or poor_arr.size == 0:
        return (float("nan"), float("nan"), "insufficient")

    combined = np.concatenate([good_arr, poor_arr])
    has_ties = np.unique(combined).size < combined.size
    method = "asymptotic" if has_ties else "exact"

    result = mannwhitneyu(good_arr, poor_arr, alternative="two-sided", method=method)
    return (float(result.statistic), float(result.pvalue), method)


def rank_biserial(U: float, n_good: int, n_poor: int) -> float:
    """Rank-biserial correlation from the GOOD-group U statistic.

    ``r = 2U/(n_good * n_poor) - 1`` in [-1, 1]; positive means the endpoint is
    larger in Good. Requires that ``U`` came from ``mannwhitneyu(good, poor)``.
    """
    n_good = int(n_good)
    n_poor = int(n_poor)
    if n_good <= 0 or n_poor <= 0:
        return float("nan")
    if U is None or not np.isfinite(U):
        return float("nan")
    return float(2.0 * float(U) / (n_good * n_poor) - 1.0)


def _u_and_pairs(good: np.ndarray, poor: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """NaN-aware vectorised U statistic.

    Parameters
    ----------
    good : (..., n_good)
    poor : (..., n_poor)

    Returns
    -------
    (U, n_pairs)
        ``U`` counts ``good > poor`` plus half the ties over pairs where both
        values are finite; ``n_pairs`` is the number of such pairs. Both have
        the leading (broadcast) shape.
    """
    good_valid = np.isfinite(good)
    poor_valid = np.isfinite(poor)
    n_pairs = good_valid.sum(axis=-1).astype(np.float64) * poor_valid.sum(axis=-1).astype(np.float64)

    u_total = np.zeros(good.shape[:-1], dtype=np.float64)
    # Loop over the SHORTER axis only (<= 13 here); everything else stays
    # vectorised, so the peak allocation is (chunk, n_resamples, n_poor).
    # No NaN masking is needed inside the comparisons: NaN > x and NaN == x are
    # both False, so invalid pairs contribute nothing to U on their own. Keeping
    # the comparisons boolean (rather than np.where(..., 1.0, 0.0)) avoids two
    # float64 temporaries per iteration, which dominates the runtime at scale.
    for i in range(good.shape[-1]):
        g_i = good[..., i][..., None]
        u_total += (g_i > poor).sum(axis=-1)
        u_total += 0.5 * (g_i == poor).sum(axis=-1)
    return u_total, n_pairs


def bootstrap_ci(
    value_matrix: np.ndarray,
    good_idx: Sequence[int],
    poor_idx: Sequence[int],
    n_resamples: int,
    seed: int,
    *,
    chunk_target_elements: int = 20_000_000,
    ci: float = 95.0,
) -> dict:
    """Patient-level bootstrap CIs for median difference and rank-biserial.

    Parameters
    ----------
    value_matrix : (n_endpoints, n_patients) float array
        One row per endpoint, one column per patient. NaN marks a patient for
        whom the endpoint is structurally unobservable; NaNs are carried
        through the resample (they are *not* dropped up-front, because dropping
        would change which patient is being resampled).
    good_idx, poor_idx : column indices of the two outcome groups.
    n_resamples : bootstrap replicates.
    seed : feeds ``np.random.default_rng``.

    Method
    ------
    PATIENTS are resampled with replacement **within** outcome group; cells are
    never resampled and a patient never crosses groups. The
    ``(n_resamples, n_good)`` and ``(n_resamples, n_poor)`` index arrays are
    drawn ONCE and reused for every endpoint, so all endpoints share the same
    bootstrap world (this also makes the whole thing a handful of vectorised
    array ops). Endpoints are processed in chunks sized to keep peak memory
    near ``chunk_target_elements`` doubles -- a per-endpoint Python loop over
    1e4 endpoints x 1e4 resamples would never finish.

    Returns
    -------
    dict with per-endpoint arrays of length n_endpoints:
        ``median_difference``, ``median_difference_ci_low/high``,
        ``rank_biserial``, ``rank_biserial_ci_low/high``,
        ``n_good_evaluable``, ``n_poor_evaluable``,
        plus scalars ``n_resamples``, ``seed``, ``ci``, ``chunk_size``,
        ``good_resample_shape``, ``poor_resample_shape``.
    """
    values = np.asarray(value_matrix, dtype=float)
    if values.ndim == 1:
        values = values[None, :]
    if values.ndim != 2:
        raise ValueError("value_matrix must be 2-D (n_endpoints x n_patients)")

    good_idx = np.asarray(good_idx, dtype=int)
    poor_idx = np.asarray(poor_idx, dtype=int)
    if good_idx.size == 0 or poor_idx.size == 0:
        raise ValueError("both outcome groups need at least one patient")
    if set(good_idx.tolist()) & set(poor_idx.tolist()):
        raise ValueError("a patient appears in both outcome groups")
    n_patients = values.shape[1]
    if good_idx.max(initial=0) >= n_patients or poor_idx.max(initial=0) >= n_patients:
        raise IndexError("patient index out of range for value_matrix")

    n_endpoints = values.shape[0]
    n_good = int(good_idx.size)
    n_poor = int(poor_idx.size)
    n_resamples = int(n_resamples)

    rng = np.random.default_rng(int(seed))
    good_boot = rng.integers(0, n_good, size=(n_resamples, n_good))
    poor_boot = rng.integers(0, n_poor, size=(n_resamples, n_poor))

    good_values = values[:, good_idx]
    poor_values = values[:, poor_idx]

    # observed point estimates
    with _quiet_all_nan():
        obs_median_diff = np.nanmedian(good_values, axis=1) - np.nanmedian(poor_values, axis=1)
    obs_u, obs_pairs = _u_and_pairs(good_values, poor_values)
    with _quiet_all_nan():
        obs_rb = np.where(obs_pairs > 0, 2.0 * obs_u / np.where(obs_pairs > 0, obs_pairs, 1.0) - 1.0, np.nan)

    lo_q = (100.0 - ci) / 2.0
    hi_q = 100.0 - lo_q

    md_lo = np.full(n_endpoints, np.nan)
    md_hi = np.full(n_endpoints, np.nan)
    rb_lo = np.full(n_endpoints, np.nan)
    rb_hi = np.full(n_endpoints, np.nan)

    per_endpoint_elements = max(n_resamples * max(n_good, n_poor), 1)
    chunk_size = max(1, int(chunk_target_elements // per_endpoint_elements))

    for start in range(0, n_endpoints, chunk_size):
        stop = min(start + chunk_size, n_endpoints)
        g_chunk = good_values[start:stop]              # (C, n_good)
        p_chunk = poor_values[start:stop]              # (C, n_poor)
        g_boot = g_chunk[:, good_boot]                 # (C, B, n_good)
        p_boot = p_chunk[:, poor_boot]                 # (C, B, n_poor)

        with _quiet_all_nan():
            boot_md = np.nanmedian(g_boot, axis=2) - np.nanmedian(p_boot, axis=2)
        u_boot, pairs_boot = _u_and_pairs(g_boot, p_boot)
        with _quiet_all_nan():
            boot_rb = np.where(
                pairs_boot > 0,
                2.0 * u_boot / np.where(pairs_boot > 0, pairs_boot, 1.0) - 1.0,
                np.nan,
            )
            md_q = np.nanpercentile(boot_md, [lo_q, hi_q], axis=1)
            rb_q = np.nanpercentile(boot_rb, [lo_q, hi_q], axis=1)
        md_lo[start:stop], md_hi[start:stop] = md_q[0], md_q[1]
        rb_lo[start:stop], rb_hi[start:stop] = rb_q[0], rb_q[1]

        del g_boot, p_boot, boot_md, boot_rb, u_boot, pairs_boot

    return {
        "median_difference": obs_median_diff,
        "median_difference_ci_low": md_lo,
        "median_difference_ci_high": md_hi,
        "rank_biserial": obs_rb,
        "rank_biserial_ci_low": rb_lo,
        "rank_biserial_ci_high": rb_hi,
        "n_good_evaluable": np.isfinite(good_values).sum(axis=1).astype(int),
        "n_poor_evaluable": np.isfinite(poor_values).sum(axis=1).astype(int),
        "n_resamples": n_resamples,
        "seed": int(seed),
        "ci": float(ci),
        "chunk_size": int(chunk_size),
        "good_resample_shape": tuple(good_boot.shape),
        "poor_resample_shape": tuple(poor_boot.shape),
    }


def permutation_p(
    values: Sequence[float],
    labels: Sequence[Any],
    n_perm: int,
    seed: int,
    min_per_group: int,
    *,
    good_label: Any = "Good",
) -> tuple[float, int]:
    """Label-permutation p-value for one endpoint, honouring evaluability.

    The PATIENT labels are permuted with group sizes preserved. Because some
    patients are structurally missing (NaN) for this endpoint, the number of
    *evaluable* patients per group changes from permutation to permutation; it
    is recomputed every time and an allocation that leaves either group with
    fewer than ``min_per_group`` evaluable patients is skipped rather than
    scored. The reported p-value is therefore conditional on evaluability, and
    ``n_valid`` tells you how much of the permutation space survived.

    Statistic: absolute rank-biserial (equivalently the two-sided
    Mann-Whitney U effect), computed NaN-aware.

    Returns
    -------
    (p, n_valid)
        ``p = (extreme + 1) / (n_valid + 1)``; NaN with ``n_valid == 0`` if the
        observed data are themselves not evaluable or no permutation survived.
    """
    values = _as_1d_float(values)
    labels = np.asarray(labels)
    if values.size != labels.size:
        raise ValueError("values and labels must be the same length")

    is_good = labels == good_label
    if not is_good.any() or is_good.all():
        raise ValueError(f"labels must contain both '{good_label}' and the other group")

    finite = np.isfinite(values)
    obs_u, obs_pairs = _u_and_pairs(values[is_good & finite], values[~is_good & finite])
    n_obs_good = int((is_good & finite).sum())
    n_obs_poor = int((~is_good & finite).sum())
    if n_obs_good < min_per_group or n_obs_poor < min_per_group or obs_pairs <= 0:
        return (float("nan"), 0)
    obs_stat = abs(2.0 * float(obs_u) / float(obs_pairs) - 1.0)

    rng = np.random.default_rng(int(seed))
    n_total = values.size
    n_good = int(is_good.sum())

    extreme = 0
    n_valid = 0
    tol = 1e-12
    for _ in range(int(n_perm)):
        perm = rng.permutation(n_total)
        good_pos = perm[:n_good]
        poor_pos = perm[n_good:]
        g_vals = values[good_pos]
        p_vals = values[poor_pos]
        g_ok = np.isfinite(g_vals)
        p_ok = np.isfinite(p_vals)
        if g_ok.sum() < min_per_group or p_ok.sum() < min_per_group:
            continue
        u_perm, pairs_perm = _u_and_pairs(g_vals[g_ok], p_vals[p_ok])
        if pairs_perm <= 0:
            continue
        stat = abs(2.0 * float(u_perm) / float(pairs_perm) - 1.0)
        n_valid += 1
        if stat >= obs_stat - tol:
            extreme += 1

    if n_valid == 0:
        return (float("nan"), 0)
    return ((extreme + 1) / (n_valid + 1), n_valid)


def bh(pvals: Sequence[float], *, method: str = "fdr_bh") -> np.ndarray:
    """Benjamini-Hochberg q-values, NaN-safe.

    Positions holding NaN are treated as *not tested*: they are excluded from
    the family (so they do not inflate ``m``) and come back as NaN. Pass NaN
    for every endpoint whose ``tested`` flag is False -- that is the contract
    the notebooks rely on to guarantee BH is applied to all and only the tested
    endpoints.
    """
    from statsmodels.stats.multitest import multipletests

    arr = _as_1d_float(pvals)
    out = np.full(arr.shape, np.nan, dtype=float)
    tested = np.isfinite(arr)
    if not tested.any():
        return out
    subset = arr[tested]
    if (subset < 0).any() or (subset > 1).any():
        raise ValueError("p-values outside [0, 1]")
    _, qvals, _, _ = multipletests(subset, method=method)
    out[tested] = qvals
    return out


def fisher_availability(avail_good: Any, avail_poor: Any) -> tuple[float, float]:
    """Fisher exact test on endpoint availability (Good vs Poor).

    Accepts either
      * per-patient boolean/0-1 arrays (``avail_good[i]`` = endpoint evaluable
        in Good patient i), or
      * a 2-element ``(n_available, n_total)`` pair of integers per group.

    Returns ``(odds_ratio, p)`` for the 2x2 table
    ``[[avail_good, unavail_good], [avail_poor, unavail_poor]]``; an odds ratio
    > 1 means the endpoint is more often *available* in Good.
    """
    from scipy.stats import fisher_exact

    def _counts(value: Any, name: str) -> tuple[int, int]:
        arr = np.asarray(value)
        if arr.ndim == 1 and arr.size == 2 and arr.dtype.kind in "iu":
            n_avail, n_total = int(arr[0]), int(arr[1])
            if n_total < n_avail or n_avail < 0:
                raise ValueError(f"{name}: (n_available, n_total) is inconsistent")
            return n_avail, n_total - n_avail
        flat = arr.ravel()
        if flat.size == 0:
            raise ValueError(f"{name} is empty")
        as_bool = flat.astype(bool)
        return int(as_bool.sum()), int((~as_bool).sum())

    a_good, u_good = _counts(avail_good, "avail_good")
    a_poor, u_poor = _counts(avail_poor, "avail_poor")
    odds_ratio, pvalue = fisher_exact([[a_good, u_good], [a_poor, u_poor]])
    return (float(odds_ratio), float(pvalue))


# =============================================================================
# safety
# =============================================================================

def safe_unlink(path: str | os.PathLike, allowed_root: str | os.PathLike) -> None:
    """Delete a file only if it resolves inside ``allowed_root``.

    Guards a real hazard inherited from the legacy tree: ``plot_figures.py``
    unlinks fourteen hardcoded filenames from whatever directory it is pointed
    at, so a mistyped output directory could delete files well outside the
    analysis. Escaping the root -- via ``..``, an absolute path or a symlink --
    raises instead of deleting. A missing file is a no-op; a directory is
    refused (this function removes files only).
    """
    root = Path(allowed_root).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"allowed_root is not an existing directory: {root}")

    target = Path(path).expanduser()
    if not target.is_absolute():
        target = root / target
    resolved = target.resolve()

    if resolved == root:
        raise ValueError(f"refusing to unlink the allowed root itself: {resolved}")
    if not resolved.is_relative_to(root):
        raise ValueError(
            f"refusing to unlink {resolved}: outside the allowed root {root}"
        )
    if resolved.is_dir():
        raise IsADirectoryError(f"refusing to unlink a directory: {resolved}")
    resolved.unlink(missing_ok=True)
