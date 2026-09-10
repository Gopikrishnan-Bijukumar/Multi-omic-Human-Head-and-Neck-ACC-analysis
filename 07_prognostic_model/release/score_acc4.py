"""The one scoring function for the locked four-gene ACC signature.

Panel: DSCAM, ODC1, NCAPG, CCNB2.  All four discovery coefficients are positive, so the
score is the unweighted mean of rank-inverse-normal expression; higher = worse prognosis.
The hierarchical Bayesian procedure selected the genes and their directions - it does not
weight the deployed score.

Two transforms, and the difference matters:

    score_batch(X)              cohort-relative. Needs the whole cohort at once and a
                                patient's score changes when the cohort changes.
                                This produced the reported DK C = 0.671.

    score_frozen(x, X_ref)      each sample is mapped into a FROZEN discovery reference
                                independently of every other sample. This is what a
                                deployed assay does. DK C = 0.689, and the discovery-locked
                                threshold of -0.062172 applies to this score only.
                                Post hoc: designed after DK was unblinded.

Both are lifted verbatim from work/27_threshold.py (rint_batch, rint_vs_reference).

    from score_acc4 import PANEL, score_batch, score_frozen, THRESHOLD
"""
import json
import os

import numpy as np
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPEC = json.load(open(os.path.join(_HERE, "model_spec.json")))

PANEL = list(_SPEC["genes"])
SIGNS = np.array([_SPEC["signs"][g] for g in PANEL], dtype=float)
THRESHOLD = _SPEC["frozen_threshold"]["value"]


def _as_matrix(X):
    """Accept a DataFrame with gene columns or a plain array already in PANEL order."""
    if hasattr(X, "columns"):
        missing = [g for g in PANEL if g not in X.columns]
        if missing:
            raise KeyError(f"missing panel genes: {missing}")
        return np.asarray(X[PANEL], dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] != len(PANEL):
        raise ValueError(f"expected (n_samples, {len(PANEL)}) in PANEL order, got {X.shape}")
    return X


def rint_batch(X):
    """Cohort-relative rank-inverse-normal transform (Blom)."""
    r = stats.rankdata(X, axis=0)
    n = X.shape[0]
    return stats.norm.ppf((r - 0.375) / (n + 0.25))


def rint_vs_reference(x_new, X_ref):
    """Map each row of x_new into X_ref's rank space, using no other row of x_new."""
    ref = np.sort(X_ref, axis=0)
    n = ref.shape[0]
    out = np.empty(np.shape(x_new), dtype=float)
    for j in range(np.shape(x_new)[1]):
        pos = np.searchsorted(ref[:, j], np.asarray(x_new)[:, j], side="left")
        out[:, j] = stats.norm.ppf((pos + 0.5) / (n + 1.0))
    return out


def score_batch(X):
    """Cohort-relative score. X: (n_samples, 4) log2-normalised expression, PANEL order."""
    return (rint_batch(_as_matrix(X)) * SIGNS).mean(1)


def score_frozen(x_new, X_ref):
    """Single-sample score against a frozen discovery reference. Post hoc; see module docstring."""
    return (rint_vs_reference(_as_matrix(x_new), _as_matrix(X_ref)) * SIGNS).mean(1)


def classify_frozen(x_new, X_ref):
    """High-risk call using the discovery-locked threshold. Uses no validation-cohort information."""
    return score_frozen(x_new, X_ref) > THRESHOLD
