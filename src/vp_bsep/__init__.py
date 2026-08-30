"""BSEP / ABCB11 inhibition from a SMILES string.

Inhibition of the bile salt export pump is the canonical molecular initiating
event for cholestatic drug-induced liver injury.

    from vp_bsep import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])      # -> DataFrame[bsep_inhib, bsep_potency_um]
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from vp_bsep.target import TARGET, TARGETS, Target, all_names
from vp_bsep.target import get as get_target
from vp_core.registry import Version, VersionedPathway

__version__ = "2.0.0"

PATHWAY = "bsep"
VERSIONS_DIR = Path(__file__).resolve().parent / "versions"


def _predict_values(model: Any, smiles: list[str], version: Version) -> np.ndarray:
    """Columns for ``version``, in the order its signature declares them.

    Signature 1 emits a probability; signature 2 adds a potency.
    """
    from rdkit import Chem, RDLogger

    from vp_core import fingerprints, xgb

    RDLogger.DisableLog("rdApp.*")
    X = fingerprints.featurize(smiles, str(version.features))
    if isinstance(model, xgb.CensoredModel):
        values = np.column_stack([model.probability(X), model.potency(X)])
    else:
        values = xgb.predict_proba(model, X).reshape(-1, 1)

    # An unparseable input is a declared NaN.
    unparseable = [Chem.MolFromSmiles(s) is None for s in smiles]
    values = values.astype(np.float32)
    values[np.asarray(unparseable)] = np.nan
    return values


_pathway = VersionedPathway(PATHWAY, VERSIONS_DIR, predict_fn=_predict_values)

__all__ = [
    "PATHWAY",
    "TARGET",
    "TARGETS",
    "VERSIONS_DIR",
    "Target",
    "Version",
    "__version__",
    "all_names",
    "current_version",
    "get",
    "get_target",
    "predict",
    "signature",
    "versions",
]


def predict(smiles: list[str], *, version: str | None = None) -> pd.DataFrame:
    """P(BSEP inhibitor) for each SMILES, using ``version`` (default: newest)."""
    return _pathway.predict(smiles, version=version)


def versions() -> list[str]:
    """Released version names, oldest first."""
    return _pathway.versions()


def current_version() -> str:
    """The newest released version."""
    return _pathway.current()


def get(version: str | None = None) -> Version:
    """Load a version and its validated manifest."""
    return _pathway.get(version)


def signature(version: str | None = None) -> dict[str, Any]:
    """The output contract of a version."""
    return get(version).manifest["signature"]
