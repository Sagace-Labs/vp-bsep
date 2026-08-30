"""BSEP / ABCB11 inhibition from a SMILES string.

Inhibition of the bile salt export pump is the canonical molecular initiating
event for cholestatic drug-induced liver injury. Blocking the transporter causes bile acids to accumulate inside the hepatocyte.

    from vp_bsep import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])      # -> DataFrame[bsep_inhib]
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from vp_bsep.target import TARGET, TARGETS, Target, all_names
from vp_bsep.target import get as get_target
from vp_core.registry import Version, VersionedPathway

__version__ = "1.0.0"

PATHWAY = "bsep"
VERSIONS_DIR = Path(__file__).resolve().parent / "versions"

_pathway = VersionedPathway(PATHWAY, VERSIONS_DIR)

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
