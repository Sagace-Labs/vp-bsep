"""The BSEP / ABCB11 target identity.

Frozen so the data source is auditable.

The ChEMBL id was verified live against the target-search API on 2026-06-04:

    CHEMBL6020     Bile salt export pump   Homo sapiens    SINGLE PROTEIN
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["CUTOFF_UM", "TARGET", "TARGETS", "Target", "all_names", "get"]


@dataclass(frozen=True)
class Target:
    name: str
    gene: str
    uniprot: str
    chembl_target: str
    potency_types: tuple[str, ...] = ("IC50", "Ki")
    organism: str = "Homo sapiens"


TARGET = Target(
    name="Bile salt export pump",
    gene="ABCB11",
    uniprot="O95342",
    chembl_target="CHEMBL6020",
)

TARGETS: dict[str, Target] = {"ABCB11": TARGET}

# Inhibitor cutoff, following the 100 µM screening threshold Morgan 2010 reports
# for BSEP. Published cutoffs span 10-300 µM.
CUTOFF_UM = 100.0


def get(name: str) -> Target:
    try:
        return TARGETS[name.upper()]
    except KeyError:
        raise KeyError(f"unknown target {name!r}; known: {sorted(TARGETS)}") from None


def all_names() -> list[str]:
    return sorted(TARGETS)
