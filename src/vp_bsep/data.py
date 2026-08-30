"""The BSEP dataset: obtain, standardise, verify.

The shipped table is the *parsed* result — one row per compound, keyed on the
standardised InChIKey — and that is what the manifest hash covers. ``fetch``
rebuilds it from ChEMBL and compares hashes, so its real job is drift
detection: a differing hash means the upstream target data changed, which is a
scientific event worth a new version.

Pipeline: page the ChEMBL activity API for CHEMBL6020 restricted to IC50/Ki,
keep nanomolar potencies, treat a ``>`` relation as a censored inactive,
standardise, then collapse to one row per compound by median potency and
re-apply the 100 µM cutoff to the aggregate.

Run as a command:

    python -m vp_bsep.data fetch --verify
    python -m vp_bsep.data verify
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from vp_bsep.target import CUTOFF_UM, TARGET

__all__ = [
    "DATA_DIR",
    "EXAMPLE_PATH",
    "TABLE_PATH",
    "build_example",
    "example",
    "fetch",
    "load",
    "verify",
]

_PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = _PACKAGE_ROOT / "data"
TABLE_PATH = DATA_DIR / "bsep_chembl.parquet"
EXAMPLE_PATH = DATA_DIR / "example" / "bsep_example.parquet"

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_MIN_INTERVAL = 0.34  # ChEMBL asks for roughly one request per second


def _missing(path: Path) -> str:
    return (
        f"{path} is not present. The dataset is not included in the wheel; run "
        "`python -m vp_bsep.data fetch` to rebuild it from ChEMBL."
    )


def load() -> pd.DataFrame:
    """The full standardised dataset. Raises with guidance when absent."""
    from vp_core import dataset

    if not TABLE_PATH.exists():
        raise FileNotFoundError(_missing(TABLE_PATH))
    return dataset.read_table(TABLE_PATH)


def example() -> pd.DataFrame:
    """The committed test fixture — small, stratified, always available."""
    from vp_core import dataset

    if not EXAMPLE_PATH.exists():
        raise FileNotFoundError(_missing(EXAMPLE_PATH))
    return dataset.read_table(EXAMPLE_PATH)


# ---------------------------------------------------------------------------
# Rebuild from source
# ---------------------------------------------------------------------------


def _page_activities(page_size: int = 1000) -> pd.DataFrame:
    import requests

    only = (
        "molecule_chembl_id,canonical_smiles,standard_type,"
        "standard_relation,standard_value,standard_units"
    )
    rows: list[dict] = []
    offset = 0
    session = requests.Session()
    while True:
        url = (
            f"{CHEMBL_BASE}/activity?target_chembl_id={TARGET.chembl_target}"
            f"&standard_type__in={','.join(TARGET.potency_types)}"
            f"&format=json&limit={page_size}&offset={offset}&only={only}"
        )
        response = session.get(url, timeout=120)
        response.raise_for_status()
        payload = response.json()
        activities = payload.get("activities", [])
        rows.extend(activities)
        total = payload.get("page_meta", {}).get("total_count")
        print(f"  fetched {len(rows)} / {total}", file=sys.stderr)
        if not activities or (total is not None and len(rows) >= total):
            break
        offset += page_size
        time.sleep(_MIN_INTERVAL)

    if not rows:
        raise RuntimeError(
            f"ChEMBL returned no activities for {TARGET.chembl_target}. The target "
            "or endpoint may have changed upstream."
        )
    return pd.DataFrame(rows)


def _to_compounds(raw: pd.DataFrame, cutoff_um: float) -> pd.DataFrame:
    from vp_core.standardise import standardise_many

    df = raw[raw["standard_units"] == "nM"].copy()
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df[df["standard_value"].notna() & (df["standard_value"] > 0)]
    df = df[df["canonical_smiles"].notna() & (df["canonical_smiles"].astype(str) != "")]
    df["potency_um"] = df["standard_value"].astype(float) / 1000.0
    # A ">" relation means the compound did not reach the threshold within the
    # tested range: a censored inactive, negative regardless of the value.
    df["censored"] = df["standard_relation"].fillna("=").eq(">")

    smiles, keys = standardise_many(df["canonical_smiles"].astype(str).tolist())
    df["smiles"] = smiles
    df["inchikey"] = keys
    df = df.dropna(subset=["smiles", "inchikey"])

    records = []
    for inchikey, group in df.groupby("inchikey"):
        measured = group.loc[~group["censored"], "potency_um"]
        median = float(np.median(measured)) if len(measured) else float("inf")
        records.append(
            {
                "inchikey": inchikey,
                "smiles": group["smiles"].iloc[0],
                "label": int(median <= cutoff_um),
                "potency_um": median,
                "n_measurements": len(group),
            }
        )
    return pd.DataFrame(records).sort_values("inchikey").reset_index(drop=True)


def fetch(*, cutoff_um: float = CUTOFF_UM, write: bool = True) -> pd.DataFrame:
    """Rebuild the dataset from ChEMBL. Returns the standardised table."""
    from vp_core import dataset

    print(f"fetching {TARGET.chembl_target} activities from ChEMBL…", file=sys.stderr)
    table = _to_compounds(_page_activities(), cutoff_um)
    problems = dataset.validate_table(table)
    if problems:
        raise ValueError(f"rebuilt table is invalid: {'; '.join(problems)}")
    if write:
        dataset.write_table(table, TABLE_PATH)
    print(
        f"{len(table)} compounds, {int(table['label'].sum())} positive "
        f"({table['label'].mean():.1%}), sha256 {dataset.dataset_hash(table)}",
        file=sys.stderr,
    )
    return table


def verify(version: str | None = None) -> dict:
    """Compare the on-disk table against the hash a released version recorded."""
    import vp_bsep
    from vp_core import dataset

    resolved = vp_bsep.get(version)
    declared = resolved.manifest.get("dataset", {}).get("sha256")
    actual = dataset.dataset_hash(load())
    return {
        "version": resolved.name,
        "declared": declared,
        "actual": actual,
        "match": declared == actual,
    }


def build_example(n: int = 200, seed: int = 0) -> pd.DataFrame:
    """Regenerate the committed fixture from the full table."""
    from vp_core import dataset

    sample = dataset.stratified_example(load(), n=n, seed=seed)
    dataset.write_table(sample, EXAMPLE_PATH)
    return sample


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_bsep.data")
    sub = parser.add_subparsers(dest="command", required=True)

    fetch_cmd = sub.add_parser("fetch", help="rebuild the dataset from ChEMBL")
    fetch_cmd.add_argument(
        "--verify",
        action="store_true",
        help="after fetching, compare the hash against the current version",
    )
    sub.add_parser("verify", help="check the on-disk table against the recorded hash")
    example_cmd = sub.add_parser("build-example", help="regenerate the test fixture")
    example_cmd.add_argument("--n", type=int, default=200)

    args = parser.parse_args(argv)

    if args.command == "fetch":
        fetch()
        if not args.verify:
            return 0
    if args.command == "build-example":
        sample = build_example(n=args.n)
        print(f"wrote {len(sample)} rows to {EXAMPLE_PATH}")
        return 0

    report = verify()
    if report["match"]:
        print(f"dataset matches {report['version']}: {report['actual']}")
        return 0
    print(
        f"DRIFT: {report['version']} recorded {report['declared']}\n"
        f"       the current table hashes to {report['actual']}\n"
        "The upstream source has changed. Record it as a new version rather than "
        "overwriting the released one.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
