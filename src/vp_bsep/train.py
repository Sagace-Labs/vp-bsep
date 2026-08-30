"""Build a new released version.

    python -m vp_bsep.train --version v2 --reason "why this version exists"

Writes ``versions/<version>/manifest.toml`` and ``weights.joblib``.

The shipped weights are fit on the **whole** dataset, with a scaffold carve
held out only for early stopping.
"""

from __future__ import annotations

import argparse
import platform
import shutil
import sys
from datetime import date
from pathlib import Path

from vp_bsep import contract as bsep_contract
from vp_bsep import data as bsep_data
from vp_bsep import model as bsep_model
from vp_bsep.target import CUTOFF_UM, TARGET

__all__ = ["build_version", "main"]

VERSIONS_DIR = Path(__file__).resolve().parent / "versions"

DEFAULT_PROTOCOL = "scaffold-shuffle-5seed@1"


def _provenance() -> dict:
    """The environment that produced the version."""
    import rdkit
    import xgboost

    return {
        "python": platform.python_version(),
        "rdkit": rdkit.__version__,
        "xgboost": xgboost.__version__,
    }


def build_version(
    version: str,
    *,
    reason: str,
    protocol: str = DEFAULT_PROTOCOL,
    supersedes: str | None = None,
    seed: int = 0,
) -> Path:
    """Fit the deployment model and write the version directory."""
    from vp_core import fingerprints, protocols
    from vp_core.splits import scaffold_train_val

    protocols.get(protocol)  # fail early on an unknown protocol

    directory = VERSIONS_DIR / version
    if directory.exists():
        raise FileExistsError(
            f"{directory} already exists. Released versions are immutable — "
            "publish a new version instead of editing this one."
        )

    table = bsep_data.load()
    smiles = table["smiles"].tolist()
    y = table["label"].to_numpy(dtype=int)

    # Deployment fit: everything, with a small scaffold carve for early stopping.
    train_idx, val_idx = scaffold_train_val(smiles, val_frac=0.10, seed=seed)
    X = fingerprints.featurize(smiles, bsep_model.FEATURES)
    fitted = bsep_model.fit(X[train_idx], y[train_idx], X[val_idx], y[val_idx], seed=seed)

    directory.mkdir(parents=True)
    try:
        return _write_version(
            directory,
            fitted,
            table,
            y,
            version=version,
            reason=reason,
            protocol=protocol,
            supersedes=supersedes,
        )
    except Exception:
        shutil.rmtree(directory, ignore_errors=True)
        raise


def _write_version(
    directory: Path,
    fitted,
    table,
    y,
    *,
    version: str,
    reason: str,
    protocol: str,
    supersedes: str | None,
) -> Path:
    import joblib

    import vp_core
    from vp_core import dataset, hashing, manifest

    weights_path = directory / "weights.joblib"
    joblib.dump(fitted, weights_path)

    record = {
        "schema": manifest.SCHEMA_VERSION,
        "pathway": "bsep",
        "version": version,
        "released": date.today().isoformat(),
        "supersedes": supersedes,
        "reason": reason,
        "signature": bsep_contract.as_manifest_table(),
        "dataset": {
            "name": "ChEMBL ABCB11 (BSEP) bioactivity",
            "source": (
                f"ChEMBL target {TARGET.chembl_target} "
                f"({TARGET.organism}, UniProt {TARGET.uniprot}), "
                f"standard_type in {list(TARGET.potency_types)}, "
                f"binarised at {CUTOFF_UM:g} uM median potency"
            ),
            "url": (
                "https://www.ebi.ac.uk/chembl/api/data/activity"
                f"?target_chembl_id={TARGET.chembl_target}"
            ),
            "retrieved": "2026-06-04",
            "licence": "CC-BY-SA-3.0",
            "redistributable": True,
            "path": "data/bsep_chembl.parquet",
            "sha256": dataset.dataset_hash(table),
            "n_rows": len(table),
            "base_rate": round(float(y.mean()), 6),
            "fetch": "python -m vp_bsep.data fetch --verify",
        },
        "model": {
            "family": "xgboost",
            "features": bsep_model.FEATURES,
            "fit": (
                "full dataset minus a 10% scaffold carve used only for early stopping"
            ),
            "weights": "weights.joblib",
            "sha256": hashing.sha256_file(weights_path),
        },
        "protocol": {
            "id": protocol,
            "provider": "vp-core",
            "core_version": vp_core.__version__,
        },
        "provenance": _provenance(),
    }

    problems = manifest.validate(record, version_dir=directory)
    if problems:
        raise ValueError(f"refusing to write an invalid manifest: {'; '.join(problems)}")
    manifest.write(directory / "manifest.toml", record)

    print(
        f"wrote {directory}\n"
        f"  {len(table)} compounds, positive rate {y.mean():.3f}\n"
        f"  best iteration {getattr(fitted, 'best_iteration', None)}\n"
        f"  next: python -m vp_bsep.evaluate --version {version}",
        file=sys.stderr,
    )
    return directory


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_bsep.train")
    parser.add_argument("--version", required=True, help="new version name, e.g. v2")
    parser.add_argument("--reason", required=True, help="why this version exists")
    parser.add_argument("--protocol", default=DEFAULT_PROTOCOL)
    parser.add_argument("--supersedes", default=None)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    build_version(
        args.version,
        reason=args.reason,
        protocol=args.protocol,
        supersedes=args.supersedes,
        seed=args.seed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
