"""Measure a version under an evaluation protocol.

    python -m vp_bsep.evaluate --version v1

Writes ``versions/<version>/metrics.json`` and regenerates ``CARD.md``. The
protocol named in the version's manifest supplies the split, the fold sizes,
the seed set and the metric list. ``--protocol`` measures the same weights
under another registered protocol; ``metrics.json`` keys every run by
protocol id.

Each seed refits the model on its own training fold.

A seed that yields an empty fold raises rather than being skipped.
"""

from __future__ import annotations

import argparse
import sys
from datetime import date

import numpy as np

from vp_bsep import data as bsep_data
from vp_bsep import model as bsep_model

__all__ = ["evaluate_version", "main"]


def evaluate_version(
    version: str,
    *,
    protocol_id: str | None = None,
    use_example: bool = False,
    write: bool = True,
) -> dict:
    """Run the version's protocol and return the metrics record."""
    import vp_bsep
    from vp_core import card, dataset, fingerprints, metrics_store, protocols
    from vp_core import metrics as metrics_mod

    resolved = vp_bsep.get(version)
    protocol_id = protocol_id or resolved.protocol_id
    protocol = protocols.get(protocol_id)

    table = bsep_data.example() if use_example else bsep_data.load()
    smiles = table["smiles"].tolist()
    y = table["label"].to_numpy(dtype=int)
    X = fingerprints.featurize(smiles, bsep_model.FEATURES)

    per_seed: list[dict[str, float]] = []
    folds: list[dict[str, int]] = []
    for seed in protocol.seeds:
        train_idx, val_idx, test_idx = protocol.split_indices(smiles, seed)
        fitted = bsep_model.fit(
            X[train_idx], y[train_idx], X[val_idx], y[val_idx], seed=seed
        )
        proba = bsep_model.predict_proba(fitted, [smiles[i] for i in test_idx])
        scored = metrics_mod.binary_metrics(y[test_idx], proba)
        per_seed.append(scored)
        folds.append(
            {
                "seed": int(seed),
                "train": len(train_idx),
                "val": len(val_idx),
                "test": len(test_idx),
            }
        )
        print(
            f"  seed {seed}: test AUC {scored['auc_roc']:.4f} "
            f"(n_test={len(test_idx)})",
            file=sys.stderr,
        )

    record = {
        "pathway": "bsep",
        "version": resolved.name,
        "protocol_id": protocol_id,
        "dataset_sha256": dataset.dataset_hash(table),
        "dataset": "example fixture" if use_example else "full",
        "evaluated": date.today().isoformat(),
        "n": {
            "total": len(table),
            "train": folds[0]["train"],
            "val": folds[0]["val"],
            "test": folds[0]["test"],
        },
        "folds": folds,
        "test": metrics_mod.aggregate(per_seed, protocol.metrics),
    }

    if write:
        if use_example:
            raise ValueError(
                "refusing to record fixture metrics as a released result — "
                "--example is for smoke-checking the harness only"
            )
        path = metrics_store.write(resolved.directory, record)
        card.write_card(resolved.directory)
        auc = record["test"]["auc_roc"]
        print(
            f"wrote {path}\n"
            f"  {protocol_id}: AUC {auc['mean']:.4f} +/- {auc['std']:.4f} "
            f"over {len(protocol.seeds)} seeds",
            file=sys.stderr,
        )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_bsep.evaluate")
    parser.add_argument("--version", default=None, help="default: the newest version")
    parser.add_argument(
        "--protocol",
        default=None,
        help="default: the protocol the version's manifest declares",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="run on the committed fixture without writing (smoke check only)",
    )
    args = parser.parse_args(argv)

    import vp_bsep

    version = args.version or vp_bsep.current_version()
    record = evaluate_version(
        version,
        protocol_id=args.protocol,
        use_example=args.example,
        write=not args.example,
    )
    if args.example:
        auc = record["test"]["auc_roc"]["mean"]
        print(f"fixture smoke AUC {auc if np.isfinite(auc) else float('nan'):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
