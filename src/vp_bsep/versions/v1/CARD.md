# bsep v1

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-08-29 · signature 1

**Why this version.** first release: XGBoost on ECFP4 over the ChEMBL ABCB11 IC50/Ki set, binarised at 100 uM

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `bsep_inhib` | float32 | 0.0–1.0 | P(inhibits BSEP/ABCB11 at IC50 or Ki below 100 uM) |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance

Protocol `scaffold-shuffle-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed, so distinct seeds give distinct test sets. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-08-29 on n_train=863, n_val=115, n_test=173.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.777 | 0.030 | 0.781, 0.805, 0.725, 0.766, 0.809 |
| auprc | 0.594 | 0.148 | 0.485, 0.372, 0.615, 0.739, 0.757 |
| mcc | 0.369 | 0.068 | 0.260, 0.369, 0.336, 0.427, 0.452 |
| brier | 0.183 | 0.023 | 0.182, 0.145, 0.212, 0.198, 0.177 |

> Comparable only with metrics carrying the same protocol id.

## Data

ChEMBL ABCB11 (BSEP) bioactivity — ChEMBL target CHEMBL6020 (Homo sapiens, UniProt O95342), standard_type in ['IC50', 'Ki'], binarised at 100 uM median potency. Retrieved 2026-06-04, licensed CC-BY-SA-3.0, redistributed here.

`1151` compounds, positive rate `0.356`, table SHA-256 `77c30dcd81a261f2…`

Regenerate and check for upstream drift with `python -m vp_bsep.data fetch --verify`.

## Model

xgboost on `morgan2_2048` features. Shipped weights: full dataset minus a 10% scaffold carve used only for early stopping

`weights.joblib` SHA-256 `b3c8c342aa26c456…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
