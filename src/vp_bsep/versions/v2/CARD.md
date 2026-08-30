# bsep v2

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-08-30 · signature 2 · supersedes v1

**Why this version.** fits the measured potency interval instead of the binary label, so a compound observed only above an assay ceiling stays a bound; adds predicted potency as a second output

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `bsep_inhib` | float32 | 0.0–1.0 | P(inhibits BSEP/ABCB11 at IC50 or Ki below 100 uM) |
| `bsep_potency_um` | float32 | — | Predicted BSEP/ABCB11 IC50 or Ki in uM. Read against the exposure expected in the liver; the probability alone fixes one threshold |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance — `scaffold-shuffle-5seed@1`

Protocol `scaffold-shuffle-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed, so distinct seeds give distinct test sets. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-08-30 on n_train=863, n_val=115, n_test=173.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.868 | 0.016 | 0.865, 0.843, 0.864, 0.879, 0.891 |
| auprc | 0.699 | 0.171 | 0.510, 0.475, 0.804, 0.838, 0.870 |
| mcc | 0.464 | 0.146 | 0.494, 0.179, 0.555, 0.512, 0.583 |
| brier | 0.134 | 0.015 | 0.118, 0.114, 0.146, 0.145, 0.147 |

## Performance — `scaffold-balanced-5seed@1`

Protocol `scaffold-balanced-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed and each group placed in the fold it overfills least, so a group larger than a fold's capacity settles in train instead of starving that fold. Same fold fractions, seeds and metrics as scaffold-shuffle-5seed@1; only the packing differs. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-19 on n_train=863, n_val=115, n_test=173.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.876 | 0.031 | 0.910, 0.882, 0.844, 0.908, 0.836 |
| auprc | 0.798 | 0.058 | 0.877, 0.809, 0.721, 0.837, 0.744 |
| mcc | 0.577 | 0.057 | 0.647, 0.559, 0.526, 0.642, 0.514 |
| brier | 0.144 | 0.017 | 0.130, 0.140, 0.159, 0.123, 0.168 |

> Comparable only with metrics carrying the same protocol id.

## Data

ChEMBL ABCB11 (BSEP) bioactivity — ChEMBL target CHEMBL6020 (Homo sapiens, UniProt O95342), standard_type in ['IC50', 'Ki'], one row per compound carrying median potency, assay ceiling and a label binarised at 100 uM. Retrieved 2026-06-04, licensed CC-BY-SA-3.0, redistributed here.

`1151` compounds, positive rate `0.356`, table SHA-256 `77c30dcd81a261f2…`

Regenerate and check for upstream drift with `python -m vp_bsep.data fetch --verify`.

## Model

xgboost-aft on `morgan2c_physchem_ion` features. Shipped weights: full dataset minus a 10% scaffold carve used for early stopping and probability calibration

`weights.joblib` SHA-256 `d53e1714570effb6…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
