# vp-bsep

Predicts inhibition of the bile salt export pump (BSEP / ABCB11) from a SMILES
string — the molecular initiating event for cholestatic drug-induced liver
injury.

## Install

    pip install vp-bsep

## Use

    from vp_bsep import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])      # -> DataFrame[bsep_inhib, bsep_potency_um]

Returns one row per input and one column per declared output: a probability of
inhibition at the 100 µM cutoff, and the predicted IC50 or Ki in µM to read
against expected exposure. Unparseable SMILES come back as NaN. Pin a version
with `predict(smiles, version="v2")`; list what is available with `versions()`.

## Current version

**v2**, signature 2, measured under protocol `scaffold-shuffle-5seed@1`. The
full record — metrics per seed, dataset hash, environment — is in
[`src/vp_bsep/versions/v2/CARD.md`](src/vp_bsep/versions/v2/CARD.md).

## Data

ChEMBL target CHEMBL6020 (human ABCB11, UniProt O95342), IC50 and Ki values
reduced to one row per compound carrying median potency, assay ceiling and a
label binarised at 100 µM, retrieved 2026-06-04 and redistributed here under
CC-BY-SA-3.0. Rebuild it and check for upstream drift with
`python -m vp_bsep.data fetch --verify`; see [`data/README.md`](data/README.md)
for the expected layout.

## Retrain

    python -m vp_bsep.train --version v3 --reason "why this version exists"
    python -m vp_bsep.evaluate --version v3

`train` fits the deployment model on the whole dataset and writes a new version
directory; `evaluate` refits per seed under the protocol and records what those
held-out models scored. Reproducibility is to the recorded dataset hash and
environment, which can change.

## Licence

Code is Apache-2.0 ([`LICENSE`](LICENSE)). The bundled dataset is CC-BY-SA-3.0
([`LICENSE-DATA`](LICENSE-DATA)).

## Cite

See [`CITATION.cff`](CITATION.cff).
