# BSEP data

    data/
      bsep_chembl.parquet          standardised table (producing the hash)
      example/
        bsep_example.parquet       small stratified fixture, used by the tests

`bsep_chembl.parquet` holds one row per compound with the three contract
columns — `inchikey`, `smiles` (standardised), `label` — plus `potency_um`,
`censored_at` and `n_measurements` as provenance. Only the contract columns are
hashed, so a provenance column may be added without moving the dataset hash.

## Origin and processing

ChEMBL target CHEMBL6020 (bile salt export pump, *Homo sapiens*, UniProt
O95342), activity records of `standard_type` IC50 or Ki, retrieved through the
public REST API. Nanomolar potencies are converted to micromolar; a `>`
relation is treated as a censored inactive; structures are standardised
(normalise, largest fragment, neutralise); measurements are collapsed to one
row per InChIKey by median potency, and the binary label applies a 100 µM
cutoff to that median.

`potency_um` is the median over uncensored measurements and is infinite when
there are none. `censored_at` is the highest concentration the compound was
tested to without inhibiting, and is NaN when no measurement was censored.

## Rebuilding it

    python -m vp_bsep.data fetch --verify

This downloads, re-parses and re-hashes, then compares against the hash the
current version recorded. A mismatch in the hash can indicate upstream changes
to the source data or the processing pipeline.

## Licence

The bundled table is derived from ChEMBL and carries CC-BY-SA-3.0; see
`../LICENSE-DATA`. That licence covers these files only, not the package code
or the trained weights.
