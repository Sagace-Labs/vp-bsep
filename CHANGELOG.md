# Changelog

Released versions are immutable. A correction to a released version is a new
patch version. This file records *why* each version exists.

## v2

Fits the measured potency interval rather than the binary label. Over half the
compounds were only ever observed above an assay ceiling; as a bound each one
still informs the fit.
The probability is calibrated against the cutoff on the held-out carve.

Features gain per-substructure counts and a descriptor panel with explicit
ionisable-group counts, because standardisation neutralises a molecule before
it is featurised.

Signature 2: `bsep_potency_um` joins `bsep_inhib`.

## v1

Initial release.

Signature 1: a single column `bsep_inhib`.
