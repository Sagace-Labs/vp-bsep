"""BSEP package tests.

Everything here runs on the committed fixture and the shipped weights.
"""

from __future__ import annotations

import numpy as np
import pytest

import vp_bsep
from vp_bsep import contract, data

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
NONSENSE = "not-a-molecule"


def test_versions_are_discoverable_and_ordered():
    versions = vp_bsep.versions()
    assert versions, "no released versions found"
    assert vp_bsep.current_version() == versions[-1]


def test_declared_signature_matches_the_shipped_manifest():
    shipped = vp_bsep.signature()
    assert shipped["version"] == contract.SIGNATURE_VERSION
    assert [o["name"] for o in shipped["outputs"]] == contract.column_names()


def test_predictions_carry_the_declared_column_names():
    frame = vp_bsep.predict([ASPIRIN, CAFFEINE])
    assert list(frame.columns) == contract.column_names()
    assert len(frame) == 2


def test_probabilities_are_in_range():
    values = vp_bsep.predict([ASPIRIN, CAFFEINE])["bsep_inhib"].to_numpy()
    assert np.all((values >= 0.0) & (values <= 1.0))


def test_unparseable_input_becomes_nan_rather_than_raising():
    frame = vp_bsep.predict([ASPIRIN, NONSENSE])
    assert np.isfinite(frame["bsep_inhib"].iloc[0])
    assert np.isnan(frame["bsep_inhib"].iloc[1])


def test_a_version_can_be_pinned():
    """Every release stays loadable under its own contract, not the newest one."""
    for name in vp_bsep.versions():
        frame = vp_bsep.predict([ASPIRIN], version=name)
        declared = [o["name"] for o in vp_bsep.signature(name)["outputs"]]
        assert list(frame.columns) == declared

    assert vp_bsep.predict([ASPIRIN]).equals(
        vp_bsep.predict([ASPIRIN], version=vp_bsep.current_version())
    )


def test_predicted_potency_is_positive():
    potency = vp_bsep.predict([ASPIRIN, CAFFEINE])["bsep_potency_um"].to_numpy()
    assert np.all(np.isfinite(potency))
    assert np.all(potency > 0.0)


def test_a_single_string_is_rejected():
    with pytest.raises(TypeError, match="list of SMILES"):
        vp_bsep.predict(ASPIRIN)  # type: ignore[arg-type]


def test_unknown_version_names_the_available_ones():
    with pytest.raises(FileNotFoundError, match="v1"):
        vp_bsep.predict([ASPIRIN], version="v999")


def test_fixture_satisfies_the_dataset_contract():
    from vp_core import dataset as core_dataset

    fixture = data.example()
    assert core_dataset.validate_table(fixture) == []
    assert set(fixture["label"]) == {0, 1}


def test_target_registry_is_uniform():
    assert vp_bsep.all_names() == ["ABCB11"]
    assert vp_bsep.get_target("abcb11").chembl_target == "CHEMBL6020"
    with pytest.raises(KeyError):
        vp_bsep.get_target("nope")


def test_evaluation_harness_runs_on_the_fixture():
    """A smoke run must never be recordable as a released result."""
    from vp_bsep.evaluate import evaluate_version

    record = evaluate_version(vp_bsep.current_version(), use_example=True, write=False)
    assert record["dataset"] == "example fixture"
    with pytest.raises(ValueError, match="refusing to record"):
        evaluate_version(vp_bsep.current_version(), use_example=True, write=True)
