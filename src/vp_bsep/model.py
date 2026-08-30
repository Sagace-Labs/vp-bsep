"""The BSEP model: censored XGBoost on counted ECFP4 with a descriptor panel.

Hyperparameters and feature choice.

The target is the measured potency interval rather than the binary label. The
probability is calibrated against the cutoff on the validation fold.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from vp_bsep.target import CUTOFF_UM
from vp_core import fingerprints, xgb

__all__ = ["FEATURES", "HYPERPARAMS", "fit", "predict"]

FEATURES = "morgan2c_physchem_ion"

HYPERPARAMS: dict[str, Any] = {
    "n_estimators": 2000,
    "learning_rate": 0.05,
    "max_depth": 6,
    "min_child_weight": 1.0,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 2.0,
    "reg_alpha": 0.05,
    "gamma": 0.1,
    "early_stopping_rounds": 40,
}


def fit(
    X_train,
    interval_train: tuple[np.ndarray, np.ndarray],
    X_val,
    interval_val: tuple[np.ndarray, np.ndarray],
    y_val,
    *,
    seed: int = 0,
) -> xgb.CensoredModel:
    """Fit one censored model. The validation fold stops boosting and calibrates."""
    return xgb.fit_censored(
        X_train,
        *interval_train,
        X_val,
        *interval_val,
        y_val,
        params=dict(HYPERPARAMS),
        cutoff_um=CUTOFF_UM,
        seed=seed,
    )


def predict(model: xgb.CensoredModel, smiles: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """``(probability, potency_um)`` for arbitrary SMILES."""
    X = fingerprints.featurize(smiles, FEATURES)
    return model.probability(X), model.potency(X)
