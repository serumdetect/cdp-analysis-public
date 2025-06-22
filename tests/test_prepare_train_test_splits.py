"""Tests for workflow/scripts/prepare_tcr_data/prepare_train_test_splits.py."""

import itertools

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.prepare_tcr_data import prepare_train_test_splits


@pytest.mark.parametrize(
    ["num_cv_folds", "cv_fold", "assertion_error"],
    [
        (None, None, False),
        (10, 0, False),
        (10, 1, False),
        (10, None, True),
        (None, 5, True),
    ],
)
def test_validate_cv_params(num_cv_folds, cv_fold, assertion_error):
    """Test validate_cv_params."""
    if assertion_error:
        with pytest.raises(AssertionError):
            prepare_train_test_splits.validate_cv_params(num_cv_folds, cv_fold)
    else:
        prepare_train_test_splits.validate_cv_params(num_cv_folds, cv_fold)


@pytest.mark.parametrize("random_state", [None, 42])
def test_cross_validate_tcr_data(q42_tcrs, random_state):
    """Test cross_validate_tcr_data.

    Tests that with 10-fold CV all samples are used exactly once.
    """
    num_cv_folds = 10
    metadata = pd.DataFrame({"sample_name": q42_tcrs["sample_name"].unique()})
    metadata["is_cancer"] = np.mod(np.arange(len(metadata)), 2).astype(bool)
    cv_data = [
        prepare_train_test_splits.cross_validate_tcr_data(
            num_cv_folds,
            cv_fold,
            q42_tcrs,
            metadata,
            "is_cancer",
            random_state=random_state,
        )
        for cv_fold in range(num_cv_folds)
    ]

    # Assert that each sample appears exactly once.
    samples = set(
        itertools.chain.from_iterable(
            [x[3]["sample_name"].tolist() for x in cv_data]
        )
    )
    assert len(samples) == len(metadata)
    assert samples == set(metadata["sample_name"])
