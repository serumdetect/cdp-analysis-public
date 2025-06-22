"""Helper functions for loading test data."""

import importlib.resources

import pandas as pd
import pytest


@pytest.fixture
def q42_tcrs():
    """Return example data from Q42."""
    fname = importlib.resources.files("tests").joinpath(
        "data", "Q0042_example_data.tsv.gz"
    )
    return pd.read_table(fname, header=0)
