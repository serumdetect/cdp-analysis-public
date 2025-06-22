"""Unit tests for workflow/scripts/prepare_tcr_data/preprocess_tcr_data.py."""

import pandas as pd
import pytest

from workflow.scripts.prepare_tcr_data import preprocess_tcr_data


def test_filter_by_blacklisted_tcrs(q42_tcrs):
    """Test filtering TCR data by blacklisted TCRs."""
    usecols = ["v_gene", "j_gene", "cdr3_nt"]
    q42_tcrs = q42_tcrs.drop_duplicates(subset=usecols)
    blacklisted_tcrs = q42_tcrs.sample(10)[usecols]
    data = preprocess_tcr_data.filter_by_blacklisted_tcrs(
        q42_tcrs, blacklisted_tcrs
    )
    assert len(data) == len(q42_tcrs) - len(blacklisted_tcrs)


@pytest.mark.parametrize(
    ["assigned_seq_run", "seq_runs"],
    [
        ("SEQ1", ["SEQ2"]),
        ("SEQ1", ["SEQ2", "SEQ3"]),
        ("SEQ1", ["SEQ1"]),
        ("SEQ1", ["SEQ1", "SEQ3"]),
    ],
)
def filter_by_seq_runs(assigned_seq_run, seq_runs, q42_tcrs):
    """Test filtering TCR data by blacklisted seq runs."""
    metadata = pd.DataFrame(
        {
            "sample_name": q42_tcrs["sample_name"].unique(),
            "seq_run": assigned_seq_run,
        }
    )
    filtered_data, _ = preprocess_tcr_data.filter_by_seq_runs(
        q42_tcrs, metadata, seq_runs
    )
    if assigned_seq_run in seq_runs:
        assert len(filtered_data) == len(q42_tcrs)
    else:
        assert len(filtered_data) == 0
