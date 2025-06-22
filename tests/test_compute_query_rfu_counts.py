"""
Unit tests for
workflow.scripts.query_feature_generation.compute_query_rfu_counts.py.
"""

import pandas as pd
import polars as pl
import pytest

from workflow.scripts.query_feature_generation import compute_query_rfu_counts


@pytest.fixture
def query_tcrs():
    QUERY_TCRS = pd.DataFrame(
        [
            # Below query TCRs will have query size = 1.
            ("sample_1", "TRBV3-1*01", "CASSLSGNTIYF", 0, 1.0),
            ("sample_1", "TRBV18*01", "CASSPLDLSLNGYTF", 1, 60.0 / 62.0),
            ("sample_1", "TRBV4-1*01", "CASSFGTAQETQYF", 2, 60.0 / 62.0),
            # Below query TCRs are identical and hence will have query size = 2.
            # The weight is hence 1 since the removal of any single reference
            # sample will bring the total size above the cutoff.
            ("sample_1", "TRBV6-5*01", "CASSYSSRNTEAFF", 3, 1.0),
            ("sample_1", "TRBV6-5*01", "CASSYSSRNTEAFF", 3, 1.0),
        ],
        columns=[
            "sample_name",
            "v_gene",
            "cdr3",
            "centroid",
            "weight_at_min_size~3",
        ],
    )
    return QUERY_TCRS


MIN_SIZE = 3


def test_compute_query_tcr_weights(q42_tcrs, query_tcrs):
    """Test _compute_query_tcr_weights()

    The top q42 TCRs (md5sum: b86bff6dd35d06a11a149af6bbb843bf) by size are as
    follows. There are 62 samples and no sample has the same V gene and CDR3
    more than once. Asterisks mark the reference TCRs that are in `QUERY_TCRS`.

               v_gene              cdr3  size
    3016   TRBV3-1*01      CASSLSGNTIYF     3  *
    652     TRBV18*01   CASSPLDLSLNGYTF     2  *
    3366   TRBV4-1*01    CASSFGTAQETQYF     2  *
    4882   TRBV6-5*01    CASSYSSRNTEAFF     2
    519     TRBV15*01    CATSRDWDRDIQYF     2
    5838   TRBV7-9*01   CASSPTVAGDNEQFF     2
    3748   TRBV5-1*01   CASSLLGQGPYEQYF     2
    3955   TRBV5-4*01  CASSSKASGSTDTQYF     2
    2761    TRBV28*01      CASTLFGTEAFF     2
    1644  TRBV20-1*01   CSARGGTESNSPLHF     2

    At min_size 2, the size 3 TCR passes with the removal of any single
    reference sample. As such, the weight is 1.0. However, for size 2 TCRs,
    there are 2/62 reference samples that would bring the size below the cutoff.
    As such, the weight is 60/62.
    """
    # Sanity check that an error is raised if there are centroids in the query
    # TCRs that are not in the reference TCRs.
    NUM_REF_SAMPLES = len(q42_tcrs["sample_name"].unique())
    with pytest.raises(AssertionError):
        compute_query_rfu_counts._compute_query_tcr_weights(
            pl.from_pandas(query_tcrs),
            pl.from_pandas(q42_tcrs.assign(centroid=-1)),
            min_size=MIN_SIZE,
            num_ref_samples=NUM_REF_SAMPLES,
        )

    weighted_query_tcrs = compute_query_rfu_counts._compute_query_tcr_weights(
        pl.from_pandas(query_tcrs.assign(centroid=0)),
        pl.from_pandas(q42_tcrs.assign(centroid=0)),
        min_size=MIN_SIZE,
        num_ref_samples=NUM_REF_SAMPLES,
    ).collect()
    assert all(
        weighted_query_tcrs["tcr_weight"]
        == weighted_query_tcrs[f"weight_at_min_size~{MIN_SIZE}"]
    )


def test_filtered_query_rfu_count(query_tcrs, q42_tcrs, monkeypatch):
    """Test filtered_query_rfu_count()."""

    def mock_fn(tcrs_polars_df, *args):
        return tcrs_polars_df.rename(
            {f"weight_at_min_size~{MIN_SIZE}": "tcr_weight"}
        ).lazy()

    monkeypatch.setattr(
        compute_query_rfu_counts, "_compute_query_tcr_weights", mock_fn
    )
    result = compute_query_rfu_counts.filtered_query_rfu_count(
        pl.from_pandas(query_tcrs),
        pl.from_pandas(q42_tcrs.assign(centroid=range(len(q42_tcrs)))),
        min_size=3,
        num_ref_samples=len(q42_tcrs["sample_name"].unique()),
    )
    assert result.sort("centroid").equals(
        pl.from_pandas(query_tcrs)
        .rename({f"weight_at_min_size~{MIN_SIZE}": "tcr_weight"})
        .group_by("centroid")
        .agg(pl.col("tcr_weight").sum().alias("tcr_count"))
        .sort("centroid")
    )
