"""Snakemake script for indexing a TCR dataset."""

import joblib
import pandas as pd
import tcrnn
import tcrnn.tcrdist


def main(input_ftr, normalize_dist, threads, n_neighbors, output):
    # Load the dataset. Convert the columns that TcrRepertoire uses into string
    # dtypes because "large_string[pyarrow]" does not understand str accessors.
    dtypes = {"v_gene": "string", "cdr3": "string"}
    tcr_df = pd.read_feather(input_ftr).astype(dtypes)
    if tcr_df["v_gene"].isna().any():
        raise ValueError("v_gene column contains missing values")
    tcr_repertoire = tcrnn.TcrRepertoire(
        tcr_df,
        tcrnn.tcrdist.TcrDist(normalize=bool(normalize_dist)),
        n_jobs=threads,
        progress=300.0,
    )
    # For now, use a constant random state to ensure reproducibility.
    tcr_repertoire.index(n_neighbors=n_neighbors, random_state=0)
    joblib.dump(tcr_repertoire, output)


if __name__ == "__main__":
    main(
        snakemake.input[0],
        int(snakemake.wildcards.normalize_dist),
        snakemake.threads,
        snakemake.params.n_neighbors,
        snakemake.output[0],
    )
