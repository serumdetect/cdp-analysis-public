"""Snakemake script for computing CDP decision graph."""

import joblib
import pandas as pd
import tcrnn.cfsfdp


def main(input_sav, d_c, output_ftr):
    # Load the TCR index.
    tcr_index = joblib.load(input_sav)

    # Compute density.
    tcr_index.progress = False
    tcr_cdp = tcrnn.cfsfdp.NnCfsfdp(tcr_index, d_c=d_c)
    tcr_cdp_df = pd.concat(
        [
            tcr_cdp.density.index.to_frame(index=False),
            tcr_cdp.density.rename("density").reset_index(drop=True),
        ],
        axis=1,
    )
    tcr_cdp_df.to_feather(output_ftr)


if __name__ == "__main__":
    main(
        snakemake.input[0],
        float(snakemake.wildcards.d_c),
        snakemake.output[0],
    )
