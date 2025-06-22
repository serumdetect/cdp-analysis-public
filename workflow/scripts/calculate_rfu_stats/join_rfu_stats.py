import pandas as pd


def main(
    glm_stats_filepath,
    output_filepath,
    vj_fractions_filepath,
):
    glm_stats = pd.read_csv(glm_stats_filepath, sep="\t", index_col="rfu")
    vj_fractions = pd.read_csv(vj_fractions_filepath, sep="\t", index_col="rfu")
    rfu_stats = glm_stats.merge(
        vj_fractions, how="left", left_index=True, right_index=True
    )
    assert rfu_stats.columns.nunique() == rfu_stats.shape[1]
    rfu_stats.to_csv(output_filepath, sep="\t", index=True, compression="gzip")


if __name__ == "__main__":
    main(
        glm_stats_filepath=snakemake.input.glm_stats,
        output_filepath=snakemake.output[0],
        vj_fractions_filepath=snakemake.input.vj_fractions,
    )
