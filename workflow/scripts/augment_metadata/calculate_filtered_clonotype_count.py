import loguru
import pandas as pd


def main(
    input_filepath,
    output_filepath,
):
    loguru.logger.info("Loading TCRs...")
    tcrs = pd.read_feather(input_filepath, columns=["sample_name"])
    loguru.logger.info("Counting clonotypes...")
    clonotype_counts = (
        tcrs.groupby("sample_name", observed=True)
        .size()
        .rename("filtered_clonotype_count")
    )
    loguru.logger.info("Saving clonotype counts...")
    clonotype_counts.to_csv(output_filepath, sep="\t")
    loguru.logger.success("Done.")


if __name__ == "__main__":
    loguru.logger.remove()
    loguru.logger.add(snakemake.log[0])
    main(
        input_filepath=snakemake.input[0],
        output_filepath=snakemake.output[0],
    )
