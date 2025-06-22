import loguru
import pandas as pd


def main(
    input_filepath,
    output_filepath,
):
    loguru.logger.info("Loading TCRs...")
    tcrs = pd.read_feather(input_filepath, columns=["sample_name", "templates"])
    loguru.logger.info("Counting UMIs...")
    template_counts = (
        tcrs.groupby("sample_name")["templates"]
        .median()
        .rename("filtered_median_umi_count")
    )
    loguru.logger.info("Saving median UMI counts...")
    template_counts.to_csv(output_filepath, sep="\t")
    loguru.logger.success("Done.")


if __name__ == "__main__":
    loguru.logger.remove()
    loguru.logger.add(snakemake.log[0])
    main(
        input_filepath=snakemake.input[0],
        output_filepath=snakemake.output[0],
    )
