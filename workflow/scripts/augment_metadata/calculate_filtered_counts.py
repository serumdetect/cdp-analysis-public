import loguru
import polars as pl


def main(
    input_filepath,
    output_filepath,
):
    loguru.logger.info("Loading TCRs...")
    tcrs = pl.read_ipc(input_filepath)

    loguru.logger.info("Counting clonotypes...")
    clonotype_counts = tcrs.group_by("sample_name").agg(
        pl.len().alias("filtered_clonotype_count")
    )

    loguru.logger.info("Counting median UMI counts...")
    median_umi_counts = tcrs.group_by("sample_name").agg(
        pl.col("templates").median().alias("filtered_median_umi_count")
    )

    loguru.logger.info("Counting UMIs...")
    template_counts = tcrs.group_by("sample_name").agg(
        pl.col("templates").sum().alias("filtered_umi_count")
    )

    filtered_counts = clonotype_counts.join(
        median_umi_counts, on="sample_name", validate="1:1"
    ).join(template_counts, on="sample_name", validate="1:1")

    loguru.logger.info("Saving filtered counts...")
    filtered_counts.write_csv(output_filepath, separator="\t")
    loguru.logger.success("Done.")


if __name__ == "__main__":
    loguru.logger.remove()
    loguru.logger.add(snakemake.log[0])
    main(
        input_filepath=snakemake.input[0],
        output_filepath=snakemake.output[0],
    )
