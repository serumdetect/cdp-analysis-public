import loguru
import pandas as pd


def main(
    filtered_counts_filepath,
    metadata_filepath,
    output_filepath,
):
    loguru.logger.info("Loading clonotype counts...")
    filtered_counts = pd.read_csv(
        filtered_counts_filepath, sep="\t", index_col="sample_name"
    )
    loguru.logger.info("Loading metadata...")
    metadata = pd.read_csv(metadata_filepath, sep="\t", index_col="sample_name")
    loguru.logger.info("Augmenting metadata...")
    metadata = metadata.join(filtered_counts, how="left")
    loguru.logger.info("Saving augmented metadata...")
    metadata.to_csv(output_filepath, sep="\t")
    loguru.logger.success("Done.")


if __name__ == "__main__":
    loguru.logger.remove()
    loguru.logger.add(snakemake.log[0])
    main(
        filtered_counts_filepath=snakemake.input.filtered_counts,
        metadata_filepath=snakemake.input.metadata,
        output_filepath=snakemake.output[0],
    )
