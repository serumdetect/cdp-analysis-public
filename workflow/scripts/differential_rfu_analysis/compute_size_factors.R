#' Snakemake script for computing size factors for a given counts matrix.


library(scuttle)


main = function(counts_mat_tsv, output_tsv) {
    temp = read.table(
        counts_mat_tsv, header = T, row.names = 1, sep = "\t", check.names = F
    )
    counts_mat = as.matrix(temp)
    rownames(counts_mat) = rownames(temp)
    colnames(counts_mat) = colnames(temp)
    rm(temp)

    size_factors = setNames(
        scuttle::pooledSizeFactors(counts_mat), colnames(counts_mat)
    )

    write.table(
        data.frame(size_factors),
        output_tsv,
        col.names = F,
        row.names = T,
        sep = "\t",
        quote = F
    )
}


if (!interactive()) {
    main(
        counts_mat_tsv = snakemake@input[[1]],
        output_tsv = snakemake@output[[1]]
    )
}
