#' Unit tests for workflow/scripts/feature_generation/adjust_rfu_counts.R



library(testthat)


source("../workflow/scripts/feature_generation/adjust_rfu_counts.R")


# Prepare some data for running tests in this file.
COVARIATES = LETTERS[11:15]
FEATURES = LETTERS[1:10]
SAMPLES = letters[1:20]
RFU_COUNTS = matrix(
    0,
    nrow = length(FEATURES),
    ncol = length(SAMPLES),
    dimnames = list(FEATURES, SAMPLES)
)
QUERY_MODEL_MATRIX = REF_MODEL_MATRIX = matrix(
    0,
    nrow = length(SAMPLES),
    ncol = length(COVARIATES),
    dimnames = list(SAMPLES, COVARIATES)
)
BETA_MAT = matrix(
    0,
    nrow = length(FEATURES),
    ncol = length(COVARIATES),
    dimnames = list(FEATURES, COVARIATES)
)
SIZE_FACTORS = setNames(rep(1, length(SAMPLES)), SAMPLES)


test_that(
    "Test align_matrix_to_rows().",
    {
        # Missing values in the query matrix are replaced with the specified
        # fill value.
        result = align_matrix_to_rows(
            mat = matrix(
                1:6, nrow = 2, dimnames = list(letters[1:2], NULL)
            ),
            rownames = letters[1:4],
            0.0
        )
        expect_true(all(result[letters[3:4], ] == 0))
        expect_true(all(!is.na(result)))
    }
)


test_that(
    "Test .unify_factor_levels().",
    {
        # After running this function, target_df should have the same column
        # order and factor levels as source_df.
        source_df = data.frame(
            a = factor(c("a", "a", "a"), levels = c("a", "b", "c")),
            b = factor(c("d", "d", "d"), levels = c("d", "e", "f"))
        )
        target_df = data.frame(
            b = factor(c("e", "e", "e"), levels = c("f", "e", "d")),
            a = factor(c("b", "b", "b"), levels = c("c", "b", "a"))
        )
        result = .unify_factor_levels(target_df, source_df)
        expect_false(all(colnames(target_df) == colnames(source_df)))
        expect_true(all(colnames(result) == colnames(source_df)))
        expect_false(all(levels(target_df$a) == levels(source_df$a)))
        expect_true(all(levels(result$a) == levels(source_df$a)))
        expect_false(all(levels(target_df$b) == levels(source_df$b)))
        expect_true(all(levels(result$b) == levels(source_df$b)))
    }
)


test_that(
    "Test sanity_check_inputs()",
    {
        # First check that the example data returns no error.
        expect_no_error(
            sanity_check_inputs(
                RFU_COUNTS,
                QUERY_MODEL_MATRIX,
                REF_MODEL_MATRIX,
                BETA_MAT,
                SIZE_FACTORS
            )
        )

        # All of the reference matrix, query matrix, and beta matrix must have
        # the same columns (covariates). Make sure an error is raised otherwise.
        expect_warning(
            expect_error(
                sanity_check_inputs(
                    RFU_COUNTS,
                    QUERY_MODEL_MATRIX[, -1],
                    REF_MODEL_MATRIX,
                    BETA_MAT,
                    SIZE_FACTORS
                )
            ),
            "longer object length is not a multiple of shorter object length"
        )
        expect_warning(
            expect_error(
                sanity_check_inputs(
                    RFU_COUNTS,
                    QUERY_MODEL_MATRIX,
                    REF_MODEL_MATRIX[, -1],
                    BETA_MAT,
                    SIZE_FACTORS
                )
            ),
            "longer object length is not a multiple of shorter object length"
        )
        expect_error(
            sanity_check_inputs(
                RFU_COUNTS,
                QUERY_MODEL_MATRIX[, rev(seq_len(ncol(QUERY_MODEL_MATRIX)))],
                REF_MODEL_MATRIX,
                BETA_MAT,
                SIZE_FACTORS
            )
        )
        expect_error(
            sanity_check_inputs(
                RFU_COUNTS,
                QUERY_MODEL_MATRIX,
                REF_MODEL_MATRIX[, rev(seq_len(ncol(REF_MODEL_MATRIX)))],
                BETA_MAT,
                SIZE_FACTORS
            )
        )

        # An error should be raised if the row names RFU counts doesn't match
        # the row names of the beta matrix.
        expect_error(
            sanity_check_inputs(
                RFU_COUNTS,
                QUERY_MODEL_MATRIX,
                REF_MODEL_MATRIX,
                BETA_MAT[rev(seq_len(nrow(BETA_MAT))), ],
                SIZE_FACTORS
            ),
            "RFUs do not match between RFU counts and beta matrix."
        )

        # An error should be raised if the column names of RFU counts (samples)
        # does not match the row names of the query model matrix.
        expect_error(
            sanity_check_inputs(
                RFU_COUNTS[, rev(seq_len(ncol(RFU_COUNTS)))],
                QUERY_MODEL_MATRIX,
                REF_MODEL_MATRIX,
                BETA_MAT,
                SIZE_FACTORS
            ),
            paste0(
                "Sample names do not match between RFU counts and query model",
                " matrix."
            )
        )

        # Size factors must have the same length as RFU counts columns and the
        # sample names must match.
        expect_error(
            sanity_check_inputs(
                RFU_COUNTS,
                QUERY_MODEL_MATRIX,
                REF_MODEL_MATRIX,
                BETA_MAT,
                SIZE_FACTORS[-1]
            ),
            "Size factors do not match the number of samples."
        )
        expect_error(
            sanity_check_inputs(
                RFU_COUNTS,
                QUERY_MODEL_MATRIX,
                REF_MODEL_MATRIX,
                BETA_MAT,
                rev(SIZE_FACTORS)
            ),
            "Sample names do not match between size factors and RFU counts."
        )
    }
)


test_that(
    "Test .get_variables_to_adjust()",
    {
        expect_no_error(.get_variables_to_adjust(letters[1:5], "a"))
        expect_error(.get_variables_to_adjust(letters[1:5], "A"))
    }
)


test_that(
    "Test covariate_adjusted_counts()",
    {
        # Test that the function returns the expected output.
        unadjusted_counts = setNames(rep(1, 5), LETTERS[1:5])
        beta = setNames(c(1, 2, 3), letters[1:3])
        query_model_matrix = matrix(
            1:15,
            nrow = length(unadjusted_counts),
            dimnames = list(names(unadjusted_counts), names(beta))
        )

        # Make sure counts stay the same if nothing is to be adjusted for.
        result = covariate_adjusted_counts(
            unadjusted_counts, beta, query_model_matrix, letters[1:3]
        )
        expect_true(all(result == unadjusted_counts))

        # If covariate "b" is to be adjusted for, its coefficient is 2 and its
        # covariate values are 6:10. This gives effects of 2 * 6:10. After
        # correction by moving these fitted effects, we get adjusted counts of 1
        # / 2 ^ (2 * (6:10)).
        result = covariate_adjusted_counts(
            unadjusted_counts, beta, query_model_matrix, c("a", "c")
        )
        expect_true(all(result == 2 ^ (-2 * (6:10))))
    }
)
