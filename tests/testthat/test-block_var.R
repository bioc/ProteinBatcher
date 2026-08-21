# Unit tests for the configurable block_var parameter in
# test_limma_customized() and .tl_fit_limma().

# Helper: build an SE with enough samples for duplicateCorrelation
# to estimate a non-degenerate correlation.  Uses 8 samples (4 per
# condition) with 4 blocks, each block appearing once per condition.
.make_se_with_block <- function(block_col, block_vals) {
    set.seed(42)
    n_proteins <- 20
    X <- matrix(rnorm(n_proteins * 8, mean = 10, sd = 2),
                nrow = n_proteins, ncol = 8)
    colnames(X) <- paste0("s", 1:8)
    rownames(X) <- paste0("p", seq_len(n_proteins))

    # Inject a condition effect so the contrast is estimable
    X[, 5:8] <- X[, 5:8] + 3

    cd_list <- list(
        condition = c("A", "A", "A", "A", "B", "B", "B", "B"),
        label     = paste0("s", 1:8)
    )
    cd_list[[block_col]] <- block_vals

    cd <- do.call(S4Vectors::DataFrame, cd_list)

    rd <- S4Vectors::DataFrame(
        name  = rownames(X),
        ID    = rownames(X),
        Genes = paste0("G", seq_len(n_proteins))
    )
    rownames(rd) <- rownames(X)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    imp_map <- matrix("none", nrow = n_proteins, ncol = 8,
                      dimnames = list(rownames(X), colnames(X)))
    S4Vectors::metadata(se)$imputation_map <- imp_map
    se
}

# ---------------------------------------------------------------------------
# block_var = "donor_id" (default) — works when donor_id present
# ---------------------------------------------------------------------------
test_that("block_effect=TRUE with block_var='donor_id' runs successfully", {
    if (!requireNamespace("limma", quietly = TRUE))
        testthat::skip("limma not installed.")

    se <- .make_se_with_block("donor_id", c("d1", "d2", "d3", "d4",
                                             "d1", "d2", "d3", "d4"))

    se_out <- ProteinBatcher::test_limma_customized(
        se = se,
        test = "B_vs_A",
        test_interaction = "NA",
        design_formula = ~ 0 + condition,
        ref_condition = "A",
        paired = FALSE,
        block_effect = TRUE,
        block_var = "donor_id"
    )

    expect_s4_class(se_out, "SummarizedExperiment")
    rd <- SummarizedExperiment::rowData(se_out)
    padj_cols <- grep("p\\.?adj|adj\\.p|FDR", colnames(rd),
                      value = TRUE, ignore.case = TRUE)
    expect_true(length(padj_cols) >= 1)
})

# ---------------------------------------------------------------------------
# block_var = custom column name (e.g. "batch")
# ---------------------------------------------------------------------------
test_that("block_effect=TRUE with custom block_var='batch' runs successfully", {
    if (!requireNamespace("limma", quietly = TRUE))
        testthat::skip("limma not installed.")

    se <- .make_se_with_block("batch", c("b1", "b2", "b3", "b4",
                                          "b1", "b2", "b3", "b4"))

    se_out <- ProteinBatcher::test_limma_customized(
        se = se,
        test = "B_vs_A",
        test_interaction = "NA",
        design_formula = ~ 0 + condition,
        ref_condition = "A",
        paired = FALSE,
        block_effect = TRUE,
        block_var = "batch"
    )

    expect_s4_class(se_out, "SummarizedExperiment")
})

# ---------------------------------------------------------------------------
# block_var missing from colData — errors cleanly
# ---------------------------------------------------------------------------
test_that("block_effect=TRUE errors when block_var column is absent", {
    if (!requireNamespace("limma", quietly = TRUE))
        testthat::skip("limma not installed.")

    # SE has donor_id but we ask for a non-existent column
    se <- .make_se_with_block("donor_id", c("d1", "d2", "d3", "d4",
                                             "d1", "d2", "d3", "d4"))

    expect_error(
        ProteinBatcher::test_limma_customized(
            se = se,
            test = "B_vs_A",
            test_interaction = "NA",
            design_formula = ~ 0 + condition,
            ref_condition = "A",
            paired = FALSE,
            block_effect = TRUE,
            block_var = "nonexistent_col"
        ),
        "Block variable.*missing"
    )
})

# ---------------------------------------------------------------------------
# block_var with only one level — errors cleanly
# ---------------------------------------------------------------------------
test_that("block_effect=TRUE errors when block_var has only one level", {
    if (!requireNamespace("limma", quietly = TRUE))
        testthat::skip("limma not installed.")

    se <- .make_se_with_block("donor_id", rep("d1", 8))

    expect_error(
        ProteinBatcher::test_limma_customized(
            se = se,
            test = "B_vs_A",
            test_interaction = "NA",
            design_formula = ~ 0 + condition,
            ref_condition = "A",
            paired = FALSE,
            block_effect = TRUE,
            block_var = "donor_id"
        ),
        "only has one level"
    )
})

# ---------------------------------------------------------------------------
# block_effect=FALSE ignores block_var entirely
# ---------------------------------------------------------------------------
test_that("block_effect=FALSE runs regardless of block_var value", {
    if (!requireNamespace("limma", quietly = TRUE))
        testthat::skip("limma not installed.")

    se <- .make_se_with_block("donor_id", c("d1", "d2", "d3", "d4",
                                             "d1", "d2", "d3", "d4"))

    # block_var points to a non-existent column, but block_effect=FALSE
    # so it should be ignored
    se_out <- ProteinBatcher::test_limma_customized(
        se = se,
        test = "B_vs_A",
        test_interaction = "NA",
        design_formula = ~ 0 + condition,
        ref_condition = "A",
        paired = FALSE,
        block_effect = FALSE,
        block_var = "nonexistent_col"
    )

    expect_s4_class(se_out, "SummarizedExperiment")
})
