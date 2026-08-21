# Unit tests for optional annotation columns (batch, donor_id, replicate).
# After v0.99.7 only file, sample, sample_name, and condition are required;
# the remaining columns are optional and handled gracefully throughout.

# ---------------------------------------------------------------------------
# .pp_validate_inputs — minimal header accepted
# ---------------------------------------------------------------------------
test_that(".pp_validate_inputs accepts annotation with only 4 required columns", {
    tmp <- tempfile(fileext = ".tsv"); writeLines("dummy", tmp)

    ann_min <- tempfile(fileext = ".tsv")
    min_header <- paste(
        c("file", "sample", "sample_name", "condition"),
        collapse = "\t"
    )
    writeLines(c(min_header, "f1\ts1\tA_1\tA"), ann_min)

    out <- ProteinBatcher:::.pp_validate_inputs(
        path_pgmatrix = tmp, path_annotation = ann_min, level = "protein",
        type = "DIA", percent_missing = 50, path_output = tempdir(),
        experiment = "EXP"
    )
    expect_type(out, "list")
    expect_identical(out$experiment, "EXP")
})

test_that(".pp_validate_inputs still rejects annotation missing a required column", {
    tmp <- tempfile(fileext = ".tsv"); writeLines("dummy", tmp)

    # Missing 'condition' — only 3 of 4 required columns
    ann_bad <- tempfile(fileext = ".tsv")
    writeLines(c("file\tsample\tsample_name", "f1\ts1\tA_1"), ann_bad)

    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, ann_bad, "protein", "DIA", 50, tempdir(), "EXP"),
        "missing required columns"
    )
})

# ---------------------------------------------------------------------------
# readExpDesign — replicate placeholder + label from file
# ---------------------------------------------------------------------------
test_that("readExpDesign creates replicate placeholder when absent", {
    ann <- tempfile(fileext = ".tsv")
    df <- data.frame(
        file = c("f1", "f2"),
        sample = c("s1", "s2"),
        sample_name = c("A_1", "B_1"),
        condition = c("A", "B"),
        stringsAsFactors = FALSE
    )
    utils::write.table(df, ann, sep = "\t", row.names = FALSE, quote = FALSE)

    out <- ProteinBatcher:::readExpDesign(ann, type = "DIA")
    expect_s3_class(out, "data.frame")
    expect_true("replicate" %in% colnames(out))
    expect_true(all(is.na(out$replicate)))
})

test_that("readExpDesign always creates label from file", {
    ann <- tempfile(fileext = ".tsv")
    df <- data.frame(
        file = c("f1", "f2"),
        sample = c("s1", "s2"),
        sample_name = c("A_1", "B_1"),
        condition = c("A", "B"),
        stringsAsFactors = FALSE
    )
    utils::write.table(df, ann, sep = "\t", row.names = FALSE, quote = FALSE)

    out <- ProteinBatcher:::readExpDesign(ann, type = "DIA")
    expect_true("label" %in% colnames(out))
    expect_identical(out$label, out$file)
})

test_that("readExpDesign preserves replicate when present", {
    ann <- tempfile(fileext = ".tsv")
    df <- data.frame(
        file = c("f1", "f2"),
        sample = c("s1", "s2"),
        sample_name = c("A_1", "B_1"),
        condition = c("A", "B"),
        replicate = c(1L, 2L),
        stringsAsFactors = FALSE
    )
    utils::write.table(df, ann, sep = "\t", row.names = FALSE, quote = FALSE)

    out <- ProteinBatcher:::readExpDesign(ann, type = "DIA")
    expect_true("replicate" %in% colnames(out))
    expect_identical(as.character(out$replicate), c("1", "2"))
})

# ---------------------------------------------------------------------------
# make_se_customized — works without replicate column
# ---------------------------------------------------------------------------
test_that("make_se_customized works with only label, condition, and sample_name", {
    proteins <- data.frame(
        name = c("P1", "P2"),
        ID = c("P1", "P2"),
        s1 = c(10, 20),
        s2 = c(11, 21),
        stringsAsFactors = FALSE
    )
    cols <- which(colnames(proteins) %in% c("s1", "s2"))

    expdesign <- data.frame(
        label = c("s1", "s2"),
        sample_name = c("A_1", "B_1"),
        condition = c("A", "B"),
        # replicate intentionally absent
        stringsAsFactors = FALSE
    )

    se <- ProteinBatcher:::make_se_customized(
        proteins_unique = proteins,
        columns = cols,
        expdesign = expdesign,
        log2transform = TRUE
    )

    expect_s4_class(se, "SummarizedExperiment")
    cd <- SummarizedExperiment::colData(se)
    expect_true("replicate" %in% colnames(cd))
    expect_true("condition" %in% colnames(cd))
})

# ---------------------------------------------------------------------------
# .tl_check_inputs — relaxed colData check
# ---------------------------------------------------------------------------
test_that(".tl_check_inputs accepts colData with only label and condition", {
    X <- matrix(rnorm(2 * 4), nrow = 2)
    colnames(X) <- paste0("s", 1:4)
    rownames(X) <- paste0("p", 1:2)

    cd <- S4Vectors::DataFrame(
        label     = paste0("s", 1:4),
        condition = c("A", "A", "B", "B")
        # batch, donor_id, replicate intentionally absent
    )

    rd <- S4Vectors::DataFrame(name = c("p1", "p2"))
    rownames(rd) <- rownames(X)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    imp_map <- matrix("none", nrow = 2, ncol = 4,
                      dimnames = list(rownames(X), colnames(X)))
    S4Vectors::metadata(se)$imputation_map <- imp_map

    # Should not error — only label and condition are required now
    expect_invisible(
        ProteinBatcher:::.tl_check_inputs(
            se, ~ 0 + condition, "A", "B_vs_A", "NA"
        )
    )
})
