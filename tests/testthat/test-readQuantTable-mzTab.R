mztab_fixture <- system.file(
    "extdata", "example_pg_matrix.mzTab",
    package = "ProteinBatcher"
)
# When running tests from a source checkout (not an installed package),
# system.file() on the package itself returns "", so fall back to the path
# relative to the testthat working directory.
if (!nzchar(mztab_fixture)) {
    mztab_fixture <- testthat::test_path("..", "..", "inst", "extdata",
                                         "example_pg_matrix.mzTab")
}

test_that(".pp_is_mztab() detects mzTab files by extension and by content", {
    expect_true(.pp_is_mztab("foo.mzTab"))
    expect_true(.pp_is_mztab("foo.MZTAB"))
    expect_false(.pp_is_mztab("foo.pg_matrix.tsv"))
    expect_true(.pp_is_mztab(mztab_fixture))
})

test_that(".pp_read_mztab_pg_matrix() parses the PRH/PRT protein section", {
    skip_if_not(file.exists(mztab_fixture), "mzTab fixture not found")
    out <- .pp_read_mztab_pg_matrix(mztab_fixture)

    expect_s3_class(out, "data.frame")
    expect_true(all(c("Protein.Group", "Protein.Names", "Genes") %in%
                        colnames(out)))
    expect_equal(nrow(out), 3L)
    expect_equal(out$Protein.Group, c("P05412", "Q13835", "P05413"))
    # gene symbols recovered from the opt_global_gene_name column
    expect_equal(out$Genes, c("JUN", "PCDH1", "FABP3"))
})

test_that(".pp_read_mztab_pg_matrix() renames assay columns via ms_run location", {
    skip_if_not(file.exists(mztab_fixture), "mzTab fixture not found")
    out <- .pp_read_mztab_pg_matrix(mztab_fixture)

    sample_cols <- setdiff(colnames(out),
                           c("Protein.Group", "Protein.Names", "Genes"))
    expect_setequal(sample_cols,
                    c("HaCaT_IL13_rep1.raw", "HaCaT_NoTreated_rep1.raw"))
})

test_that(".pp_read_mztab_pg_matrix() maps mzTab's 'null' token to NA", {
    skip_if_not(file.exists(mztab_fixture), "mzTab fixture not found")
    out <- .pp_read_mztab_pg_matrix(mztab_fixture)

    fabp3 <- out[out$Protein.Group == "P05413", ]
    expect_true(is.na(fabp3[["HaCaT_IL13_rep1.raw"]]))
    expect_equal(fabp3[["HaCaT_NoTreated_rep1.raw"]], 9.77)
})

#' Write `lines` to a fresh temp file with the given extension and return
#' its path; avoids adding a hard dependency on the `withr` package.
#' @noRd
.write_temp_mztab <- function(lines, fileext = ".mzTab") {
    tmp <- tempfile(fileext = fileext)
    writeLines(lines, tmp)
    tmp
}

test_that(".pp_read_mztab_pg_matrix() falls back to accession without a gene column", {
    tmp <- .write_temp_mztab(c(
        "MTD\tmzTab-version\t1.0.0",
        "MTD\tms_run[1]-location\tfile:///data/run1.raw",
        "MTD\tassay[1]-ms_run_ref\tms_run[1]",
        "PRH\taccession\tdescription\tprotein_abundance_assay[1]",
        "PRT\tP12345\tSome protein\t7.5"
    ))
    out <- suppressMessages(.pp_read_mztab_pg_matrix(tmp))
    expect_equal(out$Genes, "P12345")
})

test_that(".pp_read_mztab_pg_matrix() errors without accession or abundance columns", {
    no_accession <- .write_temp_mztab(c(
        "PRH\tdescription\tprotein_abundance_assay[1]",
        "PRT\tSome protein\t7.5"
    ))
    expect_error(.pp_read_mztab_pg_matrix(no_accession), "accession")

    no_abundance <- .write_temp_mztab(c(
        "PRH\taccession\tdescription",
        "PRT\tP12345\tSome protein"
    ))
    expect_error(.pp_read_mztab_pg_matrix(no_abundance),
                 "protein_abundance_assay")
})

test_that("readQuantTable() dispatches to the mzTab reader for type = 'DIA'", {
    skip_if_not(file.exists(mztab_fixture), "mzTab fixture not found")
    out <- readQuantTable(mztab_fixture, type = "DIA")
    expect_true(all(c("Protein.Group", "Genes") %in% colnames(out)))
    expect_equal(nrow(out), 3L)
})
