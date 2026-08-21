#' Detect whether a quantitation file is in mzTab format
#'
#' Detection is based on the file extension first (\code{.mzTab}, any case),
#' falling back to sniffing the first non-empty line for the \code{MTD} or
#' \code{COM} line-prefixes mandated by the mzTab specification. This lets
#' \code{.mztab} files without a matching extension still be recognised.
#'
#' @keywords internal
#' @noRd
.pp_is_mztab <- function(path) {
    if (grepl("\\.mztab$", path, ignore.case = TRUE)) return(TRUE)
    first_line <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                           error = function(e) "")
    if (!length(first_line)) return(FALSE)
    grepl("^(MTD|COM)\\t", first_line[1])
}

#' Parse the metadata (MTD) section of an mzTab file into a named list
#'
#' Each `MTD` line has the form `MTD\\t<key>\\t<value>`. Lines that do not
#' carry a value are skipped.
#'
#' @keywords internal
#' @noRd
.pp_parse_mztab_metadata <- function(lines) {
    mtd_lines <- lines[startsWith(lines, "MTD\t")]
    if (!length(mtd_lines)) return(list())
    parts <- strsplit(mtd_lines, "\t", fixed = TRUE)
    keys <- vapply(parts, function(x) x[2], character(1))
    vals <- vapply(parts, function(x) {
        if (length(x) >= 3 && nzchar(x[3])) x[3] else NA_character_
    }, character(1))
    keep <- !is.na(keys)
    stats::setNames(as.list(vals[keep]), keys[keep])
}

#' Derive a sample/file column name for each protein_abundance_assay[n] column
#'
#' mzTab reports quantitative values per *assay*, each of which is linked to
#' an \code{ms_run} via \code{assay[n]-ms_run_ref}, and each \code{ms_run} is
#' in turn linked to a raw file via \code{ms_run[m]-location}. ProteinBatcher
#' matches sample-annotation rows to the quantitation matrix by exact column
#' name (see \code{readExpDesign()}), so we recover the original file name
#' whenever the metadata section allows it, and fall back to a generic
#' \code{assay_<n>} label otherwise.
#'
#' @keywords internal
#' @noRd
.pp_mztab_assay_sample_names <- function(abundance_cols, mtd) {
    assay_n <- sub("^protein_abundance_assay\\[(\\d+)\\]$", "\\1",
                   abundance_cols)
    vapply(assay_n, function(i) {
        run_ref <- mtd[[paste0("assay[", i, "]-ms_run_ref")]]
        if (is.null(run_ref) || is.na(run_ref) || !nzchar(run_ref))
            return(paste0("assay_", i))
        run_id <- sub("^ms_run\\[(\\d+)\\]$", "\\1", trimws(run_ref))
        loc <- mtd[[paste0("ms_run[", run_id, "]-location")]]
        if (is.null(loc) || is.na(loc) || !nzchar(loc))
            return(paste0("assay_", i))
        basename(sub("^file://", "", loc))
    }, character(1), USE.NAMES = FALSE)
}

#' Recover a gene-symbol column from an mzTab protein (PRT) table
#'
#' The mzTab protein section has no mandatory gene-symbol column. We look for
#' the conventional optional columns used by proteomics mzTab exporters
#' (\code{opt_global_gene_name}, \code{opt_global_Genes}, ...), or a
#' non-standard \code{Genes} column, before falling back to the protein
#' accession, mirroring the fallback already used by \code{make_unique()}
#' for missing gene names.
#'
#' @keywords internal
#' @noRd
.pp_mztab_gene_column <- function(prt, header) {
    gene_like <- grep("^opt_.*_(gene_?name|genes)$|^Genes$", header,
                      value = TRUE, ignore.case = TRUE)
    if (length(gene_like)) return(prt[[gene_like[1]]])
    message(
        "mzTab protein section has no gene-symbol column (looked for ",
        "'opt_*_gene_name' / 'opt_*_genes' / 'Genes'); using the protein ",
        "accession as 'Genes' instead."
    )
    prt$accession
}

#' Read the protein (PRT) section of a proteomics mzTab file
#'
#' Implements a minimal reader for the \strong{mzTab (proteomics) 1.0/2.0}
#' format as specified by HUPO-PSI (metadata section \code{MTD}, protein
#' header \code{PRH} and protein rows \code{PRT}), returning a data.frame
#' shaped like the one produced by the \code{.pg_matrix.tsv} branch of
#' \code{readQuantTable()} (\code{Protein.Group}, \code{Protein.Names},
#' \code{Genes}, followed by one column per sample).
#'
#' This targets the classical, PSI-standard proteomics mzTab dialect, not
#' mzTab-M (the metabolomics/lipidomics variant, which uses a different
#' schema). Assay-level (\code{protein_abundance_assay[n]}) abundances are
#' required; study-variable summaries alone are not supported, as
#' ProteinBatcher needs per-sample values.
#'
#' No confirmed example of a DIA-NN-produced \code{.pg_matrix.mzTab} file
#' was available at implementation time, so this parser follows the public
#' HUPO-PSI specification directly rather than a vendor-specific dialect. If
#' a real export deviates from the standard (e.g. a non-standard
#' gene-symbol column, or abundances reported only as study variables),
#' this function will need a follow-up adjustment against an actual file.
#'
#' @param path Path to a \code{.mzTab} (or \code{.pg_matrix.mzTab}) file.
#' @keywords internal
#' @noRd
.pp_read_mztab_pg_matrix <- function(path) {
    lines <- readLines(path, warn = FALSE)
    mtd <- .pp_parse_mztab_metadata(lines)

    prh_idx <- which(startsWith(lines, "PRH\t"))
    if (!length(prh_idx))
        stop("No protein section (PRH/PRT lines) found in mzTab file: ",
             path, call. = FALSE)
    header <- strsplit(lines[prh_idx[1]], "\t", fixed = TRUE)[[1]][-1]

    prt_lines <- lines[startsWith(lines, "PRT\t")]
    if (!length(prt_lines))
        stop("mzTab file has a PRH header but no PRT rows: ", path,
             call. = FALSE)

    rows <- lapply(strsplit(prt_lines, "\t", fixed = TRUE), function(x) {
        x <- x[-1]
        if (length(x) < length(header))
            x <- c(x, rep(NA_character_, length(header) - length(x)))
        x[seq_along(header)]
    })
    prt <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    colnames(prt) <- header
    prt[prt == "null"] <- NA # mzTab's own null token

    if (!"accession" %in% colnames(prt))
        stop("mzTab protein section is missing the required 'accession' ",
             "column.", call. = FALSE)

    abundance_cols <- grep("^protein_abundance_assay\\[\\d+\\]$", header,
                           value = TRUE)
    if (!length(abundance_cols))
        stop(
            "mzTab protein section has no 'protein_abundance_assay[n]' ",
            "columns. ProteinBatcher requires assay-level (per-sample) ",
            "protein abundances; study_variable-only summaries are not ",
            "currently supported.", call. = FALSE
        )

    sample_names <- .pp_mztab_assay_sample_names(abundance_cols, mtd)
    genes <- .pp_mztab_gene_column(prt, header)

    abundances <- as.data.frame(
        lapply(prt[abundance_cols], function(x) suppressWarnings(as.numeric(x))),
        stringsAsFactors = FALSE
    )
    colnames(abundances) <- sample_names

    cbind(
        data.frame(
            Protein.Group = prt$accession,
            Protein.Names = if ("description" %in% colnames(prt))
                prt$description else NA_character_,
            Genes = genes,
            stringsAsFactors = FALSE
        ),
        abundances
    )
}
