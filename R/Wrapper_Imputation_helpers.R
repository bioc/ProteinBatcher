#' Validate inputs for impute_se
#' @keywords internal
#' @noRd
.impute_validate_inputs <- function(se, threshold,
                                    ldv_source = c("global",
                                                   "per-condition")) {
    if (!methods::is(se, "SummarizedExperiment"))
        stop("`se` must be aSummarizedExperiment.")
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || threshold <= 0 || threshold > 1)
        stop("`threshold` must be a single numeric in (0, 1].")
    ldv_source <- match.arg(ldv_source)
    if (!("condition" %in% colnames(SummarizedExperiment::colData(se))))
        stop("`colData(se)$condition` is required.")
    X <- SummarizedExperiment::assay(se)
    if (!is.numeric(X)) stop("`assay(se)` must be numeric.")
    list(se = se, threshold = threshold, ldv_source = ldv_source)
}

#' Exclude samples with NA condition (message emitted)
#' @keywords internal
#' @noRd
.impute_exclude_na_condition <- function(X, cd) {
    cond <- as.character(cd$condition)
    if (anyNA(cond)) {
        keep <- !is.na(cond)
        message(sprintf("Excluding %d sample(s) with NA condition.",
                        sum(!keep)))
        return(list(X = X[, keep, drop = FALSE], cond = cond[keep],
                    dropped = TRUE))
    }
    list(X = X, cond = cond, dropped = FALSE)
}

#' Step 1: within-condition mean imputation (< threshold)
#' @keywords internal
#' @noRd
.impute_mean_within_condition <- function(X, cond, threshold, imp_map) {
    X1 <- X
    for (cn in unique(cond)) {
        idx <- which(cond == cn); if (!length(idx)) next
        sub <- X[, idx, drop = FALSE]
        prop_na <- rowMeans(is.na(sub))
        need <- prop_na > 0 & prop_na < threshold
        if (any(need)) {
            mu <- rowMeans(sub, na.rm = TRUE)
            for (i in which(need)) {
                miss <- idx[which(is.na(X[i, idx]))]
                if (length(miss)) {
                    X1[i, miss] <- mu[i]
                    imp_map[i, miss] <- "mean"
                }
            }
        }
    }
    list(X = X1, map = imp_map)
}

#' Compute LDV mean(s)
#' @keywords internal
#' @noRd
.impute_ldv_mu <- function(X1, cond, ldv_source) {
    if (ldv_source == "global") return(list(type = "global",
                                            mu = min(X1,na.rm = TRUE)))
    vals <- vapply(unique(cond), function(cn) {
        idx <- which(cond == cn); min(X1[, idx, drop = FALSE], na.rm = TRUE)
    }, numeric(1))
    names(vals) <- unique(cond)
    list(type = "per_condition", mu = vals)
}

#' Per-protein SDs with robust fallback
#' @keywords internal
#' @noRd
.impute_sd_protein <- function(X, sd_fallback = NULL) {
    sdv <- apply(X, 1L, stats::sd, na.rm = TRUE)
    if (is.null(sd_fallback)) {
        sd_fallback <- stats::median(sdv[is.finite(sdv) & sdv > 0],
                                     na.rm = TRUE)
        if (!is.finite(sd_fallback)) sd_fallback <- 0.1
    }
    sdv[!is.finite(sdv) | sdv <= 0] <- sd_fallback
    list(sd = sdv, fallback = sd_fallback)
}

#' Step 2: LDV imputation for remaining NAs
#' @keywords internal
#' @noRd
.impute_ldv_fill <- function(X1, cond, ldv, sd_prot, ldv_source, lower_bound,
                             imp_map) {
    X2 <- X1
    na_pos <- which(is.na(X2), arr.ind = TRUE)
    if (nrow(na_pos)) {
        for (k in seq_len(nrow(na_pos))) {
            i <- na_pos[k, 1]; j <- na_pos[k, 2]
            mu_j <- if (ldv$type == "global") ldv$mu else ldv$mu[cond[j]]
            draw <- stats::rnorm(1L, mean = mu_j, sd = sd_prot[i])
            if (!is.null(lower_bound)) draw <- max(lower_bound, draw)
            X2[i, j] <- draw; imp_map[i, j] <- "ldv"
        }
    }
    list(X = X2, map = imp_map)
}

#' Write back assay, expand map if needed, and annotate SE
#' @keywords internal
#' @noRd
.impute_write_back <- function(se, X0, X_imputed, imp_map, dropped) {
    se_out <- se
    if (!dropped) {
        SummarizedExperiment::assay(se_out) <- X_imputed
    } else {
        full <- X0; full[, colnames(X_imputed)] <- X_imputed
        SummarizedExperiment::assay(se_out) <- full
        full_map <- matrix("none", nrow = nrow(imp_map), ncol = ncol(X0),
                           dimnames = list(rownames(imp_map), colnames(X0)))
        full_map[, colnames(X_imputed)] <- imp_map
        imp_map <- full_map
    }
    S4Vectors::metadata(se_out)$imputation_map <- imp_map
    SummarizedExperiment::rowData(se_out)$imputed <-
        apply(imp_map, 1L, function(x) any(x == "mean" | x == "ldv"))
    se_out
}

#' @keywords internal
#' @noRd
.pp_validate_inputs <- function(path_pgmatrix, path_annotation,
                                level, type, percent_missing,
                                path_output, experiment,
                                annotation_format = c("standard", "sdrf"),
                                sdrf_map = NULL)
{
    annotation_format <- match.arg(annotation_format)

    if (!file.exists(path_pgmatrix))
        stop("path_pgmatrix does not exist.")
    if (!file.exists(path_annotation))
        stop("path_annotation does not exist.")
    if (!is.numeric(percent_missing) || percent_missing <= 0 ||
        percent_missing > 100)
        stop("percent_missing must be 1-100.")
    if (!dir.exists(path_output))
        stop("path_output folder does not exist.")
    if (!is.character(experiment) || experiment == "")
        stop("experiment must be a non-empty string.")

    # Read ONLY the header (fast) and try common delimiters
    # 1) tab (typical for .tsv/.txt)
    ann_head <- tryCatch(
        utils::read.delim(path_annotation, sep = "\t", header = TRUE, nrows = 1,
                          quote = "", comment.char = "", check.names = FALSE,
                          stringsAsFactors = FALSE),
        error = function(e) NULL
    )
    # If it failed or produced a single column, try comma
    if (is.null(ann_head) || ncol(ann_head) <= 1) {
        ann_head <- tryCatch(
            utils::read.csv(path_annotation, sep = ",",header = TRUE, nrows = 1,
                            quote = "", comment.char = "", check.names = FALSE,
                            stringsAsFactors = FALSE),
            error = function(e) NULL
        )
    }
    # Last attempt: semicolon (common in European Excel exports)
    if (is.null(ann_head) || ncol(ann_head) <= 1) {
        ann_head <- tryCatch(
            utils::read.csv(path_annotation, sep = ";",header = TRUE, nrows = 1,
                            quote = "", comment.char = "", check.names = FALSE,
                            stringsAsFactors = FALSE),
            error = function(e) NULL
        )
    }
    if (is.null(ann_head))
        stop("Could not read path_annotation header
                (unknown delimiter or malformed file).")

    # Normalize column names (trim whitespace)
    coln <- trimws(colnames(ann_head))

    # Detect duplicate column names
    if (anyDuplicated(coln)) {
        dups <- unique(coln[duplicated(coln)])
        stop("path_annotation has duplicated column names: ",
             paste(dups, collapse = ", "))
    }

    if (annotation_format == "standard") {
        required_cols <- c("file","sample","sample_name","condition")
        missing_cols <- setdiff(required_cols, coln)
        if (length(missing_cols) > 0) {
            stop(
                "path_annotation is missing required columns: ",
                paste(missing_cols, collapse = ", "),
                ". Found columns: ",
                paste(coln, collapse = ", ")
            )
        }
    } else {
        if (!"comment[data file]" %in% coln) {
            stop("SDRF annotation is missing the 'comment[data file]' column.",
                 " Found columns: ", paste(coln, collapse = ", "))
        }
        if (!any(c("source name", "assay name") %in% coln)) {
            stop("SDRF annotation is missing 'source name' or 'assay name'.",
                 " Found columns: ", paste(coln, collapse = ", "))
        }
        if (is.null(sdrf_map) || is.null(sdrf_map$condition)) {
            stop("For annotation_format = 'sdrf' you must supply ",
                 "sdrf_map = list(condition = \"factor value[...]\", ...) ",
                 "indicating which SDRF column plays the role of 'condition'.")
        }
    }
    list(path_pgmatrix = path_pgmatrix,
         path_annotation = path_annotation,
         level = level,
         type = type,
         percent_missing = percent_missing,
         path_output = path_output,
         experiment = experiment,
         annotation_format = annotation_format,
         sdrf_map = sdrf_map)
}

#' @keywords internal
#' @noRd
.pp_import_se <- function(path_pgmatrix, path_annotation, level, type,
                          annotation_format = "standard", sdrf_map = NULL) {
    se <- make_se_from_files(
        path_pgmatrix,
        path_annotation,
        level = level,
        type = type,
        annotation_format = annotation_format,
        sdrf_map = sdrf_map
    )
    if (!methods::is(se, "SummarizedExperiment"))
        stop("make_se_from_files did not return a SummarizedExperiment.")
    se
}

#' @keywords internal
#' @noRd
.pp_write_filtered <- function(removed_df, path_output, experiment){
    out <- file.path(path_output,
                     paste0(experiment, "_filtered_out_proteins.csv"))
    utils::write.csv(removed_df, out, row.names = FALSE)
    invisible(out)
}

#' @keywords internal
#' @noRd
.pp_write_before_after <- function(se_before, se_after,
                                   path_output, experiment){
    before <- cbind(gene = SummarizedExperiment::rowData(se_before)$Genes,
                    as.data.frame(SummarizedExperiment::assay(se_before)))

    after  <- cbind(gene = SummarizedExperiment::rowData(se_after)$Genes,
                    as.data.frame(SummarizedExperiment::assay(se_after)))

    out_before <- file.path(path_output,
                            paste0(experiment, "_Imputation_before.csv"))
    out_after  <- file.path(path_output,
                            paste0(experiment, "_Imputation_after.csv"))

    utils::write.csv(before, out_before, row.names = FALSE)
    utils::write.csv(after,  out_after,  row.names = FALSE)
    invisible(c(before = out_before, after = out_after))
}

#' @keywords internal
#' @noRd
readQuantTable <- function (quant_table_path, type = "TMT", level = NULL,
                            log2transform = FALSE, exp_type = NULL,
                            additional_cols = NULL)
{
    # Dispatch to the mzTab reader when the file is in mzTab format.
    # The mzTab parser returns a data.frame with the same shape as the
    # DIA pg_matrix branch (Protein.Group, Protein.Names, Genes, + samples).
    if (type == "DIA" && .pp_is_mztab(quant_table_path)) {
        return(.pp_read_mztab_pg_matrix(quant_table_path))
    }
    temp_data <- utils::read.table(quant_table_path, header = TRUE,
                                   fill = TRUE, sep = "\t", quote = "",
                                   comment.char = "", blank.lines.skip = FALSE,
                                   check.names = FALSE)
    colnames(temp_data) <- unique(colnames(temp_data), "_")
    if(type == "DIA"){
        temp <- data.table::melt.data.table(
            data.table::setDT(
                temp_data[,!colnames(temp_data) %in%
                              c("Proteotypic","Precursor.Charge","Precursor.Id",
                                "Modified.Sequence","First.Protein.Description",
                                "All Mapped Proteins", "All Mapped Genes")]),
            id.vars = c("Protein.Group", "Protein.Names", "Genes"),
            variable.name = "File.Name")
        temp_data <- as.data.frame(
            data.table::dcast.data.table(
                temp, Protein.Group + Protein.Names + Genes ~ File.Name,
                value.var = "value",
                fun.aggregate = function(x){if (all(is.na(x))){NA_real_
                }else max(x, na.rm = TRUE)}))
    }else{stop("This current type is not supported.")}
    return(temp_data)
}

#' @keywords internal
#' @noRd
.map_sdrf_columns <- function(df, sdrf_map = NULL) {
    cn <- colnames(df)
    if ("comment[data file]" %in% cn) {
        colnames(df)[cn == "comment[data file]"] <- "file"
        cn <- colnames(df)
    }
    anchor <- intersect(c("source name", "assay name"), cn)
    if (!"sample_name" %in% cn) {
        if (length(anchor) == 0) {
            stop("SDRF file must contain 'source name' or 'assay name'.",
                 call. = FALSE)
        }
        colnames(df)[cn == anchor[1]] <- "sample_name"
        cn <- colnames(df)
    }
    if (!"sample" %in% cn) { df$sample <- df$sample_name; cn <- colnames(df) }

    if (!is.null(sdrf_map)) {
        for (target in names(sdrf_map)) {
            src <- sdrf_map[[target]]
            if (!src %in% cn) {
                stop("SDRF column '", src, "' (mapped to '", target,
                     "') was not found.", call. = FALSE)
            }
            colnames(df)[cn == src] <- target
            cn <- colnames(df)
        }
    }
    if (!"condition" %in% colnames(df)) {
        stop("'condition' could not be resolved. Provide sdrf_map = ",
             "list(condition = \"factor value[...]\").", call. = FALSE)
    }
    df
}

#' @keywords internal
#' @noRd
readExpDesign <- function (exp_anno_path, type = "TMT", lfq_type = "Intensity",
                           lowercase = FALSE, format = c("standard", "sdrf"),
                           sdrf_map = NULL)
{
    format <- match.arg(format)
    temp_df <- utils::read.table(exp_anno_path, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE, check.names = FALSE)
    if (format == "sdrf") temp_df <- .map_sdrf_columns(temp_df, sdrf_map)
    if (type == "DIA") {
        if (lowercase) {
            colnames(temp_df) <- tolower(colnames(temp_df))
        }
        if (!"file" %in% colnames(temp_df)) {
            stop("'file' column is not present in the experiment annotation
                    file.")
        }
        if (!"sample_name" %in% colnames(temp_df)) {
            stop("'sample_name' column is not present in the experiment
                    annotation file.")
        }
        if (!"sample" %in% colnames(temp_df)) {
            stop("'sample' column is not present in the experiment
                    annotation file.")
        }
        if (!is.character(temp_df$sample_name)) {
            temp_df$sample_name <- as.character(temp_df$sample_name)
        }
        if (anyDuplicated(temp_df$sample_name)) {
            stop("Duplicated 'sample_name' detected in the experiment
                    annotation file.")
        }
        if (anyDuplicated(temp_df$sample)) {
            stop("Duplicated 'sample' detected in the experiment
                    annotation file.")
        }
        temp_df$condition <- make.names(temp_df$condition)
        # 'replicate' is optional — create a placeholder if absent so
        # downstream code (make_se_customized, test_limma_customized) does
        # not crash on colData access.
        if (!"replicate" %in% colnames(temp_df)) {
            temp_df$replicate <- NA_character_
        }
        # 'label' is required by make_se_customized for row matching.
        # Always create it from 'file'.
        temp_df$label <- temp_df$file
    }else{
        stop("This type is currently not supported.")
    }
    return(temp_df)
}

#' @keywords internal
#' @noRd
make_unique <- function (proteins, names, ids, delim = ";")
{
    col_names <- colnames(proteins)
    if (!names %in% col_names) {
        stop("'", names,"' is not a column in '", deparse(substitute(proteins)),
             "'", call. = FALSE)
    }
    if (!ids %in% col_names) {
        stop("'", ids, "' is not a column in '", deparse(substitute(proteins)),
             "'", call. = FALSE)
    }
    if (inherits(proteins, "tbl_df")) {
        proteins <- as.data.frame(proteins)
    }
    double_NAs <- apply(proteins[, c(names, ids)], 1, function(x) all(is.na(x)))
    if (any(double_NAs)) {
        stop("NAs in both the 'names' and 'ids' columns")
    }
    proteins_unique <- proteins
    proteins_unique$name <- proteins_unique[[names]]
    proteins_unique$ID <- proteins_unique[[ids]]
    nm <- proteins_unique$name
    use_id <- is.na(nm) | nm == ""
    nm[use_id] <- proteins_unique$ID[use_id]
    proteins_unique$name <- make.unique(nm)
    return(proteins_unique)
}

#' @keywords internal
#' @noRd
make_se_customized <- function (proteins_unique, columns, expdesign,
                                log2transform = FALSE, exp = "LFQ",
                                lfq_type = NULL, level = NULL, exp_type = NULL)
{
    if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
        stop("'name' and/or 'ID' columns are not present in '",
             deparse(substitute(proteins_unique)), "'.\nRun make_unique() to
                obtain the required columns",
             call. = FALSE)
    }
    if (any(!c("label", "condition") %in% colnames(expdesign))) {
        stop("'label' and/or 'condition' columns",
             "are not present in the experimental design", call. = FALSE)
    }
    # 'replicate' is optional — create a placeholder if absent
    if (!"replicate" %in% colnames(expdesign)) {
        expdesign$replicate <- NA_character_
    }
    if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
        stop("specified 'columns' should be numeric", "\nRun
                make_se_customized() with the appropriate columns as argument",
             call. = FALSE)
    }
    if (inherits(proteins_unique, "tbl_df")) {
        proteins_unique <- as.data.frame(proteins_unique)
    }
    if (inherits(expdesign, "tbl_df")) {
        expdesign <- as.data.frame(expdesign)
    }
    rownames(proteins_unique) <- proteins_unique$ID
    raw <- proteins_unique[, columns]
    raw[raw == 0] <- NA
    if (log2transform) {
        raw <- log2(raw)
    }
    rownames(expdesign) <- expdesign$label
    matched <- match(make.names(expdesign$label), make.names(colnames(raw)))
    if (any(is.na(matched))) {
        print(make.names(expdesign$label))
        print(make.names(colnames(raw)))
        stop("None of the labels in the experimental design match ",
             "with column names in 'proteins_unique'", "\nRun make_se() with
                the correct labels in the experimental design", "and/or correct
                columns specification")
    }
    rownames(expdesign) <- expdesign$sample_name
    colnames(raw)[matched] <- expdesign$sample_name
    raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
    row_data <- proteins_unique[, -columns]
    rownames(row_data) <- row_data$ID
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensity = as.matrix(raw)),
        colData = expdesign, rowData = row_data,
        metadata = list(log2transform = log2transform, exp = exp,
                        lfq_type = lfq_type, exp_type = exp_type,
                        level = level))
    if (exp == "DIA" & !is.null(level)) {
        if (level == "protein") {
            SummarizedExperiment::rowData(se)$Index <- rowData(se)$Protein.Group
        }
    }
    return(se)
}

#' @keywords internal
#' @noRd
make_se_from_files <- function (quant_table_path, exp_anno_path, type = "TMT",
                                level = NULL, exp_type = NULL,
                                log2transform = NULL, lfq_type = "Intensity",
                                gencode = FALSE, additional_cols = NULL,
                                annotation_format = c("standard", "sdrf"),
                                sdrf_map = NULL) {
    annotation_format <- match.arg(annotation_format)
    if (type == "DIA" & is.null(log2transform)) {
        log2transform <- TRUE
    }else if (is.null(level)) {
        log2transform <- FALSE
    }
    if (!level %in% c("gene", "protein", "peptide", "site", "glycan")) {
        cat(paste0("The specified level: ", level, " is not a valid level.
                    Available level is protein.\n"))
        return(NULL)
    }
    quant_table <- readQuantTable(quant_table_path, type = type,
                                  level = level, exp_type = exp_type,
                                  additional_cols = additional_cols)
    exp_design <- readExpDesign(exp_anno_path, type = type, lfq_type = lfq_type,
                                format = annotation_format, sdrf_map = sdrf_map)
    if (type == "DIA") {
        if (level == "protein") {
            if (gencode) {
                quant_table <- quant_table[grepl("^ENS",
                                                 quant_table$Protein.Group),
                ]
            }
            data_unique <- make_unique(quant_table, "Genes", "Protein.Group")
            cols <- colnames(data_unique)
            if (is.null(additional_cols)) {
                selected_cols <- which(!(cols %in%
                                             c("Protein.Group", "Protein.Ids",
                                               "Protein.Names", "Genes",
                                               "First.Protein.Description",
                                               "ID", "name")))
            }else {
                selected_cols <- which(!(cols %in%
                                             c("Protein.Group", "Protein.Ids",
                                               "Protein.Names", "Genes",
                                               "First.Protein.Description",
                                               "ID", "name", additional_cols)))
            }
            data_se <- make_se_customized(
                data_unique, selected_cols, exp_design,
                log2transform = log2transform, exp = "DIA", level = "protein")
            dimnames(data_se) <-
                list(dimnames(data_se)[[1]],
                     SummarizedExperiment::colData(data_se)$sample_name)
            SummarizedExperiment::colData(data_se)$label <-
                SummarizedExperiment::colData(data_se)$sample_name
        }else{
            stop("Current level is not supported.")
        }
        supported_types <- c("global", "phospho", "glyco", "acetyl",
                             "ubiquit")
        if (!is.null(exp_type)) {
            if (!(exp_type %in% supported_types)) {
                cat(paste0(exp_type), " is not a customized type.
                might not work.")
            }
        }
        S4Vectors::metadata(data_se)$exp_type <- exp_type
        S4Vectors::metadata(data_se)$level <- level
        return(data_se)
    }
}
