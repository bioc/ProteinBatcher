#' Validate input to limma
#' @keywords internal
#' @noRd
.tl_check_inputs <- function(se, design_formula, ref_condition, test, test_interaction){
    ok <- inherits(se, "SummarizedExperiment") &&
        methods::is(design_formula, "formula") &&
        is.character(ref_condition) && length(ref_condition) == 1L &&
        is.character(test) && length(test) >= 1L &&
        is.character(test_interaction)
    if (!ok) {
        stop("Invalid inputs to test_limma_customized(): check that `se` is a ",
             "SummarizedExperiment, `design_formula` is a formula, ",
             "`ref_condition` is a single string, `test` is a non-empty ",
             "character vector and `test_interaction` is character.",
             call. = FALSE)
    }
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    # Only 'label' and 'condition' are strictly required in colData.
    # 'replicate', 'batch', 'donor_id' are optional and only needed when
    # the user references them in the design formula or block_var.
    need_cd <- c("label", "condition")
    if (!all(need_cd %in% colnames(cd))) {
        stop("Missing in colData(se): ", paste(setdiff(need_cd, colnames(cd)), collapse=", "), call. = FALSE)
    }
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    if (!("name" %in% colnames(rd))) stop("rowData(se) must contain column 'name'.", call. = FALSE)
    imp <- S4Vectors::metadata(se)$imputation_map
    if (is.null(imp)) stop("metadata(se)$imputation_map is missing.", call. = FALSE)
}

#' Aggregation of imputation map by condition
#' @keywords internal
#' @noRd
.tl_condition_imputation_map <- function(se, cd){
    imp_map <- S4Vectors::metadata(se)$imputation_map
    conds <- levels(factor(cd$condition))
    out <- lapply(conds, function(cond){
        sub <- imp_map[, cd$condition == cond, drop = FALSE]
        apply(sub, 1, function(x){
            if (any(x == "ldv")) "ldv" else if (any(x == "mean")) "mean" else "none"
        })
    })
    names(out) <- conds
    out
}

#' Normalize the formuala for model matrix
#' @keywords internal
#' @noRd
.tl_normalize_formula_for_model_matrix <- function(f) {
    f_chr <- as.character(f)

    # case update(): . ~ RHS
    if (length(f_chr) == 3 && f_chr[2] == ".") {
        return(stats::as.formula(paste("~", f_chr[3])))
    }

    # case normal: ~ RHS
    if (length(f_chr) == 2 && f_chr[1] == "~") {
        return(f)
    }

    stop("Unsupported design_formula structure: ", deparse(f), call. = FALSE)
}

#' Robust design matrix
#' @keywords internal
#' @noRd
.tl_make_design <- function(design_formula, cd){
    design_formula <- .tl_normalize_formula_for_model_matrix(design_formula)
    design <- stats::model.matrix(design_formula, data = cd)
    colnames(design) <- make.names(gsub("^condition", "", colnames(design)))
    design
}

#' Fit limma with or without "block" (duplicateCorrelation)
#' @keywords internal
#' @noRd
.tl_fit_limma <- function(raw, design, cd, block_effect, block_var = "donor_id"){
    if (!block_effect) {
        message("Fitting limma model without block")
        return(limma::lmFit(raw, design = design))
    }
    if (!(block_var %in% colnames(cd)))
        stop("Block variable '", block_var, "' is missing in colData(se).",
             call. = FALSE)
    if (length(unique(cd[[block_var]])) <= 1)
        stop("Block variable '", block_var, "' only has one level; ",
             "duplicateCorrelation requires 2 or more levels.",
             call. = FALSE)
    message("Fitting limma model with block variable '", block_var, "'")
    corfit <- limma::duplicateCorrelation(raw, design, block = cd[[block_var]])
    limma::lmFit(raw, design, block = cd[[block_var]],
                 correlation = corfit$consensus)
}

#' Execute contrasts + BH "recomputed" after filtering genes ldv/ldv
#' @keywords internal
#' @noRd
.tl_run_contrasts <- function(fit0, design, tests, cond_imp_map, map_mn, ref_condition){
    cn_expr <- gsub("_vs_", " - ", tests)
    cn_name <- gsub(" - ", "_vs_", cn_expr)

    cm <- limma::makeContrasts(contrasts = cn_expr, levels = design)
    colnames(cm) <- cn_name
    fit <- limma::eBayes(limma::contrasts.fit(fit0, cm))

    cond_cols <- intersect(colnames(design), names(map_mn))

    one <- function(comp_name, comp_expr){
        tt <- limma::topTable(fit, coef = comp_name, number = Inf, sort.by = "none",
                              adjust.method = "none", confint = TRUE)
        gene <- rownames(tt)

        present <- cond_cols[vapply(cond_cols,function(cc) grepl(paste0("\\b", cc, "\\b"), comp_expr), logical(1))]
        c1 <- if (length(present) >= 1) map_mn[[present[1]]] else ref_condition
        c2 <- if (length(present) >= 2) map_mn[[present[2]]] else ref_condition

        m1 <- cond_imp_map[[c1]][gene]; m2 <- cond_imp_map[[c2]][gene]
        keep <- !(m1 == "ldv" & m2 == "ldv")

        adj <- rep(NA_real_, length(gene))
        adj[keep] <- stats::p.adjust(tt$P.Value[keep], method = "BH")
        tt$adj.P.Val <- adj

        # Wide block for this contrast, named "<stat>_<comparison>" to match the
        # downstream naming convention consumed by plotting/visualization code.
        block <- data.frame(
            log2FC = tt$logFC, CI.L = tt$CI.L, CI.R = tt$CI.R,
            p.val = tt$P.Value, p.adj = tt$adj.P.Val,
            row.names = gene, check.names = FALSE, stringsAsFactors = FALSE
        )
        colnames(block) <- paste0(colnames(block), "_", comp_name)
        block
    }

    blocks <- Map(one, cn_name, cn_expr)
    # All blocks share identical row order (sort.by = "none" on the same fit),
    # so a column bind reproduces tidyr::pivot_wider() exactly.
    wide <- do.call(cbind, unname(blocks))
    wide <- as.data.frame(wide, check.names = FALSE, stringsAsFactors = FALSE)
    wide$ID <- rownames(wide)
    wide
}

#' Clean merge to rowData(se)
#' @keywords internal
#' @noRd
.tl_merge_rowdata <- function(se, wide_tbl){
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    rd$ID <- rownames(rd)
    wide_tbl <- as.data.frame(wide_tbl, check.names = FALSE,
                              stringsAsFactors = FALSE)
    idx <- match(rd$ID, wide_tbl$ID)
    add_cols <- setdiff(colnames(wide_tbl), "ID")
    rd[add_cols] <- wide_tbl[idx, add_cols, drop = FALSE]
    rownames(rd) <- rd$ID
    rd$ID <- NULL
    rd
}

#' Full contrast of interaction effects and base effects
#' @keywords internal
#' @noRd
.tl_build_full_contrasts <- function(test, test_interaction){
    main <- gsub("_vs_", " - ", test)
    inter <- gsub("_vs_", " - ", test_interaction)
    out <- paste0("(", main, ") + ", inter)
    gsub(" - ", "_vs_", out)
}
