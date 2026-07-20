#' Run a the downstream proteomics pipeline (filter, impute, limma, reporting)
#'
#' This workflow performs a reproducible downstream analysis for quantitative
#' proteomics: (i) imports quantification and sample annotation into a
#' SummarizedExperiment, (ii) filters features by missingness,
#' (iii) imputes remaining missing values using a mean vs LDV (left-censored)
#' strategy, (iv) runs differential abundance testing with limma (optionally
#' including an interaction term), (v) organizes results into plot-ready
#' effect objects, and (vi) optionally exports volcano plots, deregulograms and
#' results to disk.
#'
#' @section Input files:
#' \describe{
#'   \item{\code{path_pgmatrix}}{
#'     Path to a quantitative matrix file produced by the upstream
#'     quantification tool (e.g., a TSV). Rows represent quantified features
#'     (e.g. protein groups) and columns represent samples.
#'   }
#'   \item{\code{path_annotation}}{
#'     Path to a tab-delimited sample annotation file describing the
#'     experimental design. At minimum it must encode sample-to-condition
#'     mapping. Additional columns are required
#'     (\code{batch}, \code{replicate}, \code{block}).
#'     Accepts either a standard annotation TSV (default,
#'     \code{annotation_format = "standard"}) or an SDRF-Proteomics file
#'     (\code{annotation_format = "sdrf"}, with column roles supplied via
#'     \code{sdrf_map}).
#'   }
#' }
#'
#' @section Main steps:
#' \enumerate{
#'   \item \strong{Import}: creates a raw SummarizedExperiment
#'   (\code{se_raw}).
#'   \item \strong{Missingness filter}: removes features not sufficiently
#'   observed across samples/conditions according to \code{percent_missing}.
#'   Filtered-out features are returned and may be written to disk.
#'   \item \strong{Imputation}: imputes remaining missing values using a
#'   two-step strategy: within-condition mean imputation when missingness is
#'   below a threshold and LDV (left-censored) imputation when missingness is
#'   above the threshold. The imputation method is controlled by
#'   \code{ldv_source} and \code{threshold}.
#'   \item \strong{Differential testing}: calls
#'   \code{\link{test_limma_customized}} to fit a limma model and evaluate
#'   user-specified contrasts. Interaction testing is optional.
#'   \item \strong{Organize outputs}: splits the limma results into plot-ready
#'   effect objects using \code{.pp_organize_effects}:
#'         \itemize{
#'           \item \code{all_common_effect}: main-effect statistics for all
#'           proteins (no interaction filtering).
#'           \item \code{common_effect}: proteins not significant for the
#'           interaction term (main effect stable).
#'           \item \code{interaction_effect}: proteins significant for the
#'           interaction term.
#'         }
#'   \item \strong{Optional reporting} (\code{plots = TRUE}): exports volcano
#'   plots and tables for the three effect types above, and generates
#'   deregulograms when an interaction factor with exactly two levels is
#'   detected (e.g., female vs male, day1 vs day2, ...).
#' }
#'
#' @section Differential testing inputs:
#' The limma step requires the following explicit parameters:
#' \describe{
#'   \item{\code{tests}}{
#'     Character vector of main contrasts to test, using a consistent naming
#'     scheme (e.g. \code{"A_vs_B"}). These define the primary condition
#'     effects.
#'   }
#'   \item{\code{tests_interaction}}{
#'     Character vector aligned with \code{tests} defining interaction
#'     contrasts or interaction coefficients (or \code{"NA"} to disable
#'     interaction testing). These represent factor-dependent differences in the
#'     condition effect (e.g. condition-by-batch).
#'   }
#'   \item{\code{formula}}{
#'     A model formula describing the design matrix (e.g.
#'     \code{~ 0 + condition + batch + condition:batch}). Variables referenced
#'     in the formula must be present in \code{colData(se)}.
#'   }
#'   \item{\code{reference_condition}}{
#'     Character scalar indicating the reference level for \code{condition}.
#'     The function will relevel \code{colData(se)$condition} to this value
#'     before model fitting.
#'   }
#' }
#'
#' @section Paired designs and blocking:
#' \describe{
#'   \item{\code{paired}}{
#'     Logical (forwarded to \code{test_limma_customized}). If \code{TRUE},
#'     the model is treated as a paired / repeated-measures design by ensuring
#'     \code{replicate} is included in the design formula. In the current
#'     implementation, if \code{replicate} is not already present in
#'     \code{formula}, it is added automatically. This preserves any additional
#'     covariates already present in the formula.
#'   }
#'   \item{\code{block_effect}}{
#'     Logical (forwarded to \code{test_limma_customized}). If \code{TRUE},
#'     correlation between repeated observations is modeled using
#'     \code{limma::duplicateCorrelation} and
#'     \code{lmFit(..., block=..., correlation=...)}.
#'     This requires a \code{block} column in \code{colData(se)}. Use this when
#'     there is a known blocking factor (e.g., donor/subject) inducing
#'     correlation across samples.
#'   }
#' }
#'
#' @section Volcano vs deregulogram (interpretation note):
#' Interaction volcano plots highlight proteins with significant interaction
#' coefficients (i.e., where the condition effect depends on a second factor).
#' Deregulograms are more restrictive: they visualize full effects across the
#' two levels of the interaction factor and emphasize proteins where interaction
#' occurs in the context of a relevant condition effect (effect-size and FDR
#' thresholds). Therefore, the labeled proteins in the deregulogram are not
#' expected to match one-to-one with those in the interaction volcano plot.
#'
#' @param path_pgmatrix Character scalar. Path to the quantitative matrix file.
#' @param path_annotation Character scalar. Path to the sample annotation TSV.
#' @param path_output Character scalar. Output directory where files may be
#' written.
#' @param level Character. Quantification level (e.g. \code{"protein"}).
#' @param type Character. Quantification type (e.g. \code{"DIA"}).
#' @param experiment Character scalar. Experiment identifier used in output
#' filenames.
#' @param percent_missing Numeric. Missingness threshold used for filtering
#' (percentage, e.g. 50).
#' @param ldv_source Character. One of \code{"global"} or \code{"per-condition"}
#' controlling LDV definition.
#' @param threshold Numeric in (0, 1]. Forwarded to \code{impute_se()}.
#' Missingness proportion below which within-condition mean imputation is
#' used; at or above it, LDV (left-censored) imputation is used instead.
#' @param paired Logical. Forwarded to \code{test_limma_customized()}. If
#' \code{TRUE}, treats the design as paired/repeated-measures (see
#' "Paired designs and blocking" section below).
#' @param block_effect Logical. Forwarded to \code{test_limma_customized()}.
#' If \code{TRUE}, models correlation between repeated observations via
#' \code{limma::duplicateCorrelation} (see "Paired designs and blocking"
#' section below); requires a \code{block} column in \code{colData(se)}.
#' @param annotation_format Character. One of \code{"standard"} (default) or
#' \code{"sdrf"}. Use \code{"sdrf"} when \code{path_annotation} is an
#' SDRF-Proteomics file rather than a plain annotation TSV with
#' \code{file}/\code{sample_name}/\code{condition} columns; see
#' \code{sdrf_map}.
#' @param sdrf_map Named list, only used when
#' \code{annotation_format = "sdrf"}. Maps SDRF column names (which are not
#' fixed by the standard) to the roles the pipeline needs, e.g.
#' \code{list(condition = "factor value[treatment]", batch = "comment[batch]",
#' donor_id = "characteristics[individual]")}. At minimum \code{condition}
#' must be supplied. Ignored when \code{annotation_format = "standard"}.
#' @param plots Logical. If \code{TRUE}, exports plots/tables to
#' \code{path_output}.
#'
#' @param tests Character vector of main contrasts to test.
#'   Each element must correspond to a valid condition contrast present
#'   in the design matrix (e.g. \code{"IL13_vs_NoTreated"}).
#'
#' @param tests_interaction Character vector defining interaction contrasts
#'   or interaction coefficients aligned with \code{tests}. Use \code{"NA"}
#'   to disable interaction testing for a given contrast. Typical values
#'   encode condition-by-factor interactions (e.g. \code{"IL13.batchday2"}).
#'
#' @param formula A model formula specifying the design matrix for limma.
#'   Common examples include \code{~ 0 + condition + batch + condition:batch}.
#'   All variables referenced in the formula must be present in
#'   \code{colData()}.
#'
#' @param reference_condition Character scalar specifying the reference
#'   level for the \code{condition} variable. This level is used to relevel
#'   the design prior to fitting the limma model.
#'
#' @param ... Additional arguments forwarded only to
#' \code{test_limma_customized()} (beyond \code{test}, \code{test_interaction},
#' \code{design_formula}, \code{ref_condition}, \code{paired} and
#' \code{block_effect}, which are already explicit parameters above).
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{se_raw}}{Raw SummarizedExperiment
#'   imported from files.}
#'   \item{\code{se_filt}}{Filtered SummarizedExperiment.}
#'   \item{\code{removed}}{A \code{data.frame} describing filtered-out features.}
#'   \item{\code{se_imp}}{Imputed SummarizedExperiment.}
#'   \item{\code{effects}}{
#'     A named list keyed by main contrast (the values of \code{tests}). Each
#'     element is itself a list of three \code{SummarizedExperiment} objects:
#'     \code{all_common_effect}, \code{common_effect} and
#'     \code{interaction_effect}. Per-contrast statistics are stored in their
#'     \code{rowData()} using the \code{"<stat>_<contrast>"} naming convention
#'     (e.g. \code{log2FC_IL13_vs_NoTreated}). This element is always returned,
#'     independently of \code{plots}, so a contrast can be accessed with
#'     \code{$}, e.g. \code{res$effects$IL13_vs_NoTreated$all_common_effect}.
#'   }
#' }
#'
#' @seealso \code{\link{test_limma_customized}},
#' \code{\link{plot_volcano_customized}},
#'   \code{\link{plot_deregulogram}}
#'
#' @examples
#' ## Load package
#' library(ProteinBatcher)
#'
#' ## ---------------------------------------------------------------
#' ## Input files shipped with the package
#' ## ---------------------------------------------------------------
#'
#' path_annotation <- system.file(
#'   "extdata", "annotation_HaCaT.tsv",
#'   package = "ProteinBatcher"
#' )
#'
#' path_pgmatrix <- system.file(
#'   "extdata",
#'   "2024MK017_HaCaT_Stimulation_1to24_Astral_report.pg_matrix.tsv",
#'   package = "ProteinBatcher"
#' )
#'
#' ## Output directory
#' path_output <- tempdir()
#'
#' ## ---------------------------------------------------------------
#' ## Workflow parameters (single contrast example)
#' ## ---------------------------------------------------------------
#'
#' experiment <- "HaCaT"
#' percent_missing <- 50
#'
#' formula <- ~ 0 + condition + batch + condition:batch
#'
#' tests <- c("IL13_vs_NoTreated")
#' tests_interaction <- c("IL13.batchday2")
#'
#' reference_condition <- "NoTreated"
#'
#' ## ---------------------------------------------------------------
#' ## Run downstream proteomics workflow
#' ## ---------------------------------------------------------------
#'
#' res <- run_proteomics_pipeline(
#'   path_pgmatrix        = path_pgmatrix,
#'   path_annotation      = path_annotation,
#'   path_output          = path_output,
#'   tests                = tests,
#'   tests_interaction    = tests_interaction,
#'   formula              = formula,
#'   reference_condition  = reference_condition,
#'   percent_missing      = percent_missing,
#'   ldv_source           = "per-condition",
#'   experiment           = experiment,
#'   plots                = FALSE
#' )
#'
#' ## Inspect outputs
#' names(res)
#' res$se_imp
#'
#' ## `effects` is keyed by contrast, so it can be accessed with `$`:
#' names(res$effects)
#' res$effects$IL13_vs_NoTreated$all_common_effect
#' @export
run_proteomics_pipeline <- function(
        path_pgmatrix, path_annotation, path_output, level = "protein",
        type = "DIA", experiment, percent_missing,
        ldv_source = c("global", "per-condition"), threshold = 0.3,
        tests, tests_interaction, formula, reference_condition,
        paired = FALSE, block_effect = FALSE,
        annotation_format = c("standard", "sdrf"), sdrf_map = NULL,
        plots = FALSE, ...
){
    args <- .pp_validate_inputs(path_pgmatrix, path_annotation, level, type,
                                percent_missing, path_output, experiment,
                                annotation_format = annotation_format,
                                sdrf_map = sdrf_map)
    se0 <- .pp_import_se(args$path_pgmatrix, args$path_annotation,
                         args$level, args$type,
                         annotation_format = args$annotation_format,
                         sdrf_map = args$sdrf_map)
    # 1) Filter and imputation
    filt <- filter_se_missing(se0, percentage = args$percent_missing)
    .pp_write_filtered(filt$removed, args$path_output, args$experiment)
    # `threshold` and `ldv_source` are impute_se()'s only tunable arguments;
    # they are now explicit parameters of this wrapper instead of being
    # bundled into `...`, so accidentally passing a test_limma_customized()
    # argument here (e.g. block_effect) no longer breaks impute_se().
    se_imp <- impute_se(filt$se_filt, threshold = threshold,
                        ldv_source = ldv_source)
    .pp_write_before_after(filt$se_filt, se_imp,
                           args$path_output, args$experiment)
    # 2) Differential testing
    # `paired` and `block_effect` are explicit parameters too; any remaining
    # `...` is forwarded only to test_limma_customized().
    se_limma <- test_limma_customized(
        se_imp, type = "manual", test = tests,
        test_interaction = tests_interaction, design_formula = formula,
        ref_condition  = reference_condition,
        paired = paired, block_effect = block_effect, ...
    )
    # 3) Organize outputs (plot-ready). `effects` is a named list keyed by main
    #    contrast; each element holds the all_common_effect / common_effect /
    #    interaction_effect SummarizedExperiment objects. Per-contrast statistics
    #    are stored in rowData() using the "<stat>_<contrast>" naming convention
    #    (e.g. log2FC_IL13_vs_NoTreated), so a contrast can be retrieved with `$`
    #    (e.g. effects$IL13_vs_NoTreated$all_common_effect).
    effects <- .pp_organize_effects(
        se_limma = se_limma, tests = tests,
        tests_interaction = tests_interaction, alpha = 0.05
    )
    if (plots) {
        # 4) Volcano plots (written to disk as a side effect).
        # Main effect: all proteins, ignoring interaction.
        .pp_export_volcano_tables(
            effects = effects, tests = tests, path_output = path_output,
            experiment = experiment, effect_slot = "all_common_effect",
            file_tag = "MainEffect_AllProteins", name_col = "Genes",
            name_imputed = "imputed"
        )
        # Common effect: significant for condition and NOT interacting.
        .pp_export_volcano_tables(
            effects = effects, tests = tests, path_output = path_output,
            experiment = experiment, effect_slot = "common_effect",
            file_tag = "CommonEffect", name_col = "Genes",
            name_imputed = "imputed"
        )
        # 5) Interaction volcano + deregulogram (only when interaction tested
        #    and the interaction factor has exactly two levels).
        if (!("NA" %in% tests_interaction)) {
            .pp_export_volcano_tables(
                effects = effects, tests = tests, path_output = path_output,
                experiment = experiment, effect_slot = "interaction_effect",
                file_tag = "InteractionEffect",
                plot_contrasts = tests_interaction,
                name_col = "Genes", name_imputed = "imputed"
            )
            info <- .pp_infer_interaction_factor(se_limma, tests_interaction)
            if (!is.null(info) && length(info$levels) == 2) {
                plot_deregulogram(se_limma = se_limma, tests = tests,
                                  tests_interaction = tests_interaction,
                                  factor_name = info$name,
                                  factor_levels = info$levels,
                                  path_output = path_output,
                                  experiment = experiment, alpha = 0.05,
                                  lfc = 1, label_col = "Genes"
                )
            }
        }
    }
    # 6) Single, consistent return value. `effects` is ALWAYS present and keyed
    #    by contrast, independent of `plots`, so downstream code can rely on
    #    res$effects$<contrast>$all_common_effect.
    list(
        se_raw  = se0,
        se_filt = filt$se_filt,
        removed = filt$removed,
        se_imp  = se_imp,
        effects = effects
    )
}
