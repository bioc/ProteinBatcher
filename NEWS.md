# ProteinBatcher 0.99.0

* Initial Bioconductor submission.

# ProteinBatcher 0.99.7

* Optional annotation columns: `batch`, `donor_id`, and `replicate` are no
  longer required in the standard annotation TSV. The pipeline runs with only
  `file`, `sample`, `sample_name`, and `condition`. Missing optional columns
  are handled gracefully throughout the workflow.
* Configurable block variable: `test_limma_customized()` and
  `run_proteomics_pipeline()` now accept a `block_var` parameter (default:
  `"donor_id"`) specifying which colData column to use as the blocking
  variable for `limma::duplicateCorrelation()` when `block_effect = TRUE`.
  Previously the code hardcoded a `block` column that did not exist in the
  annotation schema.
* Bug fix: `readQuantTable()` now dispatches to the mzTab reader
  (`.pp_read_mztab_pg_matrix()`) when the input file is in mzTab format,
  instead of always treating it as a DIA pg_matrix TSV.
* Bug fix: Contrast column names in `.tl_run_contrasts()` were not preserved
  by `limma::makeContrasts(contrasts = ...)`, causing `topTable()` to fail
  with a subscript-out-of-bounds error. Contrast names are now set explicitly
  after `makeContrasts()`.
* Removed unused dependencies: `readxl` (Imports) and `htmltools` (Suggests)
  to resolve DESCRIPTION/NAMESPACE mismatches reported by `R CMD check`.
