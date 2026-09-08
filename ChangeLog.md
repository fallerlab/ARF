## v2.2.1 — 2026-09-08

### Fixed
  - `dripARF.R`: `RPs_toreport` now wrapped in `na.omit()` at all four construction sites,
    preventing `NA` RP names from propagating into `endsWith()` and generating all-NA
    logical subscripts that crashed `NES_rand_zscore` assignment.
  - `dripARF.R`: `NES_rand_zscore` column now initialised with `rep(NA_real_, nrow(...))` instead
    of bare `NA`, fixing a zero-row data-frame replacement error when GSEA returns no
    qualifying pathways.
  - `dripARF.R`, `driftARF.R`: `tochange[is.na(tochange)] <- FALSE` guard added after
    every `endsWith()` call to handle any residual `NA` entries in the logical index.
  - `dripARF.R`, `dricARF.R`: `.check_targetDir()` helper added; called upfront in all
    functions that write files, giving a clear error if the target directory does not
    exist or is not writable (replaces cryptic `cannot open the connection` failure).
  - `dripARF.R`: output CSV filenames sanitised with `gsub("[^A-Za-z0-9._-]", "_", ...)`
    and assembled with `file.path()` to prevent invalid paths from special characters in
    group/comparison names.
  - `DESCRIPTION`: minimum R version declared as `R (>= 4.1.0)`.

## v2.2 — 2026-05-01

### Added
  - `ssRPSEA.R`: internal ssGSEA-based weighting of RPSEA scores. Per-sample RP activity
    is estimated with `GSVA::gsva` (ssGSEA mode), differential activity between conditions
    is modelled with `limma`, and the resulting weight is applied to ES2 (NES_randZ) to
    produce `weighted.RPSEA.NES_randZ` in all output CSVs. Optional interaction scatter
    plot saved when `ssRPSEAplots = TRUE`.
  - `docs/dricARF_use_cases_documentation.md`
  - `DESCRIPTION`: `GSVA`, `limma`, and `tidyr` added to `Imports` (required by ssRPSEA);
    `dplyr`, `magrittr`, `readr`, `RColorBrewer`, `circlize`, `wesanderson` added to
    replace removed blanket dependencies.

### Changed
  - `DESCRIPTION`: licence field corrected from `use_gpl3_license()` to `GPL (>= 3)`;
    `bedr` and `tidyverse` removed from `Imports`; `HelloRanges` moved to `Suggests`.
  - `dricARF_result_scatterplot`: `highlightRPs` argument now functional — labels chosen
    RP points in both panels of the combined plot.
  - All roxygen2 documentation overhauled across `ARF_platform.R`, `dripARF.R`, and
    `dricARF.R`: corrected `@param` names, completed `@return` tags, fixed `@examples`,
    resolved typos, and corrected wrong titles/keywords.
  - `README.md` comprehensively updated.

### Removed
  - `GSEAplots` parameter removed from `dripARF_predict_heterogenity`, `dripARF`, and
    `dricARF`; the `enrichplot::gseaplot2` per-RP GSEA plot block has been dropped.
  - `enrichplot` removed from `Imports`.

## v2.1
  - All species support for both dricARF and dripARF.
  - New functions for PDB parsing and lift-overing; ARF_parse_PDB_ribosome(), ARF_convert_Ribo3D_pos(), 
    dripARF_get_RP_proximity_sets(), and dricARF_liftover_collision_sets().
  - dricARF output figure updated.

## v2.0
  - Integrating 3D ribosome structure analysis to ARF
  - All relevant functions now require the rRNA fasta file that is used for mapping.
  - First release of dricARF - Differential Ribosome Collision Prediction pipeline
  - Bug fixes (closing issues until now)

## v1.0
  - Publication release
  - A few bug fixes

## v0.9.1 - 2022-04-14
  - Paper revision release
  - rRNA-RP contact point distance threshold changed
  - New experimental function that allows altering the Contact Point threshold and possible directionality in predictions
  - A few technical updates
  - RP-rRNA distance matrices are added to the repo

## v0.9 - 2021-11-26
  - First release
  - Paper submission version
