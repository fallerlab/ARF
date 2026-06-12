## v2.3 — 2026-06-12

### Added
  - `driftARF` — **D**ynamic **RI**bosomal **F**eature **T**racking: a continuous-feature
    subtool alongside dripARF/dricARF. Instead of contrasting discrete groups, it correlates
    a numeric per-sample feature (a column in the samples file) with position-specific rRNA
    fragment abundances (Pearson `r`/`p` on the VST-normalised matrix) and runs the existing
    RPSEA + ORA engine to nominate RP heterogeneity sets, ribosome collision sets, or
    user-given sets whose fragmentation tracks the feature. New file `R/driftARF.R`, exporting
    `driftARF()` and `driftARF_predict_progression()`. Loops over multiple `features`; supports
    `measureID` (`abs_cor_measure` default, `cor_measure` signed, `abs_r`, `r`),
    `include_collision`, and reuses `dripARF_simplify_results` /
    `dripARF_result_scatterplot` / `dripARF_result_heatmap` unchanged (the output `comp`
    column holds the feature name).
  - `run_limma_cor_analysis` (`R/ssRPSEA.R`): continuous-covariate analogue of
    `run_limma_DE_analysis` for driftARF. Reuses the same per-sample ssGSEA scores
    (`run_DESeq2_norm`) but fits set activity against the z-scored feature (`~ feature_z`)
    instead of a group contrast, producing the driftARF `ssRPSEA.weight` and
    `weighted.RPSEA.NES_randZ` with the identical `(1 - adj.P.Val) / (1 + |logFC|)` weight.
    Per-feature interaction plot saved when `ssRPSEAplots = TRUE`.

### Changed
  - `DESCRIPTION`: version bumped to 2.3.

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
