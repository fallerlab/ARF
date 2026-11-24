## v3.0 — 2026-10-14

### Added
  - Collision interface pipeline: `collision_set_builder.R`, `python_backend.R`,
    `structure_io.R`, `structure_sasa.R`, `structure_validation.R`. New public
    functions: `ARF_download_structure_files`, `ARF_run_SASA`,
    `ARF_build_SASA_candidates`, `ARF_build_ColInt_set`, `ARF_parse_structure_fasta`,
    `ARF_build_collision_sets`, `dricARF_build_SASA_set`, `configure_sasa_python`,
    `check_sasa_python`. Requires Python ≥ 3.10 with `gemmi` and `freesasa`.
  - `ARF_platform.R`: `auto_pair_rRNAs()` — pairs source and target rRNAs by
    global alignment score. Used as the default `rRNA_pairs = NULL` behaviour
    across `ARF_build_collision_sets()`, `dricARF_build_SASA_set()`,
    `ARF_build_ColInt_set()`, and `ARF_convert_Ribo3D_pos()`.
  - `dricARF.R` / `dripARF.R`: `rRNA_pairs` parameter added to `dricARF()`,
    `dripARF()`, `dripARF_predict_heterogenity()`, and
    `dripARF_report_RPset_group_counts()` for explicit source→target rRNA pairing.
  - `dricARF.R`: `RP_proximity_df` parameter added to
    `dricARF_liftover_collision_sets()` with `parent.frame()` fallback.
  - `docs/ARF_collision_sets.Rmd`: user guide for the collision interface
    pipeline with Quick Start (single and multiple structures) and full
    `dricARF()` prediction examples.
  - `docs/ARF_use_cases.Rmd`: converted and expanded from
    `dricARF_use_cases_documentation.md` (v2.2).
  - `inst/extdata/EscherichiaColi.fa`: E. coli rRNA FASTA added to package data.
  - `inst/python/structure_sasa.py`: Python SASA backend.
  - `tests/testthat/`: testthat suite added (19 tests across 4 new test files).
    Pipeline-level functions that require a CIF structure or Python environment
    (`ARF_build_collision_sets`, `dricARF_build_SASA_set`, `ARF_build_ColInt_set`,
    `configure_sasa_python`) are not unit-tested; coverage is limited to functions
    that can be exercised without external file fixtures.
  - `DESCRIPTION`: `reticulate` added to `Imports`; `testthat` added to
    `Suggests`; `SystemRequirements` updated; version bumped to 3.0.0.

### Known limitations
  - `dricARF_build_SASA_set`: liftover indexes `pos_map` by `auth_seq_id`; correct
    for all standard PDB rRNA structures (1-based sequential numbering) but would
    silently mis-map residues in structures with non-1-based `auth_seq_id`. Proper
    fix: propagate `label_seq_id` through `ARF_build_SASA_candidates` and use it
    as the positional index instead.

### Fixed
  - `dricARF.R`: residue positions now matched by character key instead of
    integer index, preventing wrong-residue selection on non-contiguous vectors.
  - `dricARF.R`, `ARF_platform.R`: T→U and X→N conversion applied before
    alignment to handle structure FASTAs using DNA-alphabet characters.
  - `dripARF.R`: `excludedCols` made available throughout the per-RP loop
    regardless of whether thresholds are supplied.
  - `collision_set_builder.R`: several fixes to new functions — `fasta_file`
    return value, `existsFunction()` replaced with `exists()`, `n_threads`
    wired through to Python backend, `sasa_A2_diff > 0` filter corrected,
    `!is.na(sasa_A2)` pre-filter added, ambiguous-base conversion extended to
    target rRNAs.
  - `structure_sasa.R`: stale option name `rSASA.python.envname` corrected.
  - `structure_validation.R`: `n_threads > 1` with `LeeRichards` now warns and
    clamps to 1 to prevent a FreeSASA segfault.
  - `python_backend.R` (`configure_sasa_python`): direct Python path support
    via `system2()` bypassing reticulate conda detection; `env_type` option now
    set on all exit paths; confirmation prompt added for `force = TRUE`.

### Changed
  - `ARF_build_SASA_candidates`, `dricARF_build_SASA_set`: Rib.Col. threshold
    changed to z > 1.64 on the raw differential SASA.
  - `dricARF()`, `dripARF()` and related functions: when a preset organism is
    used with a custom rRNA FASTA whose sequence names differ from the preset,
    ARF now auto-lifts coordinates via ClustalW alignment and emits an
    informational message. Previously this mismatch silently produced empty results.
  - `ARF_platform.R` (`.remap_rp_proximity_to_fasta`): row names now follow the
    `rRNA_resno` convention; gap positions dropped from `gsea_sets_RP`;
    `rRNA_pairs_used` returned to avoid redundant re-alignment downstream.

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
