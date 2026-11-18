## v2.2.2 — 2026-10-14

### Added
  - `collision_set_builder.R`, `python_backend.R`, `structure_io.R`,
    `structure_sasa.R`, `structure_validation.R`: new SASA-based collision
    interface analysis functions. Public API: `ARF_download_structure_files`,
    `ARF_run_SASA`, `ARF_build_SASA_candidates`, `ARF_build_ColInt_set`,
    `dricARF_build_SASA_set`, `configure_sasa_python`, `check_sasa_python`.
    Internal helpers (not exported): `ARF_build_ColInt_proximity`,
    `dripARF_get_collision_proximity_sets`, `calculate_structure_sasa`,
    `read_structure_annotations`.
    Requires Python ≥ 3.10 with `gemmi` and `freesasa` (via reticulate).
  - `collision_set_builder.R`: `ARF_parse_structure_fasta()` — parses a RCSB
    PDB FASTA into a conversion table (chain → RP name) and a species-prefixed
    rRNA `BStringSet`; optionally writes the rRNA sequences to a FASTA file.
  - `collision_set_builder.R`: `ARF_build_collision_sets()` — one-call wrapper
    that runs the entire collision-interface pipeline (download → SASA → FASTA
    parsing → SASA candidates → Rib.Col. sets → Col.Int. sets) and returns a
    combined `gsea_sets_Collision` data.frame ready for `dricARF()`.
  - `ARF_platform.R`: `auto_pair_rRNAs(source_fasta, target_fasta)` — pairs
    each rRNA in `source_fasta` with its best-matching rRNA in `target_fasta`
    via global pairwise alignment (`Biostrings::pairwiseAlignment`). Used
    internally as the default `rRNA_pairs = NULL` behaviour in
    `ARF_build_collision_sets()`, `dricARF_build_SASA_set()`,
    `ARF_build_ColInt_set()`, and `ARF_convert_Ribo3D_pos()`; can also be
    called directly.
  - `docs/ARF_collision_sets.Rmd`: user guide for the collision interface
    pipeline; includes a Quick Start section (`ARF_build_collision_sets()`),
    a step-by-step walkthrough, and a prediction section using `dricARF()`.
  - `docs/ARF_use_cases.Rmd`: converted and expanded from
    `dricARF_use_cases_documentation.md` (v2.2).
  - `inst/extdata/EscherichiaColi.fa`: E. coli rRNA FASTA added to package
    data (accessible via `system.file()`).
  - `inst/python/structure_sasa.py`: Python SASA backend.
  - `DESCRIPTION`: `reticulate` added to `Imports`;
    `SystemRequirements: Python (>= 3.10), Python packages: gemmi, freesasa` added.
  - `tests/testthat/`: testthat suite added covering `ARF_check_organism`,
    `read_ARF_samples_file`, `dripARF_simplify_results`,
    `dripARF_get_RP_proximity_sets`, and `dripARF_read_rRNA_fragments`.
  - `DESCRIPTION`: `testthat (>= 3.0.0)` added to `Suggests`;
    `Config/testthat/edition: 3` added.

### Fixed
  - `dricARF.R`: `proxpos` now resolved via character-matching
    (`paste(rRNA, resno) %in% target_positions[[RP]]`) instead of direct integer
    indexing, preventing wrong-residue selection when position vectors are not
    contiguous. Added `if (length(proxpos) == 0) return(NULL)` null-guard.
  - `dricARF.R`, `ARF_platform.R`: T→U and X→N nucleotide conversion
    (`chartr("TtXx", "UuNn", ...)`) applied to `yeast_rRNAs`, `human_rRNAs`,
    and `target_rRNAs` before alignment, preventing mismatches when structure
    FASTAs use DNA-alphabet characters.
  - `dripARF.R`: `excludedCols` computation moved outside the
    `if (is.null(thresholds))` block so it is always available inside the
    per-RP loop regardless of whether thresholds are supplied by the caller.
  - `structure_sasa.R`: stale option name `rSASA.python.envname` corrected to
    `ARF.python.env_type`.
  - `collision_set_builder.R` (`ARF_download_structure_files`): `fasta_file`
    was never reassigned after download (local variable `rRNAs_file` was set
    instead), causing the return list to carry `fasta_file = NULL`.
  - `collision_set_builder.R` (`ARF_build_SASA_candidates`): `existsFunction()`
    replaced with `exists(..., where = asNamespace("ARF"), mode = "function")` —
    `existsFunction` is not a base R function and would error at runtime.
  - `collision_set_builder.R` (`ARF_run_SASA`): `n_threads` added to function
    signature and passed through to `calculate_structure_sasa` instead of being
    hardcoded to `1L`.
  - `collision_set_builder.R` (`ARF_build_SASA_candidates`): `sasa_A2_diff == 0`
    filter replaced with `sasa_A2_diff > 0` — collision interface residues are
    differentially exposed across chains (buried in one ribosome copy, exposed in
    the other); the previous filter selected the complement and returned 0 rows
    for disome structures.
  - `collision_set_builder.R` (`ARF_build_SASA_candidates`): added
    `!is.na(sasa_A2)` pre-filter before `group_by`/`reframe` to prevent
    `max()`/`min()` warnings on all-NA groups.
  - `structure_validation.R`: `n_threads > 1` with `algorithm = "LeeRichards"`
    now triggers a warning and clamps to `1L` at the R level, preventing the
    FreeSASA segfault that occurred when the constraint was violated.
  - `docs/ARF_collision_sets.Rmd`: `rRNAs` in Step 2 now strips the species
    prefix before passing to `ARF_build_SASA_candidates`, matching the
    unprefixed `entity_description` column in the SASA TSV.
  - `collision_set_builder.R` (`dricARF_build_SASA_set`): `chartr("TtXx",
    "UuNn", ...)` now applied to `target_rRNAs` as well as `structure_rRNAs`
    before alignment, preventing `'X' not in lookup table` errors when the
    target FASTA contains ambiguous bases.
  - `docs/ARF_collision_sets.Rmd` (prediction section): now uses
    `6XZA.rRNAs.fasta` as the reference throughout (`rRNA_pairs`, RP proximity
    liftover, and `dricARF()`) to match the chromosome names in the E. coli
    Ribo-seq bedGraph files; the previous `EscherichiaColi.fa` target caused a
    gene-ID mismatch in the GSVA step.
  - `docs/ARF_collision_sets.Rmd` (prediction section): example updated to use
    `Escherichia_coli_PRJNA1126309` test data (WT / hrpA / hrpA_smrB groups).
  - `collision_set_builder.R`, `ARF_platform.R`: `rRNA_pairs` argument now
    defaults to `NULL` in `ARF_build_collision_sets()`, `dricARF_build_SASA_set()`,
    `ARF_build_ColInt_set()`, and `ARF_convert_Ribo3D_pos()`; auto-pairing via
    `auto_pair_rRNAs()` runs when no explicit pairs are supplied.
  - `docs/ARF_collision_sets.Rmd`, `docs/ARF_use_cases.Rmd`, `README.md`:
    `rRNA_pairs` removed from all running code examples; function-reference
    signatures and step-by-step sections updated to reflect the optional
    argument and document when manual override is needed.
  - `python_backend.R` (`configure_sasa_python`): when a direct Python path is
    supplied via `python = ...`, packages are checked via `system2()` (no
    reticulate touch); if already present, the path is stored in
    `options(ARF.python.path)` and no virtualenv is created; if absent,
    `python -m venv` is used to create the virtualenv, bypassing reticulate's
    conda-binary detection (required for micromamba environments).
  - `python_backend.R` (`configure_sasa_python`): `options(ARF.python.env_type)`
    now set on both the early-exit (env already exists) and successful-creation
    paths, so `check_sasa_python()` picks up the correct backend type
    automatically after `configure_sasa_python()` returns.
  - `python_backend.R` (`configure_sasa_python`): interactive confirmation
    prompt added when `force = TRUE` and an existing environment is found,
    preventing accidental removal in interactive sessions.
  - `docs/ARF_collision_sets.Rmd`: `python-config` chunk added before the Quick
    Start; auto-detects whether `sasa-python` exists as a conda or virtualenv
    and calls `configure_sasa_python()` accordingly, so the document knits
    without requiring manual Python setup in the session.
  - `DESCRIPTION`: version bumped to 2.2.2; typos in `Description` field
    corrected ("consists of", "heterogeneity", "Both"); `LazyData: true` and
    `RoxygenNote: 7.3.2` added; `Config/roxygen2/version` removed (generated
    field).
  - `ARF_platform.R`: four internal helpers for auto-liftover when a preset
    organism is used with a custom rRNA FASTA whose names differ from the preset
    coordinate system: `.preset_rRNA_to_extdata()`, `.extdata_to_standard_id()`,
    `.detect_rRNA_mismatch()`, `.remap_rp_proximity_to_fasta()`. Not exported.
  - `dricARF.R` / `dripARF.R`: `rRNA_pairs` parameter added to `dricARF()`,
    `dripARF()`, `dripARF_predict_heterogenity()`, and
    `dripARF_report_RPset_group_counts()`. Allows callers to supply explicit
    source→target rRNA pairings for the auto-liftover step; when `NULL` and a
    mismatch is detected, `auto_pair_rRNAs()` is used automatically.
  - `dricARF.R`: `RP_proximity_df` parameter added to
    `dricARF_liftover_collision_sets()` with a `parent.frame()` fallback for
    backward compatibility.
  - `tests/testthat/test-auto_pair_rRNAs.R` (5 tests),
    `tests/testthat/test-ARF_parse_structure_fasta.R` (6 tests),
    `tests/testthat/test-ARF_build_SASA_candidates.R` (5 tests), and
    `tests/testthat/test-dricARF_liftover_collision_sets.R` (3 tests) added.

### Changed
  - `dricARF()`, `dripARF()`, `dripARF_predict_heterogenity()`,
    `dripARF_report_RPset_group_counts()`: when a preset organism (`hs`, `mm`,
    or `sc`) is specified but the rRNA sequence names in the user-supplied FASTA
    do not match those in the preset `RP_proximity_df`, ARF now automatically
    lifts over the preset coordinate data to the user FASTA via ClustalW
    alignment and emits an informational `message()`. Previously, this mismatch
    would silently produce zero matches and empty results. No change in behaviour
    when the FASTA names already match the preset.
  - `ARF_platform.R` (`.remap_rp_proximity_to_fasta`): row names of the
    remapped `RP_proximity_df` are now set as `paste(rRNA, resno, sep="_")`,
    consistent with the convention used throughout the package. Previously they
    were cleared to `NULL`, which caused row-name indexing in `dricARF()` to
    silently append rows instead of updating them.
  - `ARF_platform.R` (`.remap_rp_proximity_to_fasta`): rRNA-position gene
    entries that could not be mapped to an aligned target position (gap positions)
    are now dropped from `gsea_sets_RP` rather than kept with their stale preset
    name, preventing silent set-size inflation.
  - `ARF_platform.R` (`.remap_rp_proximity_to_fasta`): now returns
    `rRNA_pairs_used` alongside the remapped objects so callers can pass the
    already-computed pairs onward without re-running alignment.
  - `dricARF()`: passes `remapped$rRNA_pairs_used` to
    `dricARF_liftover_collision_sets()` after auto-liftover, avoiding a redundant
    second ClustalW alignment run.

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
