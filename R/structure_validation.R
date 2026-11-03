# structure_validation.R
#
# Input validation helpers for rSASA.
#
# validate_structure_file()   — path checks, returns normalised path
# validate_output_file()      — output path checks, creates parent dir
# validate_sasa_parameters()  — full parameter validation, returns named list
#                               of coerced values ready for py$calculate_rna_sasa()


# Structure file

#' Validate a structure file path
#'
#' Checks that the path is a single, non-empty, readable file with a
#' recognised structural-biology extension.  Returns the canonical
#' (normalised) path so callers never have to call \code{normalizePath()}
#' themselves.
#'
#' @param structure_file Character scalar.  Path to a PDB or mmCIF file.
#'
#' @return Character scalar — the normalised absolute path.
#'
#' @keywords internal
validate_structure_file <- function(structure_file) {

    if (!is.character(structure_file) ||
        length(structure_file) != 1L ||
        is.na(structure_file) ||
        !nzchar(structure_file)) {
        stop(
            "'structure_file' must be a single non-empty character string.",
            call. = FALSE
        )
    }

    structure_file <- path.expand(structure_file)

    if (!file.exists(structure_file)) {
        stop(
            "Structure file does not exist: ", structure_file,
            call. = FALSE
        )
    }

    if (file.access(structure_file, 4L) != 0L) {
        stop(
            "Structure file is not readable: ", structure_file,
            call. = FALSE
        )
    }

    ext <- tolower(tools::file_ext(structure_file))

    if (!ext %in% c("pdb", "ent", "cif", "mmcif")) {
        warning(
            "Unusual structure-file extension: .", ext,
            ". Gemmi may still be able to determine the format.",
            call. = FALSE
        )
    }

    normalizePath(structure_file, mustWork = TRUE)
}


# Output file

#' Validate an output file path
#'
#' Checks that the path is writable and, optionally, does not already
#' exist.  Creates the parent directory if it is missing.
#'
#' @param output_file Character scalar.  Desired output path.
#' @param overwrite   Logical scalar.  Whether to allow overwriting an
#'   existing file.  Default \code{FALSE}.
#'
#' @return The normalised absolute path (character scalar).
#'
#' @keywords internal
validate_output_file <- function(output_file, overwrite = FALSE) {

    if (!is.character(output_file) ||
        length(output_file) != 1L ||
        is.na(output_file) ||
        !nzchar(output_file)) {
        stop(
            "'output_file' must be a single non-empty character string.",
            call. = FALSE
        )
    }

    output_file <- path.expand(output_file)

    if (file.exists(output_file)) {
        if (!isTRUE(overwrite)) {
            stop(
                "Output file already exists: ", output_file, "\n",
                "Set overwrite = TRUE to replace it.",
                call. = FALSE
            )
        }
        if (file.access(output_file, 2L) != 0L) {
            stop(
                "Output file is not writable: ", output_file,
                call. = FALSE
            )
        }
    }

    output_dir <- dirname(output_file)

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        if (!dir.exists(output_dir)) {
            stop(
                "Could not create output directory: ", output_dir,
                call. = FALSE
            )
        }
    }

    output_file
}


# SASA parameters

#' Validate and coerce SASA calculation parameters
#'
#' Validates every parameter accepted by \code{calculate_structure_sasa()},
#' coerces numeric-integer arguments to \code{integer}, and returns all
#' values as a named list.  Callers should use the returned list rather
#' than the original arguments so that coercions (e.g. \code{n_slices}
#' becoming an \code{integer}) are applied consistently.
#'
#' @param algorithm      Character.  \code{"LeeRichards"} or
#'   \code{"ShrakeRupley"}.
#' @param probe_radius   Positive numeric.  Probe radius in Ångström.
#' @param n_slices       Positive integer.  Lee-Richards slices.
#' @param n_points       Positive integer.  Shrake-Rupley points.
#' @param n_threads      Positive integer.  FreeSASA threads.
#' @param polymer_type   Character.  \code{"RNA"}, \code{"DNA"}, or
#'   \code{"ALL"}.
#' @param model          \code{NULL} or a single non-\code{NA} scalar
#'   used to select a specific model from the structure.
#' @param include_hydrogen Logical scalar.
#' @param zero_tolerance   Non-negative numeric scalar.
#' @param precision        Positive integer.  Decimal places in output.
#'
#' @return A named list with elements \code{algorithm},
#'   \code{probe_radius}, \code{n_slices}, \code{n_points},
#'   \code{n_threads}, \code{polymer_type}, \code{model},
#'   \code{include_hydrogen}, \code{zero_tolerance}, \code{precision}.
#'   Integer arguments are guaranteed to be of class \code{integer}.
#'
#' @keywords internal
validate_sasa_parameters <- function(
    algorithm,
    probe_radius,
    n_slices,
    n_points,
    n_threads,
    polymer_type,
    model,
    include_hydrogen,
    zero_tolerance,
    precision
) {

    # -- algorithm 
    algorithm <- match.arg(algorithm, c("LeeRichards", "ShrakeRupley"))

    # -- probe_radius
    if (!is.numeric(probe_radius) ||
        length(probe_radius) != 1L ||
        !is.finite(probe_radius) ||
        probe_radius <= 0) {
        stop("'probe_radius' must be a finite positive number.", call. = FALSE)
    }
    probe_radius <- as.double(probe_radius)

    # -- n_slices
    n_slices_int <- suppressWarnings(as.integer(n_slices))
    if (!is.numeric(n_slices) ||
        length(n_slices) != 1L ||
        is.na(n_slices_int) ||
        n_slices_int < 1L) {
        stop("'n_slices' must be a positive integer.", call. = FALSE)
    }
    n_slices <- n_slices_int

    # -- n_points
    n_points_int <- suppressWarnings(as.integer(n_points))
    if (!is.numeric(n_points) ||
        length(n_points) != 1L ||
        is.na(n_points_int) ||
        n_points_int < 1L) {
        stop("'n_points' must be a positive integer.", call. = FALSE)
    }
    n_points <- n_points_int

    # -- n_threads
    n_threads_int <- suppressWarnings(as.integer(n_threads))
    if (!is.numeric(n_threads) ||
        length(n_threads) != 1L ||
        is.na(n_threads_int) ||
        n_threads_int < 1L) {
        stop("'n_threads' must be a positive integer.", call. = FALSE)
    }
    n_threads <- n_threads_int
    # FreeSASA's Lee-Richards algorithm is single-threaded only.
    if (algorithm == "LeeRichards" && n_threads > 1L) {
        warning(
            "'n_threads' > 1 is not supported by the LeeRichards algorithm; ",
            "clamping to 1. Use algorithm = \"ShrakeRupley\" for multi-threading.",
            call. = FALSE
        )
        n_threads <- 1L
    }

    # -- polymer_type
    polymer_type <- match.arg(polymer_type, c("RNA", "DNA", "ALL"))

    # -- model
    if (!is.null(model)) {
        if (length(model) != 1L || is.na(model)) {
            stop(
                "'model' must be NULL or a single non-NA model identifier.",
                call. = FALSE
            )
        }
        # Accept numeric (coerce to integer) or character model identifiers.
        if (is.numeric(model)) {
            model_int <- suppressWarnings(as.integer(model))
            if (is.na(model_int) || model_int < 1L) {
                stop(
                    "'model' must be a positive integer or a character string.",
                    call. = FALSE
                )
            }
            model <- model_int
        }
    }

    # -- include_hydrogen
    if (!is.logical(include_hydrogen) ||
        length(include_hydrogen) != 1L ||
        is.na(include_hydrogen)) {
        stop("'include_hydrogen' must be TRUE or FALSE.", call. = FALSE)
    }

    # -- zero_tolerance
    if (!is.numeric(zero_tolerance) ||
        length(zero_tolerance) != 1L ||
        !is.finite(zero_tolerance) ||
        zero_tolerance < 0) {
        stop(
            "'zero_tolerance' must be a finite non-negative number.",
            call. = FALSE
        )
    }
    zero_tolerance <- as.double(zero_tolerance)

    # -- precision
    precision_int <- suppressWarnings(as.integer(precision))
    if (!is.numeric(precision) ||
        length(precision) != 1L ||
        is.na(precision_int) ||
        precision_int < 0L) {
        stop(
            "'precision' must be a non-negative integer.",
            call. = FALSE
        )
    }
    precision <- precision_int

    # Return all validated, coerced values as a list.
    list(
        algorithm        = algorithm,
        probe_radius     = probe_radius,
        n_slices         = n_slices,
        n_points         = n_points,
        n_threads        = n_threads,
        polymer_type     = polymer_type,
        model            = model,
        include_hydrogen = include_hydrogen,
        zero_tolerance   = zero_tolerance,
        precision        = precision
    )
}