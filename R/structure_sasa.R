# structure_sasa.R
#
# Main R entry point for residue-level SASA calculation.
#
# calculate_structure_sasa()  — validates inputs, calls the Python backend,
#                               returns a data.frame of residue annotations.


#' Calculate residue-level SASA from a PDB or mmCIF structure
#'
#' Runs the Python Gemmi + FreeSASA backend and returns a
#' \code{data.frame} of residue-level structure annotations.  The Python
#' environment is located automatically (virtualenv then conda) by
#' \code{\link{load_sasa_python}}.
#'
#' @param structure_file  Character.  Path to a PDB or mmCIF structure
#'   file.
#' @param output_file     Character or \code{NULL}.  Path for the
#'   intermediate TSV written by the Python backend.  When \code{NULL}
#'   (the default) a temporary file is used and deleted on exit.
#' @param overwrite       Logical.  Whether to overwrite \code{output_file}
#'   if it already exists.  Ignored when \code{output_file} is
#'   \code{NULL}.  Default \code{FALSE}.
#' @param algorithm       Character.  SASA algorithm:
#'   \code{"LeeRichards"} (default) or \code{"ShrakeRupley"}.
#' @param probe_radius    Positive numeric.  Probe radius in Ångström.
#'   Default \code{1.4}.
#' @param n_slices        Positive integer.  Number of Lee-Richards
#'   slices.  Default \code{20}.
#' @param n_points        Positive integer.  Number of Shrake-Rupley
#'   points.  Default \code{1000}.
#' @param n_threads       Positive integer.  Number of FreeSASA threads.
#'   Default \code{1}.
#' @param model           \code{NULL} (use the first model) or a single
#'   model identifier (integer index or character name).
#' @param polymer_type    Character.  Residue filter: \code{"RNA"}
#'   (default), \code{"DNA"}, or \code{"ALL"}.
#' @param include_hydrogen Logical.  Include hydrogen atoms.
#'   Default \code{FALSE}.
#' @param zero_tolerance  Non-negative numeric.  SASA values below this
#'   threshold are treated as zero.  Default \code{1e-8}.
#' @param precision       Non-negative integer.  Decimal places in the
#'   TSV output.  Default \code{4}.
#' @param envname         Character.  Name of the Python environment to
#'   use.  Defaults to
#'   \code{getOption("rSASA.python.envname", "sasa-python")}.
#'
#' @return A \code{data.frame} of residue-level structure annotations.
#'   The following attributes are set on the return value:
#'   \describe{
#'     \item{\code{structure_file}}{Absolute path to the structure file.}
#'     \item{\code{python_environment}}{Name of the Python environment.}
#'     \item{\code{sasa_parameters}}{Named list of validated parameter
#'       values actually passed to the backend.}
#'   }
#'
#' @seealso \code{\link{configure_sasa_python}}, \code{\link{check_sasa_python}}
#'
#' @noRd
calculate_structure_sasa <- function(
    structure_file,
    output_file       = NULL,
    overwrite         = FALSE,
    algorithm         = "LeeRichards",
    probe_radius      = 1.4,
    n_slices          = 20L,
    n_points          = 1000L,
    n_threads         = 1L,
    model             = NULL,
    polymer_type      = "RNA",
    include_hydrogen  = FALSE,
    zero_tolerance    = 1e-8,
    precision         = 4L,
    envname           = getOption("ARF.python.envname", "sasa-python")
) {

    # -- Input validation 

    structure_file <- validate_structure_file(structure_file)

    params <- validate_sasa_parameters(
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

    # -- Output file handling -

    temporary_output <- is.null(output_file)

    if (temporary_output) {

        output_file <- tempfile(
            pattern = "rSASA_sasa_",
            fileext = ".tsv"
        )
        on.exit(unlink(output_file), add = TRUE)

    } else {

        output_file <- validate_output_file(output_file, overwrite = overwrite)

    }

    # -- Python backend --

    direct_python <- getOption("ARF.python.path")

    if (!is.null(direct_python)) {
        # Call the CLI directly — no reticulate initialisation, no conda-binary
        # detection. Required when ARF.python.path points to a micromamba or
        # other conda-layout Python that reticulate cannot load via use_python().
        python_dir <- system.file("python", package = .ARF_PKG, mustWork = FALSE)
        if (!nzchar(python_dir) || !dir.exists(python_dir)) {
            candidates <- c(
                file.path(getwd(), "inst", "python"),
                file.path(dirname(getwd()), "inst", "python")
            )
            found <- Find(dir.exists, candidates)
            python_dir <- if (!is.null(found)) found else ""
        }
        if (!nzchar(python_dir))
            stop("Could not locate the ARF Python backend (inst/python).",
                 call. = FALSE)

        python_script <- file.path(python_dir, "structure_sasa.py")
        if (!file.exists(python_script))
            stop("Python backend file not found: ", python_script, call. = FALSE)

        cli_args <- c(
            python_script,
            structure_file,
            "--output",         output_file,
            "--algorithm",      params$algorithm,
            "--probe-radius",   params$probe_radius,
            "--n-slices",       params$n_slices,
            "--n-points",       params$n_points,
            "--n-threads",      params$n_threads,
            "--polymer",        params$polymer_type,
            "--zero-tolerance", params$zero_tolerance,
            "--precision",      params$precision,
            "--overwrite"
        )
        if (!is.null(params$model))
            cli_args <- c(cli_args, "--model", params$model)
        if (isTRUE(params$include_hydrogen))
            cli_args <- c(cli_args, "--include-hydrogen")

        stderr_out <- tempfile()
        on.exit(unlink(stderr_out), add = TRUE)
        ret <- system2(direct_python, args = cli_args,
                       stdout = FALSE, stderr = stderr_out)
        if (ret != 0L) {
            err_msg <- tryCatch(
                paste(readLines(stderr_out), collapse = "\n"),
                error = function(e) "(no stderr captured)"
            )
            stop(
                "The Python SASA backend exited with code ", ret, ":\n",
                err_msg, "\n\n",
                "Tip: run check_sasa_python() to verify the environment.",
                call. = FALSE
            )
        }

    } else {
        py <- load_sasa_python(envname = envname)

        tryCatch(
            py$calculate_rna_sasa(
                structure_file   = structure_file,
                output_file      = output_file,
                algorithm        = params$algorithm,
                probe_radius     = params$probe_radius,
                n_slices         = params$n_slices,
                n_points         = params$n_points,
                n_threads        = params$n_threads,
                model            = params$model,
                polymer_type     = params$polymer_type,
                include_hydrogen = params$include_hydrogen,
                zero_tolerance   = params$zero_tolerance,
                precision        = params$precision
            ),
            error = function(e) {
                stop(
                    "The Python SASA backend raised an error:\n  ",
                    conditionMessage(e), "\n\n",
                    "Tip: run check_sasa_python() to verify the environment, ",
                    "or configure_sasa_python(force = TRUE) to reinstall it.",
                    call. = FALSE
                )
            }
        )
    }

    if (!file.exists(output_file)) {
        stop(
            "The Python SASA backend completed without producing the ",
            "expected output file: ", output_file,
            call. = FALSE
        )
    }

    # -- Read result 

    result <- read_structure_annotations(output_file)

    attr(result, "structure_file")    <- structure_file
    attr(result, "python_environment") <- envname
    attr(result, "sasa_parameters")   <- params

    result
}