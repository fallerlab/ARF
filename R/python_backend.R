# python_backend.R
#
# Manages the Python environment used by ARF.
#
# Supports both virtualenv (default, no conda required) and conda.
# Choose once per session by setting:
#
#   options(ARF.python.env_type = "conda")    # in .Rprofile or at session start
#
# All public functions default to that option, falling back to "virtualenv".


# Package name used in system.file() calls
.ARF_PKG <- "ARF"


# Internal helpers

#' Read the preferred env_type from the package option
#' @noRd
.env_type_default <- function() {
    getOption("ARF.python.env_type", default = "virtualenv")
}

#' Check whether a named Python environment exists
#'
#' Wraps conda_list() in tryCatch so that calling this on a system
#' without conda does not raise an error.
#' @noRd
.env_exists <- function(envname, env_type) {
    if (env_type == "virtualenv") {
        envname %in% reticulate::virtualenv_list()
    } else {
        tryCatch(
            envname %in% reticulate::conda_list()$name,
            error = function(e) FALSE
        )
    }
}

#' Check whether required packages are importable from a Python binary
#'
#' Uses system2() so reticulate is never touched and no conda detection fires.
#' @noRd
.python_has_packages <- function(python, packages = c("gemmi", "freesasa")) {
    # Write imports to a temp file — avoids -c and shell semicolon splitting.
    tmp <- tempfile(fileext = ".py")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(paste0("import ", packages), tmp)
    system2(python, args = tmp, stdout = FALSE, stderr = FALSE) == 0L
}

#' Activate a named Python environment
#' @noRd
.use_env <- function(envname, env_type) {
    if (env_type == "virtualenv") {
        reticulate::use_virtualenv(envname, required = TRUE)
    } else {
        reticulate::use_condaenv(envname, required = TRUE)
    }
}

#' Remove a named Python environment
#' @noRd
.remove_env <- function(envname, env_type) {
    if (env_type == "virtualenv") {
        reticulate::virtualenv_remove(envname, confirm = FALSE)
    } else {
        # reticulate does not expose a conda_remove() wrapper in all versions;
        # call the conda binary directly.
        conda_bin <- tryCatch(
            reticulate::conda_binary(),
            error = function(e) {
                stop(
                    "Could not find the conda binary. ",
                    "Is conda installed and on PATH?",
                    call. = FALSE
                )
            }
        )
        system2(
            conda_bin,
            args   = c("env", "remove", "--name", envname, "--yes", "--quiet"),
            stdout = FALSE,
            stderr = FALSE
        )
    }
}


# Public API

#' Configure the Python backend used by ARF
#'
#' Creates a dedicated Python virtual environment or conda environment
#' containing Gemmi and FreeSASA.
#'
#' @param envname        Character. Environment name. Default: \code{"sasa-python"}.
#' @param python_version Character. Python version to install via pyenv when
#'   \code{python} is \code{NULL}. Default: \code{"3.11"}.
#' @param python         Character or \code{NULL}. Path to an existing Python
#'   binary (e.g. from a micromamba/conda/system environment). When supplied:
#'   if \code{gemmi} and \code{freesasa} are already importable from that
#'   Python, ARF registers it directly and skips virtualenv creation entirely;
#'   otherwise a virtualenv is created from it and the packages are pip-installed.
#'   Default: \code{NULL}.
#' @param env_type       Character. \code{"virtualenv"} (default) or \code{"conda"}.
#' @param force          Logical. Ignore cached configuration and re-run setup.
#'   Default: \code{FALSE}.
#'
#' @return The Python path or environment name, invisibly.
#'
#' @examples
#' \dontrun{
#' # Option A — fully automatic (downloads Python via pyenv)
#' configure_sasa_python()
#'
#' # Option B — use an existing Python (packages already installed)
#' configure_sasa_python(python = "/path/to/python3")
#'
#' # Option B — use an existing Python, create a venv from it if needed
#' configure_sasa_python(python = "/path/to/python3", envname = "sasa-python")
#' }
#'
#' @export
configure_sasa_python <- function(
    envname        = "sasa-python",
    python_version = "3.11",
    python         = NULL,
    env_type       = .env_type_default(),
    force          = FALSE
) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package 'reticulate' is required. Install it with:\n",
            "  install.packages('reticulate')",
            call. = FALSE
        )
    }

    # -- Direct Python path supplied
    if (!is.null(python)) {
        if (!file.exists(python))
            stop("'python' path does not exist: ", python, call. = FALSE)

        already_configured <- identical(getOption("ARF.python.path"), python)
        if (already_configured && !force) {
            message("Python '", python, "' is already registered. ",
                    "Use force = TRUE to re-check.")
            return(invisible(python))
        }

        # Check whether gemmi + freesasa are already importable — via system2()
        # so reticulate is never touched and no conda-binary detection fires.
        message("Checking for gemmi and freesasa in: ", python)
        if (.python_has_packages(python)) {
            message("Packages found — registering Python directly (no venv needed).")
            options(ARF.python.path = python)
            return(invisible(python))
        }

        # Packages not present — create a virtualenv from this Python via the
        # stdlib venv module, bypassing reticulate's conda-meta detection.
        message("Packages not found — creating virtualenv '", envname,
                "' from: ", python)
        venv_path <- file.path(reticulate::virtualenv_root(), envname)
        ret <- system2(python, args = c("-m", "venv", venv_path),
                       stdout = FALSE, stderr = FALSE)
        if (ret != 0L)
            stop("Failed to create virtualenv at: ", venv_path, call. = FALSE)
        reticulate::virtualenv_install(
            envname          = envname,
            packages         = c("gemmi", "freesasa"),
            ignore_installed = FALSE
        )
        message("Environment '", envname, "' is ready.")
        return(invisible(envname))
    }

    # -- No Python path — use reticulate's managed virtualenv or conda
    env_type   <- match.arg(env_type, c("virtualenv", "conda"))
    env_exists <- .env_exists(envname, env_type)

    if (env_exists && !force) {
        options(ARF.python.env_type = env_type)
        message("Environment '", envname, "' (", env_type, ") already exists. ",
                "Use force = TRUE to reinstall.")
        return(invisible(envname))
    }

    if (force && env_exists) {
        if (interactive()) {
            answer <- readline(sprintf(
                "Remove existing '%s' environment '%s'? [y/N] ", env_type, envname
            ))
            if (!tolower(trimws(answer)) %in% c("y", "yes")) {
                message("Aborted — environment left unchanged.")
                return(invisible(envname))
            }
        }
        message("Removing existing environment: '", envname, "'")
        .remove_env(envname, env_type)
    }

    message("Creating ", env_type, " environment '", envname,
            "' (Python ", python_version, ") ...")

    if (env_type == "virtualenv") {
        reticulate::install_python(version = python_version)
        reticulate::virtualenv_create(envname = envname, python = python_version)
        reticulate::virtualenv_install(
            envname          = envname,
            packages         = c("gemmi", "freesasa"),
            ignore_installed = FALSE
        )
    } else {
        reticulate::conda_create(
            envname        = envname,
            python_version = python_version
        )
        reticulate::conda_install(
            envname  = envname,
            packages = c("gemmi", "freesasa"),
            channel  = "conda-forge",
            pip      = FALSE
        )
    }

    options(ARF.python.env_type = env_type)
    message("Environment '", envname, "' is ready.")
    invisible(envname)
}


#' Check the ARF Python backend
#'
#' Verifies that the Python environment exists and that Gemmi and FreeSASA
#' are importable. Prints a one-line summary and returns the status invisibly.
#'
#' @param envname  Character. Environment name. Default: \code{"sasa-python"}.
#' @param env_type Character. \code{"virtualenv"} or \code{"conda"}.
#'   Defaults to \code{getOption("ARF.python.env_type", "virtualenv")}.
#'
#' @return A named list with elements \code{environment}, \code{env_type},
#'   \code{available}, \code{python}, \code{gemmi}, \code{freesasa}.
#'   Returned invisibly so the function can be used for its printed output.
#'
#' @export
check_sasa_python <- function(
    envname  = "sasa-python",
    env_type = .env_type_default()
) {
    ok <- function(x) if (isTRUE(x)) "[OK]" else "[MISSING]"

    # When a direct Python path is registered (e.g. micromamba), check via
    # system2() — no reticulate initialisation, no conda-binary detection.
    direct_python <- getOption("ARF.python.path")
    if (!is.null(direct_python)) {
        gemmi_ok    <- .python_has_packages(direct_python, "gemmi")
        freesasa_ok <- .python_has_packages(direct_python, "freesasa")
        result <- list(
            environment = direct_python,
            env_type    = "direct",
            available   = gemmi_ok && freesasa_ok,
            python      = direct_python,
            gemmi       = gemmi_ok,
            freesasa    = freesasa_ok
        )
        message(
            "Python      : ", direct_python, " (registered via ARF.python.path)\n",
            "gemmi       : ", ok(gemmi_ok), "\n",
            "freesasa    : ", ok(freesasa_ok)
        )
        return(invisible(result))
    }

    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package 'reticulate' is required. Install it with:\n",
            "  install.packages('reticulate')",
            call. = FALSE
        )
    }

    env_type <- match.arg(env_type, c("virtualenv", "conda"))

    result <- list(
        environment = envname,
        env_type    = env_type,
        available   = FALSE,
        python      = NULL,
        gemmi       = FALSE,
        freesasa    = FALSE
    )

    if (!.env_exists(envname, env_type)) {
        msg <- sprintf(
            "Environment '%s' (%s) not found. Run configure_sasa_python().",
            envname, env_type
        )
        message(msg)
        result$reason <- msg
        return(invisible(result))
    }

    tryCatch({
        .use_env(envname, env_type)
        result$available <- TRUE
        result$python    <- reticulate::py_config()$python
        result$gemmi     <- reticulate::py_module_available("gemmi")
        result$freesasa  <- reticulate::py_module_available("freesasa")
    }, error = function(e) {
        result$reason <<- conditionMessage(e)
        message("Failed to activate environment: ", conditionMessage(e))
    })

    message(
        "Environment : ", envname, " (", env_type, ")\n",
        "Python      : ", if (is.null(result$python)) "(unavailable)" else result$python, "\n",
        "gemmi       : ", ok(result$gemmi), "\n",
        "freesasa    : ", ok(result$freesasa)
    )

    invisible(result)
}


#' Load the ARF Python structure module
#'
#' Activates the Python environment and returns the imported
#' \code{structure_sasa} module. The environment type (virtualenv or conda)
#' is detected automatically from what was created, so this function does
#' not need an \code{env_type} argument.
#'
#' When \code{options(ARF.python.path)} is set (direct micromamba/conda Python),
#' \code{calculate_structure_sasa()} calls the Python CLI via \code{system2()}
#' and never calls this function; this function is only reached when a managed
#' virtualenv or conda environment is used.
#'
#' Note: Python can only be initialised once per R session. If a different
#' Python environment was already loaded, restart R and try again.
#'
#' @param envname Character. Environment name. Default: \code{"sasa-python"}.
#'
#' @return The imported \code{structure_sasa} Python module.
#'
#' @keywords internal
load_sasa_python <- function(envname = "sasa-python") {

    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package 'reticulate' is required.",
            call. = FALSE
        )
    }

    # Auto-detect environment type so users don't have to specify it here.
    # Note: when ARF.python.path is set (direct micromamba/conda Python),
    # calculate_structure_sasa() calls the CLI via system2() and never reaches
    # this function, so no conda-binary detection fires in that path.
    in_virtualenv <- envname %in% reticulate::virtualenv_list()
    in_conda      <- tryCatch(
        envname %in% reticulate::conda_list()$name,
        error = function(e) FALSE
    )

    if (in_virtualenv) {
        reticulate::use_virtualenv(envname, required = TRUE)
    } else if (in_conda) {
        reticulate::use_condaenv(envname, required = TRUE)
    } else {
        stop(
            "Python environment '", envname, "' not found.\n",
            "Run configure_sasa_python() to create it.",
            call. = FALSE
        )
    }

    # Verify that both required modules are importable before proceeding.
    missing_modules <- Filter(
        Negate(reticulate::py_module_available),
        c("gemmi", "freesasa")
    )

    if (length(missing_modules) > 0) {
        stop(
            "Required Python module(s) not available in environment '", envname, "': ",
            paste(missing_modules, collapse = ", "), ".\n",
            "Run configure_sasa_python(force = TRUE) to reinstall.",
            call. = FALSE
        )
    }

    # Locate inst/python — works both from an installed package and during
    # package development (devtools::load_all()).
    python_dir <- system.file("python", package = .ARF_PKG, mustWork = FALSE)

    if (!nzchar(python_dir) || !dir.exists(python_dir)) {
        # Development fallback: search for inst/python relative to getwd()
        # (works when the working directory is the package root or a subdirectory).
        candidates <- c(
            file.path(getwd(), "inst", "python"),
            file.path(dirname(getwd()), "inst", "python")
        )
        found <- Find(dir.exists, candidates)
        python_dir <- if (!is.null(found)) found else ""
    }

    if (!nzchar(python_dir)) {
        stop(
            "Could not locate the ARF Python backend.\n",
            "Expected file: inst/python/structure_sasa.py",
            call. = FALSE
        )
    }

    python_file <- file.path(python_dir, "structure_sasa.py")

    if (!file.exists(python_file)) {
        stop(
            "Python backend file not found: ", python_file, "\n",
            "The package installation may be incomplete.",
            call. = FALSE
        )
    }

    reticulate::import_from_path("structure_sasa", path = python_dir)
}