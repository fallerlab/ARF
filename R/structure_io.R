# structure_io.R
#
# I/O helpers for rSASA structure annotation files.
#
# read_structure_annotations()  — reads a TSV produced by the Python backend,
#                                 validates columns, coerces types, and runs
#                                 basic sanity checks on the data.


# Constants

#' Required columns and their expected R types
#'
#' Used by \code{read_structure_annotations()} to validate and coerce the
#' TSV produced by the Python SASA backend.  Any column listed here that
#' is absent from the file raises an error.  Columns present in the file
#' but absent from this list are kept as-is.
#'
#' @noRd
.ANNOTATION_COLTYPES <- list(
    structure_id       = as.character,
    model              = as.character,
    entity_id          = as.character,
    entity_type        = as.character,
    entity_description = as.character,
    auth_chain         = as.character,
    auth_comp_id       = as.character,
    auth_seq_id        = as.character,   # keep as character: may contain ins. codes
    sasa_A2            = as.double
)

#' Names of columns that must be present in a valid annotation file
#' @noRd
.REQUIRED_COLUMNS <- names(.ANNOTATION_COLTYPES)


# Public function

#' Read structure annotation output
#'
#' Reads the residue-level TSV file produced by the rSASA Python backend,
#' validates that all required columns are present, coerces each column to
#' its expected type, and performs basic sanity checks on the result.
#'
#' The function prefers \pkg{data.table} (\code{\link[data.table]{fread}})
#' when the package is available, falling back to
#' \code{\link[utils]{read.delim}} otherwise.
#'
#' @param file Character scalar.  Path to the annotation TSV.
#'
#' @return A \code{data.frame} with at least the columns listed in
#'   \code{.REQUIRED_COLUMNS}.  Additional columns written by the backend
#'   are preserved as-is.
#'
#' @noRd
read_structure_annotations <- function(file) {

    # -- Path checks
    if (!is.character(file) ||
        length(file) != 1L ||
        is.na(file) ||
        !nzchar(file)) {
        stop(
            "'file' must be a single non-empty character string.",
            call. = FALSE
        )
    }

    if (!file.exists(file)) {
        stop(
            "Annotation file does not exist: ", file,
            call. = FALSE
        )
    }

    # Catch empty files before trying to parse them.
    if (file.info(file)$size == 0L) {
        stop(
            "Annotation file is empty: ", file,
            call. = FALSE
        )
    }

    # -- Read
    na_strings <- c("", "NA", "N/A", ".")

    if (requireNamespace("data.table", quietly = TRUE)) {
        x <- data.table::fread(
            file,
            sep         = "\t",
            header      = TRUE,
            na.strings  = na_strings,
            data.table  = FALSE,
            showProgress = FALSE
        )
    } else {
        x <- utils::read.delim(
            file,
            sep             = "\t",
            header          = TRUE,
            na.strings      = na_strings,
            stringsAsFactors = FALSE,
            check.names     = FALSE
        )
    }

    # -- Column presence
    missing_cols <- setdiff(.REQUIRED_COLUMNS, colnames(x))

    if (length(missing_cols) > 0L) {
        stop(
            "Annotation file is missing required columns: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }

    # -- Type coercions
    for (col in .REQUIRED_COLUMNS) {
        coerce_fn <- .ANNOTATION_COLTYPES[[col]]
        x[[col]]  <- coerce_fn(x[[col]])
    }

    # -- Sanity checks
    if (nrow(x) == 0L) {
        warning(
            "Annotation file contains no data rows: ", file,
            call. = FALSE
        )
        return(x)
    }

    n_negative_sasa <- sum(x[["sasa_A2"]] < 0, na.rm = TRUE)

    if (n_negative_sasa > 0L) {
        warning(
            n_negative_sasa,
            " residue(s) have negative sasa_A2 values. ",
            "This may indicate a problem with the structure file or SASA calculation.",
            call. = FALSE
        )
    }

    n_na_sasa <- sum(is.na(x[["sasa_A2"]]))

    if (n_na_sasa > 0L) {
        warning(
            n_na_sasa,
            " residue(s) have missing (NA) sasa_A2 values.",
            call. = FALSE
        )
    }

    x
}