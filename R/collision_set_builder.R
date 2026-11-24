# Copyright (C) 2026  Ferhat Alkan & Edwin Sakyi Kyei-Baffour
#
#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================
# ARF — Collision Interface Analysis
# =============================================================================
#
# Computes SASA-based and 3D proximity-based collision-interface gene sets for
# GSEA analysis of ribosomal protein binding to rRNA.
#
# Main entry point: ARF_build_collision_sets()
#
# Required packages: SASA, dplyr, tidyr, Biostrings, bio3d, msa, readr,
#                    reshape2
# =============================================================================

# Suppress R CMD CHECK notes arising from dplyr/tidyr column-name pronoun use.
utils::globalVariables(c(
    "entity_description", "auth_seq_id", "auth_chain", "rRNA",
    "sasa_A2", "sasa_A2_diff", "n", "x", "y", "z", "RP",
    "protein_chain", "resno", "RNA_molecule", "distance",
    "Col.Int.", "focus_pos", "gene", "ont", "sasa_A2_diff_z", "."
))


# Internal helpers

#' Vectorised pairwise Euclidean distance matrix
#'
#' Uses the identity \eqn{||a - b||^2 = ||a||^2 + ||b||^2 - 2 a'b} so
#' the cost is one \code{tcrossprod()} call rather than N row-wise loops.
#'
#' @param A Numeric matrix, n × 3.
#' @param B Numeric matrix, m × 3.
#' @return Numeric matrix, n × m.
#' @noRd
.pairwise_euclidean <- function(A, B) {
    A       <- as.matrix(A)
    B       <- as.matrix(B)
    sq_A    <- rowSums(A^2)
    sq_B    <- rowSums(B^2)
    cross   <- tcrossprod(A, B)                 # n × m
    sq_dist <- outer(sq_A, sq_B, "+") - 2 * cross
    sqrt(pmax(sq_dist, 0))                      # pmax guards floating-point negatives
}

# Chunked variant that returns only the per-row minimum distance; avoids
# materialising the full N×M matrix when M is large (e.g. disome cross-copy).
.pairwise_euclidean_rowmin <- function(A, B, batch = 500L) {
    A     <- as.matrix(A)
    B     <- as.matrix(B)
    sq_B  <- rowSums(B^2)
    n     <- nrow(A)
    min_d <- numeric(n)
    for (i in seq(1L, n, by = batch)) {
        j      <- min(i + batch - 1L, n)
        Ai     <- A[i:j, , drop = FALSE]
        sq_Ai  <- rowSums(Ai^2)
        sq_d   <- outer(sq_Ai, sq_B, "+") - 2 * tcrossprod(Ai, B)
        min_d[i:j] <- apply(sqrt(pmax(sq_d, 0)), 1L, min)
    }
    min_d
}

#' Circular-shift random set (1-indexed, no off-by-one)
#'
#' Shifts \code{positions} by \code{shift}, wrapping within \code{n}.
#' Result is always in \code{[1, n]}.
#'
#' @noRd
.circular_shift <- function(positions, shift, n) {
    (positions - 1L + shift) %% n + 1L
}

#' Build 99 circular-shift randomisation sets for GSEA
#'
#' @param proxpos Integer vector of foreground (1-indexed) row positions.
#' @param n       Total number of rows (denominator for modular arithmetic).
#' @param ont_prefix Character; prefix for "Rand<i>_<ont_prefix>" labels.
#' @param row_labels Character vector of length \code{n}; gene labels in row order.
#' @return data.frame with columns \code{ont} and \code{gene}.
#' @noRd
.make_rand_sets <- function(proxpos, n, ont_prefix, row_labels) {
    step  <- max(1L, round(n / 100L) - 1L)
    dplyr::bind_rows(lapply(seq_len(99L), function(i) {
        randset <- .circular_shift(proxpos, i * step, n)
        data.frame(
            ont  = paste0("Rand", i, "_", ont_prefix),
            gene = row_labels[randset],
            stringsAsFactors = FALSE
        )
    }))
}

# 1a. ARF_download_files — resolve or download structure and FASTA files

#' Resolve or download a PDB/mmCIF structure file and rRNA FASTA
#'
#' If local file paths are supplied they are returned as-is (with an existence
#' check).  Otherwise the files are downloaded from RCSB into
#' \code{download_directory}.
#'
#' @param PDBid              Character. Four-letter PDB identifier.
#' @param PDB_file           Character or \code{NULL}. Path to a local
#'   structure file; if \code{NULL} the mmCIF is downloaded from RCSB.
#' @param fasta_file         Character or \code{NULL}. Path to a local
#'   FASTA file; if \code{NULL} it is downloaded from RCSB.
#' @param download_directory Character. Directory for downloaded files.
#'   Default: \code{getwd()}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{PDB_file}}{Resolved path to the structure file.}
#'     \item{\code{fasta_file}}{Resolved path to the FASTA file.}
#'   }
#'
#' @export
ARF_download_structure_files <- function(
    PDBid,
    PDB_file = NULL,
    fasta_file = NULL,
    download_directory = getwd()
) {
  stopifnot(
    is.character(PDBid), nzchar(PDBid),
    is.character(download_directory), nzchar(download_directory)
  )
  
  # -- Structure file
  if (is.null(PDB_file)) {
    PDB_file <- file.path(download_directory, paste0(PDBid, ".cif"))
    message("Downloading structure: ", PDBid, ".cif")
    utils::download.file(
      url = paste0("https://files.rcsb.org/download/", PDBid, ".cif"),
      destfile = PDB_file,
      quiet = TRUE
    )
  } else {
    if (!file.exists(PDB_file))
      stop("PDB_file does not exist: ", PDB_file)
  }
  
  # -- rRNA FASTA
  if (is.null(fasta_file)) {
    fasta_file <- file.path(download_directory, paste0(PDBid, ".fasta"))
    message("Downloading rRNA FASTA: ", PDBid, ".fasta")
    utils::download.file(
      url = paste0("https://www.rcsb.org/fasta/entry/", PDBid),
      destfile = fasta_file,
      quiet = TRUE
    )
  } else {
    if (!file.exists(fasta_file))
      stop("fasta_file does not exist: ", fasta_file)
  }

  list(PDB_file = PDB_file, fasta_file = fasta_file)
}



# 1b. ARF_run_SASA — configure Python backend and compute residue-level SASA

#' Configure the Python backend and compute residue-level SASA
#'
#' Calls \code{\link{ARF_download_files}} to resolve the structure and FASTA
#' paths, configures the Python environment (no-op if already done), then
#' delegates to \code{\link[ARF]{calculate_structure_sasa}}.
#'
#' @inheritParams ARF_download_structure_files
#' @param output_file        Character or \code{NULL}. Path for the SASA TSV
#'   written by the Python backend.  Must be non-\code{NULL} when
#'   \code{ARF_build_SASA_candidates()} will read it afterwards.
#' @param overwrite          Logical. Allow overwriting \code{output_file}.
#' @param algorithm          SASA algorithm: \code{"LeeRichards"} (default) or
#'   \code{"ShrakeRupley"}.
#' @param probe_radius       Probe radius in Å. Default \code{10} (coarse
#'   interface analysis).
#' @param n_slices,n_points,n_threads  FreeSASA tuning parameters.
#' @param model              \code{NULL} or a model identifier.
#' @param polymer_type       \code{"RNA"}, \code{"DNA"}, or \code{"ALL"}
#'   (default) — \code{"ALL"} captures the full complex SASA.
#' @param include_hydrogen   Logical.
#' @param zero_tolerance     Numeric. Values below this are treated as zero.
#' @param precision          Integer. Decimal places in the TSV output.
#' @param envname            Python environment name. Defaults to
#'   \code{getOption("ARF.python.envname", "sasa-python")}.
#' @param env_type           \code{"virtualenv"} or \code{"conda"}.
#'   Defaults to \code{getOption("ARF.python.env_type", "virtualenv")}.
#' @param python_version     Python version string for environment creation.
#' @param force              Logical. Force-reinstall the Python environment.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{annotations}}{data.frame of residue-level SASA values.}
#'     \item{\code{PDB_file}}{Resolved path to the structure file.}
#'     \item{\code{fasta_file}}{Resolved path to the FASTA file.}
#'   }
#'
#' @export
ARF_run_SASA <- function(
    PDBid,
    PDB_file = NULL,
    rRNAs_file = NULL,
    download_directory = getwd(),
    output_file = NULL,
    overwrite = FALSE,
    algorithm = "LeeRichards",
    probe_radius = 10,
    n_slices = 20L,
    n_points = 1000L,
    n_threads = 1L,
    model = NULL,
    polymer_type = "ALL",
    include_hydrogen = FALSE,
    zero_tolerance = 1e-8,
    precision = 4L,
    envname = getOption("ARF.python.envname", "sasa-python"),
    env_type = getOption("ARF.python.env_type", "virtualenv"),
    python_version = "3.11",
    force = FALSE
) {
    # -- Resolve / download files
    files <- ARF_download_structure_files(PDBid, PDB_file, rRNAs_file, download_directory)
    PDB_file <- files$PDB_file
    fasta_file <- files$fasta_file
    
    # -- Configure Python environment (skip if a direct Python is already registered)
    if (is.null(getOption("ARF.python.path"))) {
        ARF::configure_sasa_python(
            envname        = envname,
            python_version = python_version,
            env_type       = env_type,
            force          = force
        )
    }

    # -- Run SASA
    message("Running SASA calculation on: ", basename(PDB_file))
    annotations <- calculate_structure_sasa(
        structure_file   = PDB_file,
        output_file      = output_file,
        overwrite        = overwrite,
        algorithm        = algorithm,
        probe_radius     = probe_radius,
        n_slices         = n_slices,
        n_points         = n_points,
        n_threads        = n_threads,
        model            = model,
        polymer_type     = polymer_type,
        include_hydrogen = include_hydrogen,
        zero_tolerance   = zero_tolerance,
        precision        = precision,
        envname          = envname
    )

    invisible(list(
        annotations = annotations,
        PDB_file    = PDB_file,
        fasta_file  = fasta_file
    ))
}


# 2. ARF_build_SASA_candidates — filter SASA output for rRNA residues

#' Build SASA candidate residue set from SASA annotations
#'
#' Reads the TSV produced by \code{\link{ARF_run_SASA}}, labels each residue
#' with its rRNA identity, and returns residues where the z-transformed
#' differential SASA (max - min across chains) exceeds 1.64 standard deviations (z > 1.64 on raw diff).
#'
#' @param output_file   Path to the SASA TSV (written by \code{ARF_run_SASA}).
#' @param rRNAs         Character vector of rRNA entity description names.
#' @param target_species Character. Species prefix for the rRNA labels.
#'
#' @return data.frame with columns: \code{rRNA} (rRNA identity prefixed with
#'   \code{target_species}), \code{auth_seq_id} (PDB residue number),
#'   \code{sasa_A2_diff} (max − min SASA across chains, Å²), \code{n} (number
#'   of chains contributing to that residue's range), \code{sasa_A2_diff_z}
#'   (z-score of \code{sasa_A2_diff} within the residue's rRNA; rows with
#'   z > 1.64 are retained as Rib.Col. candidates).
#'
#' @export
ARF_build_SASA_candidates <- function(output_file, rRNAs, target_species) {

    if (!file.exists(output_file)) {
        stop(
            "SASA output file not found: ", output_file, "\n",
            "Run ARF_run_SASA() first.",
            call. = FALSE
        )
    }

    # Read annotations (uses SASA's validated reader when available)
    sasa_df <- if (exists("read_structure_annotations",
                         where = asNamespace("ARF"), mode = "function")) {
        read_structure_annotations(output_file)
    } else {
        utils::read.delim(output_file, sep = "\t", stringsAsFactors = FALSE,
                          check.names = FALSE)
    }

    # -- Normalise entity descriptions and label rRNA in one pass 
    sasa_df <- sasa_df |>
        dplyr::mutate(
            entity_description = gsub(
                " ", "_", gsub("^ | $|'", "", entity_description)
            ),
            rRNA = dplyr::case_when(
                entity_description %in% rRNAs ~ entity_description,
                TRUE                          ~ NA_character_
            )
        )

    # -- Compute per-residue SASA range across chains
    rRNA_rows <- sasa_df |>
        dplyr::select(entity_description, auth_seq_id, auth_chain, rRNA, sasa_A2) |>
        dplyr::filter(!is.na(rRNA), !is.na(sasa_A2))

    if (nrow(rRNA_rows) == 0L) {
        return(data.frame(rRNA = character(), auth_seq_id = character(),
                          sasa_A2_diff = numeric(), stringsAsFactors = FALSE))
    }

    temp_sasa_df <- rRNA_rows |>
        dplyr::group_by(rRNA, auth_seq_id) |>
        dplyr::reframe(
            sasa_A2_diff = max(sasa_A2) - min(sasa_A2),
            n            = dplyr::n()
        ) |>
        as.data.frame()

    # Residues with sasa_A2_diff > 0 are differentially exposed across chains —
    # in a collision disome one copy is buried at the interface while the other
    # is not. z > 1.64 on the raw differential SASA yields perfect precision
    # against validated reference sets across disome structures.
    sasa_set <- temp_sasa_df |>
        dplyr::filter(!is.na(sasa_A2_diff), sasa_A2_diff > 0) |>
        dplyr::mutate(sasa_A2_diff_z = as.numeric(scale(sasa_A2_diff))) |>
        dplyr::filter(sasa_A2_diff_z > 1.64)

    sasa_set |>
        dplyr::mutate(rRNA = paste0(target_species, "_", rRNA))
}


# 3. ARF_build_ColInt_proximity — 3-D RP–rRNA distance matrix

#' Compute RP–rRNA 3D proximity distances from a structure file
#'
#' For every RNA N1 atom, finds the closest ribosomal-protein N1 atom per
#' chain and summarises by \code{duplication_resolve} (default: \code{min}).
#' Distance computation uses a vectorised n × m matrix operation rather than
#' row-by-row dplyr loops, which is substantially faster for large structures.
#'
#' @param rRNAs_file          Path to the rRNA FASTA file.
#' @param PDB_file            Path to the mmCIF structure file.
#' @param conversion_table    data.frame with columns \code{RP_new} and
#'   \code{chainID} (comma-separated).
#' @param PDBid               Four-letter PDB identifier (used in output filenames).
#' @param out_prefix          File-name prefix for CSV outputs.
#' @param target_species      Character. Species prefix for RP labels.
#' @param duplication_resolve Function to resolve multiple chain observations.
#'   Default: \code{min}.
#' @param n_copies Integer. \code{2L} (default) treats the structure as a
#'   disome: chains are automatically partitioned into two ribosome copies via
#'   k-means clustering of chain centroids (k = 2), and cross-copy distances
#'   (rRNA N1 to RP CA + rRNA N1 of the other copy) are stored in a single
#'   \code{"Col.Int."} column.  Set to \code{1L} to use the monosome
#'   RP-proximity method instead.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{RP_proximity_df}}{Wide data.frame: rows = rRNA residues,
#'       columns = \code{rRNA}, \code{resno}, one column per RP (or one
#'       \code{"Col.Int."} column for \code{n_copies = 2}).}
#'     \item{\code{output_file}}{Path to the long-format distance CSV.}
#'   }
#'
#' @noRd
ARF_build_ColInt_proximity <- function(
    rRNAs_file,
    PDB_file,
    conversion_table,
    PDBid,
    out_prefix,
    target_species,
    duplication_resolve = min,
    n_copies            = 2L
) {
    message("Starting 3D analysis ...")

    rRNA_seq_set <- Biostrings::readBStringSet(rRNAs_file)
    chainData    <- conversion_table

    # -- Expand comma-separated chain IDs into one row per chain
    chainData_resolved <- chainData |>
        tidyr::separate_rows(chainID, sep = ",") |>
        dplyr::mutate(chainID = trimws(chainID)) |>
        dplyr::rename(chain = chainID, RP = RP_new) |>
        dplyr::mutate(RP = paste0(target_species, "_", gsub(" ", "_", RP)))

    # -- Read and annotate the structure
    PDB_data      <- bio3d::read.cif(PDB_file)
    PDB_dataframe <- as.data.frame(PDB_data[["atom"]], stringsAsFactors = FALSE) |>
        dplyr::left_join(chainData_resolved, by = "chain") |>
        dplyr::mutate(x = as.numeric(x), y = as.numeric(y), z = as.numeric(z))

    RNA_atom <- "N1"

    # Protein atoms: one N1 per residue acts as the "contact atom"
    alphaCarbons <- PDB_dataframe |>
        dplyr::filter(elety == RNA_atom, !is.na(RP), RP != "") |>
        dplyr::select(x, y, z, RP, protein_chain = chain)

    # RNA atoms: one N1 per nucleotide, restricted to annotated rRNA chains
    rRNA_data <- PDB_dataframe |>
        dplyr::filter(
            elety == RNA_atom,
            RP %in% names(rRNA_seq_set)
        ) |>
        dplyr::select(x, y, z, RP, resno, RNA_chain = chain)

    if (nrow(alphaCarbons) == 0L) {
        stop("No protein N1 atoms found after chain mapping. ",
             "Check 'conversion_table' chain IDs.", call. = FALSE)
    }
    if (nrow(rRNA_data) == 0L) {
        stop("No RNA N1 atoms found matching the rRNA FASTA names. ",
             "Check that 'rRNAs_file' names match entity descriptions in PDB.", call. = FALSE)
    }

    # ---- Disome path (n_copies == 2) ----------------------------------------
    if (n_copies == 2L) {
        # Cluster on PROTEIN chains only: rRNA copies face each other at the
        # collision interface and their centroids can be nearly co-located,
        # making k-means unreliable when rRNA is included. Protein chains are
        # distributed all around each copy and cluster cleanly. rRNA chains are
        # then assigned to the nearest protein-cluster centre.
        prot_chains <- PDB_dataframe |>
            dplyr::filter(elety == "CA",
                          !is.na(RP), RP != "",
                          !RP %in% names(rRNA_seq_set)) |>
            dplyr::pull(chain) |>
            unique()

        prot_cents <- PDB_dataframe |>
            dplyr::filter(chain %in% prot_chains) |>
            dplyr::mutate(x = as.numeric(x), y = as.numeric(y), z = as.numeric(z)) |>
            dplyr::group_by(chain) |>
            dplyr::summarise(cx = mean(x, na.rm = TRUE),
                             cy = mean(y, na.rm = TRUE),
                             cz = mean(z, na.rm = TRUE), .groups = "drop") |>
            dplyr::filter(is.finite(cx), is.finite(cy), is.finite(cz))

        if (nrow(prot_cents) < 2L)
            stop("Fewer than 2 protein chains with known atoms found. ",
                 "Cannot cluster into 2 copies.", call. = FALSE)

        set.seed(42L)
        km       <- stats::kmeans(as.matrix(prot_cents[, c("cx", "cy", "cz")]),
                                  centers = 2L, nstart = 20L)
        prot_copy_map <- stats::setNames(km$cluster, prot_cents$chain)

        # Assign rRNA chains by nearest protein cluster centre
        rRNA_chains <- unique(rRNA_data$RNA_chain)
        rRNA_cents  <- PDB_dataframe |>
            dplyr::filter(chain %in% rRNA_chains) |>
            dplyr::mutate(x = as.numeric(x), y = as.numeric(y), z = as.numeric(z)) |>
            dplyr::group_by(chain) |>
            dplyr::summarise(cx = mean(x, na.rm = TRUE),
                             cy = mean(y, na.rm = TRUE),
                             cz = mean(z, na.rm = TRUE), .groups = "drop") |>
            dplyr::filter(is.finite(cx), is.finite(cy), is.finite(cz))

        ctrs <- km$centers  # 2 × 3 matrix of protein cluster centres
        rRNA_copy <- apply(as.matrix(rRNA_cents[, c("cx", "cy", "cz")]), 1L,
            function(p) which.min(rowSums((ctrs - matrix(p, 2L, 3L, byrow = TRUE))^2))
        )
        rRNA_copy_map <- stats::setNames(rRNA_copy, rRNA_cents$chain)

        copy_map <- c(prot_copy_map, rRNA_copy_map)
        message(sprintf(
            "Chain clustering (protein-guided): copy 1 = %d chains, copy 2 = %d chains",
            sum(copy_map == 1L), sum(copy_map == 2L)
        ))

        rRNA_data$copy <- copy_map[rRNA_data$RNA_chain]

        prot_CA <- PDB_dataframe |>
            dplyr::filter(elety == "CA",
                          !is.na(RP), RP != "",
                          !RP %in% names(rRNA_seq_set)) |>
            dplyr::mutate(copy = copy_map[chain]) |>
            dplyr::filter(!is.na(copy)) |>
            dplyr::select(x, y, z, copy)

        if (nrow(prot_CA) == 0L)
            stop("No protein CA atoms found for n_copies = 2. ",
                 "Check 'conversion_table'.", call. = FALSE)

        # Same-copy threshold: 5th percentile of rRNA N1 → RP CA distances
        same_copy_min_dists <- unlist(lapply(c(1L, 2L), function(cp) {
            rna_cp <- rRNA_data[!is.na(rRNA_data$copy) & rRNA_data$copy == cp, ]
            rp_cp  <- prot_CA[prot_CA$copy == cp, ]
            if (nrow(rna_cp) == 0L || nrow(rp_cp) == 0L) return(numeric(0L))
            .pairwise_euclidean_rowmin(rna_cp[, c("x", "y", "z")],
                                       rp_cp[, c("x", "y", "z")])
        }))
        dist_thr_same <- stats::quantile(same_copy_min_dists, 0.05, na.rm = TRUE)
        message(sprintf("Same-copy RP 5th-percentile threshold: %.2f Angstrom",
                        dist_thr_same))

        # Cross-copy atoms: RP CA + rRNA N1 of the other copy
        rRNA_data_copy_col <- rRNA_data |>
            dplyr::mutate(copy = copy_map[RNA_chain]) |>
            dplyr::filter(!is.na(copy)) |>
            dplyr::select(x, y, z, copy)
        cross_atoms <- dplyr::bind_rows(
            prot_CA |> dplyr::select(x, y, z, copy),
            rRNA_data_copy_col
        )

        # For each rRNA N1 atom, min distance to cross-copy atoms
        FinalData_cross <- dplyr::bind_rows(lapply(c(1L, 2L), function(cp) {
            rna_cp   <- rRNA_data[!is.na(rRNA_data$copy) & rRNA_data$copy == cp, ]
            other_cp <- cross_atoms[cross_atoms$copy != cp, ]
            if (nrow(rna_cp) == 0L || nrow(other_cp) == 0L) return(NULL)
            message(sprintf("  Copy %d: %d rRNA atoms vs %d cross-copy atoms",
                            cp, nrow(rna_cp), nrow(other_cp)))
            data.frame(
                RP           = "Col.Int.",
                RNA_molecule = rna_cp$RP,
                RNA_residue  = rna_cp$resno,
                distance     = .pairwise_euclidean_rowmin(
                    rna_cp[, c("x", "y", "z")],
                    other_cp[, c("x", "y", "z")]
                ),
                stringsAsFactors = FALSE
            )
        })) |>
            dplyr::group_by(RNA_molecule, RNA_residue) |>
            dplyr::summarise(distance = duplication_resolve(distance), .groups = "drop") |>
            dplyr::mutate(RP = "Col.Int.") |>
            dplyr::arrange(RNA_molecule, RNA_residue)

        long_csv <- paste0(out_prefix, ".", PDBid, ".Col.ARF.minimum_distances.csv")
        readr::write_csv(FinalData_cross, long_csv)

        RP_proximity_df <- reshape2::dcast(
            FinalData_cross,
            RNA_molecule + RNA_residue ~ RP,
            value.var     = "distance",
            fun.aggregate = duplication_resolve
        )
        colnames(RP_proximity_df)[1:2] <- c("rRNA", "resno")
        RP_proximity_df[RP_proximity_df == -Inf] <- NA_real_
        RP_proximity_df[RP_proximity_df ==  Inf] <- NA_real_
        rownames(RP_proximity_df) <- paste(
            RP_proximity_df$rRNA, RP_proximity_df$resno, sep = "_"
        )
        wide_csv <- paste0(out_prefix, ".", PDBid, ".Col.ARF.minDist_matrix.csv")
        readr::write_csv(RP_proximity_df, wide_csv)
        return(list(RP_proximity_df = RP_proximity_df, output_file = long_csv))
    }
    # ---- End disome path -------------------------------------------------------

    message(sprintf(
        "Computing distance matrix: %d RNA residues × %d protein atoms ...",
        nrow(rRNA_data), nrow(alphaCarbons)
    ))

    # -- Vectorised pairwise distance computation -
    dist_mat    <- .pairwise_euclidean(
        rRNA_data[, c("x", "y", "z")],
        alphaCarbons[, c("x", "y", "z")]
    )
    unique_RPs  <- unique(alphaCarbons$RP)

    FinalData_clean <- dplyr::bind_rows(lapply(unique_RPs, function(rp) {
        prot_idx <- which(alphaCarbons$RP == rp)
        min_dist <- apply(dist_mat[, prot_idx, drop = FALSE], 1L, min)
        data.frame(
            RP           = rp,
            RNA_molecule = rRNA_data$RP,
            RNA_residue  = rRNA_data$resno,
            distance     = min_dist,
            stringsAsFactors = FALSE
        )
    })) |>
        dplyr::filter(RP != "", RNA_molecule != "") |>
        dplyr::group_by(RP, RNA_molecule, RNA_residue) |>
        dplyr::summarise(distance = duplication_resolve(distance), .groups = "drop") |>
        dplyr::arrange(RP, RNA_molecule, RNA_residue)

    # -- Write long-format CSV 
    long_csv <- paste0(out_prefix, ".", PDBid, ".Col.ARF.minimum_distances.csv")
    readr::write_csv(FinalData_clean, long_csv)

    # -- Pivot to wide matrix (rows = residues, columns = RPs)
    RP_proximity_df <- reshape2::dcast(
        FinalData_clean,
        RNA_molecule + RNA_residue ~ RP,
        value.var    = "distance",
        fun.aggregate = duplication_resolve
    )
    colnames(RP_proximity_df)[1:2] <- c("rRNA", "resno")

    # Replace ±Inf introduced by fun.aggregate on empty groups with NA
    RP_proximity_df[RP_proximity_df == -Inf] <- NA_real_
    RP_proximity_df[RP_proximity_df ==  Inf] <- NA_real_

    rownames(RP_proximity_df) <- paste(
        RP_proximity_df$rRNA, RP_proximity_df$resno, sep = "_"
    )

    wide_csv <- paste0(out_prefix, ".", PDBid, ".Col.ARF.minDist_matrix.csv")
    readr::write_csv(RP_proximity_df, wide_csv)

    list(
        RP_proximity_df = RP_proximity_df,
        output_file     = long_csv
    )
}


# 4. dripARF_get_collision_proximity_sets — build RP proximity GSEA sets

#' Build RP-proximity GSEA sets with circular-shift randomisations
#'
#' For each ribosomal protein column in \code{RP_proximity_df}, identifies the
#' 5 % closest RNA residue positions (or a custom threshold) and generates 99
#' circular-shift randomisation sets for background estimation in GSEA.
#'
#' @param RP_proximity_df  Wide data.frame: columns \code{rRNA} and
#'   \code{resno} followed by one column per RP (numeric distances).
#' @param additional_RPcols Integer vector of additional RP column indices to
#'   optionally exclude from automatic capping.
#' @param rRNAs_fasta     Optional path to a target-species FASTA; restricts
#'   residues to those in the FASTA and sets the size threshold proportionally.
#' @param thresholds      Numeric vector \code{c(dist_threshold, set_size_cap)},
#'   or \code{NULL} to compute automatically (5th-percentile distance, 5 %
#'   of total positions).
#' @param cap_added_RPcols Logical. Apply size cap to \code{additional_RPcols}
#'   as well.
#'
#' @return data.frame with columns \code{ont} (RP label or "Rand<i>_<RP>")
#'   and \code{gene} (\code{"<rRNA>_<resno>"}).
#'
#' @noRd
dripARF_get_collision_proximity_sets <- function(
    RP_proximity_df,
    additional_RPcols  = c(),
    rRNAs_fasta        = NULL,
    thresholds         = NULL,
    cap_added_RPcols   = FALSE,
    extend_nt          = 0L
) {
    message("Generating proximity sets for GSEA analysis")
    message(sprintf("Positions in RP_proximity_df: %d", nrow(RP_proximity_df)))

    RPfocus <- colnames(RP_proximity_df)[c(-1L, -2L)]

    # Excluded column indices (rRNA, resno, plus any additional RP cols).
    # Computed unconditionally so it is available inside the for-loop even
    # when thresholds are supplied by the caller.
    excluded_idx <- c(-1L, -2L)
    if (length(additional_RPcols) > 0L) {
        excluded_idx <- c(excluded_idx, -as.integer(additional_RPcols))
    }

    # -- Auto-compute thresholds if not provided
    if (is.null(thresholds)) {
        dist_thr <- stats::quantile(
            unlist(RP_proximity_df[, excluded_idx]),
            probs   = 0.05,
            na.rm   = TRUE
        )
        RPset_size_thr <- round(nrow(RP_proximity_df) * 0.05)

        if (!is.null(rRNAs_fasta)) {
            seq_set        <- Biostrings::readBStringSet(rRNAs_fasta)
            RPset_size_thr <- round(sum(Biostrings::width(seq_set)) * 0.05)
            rRNA_names     <- sapply(
                strsplit(names(seq_set), split = " "), `[`, 1L
            )
            RP_proximity_df <- RP_proximity_df[
                RP_proximity_df$rRNA %in% rRNA_names, , drop = FALSE
            ]
        }
        thresholds <- c(dist_thr, RPset_size_thr)
    }

    dist_thr        <- thresholds[1L]
    set_size_cap    <- thresholds[2L]
    standard_RPs    <- colnames(RP_proximity_df)[excluded_idx]

    # -- Collect proximity positions per RP
    RP_proxpos <- list()
    for (RP in RPfocus) {
        proxpos <- which(RP_proximity_df[[RP]] <= dist_thr)
        is_standard <- RP %in% standard_RPs

        if (length(proxpos) > set_size_cap && (is_standard || cap_added_RPcols)) {
            proxpos <- sort(
                proxpos[order(RP_proximity_df[[RP]][proxpos])[seq_len(set_size_cap)]]
            )
        }
        if (length(proxpos) > 0L) {
            if (extend_nt > 0L) {
                cur_df       <- RP_proximity_df[proxpos, c("rRNA", "resno"), drop = FALSE]
                cur_df$resno <- as.integer(cur_df$resno)
                ext_df       <- dplyr::bind_rows(lapply(seq(-extend_nt, extend_nt), function(d) {
                    data.frame(rRNA = cur_df$rRNA, resno = cur_df$resno + d)
                }))
                ext_key <- paste(ext_df$rRNA, ext_df$resno, sep = "_")
                all_key <- paste(RP_proximity_df$rRNA,
                                 as.integer(RP_proximity_df$resno), sep = "_")
                proxpos <- sort(which(all_key %in% ext_key))
            }
            RP_proxpos[[RP]] <- proxpos
        }
        message(sprintf("  %-30s %d proximity positions", RP,
                        length(RP_proxpos[[RP]])))
    }

    # -- Build GSEA sets with randomisations (moved outside the loop)
    n          <- nrow(RP_proximity_df)
    row_labels <- paste(RP_proximity_df$rRNA, RP_proximity_df$resno, sep = "_")

    gsea_sets_RP <- dplyr::bind_rows(lapply(names(RP_proxpos), function(RP) {
        proxpos <- RP_proxpos[[RP]]
        foreground <- data.frame(
            ont  = RP,
            gene = row_labels[proxpos],
            stringsAsFactors = FALSE
        )
        background <- .make_rand_sets(proxpos, n, RP, row_labels)
        rbind(foreground, background)
    }))

    gsea_sets_RP
}


# 5. ARF_build_ColInt_set — full collision-interface set pipeline

#' Build the collision-interface GSEA gene set
#'
#' Orchestrates \code{\link{ARF_build_ColInt_proximity}},
#' \code{\link[ARF]{ARF_convert_Ribo3D_pos}}, and
#' \code{\link{dripARF_get_collision_proximity_sets}} to produce the final
#' collision-interface GSEA gene set for a target species.  Selected positions
#' (5th-percentile RP distance threshold) are extended by ±1 nt in rRNA
#' coordinate space before building GSEA sets.
#'
#' @param rRNAs_file          Path to the source (structure) rRNA FASTA.
#' @param PDB_file            Path to the mmCIF structure file.
#' @param conversion_table    Chain-to-RP mapping data.frame.
#' @param PDBid               Four-letter PDB identifier.
#' @param out_prefix          Prefix for intermediate CSV files.
#' @param rRNAs               Character vector of rRNA names to retain.
#' @param duplication_resolve Function for resolving duplicate distances.
#'   Default: \code{min}.
#' @param target_species      Character. Target species name / prefix.
#' @param target_rRNAs_fasta  Path to the target-species rRNA FASTA.
#' @param rRNA_pairs          List of two-element character vectors
#'   \code{c(source_rRNA, target_rRNA)} for liftover alignment.
#' @param n_copies            Integer. \code{2L} (default) treats the structure
#'   as a disome and computes cross-ribosome distances; chains are assigned to
#'   copies automatically via k-means clustering.  Set to \code{1L} to use the
#'   monosome RP-proximity method instead.
#'
#' @return data.frame with columns \code{ont} and \code{gene}.
#'
#' @export
ARF_build_ColInt_set <- function(
    rRNAs_file,
    PDB_file,
    conversion_table,
    PDBid,
    out_prefix,
    rRNAs,
    duplication_resolve  = min,
    target_species,
    target_rRNAs_fasta,
    rRNA_pairs = NULL,
    n_copies   = 2L
) {
    # -- Step A: compute RP–rRNA distances (or cross-copy distances for n_copies = 2)
    ColInt_proximity_list <- ARF_build_ColInt_proximity(
        rRNAs_file          = rRNAs_file,
        PDB_file            = PDB_file,
        conversion_table    = conversion_table,
        PDBid               = PDBid,
        out_prefix          = out_prefix,
        target_species      = target_species,
        duplication_resolve = duplication_resolve,
        n_copies            = n_copies
    )

    # -- Step B: lift over to target species rRNA coordinates
    message("Lifting 3D positions to ", target_species, " rRNA coordinates ...")

    LO.Col_proximity_df <- ARF::ARF_convert_Ribo3D_pos(
        source_distance_file = ColInt_proximity_list$output_file,
        source_rRNAs_fasta   = rRNAs_file,
        target_species       = target_species,
        target_rRNAs_fasta   = target_rRNAs_fasta,
        rRNA_pairs           = rRNA_pairs,
        source_positions     = NULL,
        source_sets          = NULL,
        type                 = "distances"
    )

    # -- Step C: reshape to a single "Col.Int." column
    # n_copies = 2: liftover produced a single "Col.Int." column; rename directly.
    # n_copies = 1: pivot multiple RP columns into one "Col.Int." column.
    if (n_copies == 2L) {
        LO.Col_proximity_temp_df <- LO.Col_proximity_df[,
            c("rRNA", "resno", "Col.Int."), drop = FALSE] |>
            dplyr::rename(distance = `Col.Int.`) |>
            dplyr::mutate(`Col.Int.` = "Col.Int.") |>
            dplyr::arrange(distance) |>
            as.data.frame()
    } else {
        LO.Col_proximity_temp_df <- LO.Col_proximity_df[, c("rRNA", "resno", rRNAs)] |>
            tidyr::pivot_longer(
                cols      = -c(rRNA, resno),
                names_to  = "Col.Int.",
                values_to = "distance"
            ) |>
            dplyr::arrange(distance) |>
            as.data.frame()
    }

    n_unique       <- nrow(dplyr::distinct(LO.Col_proximity_temp_df, rRNA, resno))
    max_set_size   <- round(n_unique * 0.1)
    dist_thr       <- stats::quantile(LO.Col_proximity_temp_df$distance, 0.05, na.rm = TRUE)
    max_dist_value <- max(LO.Col_proximity_temp_df$distance, na.rm = TRUE)

    LO.Col_proximity_df_final <- LO.Col_proximity_temp_df |>
        dplyr::filter(distance <= dist_thr, !is.na(distance)) |>
        dplyr::group_by(`Col.Int.`) |>
        dplyr::slice_head(n = max_set_size) |>
        dplyr::ungroup() |>
        dplyr::mutate(`Col.Int.` = "Col.Int.") |>
        tidyr::pivot_wider(
            names_from   = `Col.Int.`,
            values_from  = distance,
            values_fill  = max_dist_value
        ) |>
        as.data.frame()

    rownames(LO.Col_proximity_df_final) <- paste(
        LO.Col_proximity_df_final$rRNA,
        LO.Col_proximity_df_final$resno,
        sep = "_"
    )

    # -- Step D: generate GSEA sets (±1 nt extension applied to match original method)
    dripARF_get_collision_proximity_sets(
        RP_proximity_df   = LO.Col_proximity_df_final,
        additional_RPcols = c(),
        rRNAs_fasta       = target_rRNAs_fasta,
        thresholds        = NULL,
        cap_added_RPcols  = FALSE,
        extend_nt         = 1L
    )
}


# 6a. ARF_parse_structure_fasta — extract conversion table and rRNA sequences

#' Parse a PDB FASTA into a conversion table and rRNA sequences
#'
#' Reads a RCSB-format FASTA file and extracts:
#' (a) a conversion table mapping auth chain IDs to RP names, and
#' (b) a \code{BStringSet} of rRNA-only sequences prefixed by
#' \code{target_species}.
#'
#' The function handles RCSB headers of the form
#' \code{>PDBid_chain|Chain(s) auth A|protein L1|Organism}.
#' For non-standard headers, adjust the \code{separate()} call and regex
#' patterns in the function body.
#'
#' @param fasta_file         Character. Path to the PDB FASTA file (as
#'   downloaded by \code{ARF_download_structure_files}).
#' @param target_species     Character. Short species prefix, e.g. \code{"ec"}.
#' @param rRNAs_output_file  Character or \code{NULL}.  If non-\code{NULL},
#'   the rRNA sequences are written to this path as a FASTA file.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{conversion_table}}{data.frame with columns \code{names},
#'       \code{chainID} (auth chain, comma-separated for multi-chain entries),
#'       \code{RP_name}, \code{organism}, and \code{RP_new} (cleaned name,
#'       underscores for spaces).}
#'     \item{\code{rRNA_seqs}}{BStringSet of rRNA sequences with names in the
#'       form \code{<target_species>_<RP_new>}, ready for use in downstream
#'       liftover functions.}
#'   }
#'
#' @export
ARF_parse_structure_fasta <- function(
    fasta_file,
    target_species,
    rRNAs_output_file = NULL
) {
    if (!file.exists(fasta_file))
        stop("fasta_file does not exist: ", fasta_file, call. = FALSE)

    fasta_seqs <- Biostrings::readBStringSet(fasta_file)

    # Build conversion table (chain → RP name)
    conversion_table <- data.frame(names = names(fasta_seqs),
                                   stringsAsFactors = FALSE) |>
        tidyr::separate(
            names,
            into    = c("ID", "chainID", "RP_name", "organism"),
            sep     = "\\|",
            remove  = FALSE,
            extra   = "merge",
            fill    = "right"
        ) |>
        dplyr::mutate(
            RP_name = sub(".*protein", "", RP_name),
            # New RCSB format: "Chains J[auth L7], L[auth L9]" → extract all [auth <id>]
            # Old format: "Chain(s) auth AL, AZ"               → strip prefix
            chainID = ifelse(
                grepl("\\[auth ", chainID),
                vapply(chainID, function(x) {
                    m <- regmatches(x, gregexpr("(?<=\\[auth )[^\\]]+", x, perl = TRUE))[[1L]]
                    paste(m, collapse = ",")
                }, character(1L)),
                sub("Chain(s)? .*?auth ", "", sub("\\].*", "", chainID))
            )
        ) |>
        dplyr::mutate(
            RP_new = gsub("^ | $", "", RP_name),
            RP_new = gsub(" ", "_", RP_new)
        )

    # Extract rRNA entries (sequence headers containing "[0-9]S" and "RNA")
    rRNA_seqs <- fasta_seqs[intersect(
        grep("[0-9]S", names(fasta_seqs), value = TRUE),
        grep("RNA",    names(fasta_seqs), value = TRUE)
    )]

    # Rename to target_species prefix
    names(rRNA_seqs) <- paste0(
        target_species, "_",
        conversion_table$RP_new[conversion_table$names %in% names(rRNA_seqs)]
    )

    if (!is.null(rRNAs_output_file)) {
        Biostrings::writeXStringSet(rRNA_seqs, rRNAs_output_file)
        message("rRNA FASTA written to: ", rRNAs_output_file)
    }

    list(conversion_table = conversion_table, rRNA_seqs = rRNA_seqs)
}


# 6b-helper. .build_joint_sets — joint foreground + fresh circular-shift Rands

.build_joint_sets <- function(combined_df, target_rRNAs_fasta) {
    # combined_df rows have ont values such as:
    #   "ec_9gft_disome_Rib.Col."   "Rand1_ec_9gft_disome_Rib.Col."
    #   "ec_9ggr_disome_Col.Int."   "Rand1_ec_9ggr_disome_Col.Int."
    #
    # For each set type (Rib.Col., Col.Int.) the joint foreground is the union
    # of positions across all named structures.  Fresh circular-shift Rand sets
    # are generated from that union over the full rRNA position space, so the
    # null distribution is statistically correct rather than merging individual
    # Rand sets which would have the wrong size and shift structure.
    types    <- c("Rib.Col.", "Col.Int.")
    type_pat <- paste0("(", paste(gsub("\\.", "\\\\.", types), collapse = "|"), ")$")

    is_typed <- grepl(type_pat, combined_df$ont)
    sub_df   <- combined_df[is_typed, , drop = FALSE]
    if (nrow(sub_df) == 0L)
        return(data.frame(ont = character(), gene = character(),
                          stringsAsFactors = FALSE))

    sub_df$is_fg    <- !grepl("^Rand[0-9]+_", sub_df$ont)
    sub_df$set_type <- sub(paste0("^.*_(", paste(gsub("\\.", "\\\\.", types),
                                                  collapse = "|"), ")$"),
                           "\\1", sub_df$ont, perl = TRUE)

    # Full rRNA position space from the target FASTA (same space used per-structure)
    target_seqs <- Biostrings::readBStringSet(target_rRNAs_fasta)
    row_labels  <- unlist(lapply(names(target_seqs), function(nm) {
        paste0(nm, "_", seq_len(Biostrings::width(target_seqs[nm])))
    }), use.names = FALSE)
    n <- length(row_labels)

    dplyr::bind_rows(lapply(types, function(type) {
        fg_genes <- unique(sub_df$gene[sub_df$is_fg & sub_df$set_type == type])
        if (length(fg_genes) == 0L) return(NULL)
        proxpos <- which(row_labels %in% fg_genes)
        if (length(proxpos) == 0L) return(NULL)
        foreground <- data.frame(ont  = type,
                                 gene = row_labels[proxpos],
                                 stringsAsFactors = FALSE)
        background <- .make_rand_sets(proxpos, n, type, row_labels)
        rbind(foreground, background)
    }))
}

# 6b. ARF_build_collision_sets — one-call pipeline wrapper

#' Build collision-interface GSEA gene sets from a ribosome structure
#'
#' A single-function entry point that runs the entire collision-interface
#' pipeline: resolves / downloads the structure and FASTA, computes SASA,
#' parses the FASTA, builds SASA candidates, generates \strong{Rib.Col.}
#' gene sets (SASA burial), and generates \strong{Col.Int.} gene sets (3D
#' proximity).  Returns a combined \code{data.frame} ready to pass to
#' \code{\link{dricARF}} as \code{gsea_sets_Collision}.
#'
#' Use this function for a quick end-to-end run.  For exploratory work or
#' custom intermediate inspection, use the individual functions
#' (\code{\link{ARF_run_SASA}}, \code{\link{ARF_parse_structure_fasta}},
#' \code{\link{ARF_build_SASA_candidates}}, \code{\link{dricARF_build_SASA_set}},
#' \code{\link{ARF_build_ColInt_set}}) directly.
#'
#' @param PDBid              Character. Four-letter PDB identifier.
#' @param target_species     Character. Short species prefix (e.g. \code{"ec"}).
#' @param target_rRNAs_fasta Character. Path to the target-species reference
#'   rRNA FASTA used for liftover.
#' @param rRNA_pairs         List of two-element character vectors
#'   \code{c(pdb_rRNA_name, reference_rRNA_name)}.  The first element must
#'   match (after \code{target_species}_-prefixing) a name in the parsed
#'   rRNA FASTA; the second must match a header in \code{target_rRNAs_fasta}.
#' @param PDB_file           Character or \code{NULL}.  Local structure file;
#'   \code{NULL} downloads from RCSB.
#' @param fasta_file         Character or \code{NULL}.  Local PDB FASTA;
#'   \code{NULL} downloads from RCSB.
#' @param output_dir         Character. Directory for all outputs.
#'   Default: \code{getwd()}.
#' @param set_name           Character or \code{NULL}. When supplied, all
#'   \code{ont} values in the returned sets are prefixed with this string
#'   (e.g. \code{"ec_9gft_disome"} turns \code{"Rib.Col."} into
#'   \code{"ec_9gft_disome_Rib.Col."}). Useful when combining sets from
#'   multiple structures with \code{dplyr::bind_rows()}.
#' @param overwrite          Logical. Overwrite existing intermediate files.
#' @param algorithm          SASA algorithm: \code{"LeeRichards"} (default) or
#'   \code{"ShrakeRupley"}.
#' @param probe_radius       Probe radius in Å. Default \code{10}.
#' @param n_slices,n_points,n_threads FreeSASA tuning parameters.
#' @param model              \code{NULL} or a model identifier.
#' @param polymer_type       \code{"RNA"}, \code{"DNA"}, or \code{"ALL"}.
#' @param include_hydrogen   Logical.
#' @param zero_tolerance     Numeric.
#' @param precision          Integer. Decimal places in TSV output.
#' @param duplication_resolve Function for resolving duplicate distances.
#'   Default: \code{min}.
#' @param n_copies           Integer. Passed to \code{\link{ARF_build_ColInt_set}}.
#'   \code{2L} (default) treats the structure as a disome and computes true
#'   cross-ribosome contact distances; use this when \code{PDBid} refers to a
#'   disome structure.  Set to \code{1L} to use the monosome RP-proximity method.
#' @param envname            Python environment name.
#' @param env_type           \code{"virtualenv"} or \code{"conda"}.
#' @param python_version     Python version string for environment creation.
#' @param force              Logical. Force-reinstall the Python environment.
#'
#' @return A named list (returned \code{invisible}):
#'   \describe{
#'     \item{\code{gsea_sets_Collision}}{Combined data.frame with columns
#'       \code{ont} and \code{gene}; pass directly to \code{dricARF()} as
#'       \code{gsea_sets_Collision}.}
#'     \item{\code{RibCol.gsea_sets}}{data.frame for \code{Rib.Col.} sets.}
#'     \item{\code{Col.gsea_set}}{data.frame for \code{Col.Int.} sets.}
#'     \item{\code{sasa_df}}{data.frame of SASA candidates.}
#'     \item{\code{conversion_table}}{Chain-to-RP conversion table.}
#'     \item{\code{rRNA_seqs}}{BStringSet of rRNA sequences.}
#'   }
#'
#' @seealso \code{\link{ARF_run_SASA}}, \code{\link{ARF_parse_structure_fasta}},
#'   \code{\link{ARF_build_SASA_candidates}}, \code{\link{dricARF_build_SASA_set}},
#'   \code{\link{ARF_build_ColInt_set}}, \code{\link{dricARF}}
#'
#' @export
ARF_build_collision_sets <- function(
    PDBid,
    target_species,
    target_rRNAs_fasta,
    rRNA_pairs           = NULL,
    PDB_file             = NULL,
    fasta_file           = NULL,
    output_dir           = getwd(),
    set_name             = NULL,
    overwrite            = FALSE,
    algorithm            = "LeeRichards",
    probe_radius         = 10,
    n_slices             = 20L,
    n_points             = 1000L,
    n_threads            = 1L,
    model                = NULL,
    polymer_type         = "ALL",
    include_hydrogen     = FALSE,
    zero_tolerance       = 1e-8,
    precision            = 4L,
    duplication_resolve  = min,
    n_copies             = 2L,
    envname              = getOption("ARF.python.envname", "sasa-python"),
    env_type             = getOption("ARF.python.env_type", "virtualenv"),
    python_version       = "3.11",
    force                = FALSE
) {
    # --- Multi-structure dispatch -------------------------------------------------
    if (length(PDBid) > 1L) {
        n <- length(PDBid)
        if (length(output_dir) == 1L)
            output_dir <- file.path(output_dir, PDBid)
        if (length(output_dir) != n)
            stop("`output_dir` must be length 1 or the same length as `PDBid`.")
        if (!is.null(set_name) && length(set_name) != n)
            stop("`set_name` must be NULL or the same length as `PDBid`.")
        target_species <- rep_len(target_species, n)

        results <- setNames(vector("list", n), PDBid)
        for (i in seq_len(n)) {
            results[[i]] <- ARF_build_collision_sets(
                PDBid               = PDBid[i],
                target_species      = target_species[i],
                target_rRNAs_fasta  = target_rRNAs_fasta,
                rRNA_pairs          = rRNA_pairs,
                PDB_file            = if (is.null(PDB_file))   NULL else PDB_file[i],
                fasta_file          = if (is.null(fasta_file)) NULL else fasta_file[i],
                output_dir          = output_dir[i],
                set_name            = if (is.null(set_name)) NULL else set_name[i],
                overwrite           = overwrite,
                algorithm           = algorithm,
                probe_radius        = probe_radius,
                n_slices            = n_slices,
                n_points            = n_points,
                n_threads           = n_threads,
                model               = model,
                polymer_type        = polymer_type,
                include_hydrogen    = include_hydrogen,
                zero_tolerance      = zero_tolerance,
                precision           = precision,
                duplication_resolve = duplication_resolve,
                n_copies            = n_copies,
                envname             = envname,
                env_type            = env_type,
                python_version      = python_version,
                force               = force
            )
        }
        individual <- dplyr::bind_rows(lapply(results, `[[`, "gsea_sets_Collision"))
        # When set_name is supplied, add joint Rib.Col. / Col.Int. sets (union
        # across structures per Rand level) to match the built-in preset pattern.
        joint <- if (!is.null(set_name))
            .build_joint_sets(individual, target_rRNAs_fasta) else NULL
        combined <- dplyr::bind_rows(individual, joint)
        return(invisible(list(gsea_sets_Collision = combined, structures = results)))
    }

    # --- Single-structure pipeline -----------------------------------------------
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    sasa_output_file <- file.path(output_dir, paste0(PDBid, "_sasa.tsv"))

    # --- Step 1: Download structure and FASTA + compute SASA
    # When overwrite = FALSE and the SASA TSV already exists, reuse it so the
    # expensive Python step is skipped on re-runs.
    if (!overwrite && file.exists(sasa_output_file)) {
        message("Reusing existing SASA output: ", sasa_output_file)
        if (is.null(PDB_file))   PDB_file   <- file.path(output_dir, paste0(PDBid, ".cif"))
        if (is.null(fasta_file)) fasta_file <- file.path(output_dir, paste0(PDBid, ".fasta"))
        annotations <- read_structure_annotations(sasa_output_file)
        sasa_result <- list(annotations = annotations,
                            PDB_file    = PDB_file,
                            fasta_file  = fasta_file)
    } else {
        sasa_result <- ARF_run_SASA(
            PDBid              = PDBid,
            PDB_file           = PDB_file,
            rRNAs_file         = fasta_file,
            download_directory = output_dir,
            output_file        = sasa_output_file,
            overwrite          = overwrite,
            algorithm          = algorithm,
            probe_radius       = probe_radius,
            n_slices           = n_slices,
            n_points           = n_points,
            n_threads          = n_threads,
            model              = model,
            polymer_type       = polymer_type,
            include_hydrogen   = include_hydrogen,
            zero_tolerance     = zero_tolerance,
            precision          = precision,
            envname            = envname,
            env_type           = env_type,
            python_version     = python_version,
            force              = force
        )
    }
    PDB_file   <- sasa_result$PDB_file
    fasta_file <- sasa_result$fasta_file

    # --- Step 2: Parse PDB FASTA → conversion_table + rRNA_seqs
    rRNAs_file_clean <- file.path(output_dir, paste0(PDBid, "_rRNA.fasta"))
    parsed <- ARF_parse_structure_fasta(
        fasta_file        = fasta_file,
        target_species    = target_species,
        rRNAs_output_file = if (!overwrite && file.exists(rRNAs_file_clean)) NULL
                            else rRNAs_file_clean
    )
    conversion_table <- parsed$conversion_table
    rRNA_seqs        <- parsed$rRNA_seqs

    # Auto-detect rRNA pairs once the structure FASTA is on disk
    if (is.null(rRNA_pairs))
        rRNA_pairs <- auto_pair_rRNAs(rRNAs_file_clean, target_rRNAs_fasta)

    # --- Step 3: Clean SASA annotations and write input TSV for candidates
    input_sasa_file <- file.path(output_dir, "input_sasa.tsv")
    if (overwrite || !file.exists(input_sasa_file)) {
        pdb_sasa_df <- sasa_result$annotations |>
            dplyr::mutate(
                entity_description = dplyr::case_when(
                    grepl("protein", entity_description) ~
                        sub(".*protein", "", entity_description),
                    TRUE ~ entity_description
                ),
                entity_description = gsub("^'|'$|^ | $", "", entity_description),
                entity_description = gsub(" ", "_", entity_description)
            )
        utils::write.table(pdb_sasa_df, input_sasa_file,
                           row.names = FALSE, quote = FALSE, sep = "\t")
    }

    # --- Step 4: Build SASA candidates
    rRNAs   <- sub(paste0("^", target_species, "_"), "", names(rRNA_seqs))
    sasa_df <- ARF_build_SASA_candidates(
        output_file    = input_sasa_file,
        rRNAs          = rRNAs,
        target_species = target_species
    )

    # --- Step 5: Rib.Col. gene sets (SASA burial)
    RibCol.gsea_sets <- dricARF_build_SASA_set(
        sasa_df            = sasa_df,
        rRNAs_file         = rRNAs_file_clean,
        target_rRNAs_fasta = target_rRNAs_fasta,
        rRNA_pairs         = rRNA_pairs
    )

    # --- Step 6: Col.Int. gene sets (3D proximity)
    Col.gsea_set <- ARF_build_ColInt_set(
        rRNAs_file          = rRNAs_file_clean,
        PDB_file            = PDB_file,
        conversion_table    = conversion_table,
        PDBid               = sub("\\.[^.]+$", "", basename(PDB_file)),
        out_prefix          = file.path(output_dir, paste0(PDBid, ".Ribosome.3D")),
        rRNAs               = paste(target_species, rRNAs, sep = "_"),
        duplication_resolve = duplication_resolve,
        target_species      = target_species,
        target_rRNAs_fasta  = target_rRNAs_fasta,
        rRNA_pairs          = rRNA_pairs,
        n_copies            = n_copies
    )

    gsea_sets_Collision <- dplyr::bind_rows(RibCol.gsea_sets, Col.gsea_set)

    if (!is.null(set_name)) {
        gsea_sets_Collision <- gsea_sets_Collision |>
            dplyr::mutate(ont = sub("^(Rand[0-9]+_)?(.+)$",
                                    paste0("\\1", set_name, "_\\2"),
                                    ont))
        RibCol.gsea_sets <- RibCol.gsea_sets |>
            dplyr::mutate(ont = sub("^(Rand[0-9]+_)?(.+)$",
                                    paste0("\\1", set_name, "_\\2"),
                                    ont))
        Col.gsea_set <- Col.gsea_set |>
            dplyr::mutate(ont = sub("^(Rand[0-9]+_)?(.+)$",
                                    paste0("\\1", set_name, "_\\2"),
                                    ont))
    }

    invisible(list(
        gsea_sets_Collision = gsea_sets_Collision,
        RibCol.gsea_sets    = RibCol.gsea_sets,
        Col.gsea_set        = Col.gsea_set,
        sasa_df             = sasa_df,
        conversion_table    = conversion_table,
        rRNA_seqs           = rRNA_seqs
    ))
}


# 7. dricARF_build_SASA_set — liftover + GSEA sets from SASA candidates

#' Build SASA-based GSEA sets via rRNA liftover alignment
#'
#' Aligns source and target rRNA sequences with ClustalW to lift over the SASA
#' candidate positions (z-score of raw differential SASA > 1.64), extends each
#' selected position by ±1 nt in rRNA coordinate space, then builds foreground
#' and 99 circular-shift randomisation GSEA sets.
#'
#' @param sasa_df           data.frame produced by
#'   \code{\link{ARF_build_SASA_candidates}}.
#' @param rRNAs_file        Path to the source rRNA FASTA.
#' @param target_rRNAs_fasta Path to the target-species rRNA FASTA.
#' @param rRNA_pairs        List of two-element character vectors
#'   \code{c(source_rRNA, target_rRNA)} for alignment.
#'
#' @return data.frame with columns \code{ont} and \code{gene}.
#'
#' @export
dricARF_build_SASA_set <- function(
    sasa_df,
    rRNAs_file,
    target_rRNAs_fasta,
    rRNA_pairs = NULL
) {
    # -- Load and normalise sequences
    structure_rRNAs <- Biostrings::readBStringSet(rRNAs_file)
    structure_rRNAs <- Biostrings::BStringSet(
        chartr("TtXx", "UuNn", as.character(structure_rRNAs))
    )
    target_rRNAs <- Biostrings::readBStringSet(target_rRNAs_fasta)
    target_rRNAs <- Biostrings::BStringSet(
        chartr("TtXx", "UuNn", as.character(target_rRNAs))
    )

    # -- Auto-detect rRNA pairs when not supplied
    if (is.null(rRNA_pairs))
        rRNA_pairs <- auto_pair_rRNAs(rRNAs_file, target_rRNAs_fasta)

    # -- Build source → target position maps via pairwise alignment
    message("Aligning rRNA pairs for position liftover ...")
    structure_2_t <- lapply(rRNA_pairs, function(pair) {
        src_name <- pair[1L]
        tgt_name <- pair[2L]
        src_seq  <- structure_rRNAs[src_name]
        aln      <- msa::msaClustalW(
            Biostrings::RNAStringSet(append(src_seq, target_rRNAs[tgt_name]))
        )@unmasked
        aln_chars <- strsplit(as.character(aln), split = "")

        n       <- length(aln_chars[[1L]])
        map     <- integer(length(aln_chars[[1L]]))  # pre-allocate
        spos    <- 0L
        tpos    <- 0L
        out     <- rep(NA_integer_, n)               # source length, target positions
        for (i in seq_len(n)) {
            if (aln_chars[[2L]][i] != "-") tpos <- tpos + 1L
            if (aln_chars[[1L]][i] != "-") {
                spos <- spos + 1L
                if (aln_chars[[2L]][i] != "-") out[spos] <- tpos
            }
        }
        out
    })
    names(structure_2_t) <- vapply(rRNA_pairs, `[`, character(1L), 1L)

    # -- Map SASA positions to target species residues 
    RP_proxpos <- dplyr::bind_rows(lapply(rRNA_pairs, function(pair) {
        src_name <- pair[1L]
        tgt_name <- pair[2L]
        pos_map  <- structure_2_t[[src_name]]

        source_positions <- sasa_df |>
            dplyr::filter(grepl(src_name, rRNA, fixed = TRUE)) |>
            dplyr::select(auth_seq_id) |>
            dplyr::distinct() |>
            dplyr::pull() |>
            as.integer()

        data.frame(
            rRNA      = tgt_name,
            resno     = pos_map,
            focus_pos = pos_map %in% pos_map[source_positions],
            stringsAsFactors = FALSE
        )
    })) |>
        dplyr::filter(!is.na(resno))

    # -- Build GSEA sets
    n          <- nrow(RP_proxpos)
    row_labels <- paste(RP_proxpos$rRNA, RP_proxpos$resno, sep = "_")
    proxpos    <- which(RP_proxpos$focus_pos)

    if (length(proxpos) == 0L) {
        warning("No SASA foreground positions found after liftover. ",
                "Returning empty set.", call. = FALSE)
        return(data.frame(ont = character(0L), gene = character(0L)))
    }

    # ±1 nt extension in rRNA coordinate space
    focus_df       <- RP_proxpos[proxpos, , drop = FALSE]
    focus_df$resno <- as.integer(focus_df$resno)
    ext_df  <- dplyr::bind_rows(lapply(-1L:1L, function(d) {
        data.frame(rRNA = focus_df$rRNA, resno = focus_df$resno + d)
    }))
    ext_key <- paste(ext_df$rRNA, ext_df$resno, sep = "_")
    all_key <- paste(RP_proxpos$rRNA, as.integer(RP_proxpos$resno), sep = "_")
    proxpos <- which(all_key %in% ext_key)

    foreground <- data.frame(
        ont  = "Rib.Col.",
        gene = row_labels[proxpos],
        stringsAsFactors = FALSE
    )
    background <- .make_rand_sets(proxpos, n, "Rib.Col.", row_labels)

    rbind(foreground, background)
}


