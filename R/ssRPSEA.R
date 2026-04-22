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


#' Run DESeq2 normalization and compute ssGSEA scores
#'
#' This function extracts normalized counts from a DESeq2 object, applies a log
#' transformation, optionally adjusts gene identifiers ending in numeric suffixes,
#' and computes single-sample Gene Set Enrichment Analysis (ssGSEA) scores using
#' the GSVA package.
#'
#' @param norm_counts A matrix of raw or normalized counts. Note: this argument is
#'   currently overwritten internally by \code{DESeq2::counts(dds, normalized = TRUE)},
#'   so it is not used unless the function is modified.
#' @param gsea_sets_RP A data frame containing gene set definitions. Must include
#'   at least two columns:
#'   \describe{
#'     \item{gene}{Gene identifiers}
#'     \item{ont}{Gene set (ontology/pathway) names}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts normalized counts from a DESeq2 object (\code{dds})
#'   \item Applies a \code{log1p} transformation
#'   \item Detects gene identifiers with numeric suffixes and shifts indices if
#'         zero-based numbering is present
#'   \item Filters out gene sets whose names start with "Rand"
#'   \item Computes ssGSEA scores using \code{GSVA::gsva} with method = "ssgsea"
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{ssgsea_scores}{A matrix of ssGSEA enrichment scores (gene sets x samples)}
#'   \item{geneSets}{A list of gene sets used for ssGSEA}
#' }
#'
#' @examples
#' \dontrun{
#' res <- run_DESeq2_norm(norm_counts, gsea_sets_RP)
#' head(res$ssgsea_scores)
#' }
#'
#' @importFrom DESeq2 counts
#' @importFrom GSVA gsva
#' @importFrom dplyr filter
#'
#' @export
run_DESeq2_norm <- function(norm_counts, gsea_sets_RP) {
  stopifnot(
    is.matrix(norm_counts) || is.data.frame(norm_counts),
    all(c("gene", "ont") %in% colnames(gsea_sets_RP))
  )
  
  logcounts <- log1p(norm_counts)

    all_num_ids <- rownames(logcounts) |>
        strsplit('.*_') |> 
        unlist() |> unique() |> 
        as.integer()
    all_num_ids <- all_num_ids[!is.na(all_num_ids)]
    
  if(any(0 %in% all_num_ids)){
    rownames(logcounts) <- paste0(
        sub('^(.*_)(\\d+)$', '\\1', rownames(logcounts)),
        sub('^.*_(\\d+)$', '\\1', rownames(logcounts)) |>
        as.integer()+1 
    )
  }

  ## running ssGSEA
  geneSets_df <- gsea_sets_RP |>
      dplyr::filter(!grepl('^Rand', ont))
  geneSets <- split(geneSets_df$gene, geneSets_df$ont)

  ## compute ssGSEA scores
  ssgsea_scores <- GSVA::gsva(
      expr = logcounts, gset.idx.list = geneSets,
      min.sz=10,
      max.sz=Inf,
      method = "ssgsea", ssgsea.norm = FALSE)
  
  return(list(
    ssgsea_scores = ssgsea_scores,
    geneSets = geneSets
  ))
}



##### using comparisons to check for changes in pathway
## differential activity analysis using limma

#' Perform differential analysis on ssGSEA scores using limma
#'
#' This function applies linear modeling and empirical Bayes moderation (via the
#' limma package) to identify differentially enriched gene sets based on ssGSEA
#' scores across specified sample group comparisons.
#'
#' @param ssgsea_scores A numeric matrix of ssGSEA scores with gene sets (or
#'   pathways) as rows and samples as columns.
#' @param samples A data frame containing sample metadata. Must include:
#'   \describe{
#'     \item{sampleName}{Sample identifiers matching column names of \code{ssgsea_scores}}
#'     \item{group}{Grouping variable used to define experimental conditions}
#'   }
#' @param comparisons A list of character vectors specifying pairwise comparisons.
#'   Each element should be a vector of length 2 indicating the contrast, e.g.
#'   \code{list(c("Treatment", "Control"))}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Z-score normalization of ssGSEA scores using \code{scale()}
#'   \item Cleans sample names and extracts condition labels
#'   \item Constructs a design matrix without intercept (\code{~ 0 + condition})
#'   \item Fits a linear model using \code{limma::lmFit}
#'   \item Builds contrast matrix based on user-defined comparisons
#'   \item Applies contrasts and empirical Bayes moderation
#'   \item Extracts differential results using \code{limma::topTable}
#'   \item Computes a custom weight metric (\code{ssRPSEA.weight})
#' }
#'
#' The \code{ssRPSEA.weight} is defined as:
#' \deqn{(1 - adj.P.Val) / (1 + |logFC|)}
#'
#' @return A data frame containing differential analysis results for all
#'   comparisons, including:
#'   \describe{
#'     \item{RP}{Gene set (row) identifier}
#'     \item{logFC}{Log fold change between conditions}
#'     \item{adj.P.Val}{Adjusted p-value (Benjamini-Hochberg)}
#'     \item{comparison}{Comparison label}
#'     \item{formula}{Contrast formula (e.g., "A - B")}
#'     \item{ssRPSEA.weight}{Custom weight metric}
#'   }
#'
#' @examples
#' \dontrun{
#' comps <- list(c("Treatment", "Control"))
#' res <- run_limma_DE_analysis(ssgsea_scores, samples, comps)
#' head(res)
#' }
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom tidyr separate
#'
#' @export
run_limma_DE_analysis <- function(
  ssgsea_scores,
  samples,
  comparisons
){
  if(!(all(c("sampleName", "group") %in% colnames(samples)))) {
    stop(
      error = "samples dataframe must contain 'sampleName' and 'group' columns. Please check the dataframe structure."
    )
  }
  ## check if samples in ssgsea_scores match those in samples dataframe and the comparisons are valid
  if(!(all(colnames(ssgsea_scores) %in% make.names(samples$sampleName)) &
      all(sapply(comparisons, function(comp) all(comp %in% samples$group))))) {
    stop(
      error = "Mismatch between sample names in ssgsea_scores and samples dataframe, or invalid comparisons. Please check the dataframe structure and comparisons."
    )
  }

  ssgsea_z <- scale(ssgsea_scores)

  samples$sampleName <- make.names(samples$sampleName)
  samples$condition <- factor(sub('^X', '', make.names(samples$group)))

  design <- model.matrix(~ 0 + condition, data = samples)
  colnames(design) <- levels(samples$condition)

  fit <- limma::lmFit(ssgsea_z, design)

  comparisons <- lapply(comparisons, \(x) make.names(x))
  contrast_strings <- paste0(
      sapply(comparisons, `[`, 1),
      "_vs_",
      sapply(comparisons, `[`, 2),
      " = ",
      sapply(comparisons, `[`, 1),
      " - ",
      sapply(comparisons, `[`, 2)
    )

  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_strings,
    levels = design
  )

  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  results_list <- lapply(colnames(contrast_matrix), function(comp) {
    res <- limma::topTable(
      fit2,
      coef = comp,
      number = Inf,
      adjust.method = "BH"
    )
    
    res$RP <- rownames(res)
    res$comparison <- comp
    res
  })

  lm_results <- do.call(rbind, results_list) |>
    as.data.frame() |>
    tidyr::separate(
      col = "comparison",
      into = c("comparison", "formula"),
      sep = " = ",
      remove = TRUE
    )

    lm_results$ssRPSEA.weight <- ((1 - (lm_results$adj.P.Val)) / (
        1 + abs(lm_results$logFC)))
    return(lm_results)
}


#' Create interaction scatter plots for GSEA results
#'
#' This function generates faceted scatter plots comparing weighted NES scores
#' and RPSEA NES scores from GSEA results. It highlights gene sets exceeding
#' specified thresholds and labels them for easier interpretation.
#'
#' @param GSEA_result_df A data frame containing GSEA results. Must include the
#'   following columns:
#'   \describe{
#'     \item{comparison}{Comparison labels used for faceting}
#'     \item{norm_type}{Normalization type (used for grouping)}
#'     \item{weighted.NES_randZ}{Weighted normalized enrichment scores (numeric)}
#'     \item{RPSEA.NES_randZ}{RPSEA normalized enrichment scores (numeric)}
#'     \item{RP}{(Optional but recommended) Gene set identifiers used for labeling}
#'   }
#'
#' @details
#' The function performs the following:
#' \enumerate{
#'   \item Validates presence and type of required columns
#'   \item Creates a scatter plot of \code{weighted.NES_randZ} (x-axis) versus
#'         \code{RPSEA.NES_randZ} (y-axis)
#'   \item Adds reference lines at 0 (solid) and 1 (dashed) for both axes
#'   \item Labels points where either score exceeds 1
#'   \item Facets plots by \code{comparison}
#'   \item Applies a clean theme using \code{ggplot2::theme_bw}
#' }
#'
#' Labels are added using \code{ggrepel::geom_text_repel} to reduce overlap.
#'
#' @return A \code{ggplot} object representing the interaction plot.
#'
#' @examples
#' \dontrun{
#' p <- interaction_plotter(GSEA_result_df)
#' print(p)
#' }
#'
#' @importFrom ggplot2 ggplot geom_point geom_vline geom_hline facet_grid theme_bw theme aes element_text
#' @importFrom dplyr group_by
#' @importFrom ggrepel geom_text_repel
#'
#' @export
interaction_plotter <- function(GSEA_result_df) {
    if(
      !(all(c("Description", "weighted.RPSEA.NES_randZ", "RPSEA.NES_randZ", "comparison"
          ) %in% colnames(GSEA_result_df)))) {
        stop(
          error = "Required columns not found in GSEA_result_df Please check the dataframe structure."
          )
      }
    
    if(
      !(is.numeric(GSEA_result_df$weighted.RPSEA.NES_randZ) &
        is.numeric(GSEA_result_df$RPSEA.NES_randZ))){
        stop(
          error = "weighted.RPSEA.NES_randZ and RPSEA.NES_randZ columns must be numeric. Please check the dataframe structure."
        )
      }

    interaction_plot <- GSEA_result_df |>
        ggplot2::ggplot(
        ) +
        ggplot2::geom_vline(
            xintercept = 0
        ) +
        ggplot2::geom_hline(
            yintercept = 0
        ) +
        ggplot2::geom_hline(
            yintercept = 1, linetype = 'dashed', color = 'grey'
        ) +
        ggplot2::geom_vline(
            xintercept = 1, linetype = 'dashed', color = 'grey'
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = weighted.RPSEA.NES_randZ, y = RPSEA.NES_randZ,
            # color = thresh_95
            ),
            color = 'cadetblue4', 
            alpha = 0.8, size = 1.8) +
        ggrepel::geom_text_repel(
            # data = top_rps,
            # aes(label = RP),
            ggplot2::aes(x = weighted.RPSEA.NES_randZ, y = RPSEA.NES_randZ,
            label = ifelse(
                (weighted.RPSEA.NES_randZ > 1 | RPSEA.NES_randZ > 1), Description, NA)),
            size = 3,
            max.overlaps = Inf,
            min.segment.length = 0
        )  +
        ggplot2::facet_grid(~ comparison) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
            # strip.background = element_blank(),
            strip.text = ggplot2::element_text(face = "bold"),
            legend.position = "top"
        )
  return(interaction_plot)
}
