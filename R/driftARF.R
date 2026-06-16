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


#' Per-position Pearson correlation of rRNA fragment abundance with a sample feature
#' @description Internal helper. Computes, for every rRNA position (row of the
#'   variance-stabilised abundance matrix), the Pearson correlation coefficient and p-value against
#'   a continuous sample feature, vectorised across positions. Positions with zero variance return
#'   \code{NA}.
#' @param mat Numeric matrix of positions (rows) x samples (columns), e.g. the VST-transformed
#'   abundances.
#' @param feature Numeric vector of feature values, one per column of \code{mat} and in the same
#'   order.
#' @return A data.frame with one row per position and columns \code{r} (Pearson coefficient) and
#'   \code{p} (two-sided p-value), row-named by position ID.
#' @keywords internal
driftARF_position_correlation <- function(mat, feature) {
  n <- length(feature)
  fc <- feature - mean(feature)
  mc <- mat - rowMeans(mat)
  denom <- sqrt(rowSums(mc^2) * sum(fc^2))
  r <- as.numeric((mc %*% fc) / denom)
  r[denom == 0] <- NA
  # two-sided p-value from the t distribution
  tval <- r * sqrt((n - 2) / (1 - r^2))
  p <- 2 * stats::pt(-abs(tval), df = n - 2)
  data.frame(r = r, p = p, row.names = rownames(mat))
}


#' Build the GSEA ranking statistic from per-position correlations
#' @description Internal helper. Converts per-position \code{(r, p)} correlations into a single
#'   named ranking vector for \code{clusterProfiler::GSEA}, following the same p-value flooring as
#'   dripARF (tiny p-values capped at 1e-5 so individual positions cannot dominate the ranking).
#' @param cor_df Data.frame with columns \code{r} and \code{p} as returned by
#'   \code{driftARF_position_correlation()}.
#' @param measureID One of \code{"abs_cor_measure"} (default, \code{abs(r) * -log10(p)},
#'   one-tailed), \code{"cor_measure"} (signed, \code{r * -log10(p)}, two-tailed), \code{"abs_r"}
#'   (\code{abs(r)}) or \code{"r"} (signed \code{r}).
#' @return A list with \code{measure} (named numeric ranking vector) and \code{scoreType}
#'   (\code{"pos"} or \code{"std"}) for \code{clusterProfiler::GSEA}.
#' @keywords internal
driftARF_build_measure <- function(cor_df, measureID = "abs_cor_measure") {
  p_floored <- cor_df$p
  p_floored[p_floored < 1e-5] <- 1e-5
  scoreType <- "pos"
  if (measureID == "cor_measure") {
    measure <- cor_df$r * (-log10(p_floored))
    scoreType <- "std"
  } else if (measureID == "abs_r") {
    measure <- abs(cor_df$r)
  } else if (measureID == "r") {
    measure <- cor_df$r
    scoreType <- "std"
  } else { # abs_cor_measure (default)
    measure <- abs(cor_df$r) * (-log10(p_floored))
  }
  names(measure) <- rownames(cor_df)
  list(measure = measure, scoreType = scoreType)
}


#' Predict rRNA position sets whose fragment abundance tracks a continuous feature
#' @description Core driftARF routine. For each requested feature, correlates per-position
#'   variance-stabilised rRNA fragment abundances with the (continuous, numeric) sample feature,
#'   then runs the ARF RPSEA engine (\code{clusterProfiler::GSEA} over \code{gsea_sets_RP} with
#'   z-scoring against the 99 \code{Rand*} control sets) plus an over-representation analysis (ORA)
#'   over significantly-correlated positions. This is the correlation analogue of
#'   \code{\link{dripARF_predict_heterogenity}}: only the per-position ranking statistic changes,
#'   so the same RP heterogeneity, collision and user-given sets are supported.
#' @param samples Samples data.frame from \code{read_ARF_samples_file()}. Must contain the feature
#'   column(s) named in \code{features}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism.
#' @param features Character vector of column names in \code{samples} holding continuous numeric
#'   features. driftARF loops over each feature and emits one result block (and CSV) per feature.
#' @param rRNA_counts Pre-computed rRNA count data.frame from \code{dripARF_read_rRNA_fragments()}
#'   (optional).
#' @param dripARF_dds Pre-computed DESeq2 \code{DESeqDataSet} (optional). If \code{NULL} it is built
#'   from the counts.
#' @param compare Optional column name used to build the DESeq2 model for normalisation. Default
#'   \code{NULL} fits an intercept-only design (\code{~1}), so a \code{group} column is neither
#'   required nor included in the formula. Set it to a column name only if you want that factor in
#'   the normalisation model. Either way the variance-stabilising transform is run with
#'   \code{blind=TRUE}, so the abundance scale used for correlation is not shaped by the design.
#' @param organism Organism abbreviation (\code{"hs"}, \code{"mm"}, \code{"sc"}) or \code{NULL} when
#'   supplying custom sets.
#' @param QCplot Logical, passed to \code{dripARF_read_rRNA_fragments()} (default \code{FALSE}).
#' @param targetDir Output directory for per-feature CSVs (default: working directory).
#' @param exclude Character vector of sample names to exclude.
#' @param measureID Ranking statistic; see \code{driftARF_build_measure()} (default
#'   \code{"abs_cor_measure"}).
#' @param ssRPSEAplots Logical. When \code{TRUE}, write a per-feature interaction plot of
#'   \code{weighted.RPSEA.NES_randZ} vs \code{RPSEA.NES_randZ} (default \code{FALSE}).
#' @param gsea_sets_RP Position sets in \code{ont}/\code{gene} format including the \code{Rand*}
#'   control sets (preset for \code{hs}, \code{mm}, \code{sc}; may already include merged collision
#'   sets).
#' @param RP_proximity_df RP-rRNA proximity matrix (preset for \code{hs}, \code{mm}, \code{sc}).
#' @param runID Label used in output file names (default \code{"driftARF"}).
#' @return A data.frame with one row per (set, feature), columns mirroring
#'   \code{dripARF_predict_heterogenity()} (\code{comp} = feature name, \code{Description},
#'   \code{ORA.*}, \code{RPSEA.NES}, \code{RPSEA.NES_randZ}, \code{RPSEA.padj}, ...), the ssRPSEA
#'   corroboration columns \code{ssRPSEA.weight} and \code{weighted.RPSEA.NES_randZ}, plus
#'   correlation context: \code{set.avg.r}, \code{set.avg.abs.r}, and \code{C1.avg.read.c}/\code{C2.avg.read.c}
#'   (mean normalised abundance in the below-median vs at/above-median feature halves). The
#'   \code{comp}/\code{C1}/\code{C2} columns make the result directly usable with
#'   \code{dripARF_simplify_results()}, \code{dripARF_result_scatterplot()} and
#'   \code{dripARF_result_heatmap()}.
#' @seealso \code{\link{driftARF}}, \code{\link{dripARF_predict_heterogenity}}
#' @keywords driftARF correlation RPSEA progression
#' @export
driftARF_predict_progression <- function(samples, rRNAs_fasta, features, rRNA_counts = NULL, dripARF_dds = NULL,
                                         compare = NULL, organism = NULL, QCplot = FALSE, targetDir = NA,
                                         exclude = NULL, measureID = "abs_cor_measure", ssRPSEAplots = FALSE,
                                         gsea_sets_RP = NULL, RP_proximity_df = NULL, runID = "driftARF") {

  if (is.na(targetDir)) targetDir <- getwd()

  if (!is.null(exclude))
    samples <- samples[!samples[, 1] %in% exclude, ]

  sample_col <- colnames(samples)[1]

  # Validate features are present and numeric
  missing_feats <- features[!features %in% colnames(samples)]
  if (length(missing_feats) > 0)
    stop(paste("Feature column(s) not found in samples file:", paste(missing_feats, collapse = ", ")))

  # Resolve preset sets when not supplied (mirrors dripARF_predict_heterogenity)
  if (is.null(RP_proximity_df)) {
    if (identical(organism, "hs")) {
      RP_proximity_df <- ARF:::RP_proximity_human_df
      if (is.null(gsea_sets_RP)) gsea_sets_RP <- ARF:::human_gsea_sets_RP
    } else if (identical(organism, "mm")) {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
      if (is.null(gsea_sets_RP)) gsea_sets_RP <- ARF:::mouse_gsea_sets_RP
    } else if (identical(organism, "sc")) {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
      if (is.null(gsea_sets_RP)) gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    } else {
      message(paste(c("Organism", organism, "not implemented yet! Please generate your own proximity matrix and RP-rRNA sets using ARF.\n"), collapse = " "))
      return(NULL)
    }
  }

  # Build / fetch the DESeq2 object for normalisation
  if (is.null(dripARF_dds)) {
    if (is.null(rRNA_counts))
      rRNA_counts <- dripARF_read_rRNA_fragments(samples = samples, rRNAs_fasta = rRNAs_fasta, organism = organism, QCplot = QCplot, targetDir = targetDir)
    dds <- ARF::dripARF_get_DESEQ_dds(samples = samples, rRNAs_fasta = rRNAs_fasta, rRNA_counts = rRNA_counts, compare = compare, organism = organism, exclude = exclude)
  } else {
    dds <- dripARF_dds
  }

  # Abundance scale for correlation: VST (blind so the feature is correlated against an
  # unsupervised, variance-stabilised abundance matrix), plus normalised counts for set means.
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  vsd_mat <- SummarizedExperiment::assay(vsd)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  # Sets to report (drop randomised / mito-filler control sets)
  RPs_toreport <- unique(as.character(gsea_sets_RP$ont[!substring(gsea_sets_RP$ont, 1, 3) %in% c("MRf", "FDf", "Ran")]))
  RP_pathways <- sapply(RPs_toreport, function(x) as.character(gsea_sets_RP$gene[gsea_sets_RP$ont == x]))

  # ssRPSEA (feature-correlation) corroboration: per-sample ssGSEA activity per set, then a
  # continuous-covariate limma fit of that activity against each feature. Orthogonal to the
  # position-level GSEA; yields ssRPSEA.weight -> weighted.RPSEA.NES_randZ (see run_DESeq2_norm).
  ssrpsea_weights <- NULL
  tryCatch({
    cat("\nRunning ssRPSEA (feature-correlation) normalization...\n")
    run_list <- run_DESeq2_norm(norm_counts, gsea_sets_RP)
    ssrpsea_weights <- run_limma_cor_analysis(run_list$ssgsea_scores, samples, features)
    cat("ssRPSEA weights computed.\n")
  }, error = function(e) message(paste("ssRPSEA weighting skipped:", conditionMessage(e))))

  all_results <- NULL
  for (feature in features) {
    message(paste("Running driftARF for feature:", feature, "\n"))

    fv <- suppressWarnings(as.numeric(samples[[feature]]))
    if (all(is.na(fv)))
      stop(paste("Feature column", feature, "is not numeric. driftARF accepts continuous numeric features only."))
    names(fv) <- samples[[sample_col]]

    # Align feature to the abundance matrix columns and drop samples with missing feature values
    fv <- fv[colnames(vsd_mat)]
    keep <- !is.na(fv)
    if (sum(keep) < 3) {
      message(paste("Skipping feature", feature, "- fewer than 3 samples with a numeric value.\n"))
      next
    }
    f <- fv[keep]
    submat <- vsd_mat[, keep, drop = FALSE]
    subcounts <- norm_counts[, keep, drop = FALSE]

    # Restrict to positions with non-zero coverage in the retained samples
    nonzero <- rowSums(subcounts) != 0

    cor_df <- driftARF_position_correlation(submat[nonzero, , drop = FALSE], f)
    built <- driftARF_build_measure(cor_df, measureID = measureID)
    measure <- built$measure[!is.na(built$measure)]
    geneList <- sort(measure, decreasing = TRUE)

    egmt <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = gsea_sets_RP, verbose = TRUE,
                                  minGSSize = 10, maxGSSize = 10000, pvalueCutoff = 2, scoreType = built$scoreType)
    egmt@result$NES_rand_zscore <- NA
    for (RP in RPs_toreport) {
      tochange <- endsWith(x = egmt@result$ID, suffix = RP)
      egmt@result$NES_rand_zscore[tochange] <- scale(egmt@result$NES[tochange])
    }

    # ORA over significantly-correlated positions (|r| > 0.5 & p < 0.05)
    sig_genes <- rownames(cor_df)[which(cor_df$p < 0.05 & abs(cor_df$r) > 0.5)]
    or_df <- fgsea::fora(pathways = RP_pathways, genes = sig_genes, universe = rownames(cor_df), minSize = 10)

    res_df <- data.frame(Description = or_df$pathway,
                         ORA.overlap = or_df$overlap, ORA.setSize = or_df$size, ORA.padj = or_df$padj, ORA.p = or_df$pval,
                         RPSEA.NES = NA, RPSEA.NES_randZ = NA, RPSEA.padj = NA, RPSEA.pval = NA, RPSEA.q = NA,
                         stringsAsFactors = FALSE)
    rownames(res_df) <- res_df$Description

    if (dim(egmt@result)[1] > 0) {
      res_df$RPSEA.NES        <- egmt@result[as.character(res_df$Description), "NES"]
      res_df$RPSEA.NES_randZ  <- egmt@result[as.character(res_df$Description), "NES_rand_zscore"]
      res_df$RPSEA.padj       <- p.adjust(egmt@result[as.character(res_df$Description), "pvalue"], method = "BH")
      res_df$RPSEA.pval       <- egmt@result[as.character(res_df$Description), "pvalue"]
      qcol <- intersect(c("qvalue", "qvalues"), colnames(egmt@result))
      if (length(qcol) > 0)
        res_df$RPSEA.q <- egmt@result[as.character(res_df$Description), qcol[1]]
    }

    # ssRPSEA weight for this feature -> weighted positional enrichment
    res_df$ssRPSEA.weight <- NA
    res_df$weighted.RPSEA.NES_randZ <- NA
    if (!is.null(ssrpsea_weights)) {
      fw <- ssrpsea_weights[ssrpsea_weights$comparison == feature, ]
      res_df$ssRPSEA.weight <- fw$ssRPSEA.weight[match(res_df$Description, fw$RP)]
      res_df$weighted.RPSEA.NES_randZ <- res_df$RPSEA.NES_randZ * res_df$ssRPSEA.weight
    }

    # Correlation context per set
    res_df$set.avg.r     <- sapply(res_df$Description, function(s) mean(cor_df[intersect(RP_pathways[[s]], rownames(cor_df)), "r"], na.rm = TRUE))
    res_df$set.avg.abs.r <- sapply(res_df$Description, function(s) mean(abs(cor_df[intersect(RP_pathways[[s]], rownames(cor_df)), "r"]), na.rm = TRUE))

    # Below-median vs at/above-median feature halves, as the C1/C2 abundance analogue
    med <- stats::median(f)
    low <- names(f)[f < med]
    high <- names(f)[f >= med]
    res_df[["C1.avg.read.c"]] <- sapply(res_df$Description, function(s) mean(subcounts[intersect(RP_pathways[[s]], rownames(subcounts)), low, drop = FALSE]))
    res_df[["C2.avg.read.c"]] <- sapply(res_df$Description, function(s) mean(subcounts[intersect(RP_pathways[[s]], rownames(subcounts)), high, drop = FALSE]))

    write.csv(x = res_df, file = paste(targetDir, "/", feature, "_", runID, "_", measureID, "_results.csv", sep = ""), row.names = FALSE)

    # ssRPSEA interaction plot (weighted vs raw NES_randZ) per feature
    if (ssRPSEAplots && !is.null(ssrpsea_weights)) {
      plotdf <- res_df
      plotdf$comparison <- feature
      pdf(paste(targetDir, "/", feature, "_", runID, "_", measureID, "_ssRPSEA_vs_RPSEA.NES_randZ.pdf", sep = ""), height = 5, width = 5)
      print(interaction_plotter(plotdf))
      dev.off()
    }

    # `comp` = feature name so the dripARF simplify/plot helpers work unchanged
    res_df <- data.frame(comp = feature, res_df[order(res_df$RPSEA.NES, decreasing = TRUE), ])
    all_results <- rbind(all_results, res_df)
  }

  if (!is.null(all_results))
    rownames(all_results) <- 1:(dim(all_results)[1])
  return(all_results)
}


#' driftARF - detect ribosome heterogeneity that progresses with a continuous sample feature
#' @description Entry point for the driftARF subtool (\strong{D}ynamic \strong{RI}bosomal
#'   \strong{F}eature \strong{T}racking). Given one or more continuous numeric features in the
#'   samples file, driftARF correlates each per-position rRNA fragment abundance with the feature
#'   and runs the ARF RPSEA + ORA engine to nominate position sets - ribosomal protein (RP)
#'   heterogeneity sets, ribosome collision sets, or user-given sets - whose fragmentation tracks
#'   the feature. It is the continuous-progression analogue of \code{\link{dripARF}} (RP
#'   heterogeneity) and \code{\link{dricARF}} (collisions), which instead contrast discrete groups.
#' @param samplesFile Path to the tab-separated samples file. Must contain the feature column(s)
#'   named in \code{features}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism (the same file used for alignment).
#' @param features Character vector of continuous, numeric column names in the samples file to
#'   correlate against. One result block and CSV are produced per feature.
#' @param samples_df Optional pre-read samples data.frame (skips reading \code{samplesFile}).
#' @param organism Organism abbreviation: \code{"hs"}, \code{"mm"}, \code{"sc"}, or \code{NULL} when
#'   supplying custom sets.
#' @param compare Optional column for the DESeq2 normalisation model. Default \code{NULL} uses an
#'   intercept-only design (\code{~1}); no \code{group} column is required. Provide a column name
#'   only to include that factor in the normalisation model.
#' @param QCplot Logical, whether to draw read-fragment QC plots (default \code{TRUE}).
#' @param targetDir Output directory (default: working directory).
#' @param exclude Character vector of sample names to exclude.
#' @param measureID Ranking statistic; default \code{"abs_cor_measure"}
#'   (\code{abs(r) * -log10(p)}). Use \code{"cor_measure"} for the signed, directional statistic;
#'   see \code{driftARF_build_measure()}.
#' @param ssRPSEAplots Logical. When \code{TRUE}, write a per-feature interaction plot of
#'   \code{weighted.RPSEA.NES_randZ} vs \code{RPSEA.NES_randZ} (default \code{FALSE}).
#' @param include_collision Logical. When \code{TRUE}, ribosome collision sets are merged into the
#'   RP sets (as in \code{dricARF}) so collision sets are scored alongside RP heterogeneity sets
#'   (default \code{FALSE}).
#' @param gsea_sets_RP Custom RP position sets (\code{ont}/\code{gene}) including \code{Rand*}
#'   controls; preset for \code{hs}/\code{mm}/\code{sc}.
#' @param RP_proximity_df Custom RP-rRNA proximity matrix; preset for \code{hs}/\code{mm}/\code{sc}.
#' @param gsea_sets_Collision Custom collision sets; preset for \code{hs}/\code{mm}/\code{sc}. Only
#'   used when \code{include_collision = TRUE}.
#' @return A data.frame of per-(set, feature) predictions; see
#'   \code{\link{driftARF_predict_progression}}. Per-feature CSVs are written to \code{targetDir}.
#' @seealso \code{\link{dripARF}}, \code{\link{dricARF}}, \code{\link{driftARF_predict_progression}}
#' @keywords driftARF correlation RPSEA progression feature
#' @export
#' @examples
#' \dontrun{
#' driftARF(samplesFile = "samples.tsv", rRNAs_fasta = "rRNAs/mouse_rRNAs.fa",
#'          features = c("age", "tumor_grade"), organism = "mm",
#'          targetDir = "driftARF_results/")
#' }
driftARF <- function(samplesFile, rRNAs_fasta, features, samples_df = NULL, organism = NULL, compare = NULL,
                     QCplot = TRUE, targetDir = NA, exclude = NULL, measureID = "abs_cor_measure",
                     ssRPSEAplots = FALSE, include_collision = FALSE, gsea_sets_RP = NULL, RP_proximity_df = NULL,
                     gsea_sets_Collision = NULL) {

  # Resolve preset RP proximity matrix and sets when not supplied
  if (is.null(RP_proximity_df)) {
    if (!ARF_check_organism(organism)) return(NA)
    if (organism == "hs") {
      RP_proximity_df <- ARF:::RP_proximity_human_df
    } else if (organism == "mm") {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
    } else if (organism == "sc") {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
    } else {
      message(paste(c("Organism", organism, "is not preset/implemented yet! Please generate your own sets using ARF.\n"), collapse = " "))
      return(NULL)
    }
  }
  if (is.null(gsea_sets_RP)) {
    if (!ARF_check_organism(organism)) return(NA)
    if (organism == "hs") {
      gsea_sets_RP <- ARF:::human_gsea_sets_RP
    } else if (organism == "mm") {
      gsea_sets_RP <- ARF:::mouse_gsea_sets_RP
    } else if (organism == "sc") {
      gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    } else {
      message(paste(c("Organism", organism, "is not preset/implemented yet! Please generate your own sets using ARF.\n"), collapse = " "))
      return(NULL)
    }
  }

  # Optionally merge collision sets, exactly as dricARF does
  if (include_collision) {
    if (is.null(gsea_sets_Collision)) {
      if (!ARF_check_organism(organism)) return(NA)
      if (organism == "hs") {
        gsea_sets_Collision <- ARF:::human_gsea_sets_Collision
      } else if (organism == "mm") {
        gsea_sets_Collision <- ARF:::mouse_gsea_sets_Collision
      } else if (organism == "sc") {
        gsea_sets_Collision <- ARF:::yeast_gsea_sets_Collision
      } else {
        message(paste(c("Organism", organism, "is not preset/implemented yet! Please generate your own collision sets using ARF.\n"), collapse = " "))
        return(NULL)
      }
    }
    added_sets <- unique(gsea_sets_Collision$ont)
    added_sets <- added_sets[!grepl("Rand", added_sets)]
    for (colset in added_sets) {
      RP_proximity_df[, colset] <- 100
      RP_proximity_df[gsea_sets_Collision$gene[gsea_sets_Collision$ont == colset], colset] <- 1
    }
    gsea_sets_RP <- rbind(gsea_sets_RP, gsea_sets_Collision)
  }

  if (is.na(targetDir)) targetDir <- getwd()

  if (is.null(samples_df))
    samples_df <- read_ARF_samples_file(samplesFile)
  if (!is.null(exclude))
    samples_df <- samples_df[!samples_df[, 1] %in% exclude, ]

  rRNA_counts_df <- dripARF_read_rRNA_fragments(samples = samples_df, rRNAs_fasta = rRNAs_fasta, organism = organism,
                                                QCplot = QCplot, targetDir = targetDir)

  results <- driftARF_predict_progression(samples = samples_df, rRNAs_fasta = rRNAs_fasta, features = features,
                                          rRNA_counts = rRNA_counts_df, compare = compare, organism = organism,
                                          QCplot = QCplot, targetDir = targetDir, exclude = exclude,
                                          measureID = measureID, ssRPSEAplots = ssRPSEAplots, gsea_sets_RP = gsea_sets_RP,
                                          RP_proximity_df = RP_proximity_df, runID = "driftARF")

  if (!is.null(results))
    dripARF_result_scatterplot(dripARF_results = results, targetDir = targetDir, title = "ALL driftARF predictions")

  return(results)
}
