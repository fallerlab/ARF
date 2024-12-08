% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dripARF.R
\name{dripARF_predict_heterogenity}
\alias{dripARF_predict_heterogenity}
\title{Predict Differential Ribosomal Heterogeneity candidates with rRNA count data}
\usage{
dripARF_predict_heterogenity(
  samples,
  rRNAs_fasta,
  rRNA_counts = NULL,
  dripARF_dds = NULL,
  compare = "group",
  organism = NULL,
  QCplot = FALSE,
  targetDir = NA,
  comparisons = NULL,
  exclude = NULL,
  GSEAplots = FALSE,
  gsea_sets_RP = NULL,
  RP_proximity_df = NULL,
  optimized_run = F,
  measureID = "abs_GSEA_measure_with_dynamic_p",
  runID = "dripARF"
)
}
\arguments{
\item{samples}{Samples dataframe created by read_ARF_samples_file() function.}

\item{rRNAs_fasta}{Fasta file for the 4 rRNAs of the organism.}

\item{rRNA_counts}{rRNA_counts that were read by dripARF_read_rRNA_fragments() function: (optional)}

\item{dripARF_dds}{DESEQ2 normalized rRNA_counts coming from dripARF_get_DESEQ_dds() function: (optional)}

\item{compare}{If you want to compare samples based on other grouping, choose the columnname that is given in the samplesFile (Default=group).}

\item{organism}{Organism abbreviation. Pass "hs" for human, "mm" for mouse, and "sc" for yeast.}

\item{QCplot}{TRUE or FALSE, whether to generate QC plots or not.}

\item{targetDir}{Directory to save the QCplots in.}

\item{comparisons}{List of comparisons to be included.}

\item{exclude}{List of sample names to be excluded from the analysis.}

\item{GSEAplots}{Whether to produce standard GSEA plots.}

\item{gsea_sets_RP}{RP-rRNA contact point sets to perform enrichments on.}

\item{RP_proximity_df}{RP-rRNA proximity matrix that is calculated by ARF.}

\item{optimized_run}{Run in optimized mode for time-saving.}

\item{measureID}{Alternative options for rRNA position ranking for RPSEA, default: abs_GSEA_measure_with_dynamic_p. (avail. abs_GSEA_measure abs,_GSEA_measure_with_p, abs_GSEA_measure_with_dynamic_p, S2N, GSEA_measure, GSEA_measure_with_p, GSEA_measure_with_dynamic_p, abs_w_GSEA_m, w_GSEA_m)}

\item{runID}{runID for output labeling, default: dripARF. (alt. example: dricARF).}
}
\description{
Overlap differential rRNA count data with 3d ribosome, rRNA-RP proximity data and predict heterogeneity candidates across groups.
}
\examples{
dripARF_predict_heterogenity(samples_df,  "rRNAs.fa", rRNA_counts=rRNA_counts_df, organism="hs", QCplot=TRUE)
}
\keyword{Differential}
\keyword{Heterogeneity}
\keyword{RP}
\keyword{Ribosome}
\keyword{dripARF}
\keyword{rRNA}
\keyword{ribosome}
