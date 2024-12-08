% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dripARF.R
\name{dripARF_report_RPspec_pos_results}
\alias{dripARF_report_RPspec_pos_results}
\title{Get lFC profiles of RP proximity sets}
\usage{
dripARF_report_RPspec_pos_results(
  samples,
  rRNAs_fasta,
  rRNA_counts = NULL,
  dripARF_dds = NULL,
  organism = NULL,
  compare = "group",
  comparisons = NULL,
  exclude = NULL,
  gsea_sets_RP = NULL,
  RP_proximity_df = NULL
)
}
\arguments{
\item{samples}{Samples dataframe created by read_ARF_samples_file() function.}

\item{rRNAs_fasta}{Fasta file for the 4 rRNAs of the organism.}

\item{rRNA_counts}{rRNA_counts that were read by dripARF_read_rRNA_fragments() function: (optional)}

\item{dripARF_dds}{DESEQ2 normalized rRNA_counts coming from dripARF_get_DESEQ_dds() function: (optional)}

\item{organism}{Organism abbrevation. Pass "hs" for human, "mm" for mouse, and "sc" for yeast.}

\item{compare}{If you want to compare samples based on other grouping, choose the columnname that is given in samplesFile (Default=group).}

\item{comparisons}{List of comparisons to be included.}

\item{exclude}{List of sample names to be excluded from the analysis.}

\item{gsea_sets_RP}{RP-rRNA contact point sets to perform enrichments on.}

\item{RP_proximity_df}{RP-rRNA proximity matrix that is calculated by ARF.}
}
\description{
Fish out the Differential values for RP proximity sets
}
\examples{
dripARF_report_RPspec_pos_results(samples_df, "rRNAs.fa", organism="hs")
}
\keyword{dripARF}
\keyword{for}
\keyword{logFC}
\keyword{profiles}
\keyword{rRNA}
\keyword{rRNAs}
