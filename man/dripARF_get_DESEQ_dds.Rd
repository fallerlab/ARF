% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dripARF.R
\name{dripARF_get_DESEQ_dds}
\alias{dripARF_get_DESEQ_dds}
\title{To normalize the counts with DESEQ}
\usage{
dripARF_get_DESEQ_dds(
  samples,
  rRNAs_fasta,
  rRNA_counts = NULL,
  compare = "group",
  organism = NULL,
  exclude = NULL,
  count_threshold = 100,
  QCplot = FALSE,
  targetDir = paste0(getwd(), "/")
)
}
\arguments{
\item{samples}{Samples dataframe created by read_ARF_samples_file() function.}

\item{rRNAs_fasta}{Fasta file for the 4 rRNAs of the organism.}

\item{rRNA_counts}{rRNA_counts that were read by dripARF_read_rRNA_fragments() function (optional)}

\item{compare}{If you want to compare samples based on other grouping, choose the columnname that is given in the samplesFile (Default=group).}

\item{organism}{Organism abbreviation. Pass "hs" for human, "mm" for mouse, and "sc" for yeast.}

\item{exclude}{List of sample names to be excluded from the analysis.}

\item{count_threshold}{Exclude the positions with reads less than this threshold on average. (Default: 1000)}

\item{QCplot}{TRUE or FALSE, whether to generate QC plots or not.}

\item{targetDir}{Directory to save the QC plots in. (Default: working directory, getwd() output)}
}
\description{
A function that normalizes rRNA counts with DESEQ2
}
\examples{
dripARF_get_DESEQ_dds(samples_df, "rRNAs.fa", organism="hs")
}
\keyword{DESeq2}
\keyword{abundance}
\keyword{fragment}
\keyword{normalization}
\keyword{rRNA}
