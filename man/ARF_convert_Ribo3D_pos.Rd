% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ARF_platform.R
\name{ARF_convert_Ribo3D_pos}
\alias{ARF_convert_Ribo3D_pos}
\title{Convert 3D Ribosome distance file to target organism through rRNA alignments!}
\usage{
ARF_convert_Ribo3D_pos(
  source_distance_file,
  source_rRNAs_fasta,
  target_species,
  target_rRNAs_fasta,
  rRNA_pairs = list(),
  source_positions = NULL,
  source_sets = NULL,
  type = "distances"
)
}
\arguments{
\item{source_distance_file}{XXX.ARF.minimum_distances.csv file generated with ARF_parse_PDB_ribosome() for the source 3D Ribosome}

\item{source_rRNAs_fasta}{Fasta file for the rRNAs of the source organism.}

\item{target_species}{ID for the target species, i.e. mm, sc, etc.}

\item{target_rRNAs_fasta}{Fasta file for the rRNAs of the target organism. Same file used in rRNA fragment alignment.}

\item{rRNA_pairs}{List of rRNA ID pairs matching source and target rRNAs (ncessary for lift-overing), i.e. list(c("source_28S","target_28S"), etc. ).}

\item{source_positions}{List of position vectors to convert (Alternative input)}

\item{source_sets}{GSEA sets dataframe to convert (with "ont" & "gene" columns) (Alternative input)}

\item{type}{Which conversion to perform (Currently 'distances', positionset' only)}
}
\description{
Align rRNAs from source and target organism, then, convert source species 3D ribosome data accordingly
}
\examples{
ARF_convert_ribosome3D_rRNA_pos("Ribosome.3D.4V6X.ARF.minimum_distances.csv", "4V6X.rRNAs.fasta", "dm", "drosophila_rRNAs.fa",rRNA_pairs=list(c("rRNA_28S","dm_rRNA_28S"),c("rRNA_18S","dm_rRNA_18S"),c("rRNA_5.8S","dm_rRNA_5.8S"),c("rRNA_5S","dm_rRNA_5S")))
}
\keyword{3D}
\keyword{conversion}
\keyword{coordinates}
\keyword{liftover}
\keyword{rRNA}
