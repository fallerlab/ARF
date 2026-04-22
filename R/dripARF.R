
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

#' Read rRNA quantification from .bedGraph or .bam files
#' @description Read bedgraph/bam files to create the rRNA count data.
#' @param samples Samples dataframe created by \code{read_ARF_samples_file()}. Must contain at minimum
#'   columns \code{sampleName} (column 1) and the alignment file path (column 2, bedGraph or BAM).
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism (28S, 18S, 5.8S, 5S). Sequence names
#'   are used as rRNA IDs in the output row names.
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param QCplot \code{TRUE} or \code{FALSE}, whether to generate a per-sample read-count violin plot
#'   (default: \code{FALSE}).
#' @param targetDir Directory to save the QC plots in (default: \code{getwd()}).
#' @return A data.frame with rows = rRNA positions (row names formatted as \code{rRNA_pos}, e.g.
#'   \code{human_28S_1042}) and columns = sample names, containing integer read counts at each
#'   position. Positions with any \code{NA} value across samples are removed.
#' @seealso \code{\link{dripARF_get_DESEQ_dds}}
#' @keywords Reader rRNA fragment alignment
#' @export
#' @examples
#' \dontrun{
#' dripARF_read_rRNA_fragments(samples_df, "rRNAs.fa", organism="hs")
#' dripARF_read_rRNA_fragments(samples_df, "rRNAs.fa", organism="hs", QCplot=TRUE, targetDir="./")
#' }
dripARF_read_rRNA_fragments <- function(samples, rRNAs_fasta, organism=NULL, QCplot=FALSE, targetDir=NA) {
  
  # # Check organism first
  # if (!ARF_check_organism(organism))
  #   return(NA)
  
  # Assign target directory
  if (is.na(targetDir)){
    targetDir=getwd()
  }
  
  # Read source rRNAs
  rRNAs_seq <- Biostrings::readBStringSet(file = rRNAs_fasta,use.names = T)
  names(rRNAs_seq) <- sapply(sapply(names(rRNAs_seq),strsplit,split=" ",fixed=T),"[",1)
  
  df <- NULL
  df <- setNames(data.frame(matrix(ncol = length(samples[,1]), nrow = sum(lengths(rRNAs_seq)),data=0)), 
                   samples[,1])
  rownames(df) <- unname(unlist(sapply(names(rRNAs_seq), FUN = function(x){return(paste(x,(1:lengths(rRNAs_seq[x]))-1,sep = "_"))})))
  
  for (i in 1:dim(samples)[1]) {
    sample <- samples[i,1]
    inputFile <- samples[i,2]
    message(paste("Reading bedgraph/bam file for sample",sample,"(",as.character(i),"of",as.character(dim(samples)[1]),")\n"))
    
    if(endsWith(x = inputFile, suffix=".bedGraph") | endsWith(x = inputFile, suffix=".bedgraph")) {
      temp <- read.table(inputFile, stringsAsFactors = FALSE)
      for (i in 1:dim(temp)[1]){
        if (i%%100==0){
          cat(paste("\rReading the bedgraph file %",as.character(round(100*i/(dim(temp)[1]),2))))
          flush.console() 
        }
        all <- temp$V2[i]:(temp$V3[i]-1)
        df[paste(temp$V1[i], all, sep = "_"), sample] = temp$V4[i]
      }
    } else if (endsWith(x = inputFile,suffix=".bam")) {
      tempGrange <- HelloRanges::do_bedtools_genomecov(i = inputFile, bg = TRUE)
      for (i in 1:length(tempGrange$score)){
        if (i%%100==0){
          cat(paste("\rReading the bam file %",as.character(round(100*i/length(tempGrange$score),2))))
          flush.console() 
        }
        all <- tempGrange@ranges@start[i]:(tempGrange@ranges@start[i] + tempGrange@ranges@width[i])
        df[paste(tempGrange@seqnames@values[i], all, sep = "_"), sample] = tempGrange$score[i]
      }
    } else {
      message(paste(inputFile, "is not a bam or bedGraph file!\n"))
      return(NULL)
    }
    cat("\n")
    flush.console() 
  }
  
  message("ALL bedgraphs have been read.\n")
  if (QCplot) {
    g1 <- ggplot2::ggplot(reshape2::melt(df), ggplot2::aes(x=variable, y=value)) +
      ggplot2::geom_violin()+ggplot2::scale_y_log10()+ggplot2::coord_flip()+
      ggplot2::ylab("Positional rRNA fragment counts (log-scale)")+ggplot2::xlab("Samples")
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "raw_rRNA_counts.png", sep = ""), plot = g1, limitsize = FALSE,
                    width = 8, height = 3+(floor(length(samples[,1])**0.5)*3))
  }
  
  return(df[rowSums(is.na(df))==0,])
}


#' Normalise rRNA counts with DESeq2
#' @description Normalizes rRNA counts with DESeq2, estimating size factors and dispersions using a
#'   \code{~ group} design (or a custom grouping column). Positions with mean counts below
#'   \code{count_threshold} are filtered out before fitting.
#' @param samples Samples dataframe created by \code{read_ARF_samples_file()}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism. Used to read fragments if
#'   \code{rRNA_counts} is not supplied.
#' @param rRNA_counts Pre-computed rRNA count data.frame from \code{dripARF_read_rRNA_fragments()}
#'   (optional). If \code{NULL}, counts are read from the files listed in \code{samples}.
#' @param compare Column name in the samples file to use as the grouping variable for DESeq2
#'   (default: \code{"group"}).
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param exclude Character vector of sample names to exclude from the analysis.
#' @param count_threshold Exclude positions whose mean read count across all samples is below this
#'   value (default: 100).
#' @param QCplot \code{TRUE} or \code{FALSE}, whether to generate a sample-similarity heatmap QC
#'   plot (default: \code{FALSE}).
#' @param targetDir Directory to save QC plots and the temporary \code{temp.dds.RData} file in
#'   (default: \code{getwd()}).
#' @return A DESeq2 \code{DESeqDataSet} object fitted with the specified design formula. Also saved
#'   to \code{<targetDir>/temp.dds.RData}. If \code{QCplot=TRUE}, a sample-similarity heatmap is
#'   saved as PDF and TSV.
#' @seealso \code{\link{dripARF_read_rRNA_fragments}}, \code{\link{dripARF_predict_heterogenity}}
#' @keywords DESeq2 normalization rRNA fragment abundance
#' @export
#' @examples
#' \dontrun{
#' dripARF_get_DESEQ_dds(samples_df, "rRNAs.fa", organism="hs")
#' }
dripARF_get_DESEQ_dds <- function(samples, rRNAs_fasta, rRNA_counts=NULL, compare="group", organism=NULL, exclude=NULL, count_threshold=100, 
                                  QCplot=FALSE, targetDir=paste0(getwd(),"/")){
  
  # # Check organism first
  # if (!ARF_check_organism(organism))
  #   return(NA)
  
  if(!is.null(exclude))
    samples <- samples[!samples[,1]%in%exclude,]
  
  if (is.null(rRNA_counts)) {
    rRNA_counts <- dripARF_read_rRNA_fragments(samples = samples, rRNAs_fasta=rRNAs_fasta, organism=organism, QCplot=QCplot, targetDir=targetDir)
  } else{
    rRNA_counts<-rRNA_counts[,samples[,1]]
  }
  
  s_n <- unique(samples[,compare])
  s_l <- length(s_n)
  samples$DESEQcondition <- samples[,compare]
  
  cts <- as.matrix(round(rRNA_counts, digits = 0))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = samples, design = ~DESEQcondition)
  
  keep <- rowSums(DESeq2::counts(dds)) >= count_threshold*dim(samples)[1]
  dds <- dds[keep,]
  dds <- DESeq2::DESeq(dds) # normalize with dds
  
  if (QCplot) {
    vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
    sampleDists <- dist(t(SummarizedExperiment::assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- vsd$sampleName
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
    write.table(x = sampleDistMatrix, file = paste0(targetDir,'/Sample_similarity_distance_vsdbased.tsv'), sep = "\t",quote = FALSE,
                row.names = TRUE, col.names = rownames(sampleDistMatrix))
    pdf(paste0(targetDir,'/Sample_similarity_distance_vsdbased.pdf'),
        width = 3+(floor(length(samples[,1])**0.5)*3), height = 3+(floor(length(samples[,1])**0.5)*3))
    ht <- ComplexHeatmap::pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
    ComplexHeatmap::draw(ht)
    dev.off()
  }
  
  save(file = paste0(targetDir,"/temp.dds.RData"), list = c("dds"))
  return(dds)
}



#' Get average read count of RP proximity sets
#' @description Calculate the average DESeq2-normalised read count for each RP proximity set per
#'   experimental group.
#' @param samples Samples dataframe created by \code{read_ARF_samples_file()}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism.
#' @param rRNA_counts Pre-computed rRNA count data.frame from \code{dripARF_read_rRNA_fragments()}
#'   (optional).
#' @param dripARF_dds DESeq2 \code{DESeqDataSet} from \code{dripARF_get_DESEQ_dds()} (optional).
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param compare Column name in the samples file to use as the grouping variable (default:
#'   \code{"group"}).
#' @param exclude Character vector of sample names to exclude from the analysis (optional).
#' @param gsea_sets_RP RP-rRNA contact point sets to perform enrichments on (optional; preset for
#'   \code{hs}, \code{mm}, and \code{sc}).
#' @param RP_proximity_df RP-rRNA proximity matrix (optional; preset for \code{hs}, \code{mm}, and
#'   \code{sc}; provide for custom organisms via \code{ARF_parse_PDB_ribosome()}).
#' @return A data.frame with RP proximity sets as rows and experimental groups as columns,
#'   containing DESeq2-normalised average read counts per set per group.
#' @keywords rRNA DESeq2 counts
#' @export
#' @examples
#' \dontrun{
#' dripARF_report_RPset_group_counts(samples_df, "rRNAs.fa", organism="hs")
#' }
dripARF_report_RPset_group_counts <- function(samples, rRNAs_fasta, rRNA_counts=NULL, dripARF_dds=NULL,
                                              organism=NULL, compare="group", exclude=NULL, 
                                              gsea_sets_RP=NULL, RP_proximity_df=NULL){
  
  # # Check organism first
  # if (!ARF_check_organism(organism))
  #   return(NA)
  
  # Read samples
  if(!is.null(exclude))
    samples <- samples[!samples[,1]%in%exclude,]
  
  # GEt read counts
  if (is.null(dripARF_dds)){
    if (is.null(rRNA_counts)) {
      rRNA_counts <- dripARF_read_rRNA_fragments(samples, rRNAs_fasta, organism=organism, QCplot=FALSE, targetDir=NA)
    }
    dds <- dripARF_get_DESEQ_dds(samples = samples, rRNAs_fasta=rRNAs_fasta, rRNA_counts = rRNA_counts, compare=compare, organism=organism, exclude=exclude)
  } else {
    dds <- dripARF_dds
  }
  
  
  samples$DESEQcondition <- samples[,compare]
  
  ###################################################
  if(is.null(RP_proximity_df)){
    if (organism=="hs"){
      RP_proximity_df <- ARF:::RP_proximity_human_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::human_gsea_sets_RP
    } else if (organism=="mm") {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::mouse_gsea_sets_RP
    } else if (organism=="sc") {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    } else {
      message(paste(c("Organism", organism, "Not implemented yet! Please generate your own set of proximity matrix and RP-rRNA sets using ARF.\n"),
                collapse = " "))
      return(NULL)
    }
  }
  
  RPs_toreport <- unique(as.character(gsea_sets_RP$ont[!substring(gsea_sets_RP$ont,1,3)%in%c("MRf","FDf","Ran")]))
  
  ###############################################################
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
  counts <- DESeq2::counts(dds,normalized=TRUE)
  
  results <- NULL
  for (RP in RPs_toreport){
    for (group in unique(samples$DESEQcondition)){
      if(sum(samples$DESEQcondition==group)>1){
        results<-rbind(results, data.frame(RP=RP, group=group,
          AvgCount=mean(rowMeans(counts[,samples[,1][samples$DESEQcondition==group]],na.rm = TRUE)[gsea_sets_RP$gene[gsea_sets_RP$ont==RP]],na.rm = TRUE)))
      }
    }
  }
  
  results_df <- reshape2::dcast(results, RP~group, value.var = "AvgCount")
  rownames(results_df) <- results_df$RP
  return(results_df)
}


#' Predict Differential Ribosomal Heterogeneity candidates
#' @title Predict Differential Ribosomal Heterogeneity candidates
#' @description Overlaps differential rRNA count data with 3D ribosome rRNA-RP proximity data and
#'   predicts ribosomal heterogeneity candidates across groups. For each comparison, rRNA positions
#'   are ranked by a fold-change-based enrichment measure, then Ribosomal Protein Set Enrichment
#'   Analysis (RPSEA) and Over-Representation Analysis (ORA) are performed. Note: the function name
#'   contains a historical spelling typo (\code{heterogenity} instead of \code{heterogeneity}) that
#'   is preserved for backward compatibility.
#' @param samples Samples dataframe created by \code{read_ARF_samples_file()}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism.
#' @param rRNA_counts Pre-computed rRNA count data.frame from \code{dripARF_read_rRNA_fragments()}
#'   (optional).
#' @param dripARF_dds DESeq2 \code{DESeqDataSet} from \code{dripARF_get_DESEQ_dds()} (optional).
#' @param compare Column name in the samples file to use as the grouping variable (default:
#'   \code{"group"}).
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param QCplot \code{TRUE} or \code{FALSE}, whether to generate QC plots (default: \code{FALSE}).
#' @param targetDir Directory to save result CSVs and QC plots in.
#' @param comparisons Named list of comparisons to include (each element a 2-element character
#'   vector \code{c("condition1", "condition2")}). If \code{NULL}, all pairwise comparisons are run.
#' @param exclude Character vector of sample names to exclude from the analysis.
#' @param GSEAplots \code{TRUE} or \code{FALSE}, whether to produce per-RP standard GSEA enrichment
#'   plots (saved as PDFs in \code{targetDir}, default: \code{FALSE}).
#' @param gsea_sets_RP RP-rRNA contact point sets to perform enrichments on (preset for \code{hs},
#'   \code{mm}, and \code{sc}; provide for custom organisms via \code{dripARF_get_RP_proximity_sets()}).
#' @param RP_proximity_df RP-rRNA proximity matrix (preset for \code{hs}, \code{mm}, and \code{sc};
#'   provide for custom organisms via \code{ARF_parse_PDB_ribosome()}).
#' @param optimized_run \code{TRUE} to skip redundant DESeq2 calls and reuse pre-computed objects,
#'   reducing runtime for large datasets (default: \code{FALSE}).
#' @param measureID Ranking measure used to score rRNA positions for RPSEA. Options:
#'   \code{"abs_GSEA_measure_with_dynamic_p"} (default), \code{"abs_GSEA_measure"},
#'   \code{"abs_GSEA_measure_with_p"}, \code{"S2N"}, \code{"GSEA_measure"},
#'   \code{"GSEA_measure_with_p"}, \code{"GSEA_measure_with_dynamic_p"},
#'   \code{"abs_w_GSEA_m"}, \code{"w_GSEA_m"}.
#' @param runID Label used for output file naming (default: \code{"dripARF"}).
#' @return A data.frame with one row per (RP, comparison) combination, with columns: \code{comp},
#'   \code{Description} (RP name), \code{ORA.overlap}, \code{ORA.setSize}, \code{ORA.padj},
#'   \code{ORA.p}, \code{RPSEA.NES}, \code{RPSEA.NES_randZ}, \code{RPSEA.padj}, \code{RPSEA.pval},
#'   \code{RPSEA.q}, \code{C1.avg.read.c}, \code{C2.avg.read.c}. One CSV per comparison is also
#'   written to \code{targetDir}.
#' @seealso \code{\link{dripARF_simplify_results}}, \code{\link{dripARF_result_scatterplot}},
#'   \code{\link{dripARF_result_heatmap}}
#' @keywords dripARF Differential Ribosome Heterogeneity rRNA ribosome RP
#' @export
#' @examples
#' \dontrun{
#' dripARF_predict_heterogenity(samples_df, "rRNAs.fa",
#'   rRNA_counts=rRNA_counts_df, organism="hs", QCplot=TRUE)
#' }
dripARF_predict_heterogenity <- function(samples, rRNAs_fasta, rRNA_counts=NULL, dripARF_dds=NULL,
                                         compare="group", organism=NULL, QCplot=FALSE, targetDir=NA,
                                         comparisons=NULL, exclude=NULL, ssRPSEAplots=FALSE,
                                         GSEAplots=FALSE, gsea_sets_RP=NULL, RP_proximity_df=NULL, optimized_run=F, 
                                         measureID="abs_GSEA_measure_with_dynamic_p", runID='dripARF') {
  # # Check organism first
  # if (!ARF_check_organism(organism))
  #   return(NA)
  
  # Assign target directory
  if (is.na(targetDir)){
    targetDir=getwd()
  }
  
  # Read samples
  if(!is.null(exclude))
    samples <- samples[!samples[,1]%in%exclude,]
  
  # GEt read counts
  if (is.null(dripARF_dds)){
    if (is.null(rRNA_counts)) {
      rRNA_counts <- dripARF_read_rRNA_fragments(samples = samples, rRNAs_fasta=rRNAs_fasta, organism=organism, QCplot = QCplot, targetDir=targetDir)
    }
    dds <- ARF::dripARF_get_DESEQ_dds(samples = samples, rRNAs_fasta=rRNAs_fasta, rRNA_counts = rRNA_counts, compare=compare, organism=organism, exclude=exclude)
  } else {
    dds <- dripARF_dds
  }
  
  s_n <- unique(samples[,compare])
  s_l <- length(s_n)
  samples$DESEQcondition <- samples[,compare]
  
  if(is.null(comparisons) || length(comparisons)==0) {
    comparisons <- list()
    for (i in 1:(s_l-1)) {
      for (j in (i+1):s_l) {
        comparisons[[(length(comparisons) +1)]] <- c(s_n[i],s_n[j])
        if(optimized_run==F){
          message(paste("Comparing",s_n[i],"vs",s_n[j],"\n"))
          # res <- DESeq2::results(dds, contrast=c("DESEQcondition",s_n[i],s_n[j]), lfcThreshold = 0.5, alpha = 0.05, cooksCutoff = FALSE)
          res <- DESeq2::results(dds, contrast=c("DESEQcondition",s_n[i],s_n[j]), cooksCutoff = FALSE)
          print(DESeq2::summary(res))
        }
      }
    }
  }
  
  # Read count transformations
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
  
  if (QCplot) {
    g1 <- ggplot2::ggplot(reshape2::melt(SummarizedExperiment::assay(vsd)), ggplot2::aes(x=Var2, y=value)) +
      ggplot2::geom_violin()+ggplot2::coord_flip()+
      ggplot2::xlab("Positional rRNA fragment counts (vsd) after DESeq2 normalization")+ggplot2::ylab("Samples")
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "norm_rRNA_counts.png", sep = ""), limitsize = FALSE, plot = g1,
                    width = 8, height = 3+(floor(length(samples[,1])**0.5)*3))
    
    #PCA_analysis
    pca <- prcomp(SummarizedExperiment::assay(vsd)[order(matrixStats::rowVars(SummarizedExperiment::assay(vsd)),decreasing = TRUE),])
    pca_df <- data.frame(pca$rotation)
    pca_map <- paste(colnames(pca_df),' (Var.Exp.:%',round(100*((pca$sdev^2)/sum(pca$sdev^2)),1),')', sep='')
    pca_df$sample <- rownames(pca_df)
    pca_df$group <- sapply(pca_df$sample,FUN=function(x){return(samples$DESEQcondition[samples[,1]==x])})
    
    g1 <- ggplot2::ggplot(pca_df, ggplot2::aes(x=PC1,y=PC2,color=group)) + ggplot2::geom_point() +
      ggplot2::labs(x=pca_map[1], y=pca_map[2]) + ggrepel::geom_text_repel(ggplot2::aes(label=sample), show.legend = FALSE)
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "vsd_based_PCA.png", sep = ""), plot = g1,limitsize = FALSE,
                    width = 2+(floor(length(samples[,1])**0.5)*2), height = 2+(floor(length(samples[,1])**0.5)*2))
  }
  
  ###################################################
  if(is.null(RP_proximity_df)){
    if (organism=="hs"){
      RP_proximity_df <- ARF:::RP_proximity_human_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::human_gsea_sets_RP
    } else if (organism=="mm") {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::mouse_gsea_sets_RP
    } else if (organism=="sc") {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    } else {
      message(paste(c("Organism", organism, "Not implemented yet! Please generate your own set of proximity matrix and RP-rRNA sets using ARF.\n"), 
                collapse = " "))
      return(NULL)
    }
  }
  ## Fix this later - Exclude Randomized Control sets from report
  RPs_toreport <- unique(as.character(gsea_sets_RP$ont[!substring(gsea_sets_RP$ont,1,3)%in%c("MRf","FDf","Ran")]))
  
  ########## Overrepresentation Analysis ############
  RP_pathways <- sapply(RPs_toreport,FUN=function(x){return(as.character(gsea_sets_RP$gene[gsea_sets_RP$ont==x]))})
  
  ########### Group specific means ##################
  RP_means <- dripARF_report_RPset_group_counts(samples = samples, rRNAs_fasta=rRNAs_fasta, 
                                                rRNA_counts = rRNA_counts,dripARF_dds = dds,
                                                organism = organism, compare = compare, exclude = exclude, 
                                                gsea_sets_RP=gsea_sets_RP, RP_proximity_df=RP_proximity_df)

  ### ssRPSEA hook ####
  # normalization for ssRPSEA
  cat("\nRunning ssRPSEA normalization...\n")
  counts <- DESeq2::counts(dds, normalized=TRUE)
  run_DESeq2_norm_list <- run_DESeq2_norm(counts, gsea_sets_RP)
  
  # computing ssRPSEA weights
  ssgsea_scores = run_DESeq2_norm_list$ssgsea_scores
  # geneSets = run_DESeq2_norm_list$geneSets

  # limma analysis on ssRPSEA scores
  lm_ssRPSEA_weights_df <- run_limma_DE_analysis(
    ssgsea_scores = ssgsea_scores,
    samples = samples,
    comparisons = comparisons
  )

  cat("\nssRPSEA weights computed.\n")
  ####

  ######### Gene set enrichment Analysis ############
  #library(clusterProfiler)
  #library(enrichplot)
  #library(gprofiler2)
  
  # counts <- DESeq2::counts(dds, normalized=TRUE)
  all_GSEA_results <- NULL

  # Separate for each DESEQcondition
  for (comp in comparisons) {
    message(paste("Running predictions for",comp[1],"vs",comp[2],"\n"))
    
    GSEA_result_df <- NULL
    res <- DESeq2::results(dds, contrast=c("DESEQcondition",comp[1],comp[2]), cooksCutoff = FALSE)
    
    temp_df <- res
    temp_df$weight <- scales::rescale(log10(temp_df$baseMean), to = c(0, 5))
    temp_df$padj[temp_df$padj<0.00001] = 0.00001
    temp_df$pvalue[temp_df$pvalue<0.00001] = 0.00001
    
    #measureID = "abs_GSEA_measure" #)){ #},"GSEA_measure","w_GSEA_m", "abs_w_GSEA_m")){
    
    scoreType = "pos" # for the clusterProfiler::GSEA function : "pos" "neg" one-tailed or "std" two-tailed
    used_measure <- NULL
    if (measureID=="GSEA_measure"){
      used_measure <- temp_df$log2FoldChange*(-log10(temp_df$padj))
      names(used_measure)<-rownames(temp_df)
      scoreType = "std"
    } else if (measureID=="abs_GSEA_measure"){
      used_measure <- abs(temp_df$log2FoldChange)*(-log10(temp_df$padj))
      names(used_measure)<-rownames(temp_df)
    } else if (measureID=="w_GSEA_m"){
      used_measure <- temp_df$log2FoldChange*(-log10(temp_df$padj))*temp_df$weight
      names(used_measure)<-rownames(temp_df)
      scoreType = "std"
    } else if (measureID=="abs_w_GSEA_m"){
      used_measure <- abs(temp_df$log2FoldChange)*(-log10(temp_df$padj))*temp_df$weight
      names(used_measure)<-rownames(temp_df)
    }else if (measureID=="abs_GSEA_measure_with_p"){
      used_measure <- abs(temp_df$log2FoldChange)*(-log10(temp_df$pvalue))
      names(used_measure)<-rownames(temp_df)
    }else if (measureID=="S2N"){
      # used_measure <- sapply(1:(dim(counts)[1]), FUN = function(x){
      #   a = counts[,samples$sampleName[samples$DESEQcondition==comp[1]]]
      #   b = counts[,samples$sampleName[samples$DESEQcondition==comp[2]]]
      #   return()})
      used_measure <- (matrixStats::rowMeans2(counts[,samples$sampleName[samples$DESEQcondition==comp[1]]]) -
                matrixStats::rowMeans2(counts[,samples$sampleName[samples$DESEQcondition==comp[2]]])) /
        (matrixStats::rowSds(counts[,samples$sampleName[samples$DESEQcondition==comp[1]]])+
           matrixStats::rowSds(counts[,samples$sampleName[samples$DESEQcondition==comp[2]]]))
      names(used_measure) <- rownames(counts)
    }else if (measureID=="abs_GSEA_measure_with_dynamic_p"){
      used_measure <- abs(temp_df$log2FoldChange)*(-log10(temp_df$padj))
      if (sum(duplicated(temp_df$padj))>2000) {
        used_measure <- abs(temp_df$log2FoldChange)*(-log10(temp_df$pvalue))
      }
      names(used_measure)<-rownames(temp_df)
    }else if (measureID=="GSEA_measure_with_p"){
      used_measure <- temp_df$log2FoldChange*(-log10(temp_df$pvalue))
      names(used_measure)<-rownames(temp_df)
      scoreType = "std"
    }else if (measureID=="GSEA_measure_with_dynamic_p"){
      used_measure <- temp_df$log2FoldChange*(-log10(temp_df$padj))
      if (sum(duplicated(temp_df$padj))>2000) {
        used_measure <- temp_df$log2FoldChange*(-log10(temp_df$pvalue))
      }
      names(used_measure)<-rownames(temp_df)
      scoreType = "std"
    }else{ # same as abs_GSEA_measure
      used_measure <- abs(temp_df$log2FoldChange)*(-log10(temp_df$padj))
      names(used_measure)<-rownames(temp_df)
    }
    
    used_measure <- used_measure[rowSums(DESeq2::counts(dds)[,samples$sampleName[samples$DESEQcondition%in%comp]])!=0]
    used_measure <- used_measure[!sapply(used_measure, function(x) is.na(x))]

    used_geneList <- used_measure[order(used_measure, decreasing = TRUE)]
    used_geneList_abs <- abs(used_measure)[order(abs(used_measure), decreasing = TRUE)]
    
      
    egmt_used_measure <- clusterProfiler::GSEA(geneList = used_geneList, TERM2GENE=gsea_sets_RP, verbose=TRUE,
                                               minGSSize = 10, maxGSSize = 10000,
                                               pvalueCutoff = 2, scoreType = scoreType)
    egmt_used_measure@result$NES_rand_zscore <- NA
    for (RP in RPs_toreport){
      # Change This, make it more accurate, like removing Rand_ and then check equality
      tochange <- endsWith(x = egmt_used_measure@result$ID, suffix = RP) 
      egmt_used_measure@result$NES_rand_zscore[tochange] <- scale(egmt_used_measure@result$NES[tochange])
    }
    
    ### Overrepresentation hook ####
    if(measureID=="abs_GSEA_measure_with_p"){
      or_df <- fgsea::fora(pathways = RP_pathways, genes = rownames(res)[which(res$pvalue<.05 & abs(res$log2FoldChange)>0.5)],
                         universe = rownames(res), minSize = 10)
    } else{
      or_df <- fgsea::fora(pathways = RP_pathways, genes = rownames(res)[which(res$padj<.05 & abs(res$log2FoldChange)>0.5)],
                           universe = rownames(res), minSize = 10)
    }
    
    GSEA_result_df <- data.frame(Description = or_df$pathway,
                                 ORA.overlap = or_df$overlap, ORA.setSize = or_df$size, ORA.padj = or_df$padj, ORA.p = or_df$pval,
                                 RPSEA.NES=NA, RPSEA.NES_randZ=NA, RPSEA.padj=NA, RPSEA.pval=NA, RPSEA.q=NA)
    rownames(GSEA_result_df) <- GSEA_result_df$Description
    
    if(dim(egmt_used_measure@result)[1]>0){
      if (GSEAplots){
        pdf(paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_", measureID,".pdf",sep = ""),height = 5,width = 5)
        for (i in which(egmt_used_measure@result$Description%in%RPs_toreport)){
          if(egmt_used_measure@result$NES_rand_zscore[i]>1 & egmt_used_measure@result$p.adjust[i]<0.01)
            plotTitle <- paste(egmt_used_measure$Description[i], "NES=",as.character(round(egmt_used_measure$NES[i],2)),
                  "; adjP=",as.character(round(egmt_used_measure$p.adjust[i],4)))
            if("qvalue" %in% colnames(temp)){
              plotTitle <- paste(plotTitle, "; FDR=",as.character(round(egmt_used_measure$qvalue[i],4)))
            } else if("qvalues" %in% colnames(temp)){
              plotTitle <- paste(plotTitle, "; FDR=",as.character(round(egmt_used_measure$qvalues[i],4)))
            } 
            print(enrichplot::gseaplot2(egmt_used_measure, geneSetID = i, title = plotTitle))
        }
        dev.off()
      }
      
      # Assign RPSEA scores, padjust for only valid reported sets
      GSEA_result_df$RPSEA.NES <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"NES"]
      GSEA_result_df$RPSEA.NES_randZ <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"NES_rand_zscore"]
      GSEA_result_df$RPSEA.padj <- p.adjust(p = (egmt_used_measure@result[as.character(GSEA_result_df$Description),"pvalue"]), method = "BH")
      GSEA_result_df$RPSEA.pval <- (egmt_used_measure@result[as.character(GSEA_result_df$Description),"pvalue"])
      
      #GSEA_result_df[RP,"RPSEA.q"] <- temp$qvalues
      if("qvalue" %in% colnames(egmt_used_measure@result)){
        GSEA_result_df$RPSEA.q <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"qvalue"]
      } else if("qvalues" %in% colnames(egmt_used_measure@result)){
        GSEA_result_df$RPSEA.q <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"qvalues"]
      } else{
        print("qvalue is not part of the egmt_result!!")
        GSEA_result_df[RP,"RPSEA.q"] <- NULL
      }
    }
    
    ### Add ssRPSEA weights to the GSEA result df for comparison
    ## Check for NA values in limma weights before join
    tryCatch(
      {stopifnot(!any(is.na(lm_ssRPSEA_weights_df$weight)))},
      error = function(e) print(
        "NA values found in limma weights! Check limma output.")
    )

    GSEA_result_df <- dplyr::left_join(
        GSEA_result_df |>
          dplyr::mutate(comparison = paste(comp, collapse = "_vs_")),
        lm_ssRPSEA_weights_df |>
          dplyr::select(RP, comparison, ssRPSEA.weight) |> 
          dplyr::filter(comparison == paste(comp, collapse = "_vs_")),
        by = c("Description"="RP", "comparison"="comparison")
      ) |>
      dplyr::mutate(weighted.RPSEA.NES_randZ = RPSEA.NES_randZ * ssRPSEA.weight) |>
      dplyr::arrange(desc(weighted.RPSEA.NES_randZ))

    tryCatch(
      # {!any(is.na(GSEA_result_df$weight))},
      {stopifnot(dim(GSEA_result_df)[1] >= 1)},
      error = function(e) print(
        "NA values found in GSEA result weights after join! Check limma output and GSEA result.")
    )

    ## making the interaction plot
    if (ssRPSEAplots){
      pdf(paste(targetDir,"/",
          paste(comp,collapse = "_vs_"),"_", measureID,"_ssRPSEA_vs_RPSEA.NES_randZ.pdf",sep = ""
        ),height = 5,width = 5
      )

      print(interaction_plotter(GSEA_result_df))
      dev.off()
    }
    ###

    GSEA_result_df <- GSEA_result_df |>
      dplyr::select(-comparison)
    ## shifted to the end to add the ssRPSEA weights to the output csvs as well
    GSEA_result_df[,paste0(comp[1],".avg.read.c")] <- RP_means[or_df$pathway,comp[1]]
    GSEA_result_df[,paste0(comp[2],".avg.read.c")] <- RP_means[or_df$pathway,comp[2]]

    if (!is.null(gsea_sets_RP)){
      # if(measureID=="abs_GSEA_measure") {
      #   write.csv(x =  GSEA_result_df, file = paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_dripARF_default_results.csv",sep = ""), row.names = FALSE)
      # } else {
      write.csv(x =  GSEA_result_df, file = paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_",runID,"_",measureID,"_results.csv",sep = ""), row.names = FALSE)
      # }
    }
    colnames(GSEA_result_df)[(length(colnames(GSEA_result_df))-1) : length(colnames(GSEA_result_df))] <- c("C1.avg.read.c","C2.avg.read.c")
    all_GSEA_results <- rbind(all_GSEA_results, data.frame(comp=paste(comp,collapse = "_vs_"),
                                                           GSEA_result_df[order(GSEA_result_df$RPSEA.NES,decreasing = TRUE),]))
  }

  rownames(all_GSEA_results) <- 1:(dim(all_GSEA_results)[1])
  return(all_GSEA_results)
}



#' Simplify dripARF result file
#' @description Filter dripARF results based on statistical thresholds, producing a condensed
#'   data.frame with pass/fail columns and a summary \code{top} logical flag.
#' @param dripARF_results dripARF results data.frame as returned by \code{dripARF_predict_heterogenity()}.
#' @param randZscore_thr Z-score threshold for the RPSEA NES relative to its randomised control
#'   sets. An RP passes this filter when its NES z-score is >= this value (default: 1).
#' @param ORA_adjP_thr Adjusted p-value threshold for the over-representation analysis (ORA). An RP
#'   passes this filter when ORA BH-adjusted p-value <= this value (default: 0.05).
#' @param RPSEA_adjP_thr Adjusted p-value threshold for the RPSEA enrichment score. An RP passes
#'   this filter when RPSEA BH-adjusted p-value <= this value (default: 0.05).
#' @param ORA_sig_n Minimum number of significantly differential rRNA positions (ORA overlap)
#'   required for an RP to be considered. Default: 0 in \code{dripARF_simplify_results} (i.e. any
#'   overlap is accepted); default: 1 in plotting functions.
#' @return A simplified data.frame with columns: \code{comp}, \code{Description}, \code{top}
#'   (logical, \code{TRUE} if all thresholds are passed), \code{ES1} (RPSEA NES), \code{ES2}
#'   (RPSEA NES z-score), \code{ORA} (logical), \code{RPSEA.padj}, \code{ES1.pass},
#'   \code{ES2.pass}, \code{ORA.sigN}, \code{C1.avg.read.c}, \code{C2.avg.read.c}.
#' @seealso \code{\link{dripARF_predict_heterogenity}}, \code{\link{dripARF_result_scatterplot}},
#'   \code{\link{dripARF_result_heatmap}}
#' @keywords dripARF results simplify
#' @export
#' @examples
#' \dontrun{
#' dripARF_simplify_results(dripARF_results)
#' }
dripARF_simplify_results <- function(dripARF_results, randZscore_thr=1, ORA_adjP_thr=0.05, RPSEA_adjP_thr=0.05, ORA_sig_n=0) {
  temp <- dripARF_results[, c("comp", "Description", "C1.avg.read.c", "C2.avg.read.c")]
  temp$ES1 <- dripARF_results$RPSEA.NES
  temp$ES1.pass <- (dripARF_results$RPSEA.padj<RPSEA_adjP_thr)
  temp$ES2 <- dripARF_results$RPSEA.NES_randZ
  temp$ES2.pass <- (dripARF_results$RPSEA.NES_randZ >= randZscore_thr)
  temp$RPSEA.padj <- dripARF_results$RPSEA.padj
  temp$ORA <- (dripARF_results$ORA.padj < ORA_adjP_thr)
  temp$ORA.sigN <- (dripARF_results$ORA.overlap >= ORA_sig_n)
  temp$top <- (dripARF_results$RPSEA.padj <= RPSEA_adjP_thr & dripARF_results$RPSEA.NES_randZ >= randZscore_thr &
                 dripARF_results$ORA.padj <= ORA_adjP_thr & dripARF_results$ORA.overlap >= ORA_sig_n)
  
  return(temp[,c(c("comp", "Description", "top", "ES1", "ES2", "ORA", "RPSEA.padj", "ES1.pass", "ES2.pass", "ORA.sigN", "C1.avg.read.c", "C2.avg.read.c"))])
}



#' Draw heatmap for multiple comparisons
#' @description Draw a heatmap of RPSEA NES and NES z-scores across all comparisons, filtered by one
#'   or more threshold combinations. Each element of the threshold vectors defines one threshold
#'   combination; separate PDFs are saved for each.
#' @param dripARF_results dripARF result data.frame as returned by \code{dripARF_predict_heterogenity()}.
#' @param title Title string used for the plot and output filenames.
#' @param targetDir Directory to save the PDF heatmaps in.
#' @param addedRPs Character vector of RP names to always include in the plot regardless of
#'   significance (optional).
#' @param randZscore_thr Numeric vector of z-score thresholds for the RPSEA NES relative to its
#'   randomised control sets. An RP passes this filter when its NES z-score is >= this value. Each
#'   element defines one threshold combination (default: \code{c(1, 1, 1.5, 0)}).
#' @param ORA_adjP_thr Numeric vector of adjusted p-value thresholds for the over-representation
#'   analysis (ORA). An RP passes this filter when ORA BH-adjusted p-value <= this value (default:
#'   \code{c(0.05, 2, 0.05, 2)}).
#' @param RPSEA_adjP_thr Numeric vector of adjusted p-value thresholds for the RPSEA enrichment
#'   score. An RP passes this filter when RPSEA BH-adjusted p-value <= this value (default:
#'   \code{c(0.05, 0.05, 0.05, 0.05)}).
#' @param ORA_sig_n Numeric vector specifying the minimum number of significantly differential rRNA
#'   positions (ORA overlap) required for an RP to be considered (default: \code{c(1, 1, 1, 1)}).
#' @return Invisibly returns the last \code{ComplexHeatmap} object drawn. PDFs are saved to
#'   \code{targetDir} for each threshold combination.
#' @keywords dripARF Differential Ribosome Heterogeneity Heatmap
#' @export
#' @examples
#' \dontrun{
#' dripARF_result_heatmap(dripARF_results, "Title", "/Folder/to/save/in/")
#' dripARF_result_heatmap(dripARF_results, "Title", "/Folder/to/save/in/", randZscore_thr=1)
#' }
dripARF_result_heatmap <- function(dripARF_results, title, targetDir, addedRPs=NULL,
                                   randZscore_thr=c(1,1,1.5,0), ORA_adjP_thr=c(.05,2,.05,2),
                                   RPSEA_adjP_thr=c(.05,.05,.05,.05), ORA_sig_n=c(1,1,1,1)){
  
  thrlist <- list(z=randZscore_thr,p=RPSEA_adjP_thr,o=ORA_adjP_thr,n=ORA_sig_n)
  
  palette1 = circlize::colorRamp2(c(0,1, 1.5, 2),
                                  c("white","#D9D0D3","#CCBA72","#0F0D0E"))
  
  for (i in 1:length(thrlist[[1]])) {
    message(paste("Threshold",i,"done!\n"))
    
    simplified <- dripARF_simplify_results(dripARF_results,
                                           RPSEA_adjP_thr=thrlist[["p"]][i],
                                           randZscore_thr = thrlist[["z"]][i],
                                           ORA_adjP_thr = thrlist[["o"]][i],
                                           ORA_sig_n=thrlist[["n"]][i])
    
    pdf(file = paste(targetDir,"/", title,
                     "_randZ", as.character(thrlist[["z"]][i]*100),
                     "_GSEAp", as.character(thrlist[["p"]][i]*100),
                     "_ORAp", as.character(thrlist[["o"]][i]*100),
                     "_SigPosN", as.character(thrlist[["n"]][i]),
                     ".pdf",sep=""), width = 16, height = 16)
    tryCatch(
      expr = {
        set.seed(i)
        
        allstudies_abs_sig_RPS <- unique(simplified$Description[simplified$top])
        temp <- simplified[simplified$top, ]
        
        if (length(addedRPs)>0)
          temp <- unique(rbind(temp, simplified[simplified$Description %in% addedRPs, ]))
        
        comp_all_abs_matrix <- reshape2::acast(temp, formula=Description ~ comp, value.var = "ES1")
        comp_all_abs_matrix[is.na(comp_all_abs_matrix)] <- 0
        ht <- ComplexHeatmap::Heatmap(comp_all_abs_matrix, name = "ES1", na_col = "#0F0D0E",
                                      column_title = paste0(title," ES1 - RPSEA NES ",
                                                            "\n(ORA padj<",as.character(ORA_adjP_thr[i]),
                                                            ", RPSEA padj<",as.character(RPSEA_adjP_thr[i]),
                                                            ", RPSEA_rand z-score>",as.character(randZscore_thr[i]),")"),
                                      col = palette1)
        
        last <- ht
        ComplexHeatmap::draw(ht, padding = grid::unit(c(30, 2, 2, 2), "mm"))
        
        ## NESrand
        comp_all_abs_matrix <- reshape2::acast(temp, formula=Description ~ comp, value.var = "ES2")
        comp_all_abs_matrix[is.na(comp_all_abs_matrix)] <- 0
        ht <- ComplexHeatmap::Heatmap(comp_all_abs_matrix, name = "ES2", na_col = "#0F0D0E",
                                      column_title = paste0(title," ES2 - RPSEA_rand Zscore",
                                                            "\n(ORA padj<",as.character(ORA_adjP_thr[i]),
                                                            ", RPSEA padj<",as.character(RPSEA_adjP_thr[i]),
                                                            ", RPSEA_rand z-score>",as.character(randZscore_thr[i]),")"),
                                      col = palette1)
        
        ComplexHeatmap::draw(ht, padding = grid::unit(c(30, 2, 2, 2), "mm"))
      },
      error = function(e){
        print(e)
      },
      warning = function(w){
        print(w)
      },
      finally = {
        dev.off()
      }
    )
    
  }
  
  if(exists("last"))
    return(last)
  
  return(NULL)
}


#' Draw dripARF result scatterplot for multiple comparisons
#' @description Draw a scatterplot of RPSEA NES z-score vs ORA adjusted p-value for all RPs across
#'   comparisons, with top-predicted RPs highlighted.
#' @param dripARF_results Full dripARF result data.frame as returned by
#'   \code{dripARF_predict_heterogenity()}.
#' @param targetDir Directory to save the PDF plot in.
#' @param title Title string used for the plot and output filename (default:
#'   \code{"DRH_prediction_volcanos"}).
#' @param addedRPs Character vector of RP names to always highlight regardless of significance
#'   (optional).
#' @param highlightRPs Character vector of RP names to highlight instead of the top-predicted RPs
#'   (optional).
#' @param randZscore_thr Z-score threshold for the RPSEA NES relative to its randomised control
#'   sets. An RP passes this filter when its NES z-score is >= this value (default: 1).
#' @param ORA_adjP_thr Adjusted p-value threshold for the over-representation analysis (ORA). An RP
#'   passes this filter when ORA BH-adjusted p-value <= this value (default: 0.05).
#' @param RPSEA_adjP_thr Adjusted p-value threshold for the RPSEA enrichment score. An RP passes
#'   this filter when RPSEA BH-adjusted p-value <= this value (default: 0.05).
#' @param ORA_sig_n Minimum number of significantly differential rRNA positions (ORA overlap)
#'   required for an RP to be considered (default: 1).
#' @return A \code{ggplot2} object. A PDF is also saved to \code{targetDir}.
#' @keywords dripARF Differential Ribosome Heterogeneity Scatterplot
#' @export
#' @examples
#' \dontrun{
#' dripARF_result_scatterplot(dripARF_results, "/Folder/to/save/in/")
#' }
dripARF_result_scatterplot <- function(dripARF_results, targetDir,
                                       title="DRH_prediction_volcanos", addedRPs=NULL, highlightRPs=NULL,
                                       randZscore_thr=1, ORA_adjP_thr=0.05, RPSEA_adjP_thr=0.05, ORA_sig_n=1){
  
  simplified <- dripARF_simplify_results(dripARF_results = dripARF_results,
                                         randZscore_thr = randZscore_thr, ORA_adjP_thr = ORA_adjP_thr,
                                         RPSEA_adjP_thr = RPSEA_adjP_thr, ORA_sig_n = ORA_sig_n)
  
  if (is.null(highlightRPs))
    to_highlight <- simplified[simplified$top,]
  else
    to_highlight <- simplified[simplified$Description%in%highlightRPs,]
  
  if (length(addedRPs)>0)
    to_highlight <- rbind(to_highlight, simplified[simplified$Description %in% c(addedRPs), ])
  
  cols <- c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")
  
  g1 <- ggplot2::ggplot(simplified,
                        ggplot2::aes(x=ES2,y=ES1))+
    ggplot2::geom_hline(yintercept = 0)+
    ggplot2::geom_vline(xintercept = c(1), linetype="dashed", col=cols[5])+
    ggplot2::geom_vline(xintercept = 0, col=cols[5])+
    ggplot2::facet_grid(comp~.,scales = "free")+
    ggplot2::geom_point(ggplot2::aes(col=ORA,shape=ES1.pass),alpha=.5)+
    ggplot2::geom_point(data = to_highlight, size=.1, inherit.aes = TRUE, col="black")+
    ggrepel::geom_text_repel(data = to_highlight, force = 3,
                             inherit.aes = TRUE, col="black", ggplot2::aes(label=Description),max.overlaps = 100)+
    ggplot2::scale_color_manual(values = c(cols[3],cols[2],cols[5]))+
    ggplot2::scale_shape_manual(values = c(16,17),breaks = c(TRUE,FALSE))+
    ggplot2::theme_bw()+
    ggplot2::labs(col=paste0("ORA \n padj<",as.character(ORA_adjP_thr)), shape=paste0("RPSEA\n padj<",as.character(RPSEA_adjP_thr)))+
    ggplot2::xlab("Enrichment Score 2\n(NES->Z-score within random sets)")+
    ggplot2::ylab("Enrichment Score 1 (RPSEA NES)")
  print(g1)
  ggplot2::ggsave(plot = g1, filename = paste(targetDir, "/", title, ".pdf",sep=""),
                  width = 5,
                  height = 4+(ceiling(length(unique(dripARF_results$comp)))*3),limitsize = FALSE)
  return(g1)
}


#' Get per-position DESeq2 differential abundance results for rRNA position heatmaps
#' @description Extract per-position DESeq2 differential abundance results for use in rRNA
#'   position-specific heatmaps (\code{dripARF_rRNApos_heatmaps()}). Runs DESeq2 for each pairwise
#'   comparison and returns log2 fold-change and adjusted p-values at single-nucleotide resolution.
#' @param samples Samples dataframe created by \code{read_ARF_samples_file()}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism.
#' @param rRNA_counts Pre-computed rRNA count data.frame from \code{dripARF_read_rRNA_fragments()}
#'   (optional).
#' @param dripARF_dds DESeq2 \code{DESeqDataSet} from \code{dripARF_get_DESEQ_dds()} (optional).
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param compare Column name in the samples file to use as the grouping variable (default:
#'   \code{"group"}).
#' @param comparisons Named list of comparisons to include. If \code{NULL}, all pairwise
#'   comparisons are run.
#' @param exclude Character vector of sample names to exclude from the analysis.
#' @param gsea_sets_RP RP-rRNA contact point sets (preset for \code{hs}, \code{mm}, and \code{sc}).
#' @param RP_proximity_df RP-rRNA proximity matrix (preset for \code{hs}, \code{mm}, and \code{sc}).
#' @return A data.frame with per-position DESeq2 results across all comparisons, with columns:
#'   \code{baseMean}, \code{log2FoldChange}, \code{lfcSE}, \code{stat}, \code{pvalue}, \code{padj},
#'   \code{pos} (rRNA position ID), \code{comp} (comparison label).
#' @seealso \code{\link{dripARF_rRNApos_heatmaps}}
#' @keywords dripARF rRNA logFC DESeq2
#' @export
#' @examples
#' \dontrun{
#' dripARF_report_RPspec_pos_results(samples_df, "rRNAs.fa", organism="hs")
#' }
dripARF_report_RPspec_pos_results <- function(samples, rRNAs_fasta, rRNA_counts=NULL, dripARF_dds=NULL,
                                              organism=NULL, compare="group", comparisons=NULL, exclude=NULL, 
                                              gsea_sets_RP=NULL, RP_proximity_df=NULL){
  # # Check organism first
  # if (!ARF_check_organism(organism))
  #   return(NA)
  
  # Read samples
  if(!is.null(exclude))
    samples <- samples[!samples[,1]%in%exclude,]
  
  # GEt read counts
  if (is.null(dripARF_dds)){
    if (is.null(rRNA_counts)) {
      rRNA_counts <- dripARF_read_rRNA_fragments(samples = samples, rRNAs_fasta=rRNAs_fasta, organism=organism, QCplot = FALSE)
    }
    dds <- dripARF_get_DESEQ_dds(samples = samples, rRNAs_fasta=rRNAs_fasta, rRNA_counts = rRNA_counts, compare=compare, organism=organism, exclude=exclude)
  } else {
    dds <- dripARF_dds
  }
  
  s_n <- unique(samples[,compare])
  s_l <- length(s_n)
  samples$DESEQcondition <- samples[,compare]
  
  if(is.null(comparisons) || length(comparisons)==0) {
    comparisons <- list()
    for (i in 1:(s_l-1)) {
      for (j in (i+1):s_l) {
        comparisons[[(length(comparisons) +1)]] <- c(s_n[i],s_n[j])
      }
    }
  }
  
  all_results <- NULL
  for (comp in comparisons){
    message(paste("Comparing",comp[1],"vs",comp[2],"\n"))
    res <- data.frame(DESeq2::results(dds, contrast=c("DESEQcondition",comp[1],comp[2]), cooksCutoff = FALSE ))
    res$pos <- rownames(res)
    res$comp <- paste0(comp[1],"_vs_",comp[2])
    all_results <- rbind(all_results, res)
  }
  
  rownames(all_results) <- 1:(dim(all_results)[1])
  return(all_results)
}

#' Draw Differential rRNA abundance heatmap and its overlap with RP-rRNA proximity sets
#' @description Visualises position-specific differential rRNA fragment abundance across comparisons
#'   as a heatmap, overlaid with the rRNA contact-point positions of the requested RPs. Positions
#'   that pass both the log2FC and adjusted p-value thresholds are highlighted.
#' @param dripARF_DRF Per-position DESeq2 results data.frame as returned by
#'   \code{dripARF_report_RPspec_pos_results()}.
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param RPs Character vector of RP names to include in the figure (e.g. \code{c("eL1", "uS25")}).
#' @param targetDir Directory to save the PDF heatmap in.
#' @param abs_lFC_thr Absolute log2 fold-change threshold for highlighting differential positions
#'   (default: 0.5).
#' @param adjP_thr Adjusted p-value threshold for highlighting differential positions (default:
#'   0.05).
#' @param pval_thr Nominal p-value threshold (default: 0.05). Used when adjusted p-values are not
#'   available.
#' @param title Title string and base name for the output PDF (default:
#'   \code{"rRNA_pos_spec_heatmap"}).
#' @param gsea_sets_RP RP-rRNA contact point sets (preset for \code{hs}, \code{mm}, and \code{sc}).
#' @param RP_proximity_df RP-rRNA proximity matrix (preset for \code{hs}, \code{mm}, and \code{sc}).
#' @return A \code{ComplexHeatmap} object. A PDF is also saved to \code{<targetDir>/<title>.pdf}.
#' @keywords dripARF Differential Ribosome Heterogeneity Heatmap
#' @export
#' @examples
#' \dontrun{
#' dripARF_rRNApos_heatmaps(dripARF_results, "mm", c("eL1","uS25"), "/Folder/to/save/in/")
#' }
dripARF_rRNApos_heatmaps <- function(dripARF_DRF, organism, RPs, targetDir, 
                                     abs_lFC_thr=0.5, adjP_thr=0.05, pval_thr=0.05,
                                     title="rRNA_pos_spec_heatmap", 
                                     gsea_sets_RP=NULL, RP_proximity_df=NULL){
  if(is.null(RP_proximity_df)){
    if (organism=="hs"){
      RP_proximity_df <- ARF:::RP_proximity_human_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::human_gsea_sets_RP
    } else if (organism=="mm") {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::mouse_gsea_sets_RP
    } else if (organism=="sc") {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    } else {
      message(paste(c("Organism", organism, "Not implemented yet! Please generate your own set of proximity matrix and RP-rRNA sets using ARF.\n"), 
                collapse = " "))
      return(NULL)
    }
  }
  
  RPs_toreport <- unique(as.character(gsea_sets_RP$ont[!substring(gsea_sets_RP$ont,1,3)%in%c("MRf","FDf","Ran")]))
  
  prox_col = circlize::colorRamp2(c(0,26.9999,27, 50,200,300),c("black","black","grey30","grey60","grey99","white"))
  
  profileplot <- dripARF_DRF
  profileplot$sig <- (abs(profileplot$log2FoldChange)>abs_lFC_thr) & (profileplot$padj<adjP_thr) & (profileplot$pvalue<pval_thr) 
  profileplot$rrna <- gsub(pattern = "_[0-9]*$",replacement = "", x =profileplot$pos)
  profileplot$ipos <- as.numeric(gsub(pattern = ".*_([0-9]*$)",replacement = "\\1", x = profileplot$pos))
  profileplot$tomatrix <- abs(profileplot$log2FoldChange)*profileplot$sig
  profileplot$signedmatrix <- profileplot$log2FoldChange*profileplot$sig
  
  rrnas<-unique(RP_proximity_df$rRNA)
  focused_order <- c()
  for (rrna in unique(profileplot$rrna)){
    iorder <- order(profileplot$ipos)
    focused_order<-c(focused_order,unique(profileplot$pos[iorder[profileplot$rrna[iorder]==rrna]]))
  }
  
  temp <- RP_proximity_df[focused_order,]
  rids <- c()
  extended_rids <- c()
  for (RP in RPs){
    tmp_rids <- which(temp[,RP]<27.40514)
    tmp_extended_rids <- c()
    for (i in 1:length(tmp_rids)){
      rid <- tmp_rids[i]
      rRNA <- temp$rRNA[rid]
      toadd <- max(-30+rid,min(which(temp$rRNA==rRNA))):min(rid+30,max(which(temp$rRNA==rRNA)))
      tmp_extended_rids <- unique(c(tmp_extended_rids,toadd[(!toadd%in%tmp_rids) & temp$rRNA[toadd]==rRNA]))
    }
    rids<-c(rids,tmp_rids)
    extended_rids<-c(extended_rids,tmp_extended_rids)
  }
  val_pos <- sort(unique(c(rids,extended_rids)))
  tags <- c(val_pos[1])
  for (i in 2:length(val_pos)){
    if (((temp$resno[val_pos[i-1]]+1)==temp$resno[val_pos[i]]) & (temp$rRNA[val_pos[i-1]]==temp$rRNA[val_pos[i]]))
      tags<-c(tags,tags[i-1])
    else
      tags<-c(tags,val_pos[i])
  }
  
  exclude <- c()
  for (tag in unique(tags)){
    if (sum(val_pos[tags==tag]%in%rids)==0)
      exclude<-c(exclude,-1*which(tags==tag))
    else
      tags[tags==tag] <- paste0("[",temp$resno[min(val_pos[tags==tag])],", ",temp$resno[max(val_pos[tags==tag])],"]")
  }
  
  if(length(exclude)>0) {
    get_pos <- list(val_pos=val_pos[exclude],tags=tags[exclude])
  }  else{
    get_pos <- list(val_pos=val_pos,tags=tags)
  }
  
  
  whichPos <- get_pos[[1]]
  tags <- get_pos[[2]]
  
  foo <- NULL
  for (RP in RPs){
    foo<-cbind(foo, RP_proximity_df[focused_order[whichPos],RP])
  }
  colnames(foo)<-RPs
  
  ha <- ComplexHeatmap::HeatmapAnnotation(rRNA = RP_proximity_df[focused_order,"rRNA"][whichPos],
                                          foo=foo, col=list(foo=prox_col,
                                                            rRNA=setNames(wesanderson::wes_palette("Rushmore1")[c(1,3,4,5)],rrnas), na_col = "white"),
                                          annotation_legend_param = list(foo=list(title="RP-rRNA\nproximity map\ndistance (Angstrom)",at=c(0,27,100,200,300))))
  
  comparisons <- unique(as.character(dripARF_DRF$comp))
  
  if (length(unique(profileplot$comp))>1) {
    matrixplot <- reshape2::acast(profileplot, formula=comp ~ pos, value.var = "tomatrix")[,unique(profileplot$pos)]
  } else {
    matrixplot <- reshape2::acast(profileplot, formula=comp ~ pos, value.var = "tomatrix")
  }
  
  pdf(paste0(targetDir,"/",title,".pdf"),height = 5+round((length(RPs)+length(comparisons))/4), width = 11+round(length(RPs)/10))
  if (length(unique(profileplot$comp))>1) {
    ht <- ComplexHeatmap::Heatmap(matrixplot[,focused_order[whichPos]], name = "abs(log2FC)", top_annotation = ha,
                                  cluster_columns = FALSE,show_column_names = FALSE,
                                  col=circlize::colorRamp2(c(0,abs_lFC_thr,8),c("white","cornsilk","black")),na_col = "white",
                                  column_split = factor(tags,levels = unique(tags),labels = unique(tags)),
                                  column_title_rot = 90, heatmap_legend_param = list(at=c(abs_lFC_thr,2,4,6,8)))
    ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 16, 2), "mm"))
  } else {
    tp <- t(matrixplot[,focused_order[whichPos]])
    rownames(tp)<-unique(profileplot$comp)
    ht <- ComplexHeatmap::Heatmap(tp, name = "abs(log2FC)", top_annotation = ha,
                                  cluster_columns = FALSE, show_column_names = FALSE,
                                  col=circlize::colorRamp2(c(0,abs_lFC_thr,8),c("white","cornsilk","black")),na_col = "white",
                                  column_split = factor(tags,levels = unique(tags),labels = unique(tags)),
                                  column_title_rot = 90, heatmap_legend_param = list(at=c(abs_lFC_thr,2,4,6,8)))
    ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 16, 2), "mm"))
  }
  
  dev.off()
  return(ht)
}

#' Run the complete dripARF heterogeneity prediction pipeline
#' @description Convenience wrapper that runs the full dripARF pipeline in a single call: reads
#'   sample metadata, loads rRNA alignments, normalises counts with DESeq2, performs RPSEA and ORA
#'   enrichment analyses, and saves result CSVs, a scatterplot PDF, and heatmap PDFs to
#'   \code{targetDir}.
#' @param samplesFile Path to the tab-separated samples file. Required columns: \code{sampleName},
#'   \code{bedGraphFile} (or BAM path), \code{group}. See \code{read_ARF_samples_file()} for
#'   details.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism.
#' @param samples_df Pre-loaded samples data.frame from \code{read_ARF_samples_file()} (optional;
#'   read from \code{samplesFile} if \code{NULL}).
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param compare Column name in the samples file to use as the grouping variable (default:
#'   \code{"group"}).
#' @param QCplot \code{TRUE} or \code{FALSE}, whether to generate QC plots (default: \code{TRUE}).
#' @param targetDir Directory to save all result files (CSVs, PDFs) in.
#' @param comparisons Named list of comparisons to include. If \code{NULL}, all pairwise
#'   comparisons are run.
#' @param exclude Character vector of sample names to exclude from the analysis.
#' @param GSEAplots \code{TRUE} or \code{FALSE}, whether to produce per-RP GSEA enrichment plots
#'   (default: \code{FALSE}).
#' @param gsea_sets_RP RP-rRNA contact point sets (preset for \code{hs}, \code{mm}, and \code{sc};
#'   provide for custom organisms via \code{dripARF_get_RP_proximity_sets()}).
#' @param RP_proximity_df RP-rRNA proximity matrix (preset for \code{hs}, \code{mm}, and \code{sc};
#'   provide for custom organisms via \code{ARF_parse_PDB_ribosome()}).
#' @return A data.frame -- the same as returned by \code{dripARF_predict_heterogenity()}. Scatter
#'   plot and heatmap PDFs are also saved to \code{targetDir}.
#' @seealso \code{\link{dricARF}}, \code{\link{read_ARF_samples_file}},
#'   \code{\link{dripARF_result_scatterplot}}, \code{\link{dripARF_result_heatmap}}
#' @keywords dripARF pipeline
#' @export
#' @examples
#' \dontrun{
#' dripARF("samples.txt", "rRNAs.fa", organism="mm", targetDir="/target/directory/to/save/results")
#' }
dripARF <- function(samplesFile, rRNAs_fasta, samples_df=NULL, organism=NULL, compare="group", QCplot=TRUE,  targetDir=NA,
                    comparisons=NULL, exclude=NULL, ssRPSEAplots=FALSE, GSEAplots=FALSE, gsea_sets_RP=NULL, RP_proximity_df=NULL){
  
  # # Check organism first
  # if (!ARF_check_organism(organism))
  #   return(NA)
  if(is.null(RP_proximity_df)){
    if (organism=="hs"){
      RP_proximity_df <- ARF:::RP_proximity_human_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::human_gsea_sets_RP
    } else if (organism=="mm") {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::mouse_gsea_sets_RP
    } else if (organism=="sc") {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
      if (is.null(gsea_sets_RP))
        gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    } else {
      message(paste(c("Organism", organism, "Not implemented yet! Please generate your own set of proximity matrix and RP-rRNA sets using ARF.\n"), 
                collapse = " "))
      return(NULL)
    }
  }
  
  # Assign target directory
  if (is.na(targetDir)){
    targetDir=getwd()
  }
  
  if(is.null(samples_df))
    samples_df <- read_ARF_samples_file(samplesFile)
  
  if(!is.null(exclude))
    samples_df <- samples_df[!samples_df[,1]%in%exclude,]
  
  rRNA_counts_df <- dripARF_read_rRNA_fragments(samples = samples_df, rRNAs_fasta=rRNAs_fasta, organism = organism, 
                                                QCplot = QCplot, targetDir = targetDir)
  dripARF_results <- dripARF_predict_heterogenity(samples = samples_df, rRNAs_fasta=rRNAs_fasta, rRNA_counts = rRNA_counts_df,
                                                  compare=compare, organism=organism, QCplot=QCplot, targetDir=targetDir,
                                                  comparisons = comparisons, GSEAplots=GSEAplots, ssRPSEAplots=ssRPSEAplots,
                                                  gsea_sets_RP=gsea_sets_RP, RP_proximity_df=RP_proximity_df)
  
  dripARF_result_scatterplot(dripARF_results = dripARF_results, targetDir = targetDir, title = "ALL dripARF predictions")
  dripARF_result_heatmap(dripARF_results = dripARF_results, targetDir = targetDir, title = "ALL dripARF predictions",
                         randZscore_thr = c(1), ORA_adjP_thr = c(0.05), RPSEA_adjP_thr = c(0.05), ORA_sig_n = 1)
  
  return(dripARF_results)
}



#' Generate RP-rRNA proximity sets for RPSEA
#' @description Generates GSEA-format term-to-gene sets for each ribosomal protein (RP) based on
#'   rRNA proximity. When \code{thresholds=NULL}, the distance threshold is auto-set to the 5th
#'   percentile of all RP-rRNA distances, and the set-size cap is 5\% of total rRNA positions. For
#'   each RP, 99 randomised control sets are also generated by circularly shifting the proximity
#'   positions, which are used downstream by \code{dripARF_predict_heterogenity()} to compute the
#'   RPSEA NES z-score.
#' @param RP_proximity_df RP-rRNA proximity matrix as returned by \code{ARF_parse_PDB_ribosome()} or
#'   \code{ARF_convert_Ribo3D_pos()}.
#' @param additional_RPcols Integer vector of column indexes for any custom RP columns added to
#'   \code{RP_proximity_df} beyond the standard ARF columns (needed for proper threshold and set
#'   generation, default: \code{c()}).
#' @param rRNAs_fasta FASTA file for the rRNAs of the target organism. When provided, used to
#'   calculate the set-size cap from total sequence length and to restrict positions to rRNAs
#'   present in the file (optional).
#' @param thresholds Numeric vector of length 2: \code{c(distance_threshold, max_set_size)}. If
#'   \code{NULL}, thresholds are calculated automatically (5th percentile of distances; 5\% of
#'   total rRNA length).
#' @param cap_added_RPcols Logical. Whether to also cap the set size for columns specified in
#'   \code{additional_RPcols} (default: \code{FALSE}).
#' @return A data.frame with two columns (\code{ont}, \code{gene}) suitable for use as GSEA
#'   term-to-gene sets. Includes one set per RP (rows where \code{ont} = RP name) and 99
#'   randomised control sets per RP (rows where \code{ont} = \code{RandN_RP}).
#' @seealso \code{\link{ARF_parse_PDB_ribosome}}, \code{\link{ARF_convert_Ribo3D_pos}},
#'   \code{\link{dripARF_predict_heterogenity}}
#' @keywords RPSEA RP rRNA proximity set generation
#' @export
#' @examples
#' \dontrun{
#' dripARF_get_RP_proximity_sets(RP_proximities_dataframe, rRNAs_fasta="rRNAs.fa")
#' dripARF_get_RP_proximity_sets(RP_proximities_dataframe,
#'   additional_RPcols=80:95, rRNAs_fasta="rRNAs.fa")
#' }
dripARF_get_RP_proximity_sets <- function(RP_proximity_df, additional_RPcols=c(), rRNAs_fasta=NULL, thresholds=NULL, cap_added_RPcols=F){
  gsea_sets_RP <- NULL
  RPfocus <- colnames(RP_proximity_df)[c(-1,-2)]

  if(is.null(thresholds)){
    excludedCols <- c(-1,-2)
    if(length(additional_RPcols) > 0){
      excludedCols <- c(excludedCols,-additional_RPcols)
    }
    
    dist_thr <- quantile(RP_proximity_df[,excludedCols],0.05,na.rm = T)
    RPset_size_thr <- round(dim(RP_proximity_df)[1]*0.05)
    if(!is.null(rRNAs_fasta)){
      seq_set <- Biostrings::readBStringSet(filepath = rRNAs_fasta)
      RPset_size_thr <- sum(Biostrings::width(seq_set))*0.05
      RP_proximity_df <- RP_proximity_df[RP_proximity_df$rRNA%in%sapply(strsplit(names(seq_set),split = " "),"[",1),]
    }  
    thresholds <- c(dist_thr,RPset_size_thr)
  }
  
  # Additional sets is a list of row numbers 
  RP_proxpos <- list()
  for (RP in RPfocus){
    proxpos <- which(RP_proximity_df[,RP]<thresholds[1])
    if (length(proxpos)>thresholds[2] && ((RP %in% (colnames(RP_proximity_df)[excludedCols])) || cap_added_RPcols) ){
      RP_proxpos[[RP]] <- sort(proxpos[order(RP_proximity_df[proxpos,RP])[1:thresholds[2]]])
    } else if (length(proxpos)>0){
      RP_proxpos[[RP]] <- proxpos
    }
  }
  
  gsea_sets_RP <- do.call("rbind", lapply(names(RP_proxpos), FUN = function(RP){
    tmp_df<-NULL
    proxpos <- RP_proxpos[[RP]]
    
    tmp_df <- data.frame(ont=RP,gene=paste(RP_proximity_df$rRNA[proxpos], RP_proximity_df$resno[proxpos], sep = "_"))
    rands <- c((1:100)*(round(dim(RP_proximity_df)[1]/100,digits = 0)-1))
    
    tmp_df <- rbind(tmp_df, as.data.frame(do.call("rbind", lapply(1:99, FUN = function(i){
      randset <- ((proxpos+rands[i]) %% (dim(RP_proximity_df)[1]))+1
      return(data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
                        gene=paste(RP_proximity_df$rRNA[randset], RP_proximity_df$resno[randset], sep = "_")))
    }))))
    return(tmp_df)
  }))
  
  return(gsea_sets_RP)
}


#' Test dripARF predictions across varying RP-rRNA proximity thresholds
#' @description Runs the dripARF heterogeneity prediction pipeline repeatedly for each supplied
#'   threshold combination, allowing assessment of result stability across different distance and
#'   set-size cutoffs. The output combines all runs with a \code{threshold} identifier column.
#' @param samplesFile Path to the tab-separated samples file. See \code{read_ARF_samples_file()}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism.
#' @param organism Organism abbreviation. Pass \code{"hs"} for human, \code{"mm"} for mouse, and
#'   \code{"sc"} for yeast.
#' @param compare Column name in the samples file to use as the grouping variable (default:
#'   \code{"group"}).
#' @param comparisons Named list of comparisons to include. If \code{NULL}, all pairwise
#'   comparisons are run.
#' @param exclude Character vector of sample names to exclude from the analysis.
#' @param thresholds List of numeric vectors, each defining one threshold combination. Each vector
#'   contains: 1st value = Angstrom distance threshold, 2nd value = maximum set-size threshold,
#'   optional 3rd value = directionality flag (\code{TRUE}/\code{FALSE}, default: \code{FALSE}).
#'   If \code{NULL}, auto-calculated thresholds are used (see \code{dripARF_get_RP_proximity_sets()}).
#' @param RP_proximity_df RP-rRNA proximity matrix (preset for \code{hs}, \code{mm}, and \code{sc};
#'   provide for custom organisms via \code{ARF_parse_PDB_ribosome()}).
#' @return A data.frame combining \code{dripARF_predict_heterogenity()} results across all supplied
#'   thresholds, with an additional \code{threshold} column identifying which threshold combination
#'   was used.
#' @keywords dripARF threshold alternative
#' @export
#' @examples
#' \dontrun{
#' dripARF_threshold_test("samples.txt", "rRNAs.fa", organism="mm",
#'   thresholds=list(c(20,150),c(15,100)))
#' }
dripARF_threshold_test <- function(samplesFile, rRNAs_fasta, 
                                   organism=NULL, compare="group", comparisons=NULL, exclude=NULL, thresholds=NULL,
                                   RP_proximity_df=NULL){
  
  # # Check organism first
  # if (!ARF_check_organism(organism))
  #   return(NA)
  
  if(is.null(RP_proximity_df)){
    if (organism=="hs"){
      RP_proximity_df <- ARF:::RP_proximity_human_df
    } else if (organism=="mm") {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
    } else if (organism=="sc") {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
    } else {
      message(paste(c("Organism", organism, "Not implemented yet! Please generate your own set of proximity matrix and RP-rRNA sets using ARF.\n"), 
                collapse = " "))
      return(NULL)
    }
  }
  
  samples <- read_ARF_samples_file(samplesFile)
  if(!is.null(exclude))
    samples <- samples[!samples[,1]%in%exclude,]
  
  rRNA_counts <- dripARF_read_rRNA_fragments(samples = samples, rRNAs_fasta=rRNAs_fasta, organism = organism, QCplot = FALSE)
  
  dds <- dripARF_get_DESEQ_dds(samples = samples, rRNAs_fasta=rRNAs_fasta, rRNA_counts = rRNA_counts, compare=compare, organism=organism, exclude=exclude)
  
  s_n <- unique(samples[,compare])
  s_l <- length(s_n)
  samples$DESEQcondition <- samples[,compare]
  
  # Read count transformations
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
  
  if(is.null(comparisons) || length(comparisons)==0) {
    comparisons <- list()
    for (i in 1:(s_l-1)) {
      for (j in (i+1):s_l) {
        comparisons[[(length(comparisons) +1)]] <- c(s_n[i],s_n[j])
      }
    }
  }
  
  dripARF_results_collected <- NULL
  
  
  for (thr in thresholds){
    ## first get the RPsets based on the threshold
    message(paste("Running dripARF for threshold", thr[1], thr[2], "\n"))
    
    directional = FALSE
    if (length(thr)>2) {
      if (thr[3]) {
        directional = TRUE
        message("Directional Test!\n")
      }
    }
    
    gsea_sets_RP <- NULL
    RPfocus <- colnames(RP_proximity_df)[c(-1,-2)]
    # RPfocus <- RPfocus[-1*which(startsWith(x = RPfocus,"ES") | startsWith(x = RPfocus,"YCI"))]
    
    RP_proxpos <- list()
    for (RP in RPfocus){
      proxpos <- which(RP_proximity_df[,RP]<thr[1])
      if (length(proxpos)>thr[2] & (!startsWith(RP,"ES")) & (!startsWith(RP,"YCI")) & (!startsWith(RP,"yeast"))){
        RP_proxpos[[RP]] <- sort(proxpos[order(RP_proximity_df[proxpos,RP])[1:thr[2]]])
      } else if (length(proxpos)>0){
        RP_proxpos[[RP]] <- proxpos
      }
    }
    gsea_sets_RP <- do.call("rbind", lapply(names(RP_proxpos), FUN = function(RP){
      tmp_df<-NULL
      proxpos <- RP_proxpos[[RP]]
      
      tmp_df <- data.frame(ont=RP,gene=paste(RP_proximity_df$rRNA[proxpos], RP_proximity_df$resno[proxpos], sep = "_"))
      rands <- c((1:100)*(round(dim(RP_proximity_df)[1]/100,digits = 0)-1))
      
      tmp_df <- rbind(tmp_df, as.data.frame(do.call("rbind", lapply(1:99, FUN = function(i){
        randset <- ((proxpos+rands[i]) %% (dim(RP_proximity_df)[1]))+1
        return(data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
                          gene=paste(RP_proximity_df$rRNA[randset], RP_proximity_df$resno[randset], sep = "_")))
      }))))
      return(tmp_df)
    }))
    
    ######### Gene set enrichment Analysis ############
    dripARF_results <- NULL
    ###################################################
    RPs_toreport <- unique(as.character(gsea_sets_RP$ont[!substring(gsea_sets_RP$ont,1,3)%in%c("MRf","FDf","Ran")]))
    ########## Overrepresentation Analysis ############
    RP_pathways <- sapply(RPs_toreport,FUN=function(x){return(as.character(gsea_sets_RP$gene[gsea_sets_RP$ont==x]))})
    
    # Separate for each DESEQcondition
    for (comp in comparisons) {
      message(paste("Running predictions for",comp[1],"vs",comp[2],"\n"))
      
      GSEA_result_df <- NULL
      res <- DESeq2::results(dds, contrast=c("DESEQcondition",comp[1],comp[2]), cooksCutoff = FALSE)
      
      temp_df <- res
      temp_df$weight <- scales::rescale(log10(temp_df$baseMean), to = c(0, 5))
      temp_df$padj[temp_df$padj<0.00001] = 0.00001
      
      measureID = "abs_GSEA_measure" #)){ #},"GSEA_measure","w_GSEA_m", "abs_w_GSEA_m")){
      if (directional)
        measureID = "GSEA_measure" 
      
      used_measure <- NULL
      if (measureID=="GSEA_measure"){
        used_measure <- res$log2FoldChange*(-log10(temp_df$padj))
      } else if (measureID=="abs_GSEA_measure"){
        used_measure <- abs(res$log2FoldChange)*(-log10(temp_df$padj))
      } else if (measureID=="w_GSEA_m"){
        used_measure <- res$log2FoldChange*(-log10(temp_df$padj))*temp_df$weight
      } else if (measureID=="abs_w_GSEA_m"){
        used_measure <- abs(res$log2FoldChange)*(-log10(temp_df$padj))*temp_df$weight
      }
      
      names(used_measure)<-rownames(res)
      used_measure <- used_measure[!sapply(used_measure, function(x) is.na(x))]
      used_geneList <- used_measure[order(used_measure, decreasing = TRUE)]
      used_geneList_abs <- abs(used_measure)[order(abs(used_measure), decreasing = TRUE)]
      
      if(directional) {
        egmt_used_measure <- clusterProfiler::GSEA(geneList = used_geneList, TERM2GENE=gsea_sets_RP, verbose=TRUE, minGSSize = 10, maxGSSize = 10000, pvalueCutoff = 2)
      } else {
        egmt_used_measure <- clusterProfiler::GSEA(geneList = used_geneList, TERM2GENE=gsea_sets_RP, verbose=TRUE, minGSSize = 10, maxGSSize = 10000, pvalueCutoff = 2, scoreType = "pos")
      }
      egmt_used_measure@result$NES_rand_zscore <- NA
      for (RP in RPs_toreport){
        tochange <- endsWith(x = egmt_used_measure@result$ID, suffix = RP)
        egmt_used_measure@result$NES_rand_zscore[tochange] <- scale(egmt_used_measure@result$NES[tochange])
      }
      
      ### Overrepresentation hook ####
      or_df <- fgsea::fora(pathways = RP_pathways,genes = rownames(res)[which(res$padj<.05 & abs(res$log2FoldChange)>0.5)],
                           universe = rownames(res), minSize = 10)
      
      GSEA_result_df <- data.frame(Description = or_df$pathway,
                                   ORA.overlap=or_df$overlap, ORA.setSize=or_df$size, ORA.padj=or_df$padj, ORA.p=or_df$pval,
                                   RPSEA.NES=NA, RPSEA.NES_randZ=NA, RPSEA.padj=NA, RPSEA.pval=NA, RPSEA.q=NA)
      rownames(GSEA_result_df) <- GSEA_result_df$Description
      
      
      if(dim(egmt_used_measure@result)[1]>0){
        # Assign RPSEA scores, padjust for only valid reported sets
        GSEA_result_df$RPSEA.NES <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"NES"]
        GSEA_result_df$RPSEA.NES_randZ <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"NES_rand_zscore"]
        GSEA_result_df$RPSEA.padj <- p.adjust(p = (egmt_used_measure@result[as.character(GSEA_result_df$Description),"pvalue"]), method = "BH")
        GSEA_result_df$RPSEA.pval <- (egmt_used_measure@result[as.character(GSEA_result_df$Description),"pvalue"])
        
        #GSEA_result_df[RP,"RPSEA.q"] <- temp$qvalues
        if("qvalue" %in% colnames(egmt_used_measure@result)){
          GSEA_result_df$RPSEA.q <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"qvalue"]
        } else if("qvalues" %in% colnames(egmt_used_measure@result)){
          GSEA_result_df$RPSEA.q <- egmt_used_measure@result[as.character(GSEA_result_df$Description),"qvalues"]
        } else{
          print("qvalue is not part of the egmt_result!!")
          GSEA_result_df[RP,"RPSEA.q"] <- NULL
        }
      }
      
      dripARF_results <- rbind(dripARF_results, data.frame(comp=paste(comp, collapse = "_vs_"),
                                                           GSEA_result_df[order(GSEA_result_df$RPSEA.NES,decreasing = TRUE),]))
    }
    
    dripARF_results_collected <- rbind(dripARF_results_collected, data.frame(prox_threshold=thr[1], RPset_size_threshold=thr[2], directional=directional,  
                                                                             dripARF_results))
    
    gc(verbose = T)
  }
  
  return(dripARF_results_collected)
}


visualize_geneset <- function(organism, chain_file, RP) {
  RP_proximity_df <- NULL
  if (organism=="hs"){
    RP_proximity_df <- ARF:::RP_proximity_human_df
    gsea_sets_RP <- ARF:::human_gsea_sets_RP
    chainData_raw <- read.csv(chain_file)
  } else if (organism=="sc") {
    RP_proximity_df <- ARF:::RP_proximity_yeast_df
    gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    chainData_raw <- read.csv(chain_file, header = F)
    colnames(chainData_raw) <- c("chain_lead","chain_trail","name","RP","type")
    chainData_raw$Chain_id <- sapply(1:dim(chainData_raw)[1], FUN = function(i){
      if (chainData_raw$chain_lead[i]!=chainData_raw$chain_trail[i]){
        return(paste(chainData_raw$chain_lead[i],chainData_raw$chain_trail[i],sep = "+"))
      } else {return(chainData_raw$chain_lead[i])}
    })
  } else {
    message(paste(c("Organism", organism, "Not possible structure!\n"), collapse = " "))
    return(NULL)
  }
}


#' Add pseudo-replicates to single-replicate sample groups
#' @description Creates pseudo-replicates for sample groups that contain only a single sample.
#'   DESeq2 requires at least two replicates per group to estimate dispersion; this function adds a
#'   pseudo-replicate by multiplying each position's count by a per-position Uniform(0.9, 1.1)
#'   random factor and rounding to the nearest integer. The new pseudo-replicate is appended to both
#'   the \code{samples} data.frame (with suffix \code{"_ADDED"}) and the \code{rRNA_counts}
#'   data.frame. \strong{WARNING}: this is a statistical workaround only. Biological replicates are
#'   strongly preferred and should always be used when available.
#' @param samples Samples dataframe created by \code{read_ARF_samples_file()}.
#' @param rRNA_counts rRNA count data.frame from \code{dripARF_read_rRNA_fragments()}.
#' @param QCplot \code{TRUE} or \code{FALSE}, whether to generate QC plots (default: \code{FALSE}).
#' @param targetDir Directory to save QC plots in (default: \code{getwd()}).
#' @return A list with two elements: \code{samples} (the updated samples data.frame with
#'   pseudo-replicate rows appended) and \code{rRNA_counts} (the updated counts data.frame with
#'   pseudo-replicate columns appended).
#' @keywords single sample replicate creation
#' @export
#' @examples
#' \dontrun{
#' dripARF_add_replicates(samples_df, rRNA_counts_df)
#' dripARF_add_replicates(samples_df, rRNA_counts_df, QCplot=TRUE, targetDir="./")
#' }
dripARF_add_replicates <- function(samples, rRNA_counts, QCplot=FALSE, targetDir=NA) {
  group_counts <- table(samples$group)
  for (sample in samples$sampleName[samples$group%in%names(group_counts)[group_counts==1]]){
    new_sampleName <- paste0(sample,"_ADDED")
    samples[new_sampleName,] <- data.frame(sampleName=new_sampleName, samples[samples$sampleName==sample,-1])
    samples[new_sampleName,2] <- NA
    rRNA_counts[,new_sampleName] <- round(rRNA_counts[,sample] * runif(dim(rRNA_counts)[1], min = 0.9, max = 1.1),0)
  }
  return(list(samples=samples,rRNA_counts=rRNA_counts))
}