

#' Read the input tsv file
#' @description This function allows you to read the samples data, returning a dataframe.
#' @param samplesFile "Samples.tsv" tab-seperated file that describes file locations and sample groupings.
#' @keywords Reader Samples
#' @export
#' @examples
#' read_RP4_samples_file("samples.txt")
read_RP4_samples_file <- function(samplesFile){
  samples_df <- read.csv(file = samplesFile, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  rownames(samples_df) <- samples_df$sampleName
  return(samples_df)
}

#' Organism Check
#' @description Check if the organism is valid.
#' @param organism Organism abbrevation. Pass "hs" for human and "mm" for mouse.
RP4_check_organism <- function(organism) {
  if (!(organism %in% c("hs","mm"))) {
    print("Choose a valid organism. hs(human), mm(mouse), etc.")
    return(FALSE)
  }
  return(TRUE)
}

#' Read rRNA quantification from bedGraph (maybe .bam as well) files
#' @description Read bedgraph files to create the rRNA count data.
#' @param samples Samples dataframe created by read_samples_file() function.
#' @param organism Organism abbrevation. Pass "hs" for human and "mm" for mouse.
#' @param QCplot TRUE or FALSE, whether to generate QC plots or not.
#' @param targetDir Direcory to save QC plots in.
#' @keywords Reader rRNA fragment alignment
#' @export
#' @examples
#' RP4_read_rRNA_fragments(samples_df, "hs")
RP4_read_rRNA_fragments <- function(samples, organism="hs", QCplot=FALSE, targetDir=NA) {

  # Check organism first
  if (!RP4_check_organism(organism))
    return(NA)

  # Assign target directory
  if (is.na(targetDir)){
    targetDir=getwd()
  }

  df <- NULL
  if (organism == "hs") {
    df <- setNames(data.frame(matrix(ncol = length(samples$sampleName), nrow = (5070+1869+157+121),data=0)), samples$sampleName)
    rownames(df) <- c(paste("28S_N5_NR_003287.4",0:5069,sep = "_"), paste("18S_N5_NR_003286.4",0:1868,sep = "_"),
                    paste("5.8S_N5_NR_003285.3",0:156,sep = "_"), paste("human_5S",0:120,sep = "_"))
  } else if (organism=="mm") {
    df <- setNames(data.frame(matrix(ncol = length(samples$sampleName), nrow = (4730+1870+157+121),data=0)), samples$sampleName)
    rownames(df) <- c(paste("NR_003279.1",0:4729,sep = "_"), paste("NR_003278.3",0:1869,sep = "_"),
                      paste("NR_003280.2",0:156,sep = "_"), paste("NR_030686.1",0:120,sep = "_"))
  } else {
    print(paste(c("Organism", organism, "Not implemented yet!"), collapse = " "))
    return(NULL)
  }

  for (i in 1:dim(samples)[1]) {
    sample <- samples$sampleName[i]
    temp <- read.table(samples$bedGraphFile[i], stringsAsFactors = FALSE)
    print(paste("Reading bedgraph file for sample",sample,"(",as.character(i),"of",as.character(dim(samples)[1]),")"))
    for (i in 1:dim(temp)[1]){
      if (i%%1000==0)
        print(paste(as.character(round(100*i/(dim(temp)[1]),2)),"% is done..."))
      all <- temp$V2[i]:(temp$V3[i]-1)
      df[paste(temp$V1[i], all, sep = "_"), sample] = temp$V4[i]
    }
  }

  print("ALL bedgraphs have been read.")
  if (QCplot) {
    g1 <- ggplot2::ggplot(reshape2::melt(df), ggplot2::aes(x=variable, y=value)) +
      ggplot2::geom_violin()+ggplot2::scale_y_log10()+ggplot2::coord_flip()+
      ggplot2::xlab("Positional rRNA fragment counts (log-scale)")+ggplot2::ylab("Samples")
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "raw_rRNA_counts.png", sep = ""), plot = g1, width = 8, height = 3+length(samples$sampleName))
  }

  return(df[rowSums(is.na(df))==0,])
}


#' Predict Diff.Rib.Het. candidate RPs with rRNA count data
#' @description Overlap differential rRNA count data with 3d ribosome, rRNA-RP proximity data and predict heterogeneity candidates across groups.
#' @param samples Samples dataframe created by read_samples_file() function.
#' @param organism Organism abbrevation. Pass "hs" for human and "mm" for mouse.
#' @param QCplot TRUE or FALSE, whether to generate QC plots or not.
#' @param targetDir Direcory to save QC plots in.
#' @keywords Diffferential Ribosome Heterogeneity rRNA ribosome RP
#' @export
#' @examples
#' RP4_predict_heterogenity(samples_df, rRNA_counts_df, organism="hs", QCplot=TRUE)
RP4_predict_heterogenity <- function(samples, rRNA_counts, compare="group", organism="hs", QCplot=FALSE, targetDir=NA) {
  abs_zthr <- 1
  zthr <- 1.5

  s_n <- unique(samples[,compare])
  s_l <- length(s_n)

  if (!RP4_check_organism(organism))
    return(NULL)

  # Assign target directory
  if (is.na(targetDir)){
    targetDir=getwd()
  }

  samples$DESEQcondition <- samples[,compare]
  cts <- as.matrix(round(rRNA_counts, digits = 0))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = samples, design = ~DESEQcondition)

  keep <- rowSums(DESeq2::counts(dds)) >= 1000*dim(samples)[1]
  dds <- dds[keep,]
  dds <- DESeq2::DESeq(dds)

  comparisons <- list()
  for (i in 1:(s_l-1)) {
    for (j in (i+1):s_l) {
      print(paste("Comparing",s_n[i],"vs",s_n[j]))
      comparisons[[(length(comparisons) +1)]] <- c(s_n[i],s_n[j])
      res <- DESeq2::results(dds, contrast=c("DESEQcondition",s_n[i],s_n[j]), lfcThreshold = 0.5, alpha = 0.01)
      print(DESeq2::summary(res))
    }
  }

  # Read count transformations
  vsd <- DESeq2::vst(dds, blind=FALSE)

  if (QCplot) {
    g1 <- ggplot2::ggplot(reshape2::melt(SummarizedExperiment::assay(vsd)), ggplot2::aes(x=Var2, y=value)) +
      ggplot2::geom_violin()+ggplot2::coord_flip()+
      ggplot2::xlab("Positional rRNA fragment counts (vsd) after DESeq2 normalization")+ggplot2::ylab("Samples")
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "norm_rRNA_counts.png", sep = ""), plot = g1, width = 8, height = 3+length(samples$sampleName))

    #PCA_analysis
    pca <- prcomp(SummarizedExperiment::assay(vsd)[order(matrixStats::rowVars(SummarizedExperiment::assay(vsd)),decreasing = TRUE),])
    pca_df <- data.frame(pca$rotation)
    pca_map <- paste(colnames(pca_df),' (Var.Exp.:%',round(100*((pca$sdev^2)/sum(pca$sdev^2)),1),')', sep='')
    pca_df$sample <- rownames(pca_df)
    pca_df$group <- sapply(pca_df$sample,FUN=function(x){return(samples$DESEQcondition[samples$sampleName==x])})

    g1 <- ggplot2::ggplot(pca_df, ggplot2::aes(x=PC1,y=PC2,color=group)) + ggplot2::geom_point() +
      ggplot2::labs(x=pca_map[1], y=pca_map[2]) + ggrepel::geom_text_repel(ggplot2::aes(label=sample), show.legend = FALSE)
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "vsd_based_PCA.png", sep = ""), plot = g1,
                    width = 4+(length(samples$sampleName)/3), height = 4+(length(samples$sampleName)/3))
  }

  ###################################################
  org_RP_df <- NULL
  gsea_sets_RP <- NULL
  if (organism=="hs"){
    org_RP_df <- RP4:::RP_proximity_human_df
    gsea_sets_RP <- RP4:::human_gsea_sets_RP
  } else if (organism=="mm") {
    org_RP_df <- RP4:::RP_proximity_mouse_df
    gsea_sets_RP <- RP4:::mouse_gsea_sets_RP
  } else {
    print(paste(c("Organism", organism, "Not implemented yet!"), collapse = " "))
    return(NULL)
  }

  ######### Gene set enrichment Analysis ############

  #library(clusterProfiler)
  #library(enrichplot)
  #library(gprofiler2)

  all_GSEA_results <- NULL
  # Separate for each DESEQcondition
  for (comp in comparisons) {
    print(paste("Running predictions for",comp[1],"vs",comp[2]))

    GSEA_result_df <- NULL
    gset_tag <- "RP"

    # signed GSEA measure
    res <- DESeq2::results(dds, contrast=c("DESEQcondition",comp[1],comp[2]))

    temp_df <- res
    temp_df$padj[temp_df$padj<0.00001] = 0.00001
    GSEA_measure <- res$log2FoldChange*(-log10(temp_df$padj))

    names(GSEA_measure)<-rownames(res)
    GSEA_measure <- GSEA_measure[!sapply(GSEA_measure, function(x) is.na(x))]
    #summary(GSEA_measure)

    GSEA_geneList <- GSEA_measure[order(GSEA_measure, decreasing = TRUE)]
    #summary(GSEA_geneList)
    GSEA_geneList_abs <- abs(GSEA_measure)[order(abs(GSEA_measure), decreasing = TRUE)]
    #summary(GSEA_geneList_abs)

    egmt_GSEA_measure <- clusterProfiler::GSEA(geneList = GSEA_geneList, TERM2GENE=gsea_sets_RP, verbose=TRUE, pvalueCutoff = 1)
    egmt_GSEA_measure@result$NES_rand_zscore <- NA
    for (RP in colnames(org_RP_df)[c(-1,-2)]){
      tochange <- endsWith(x = egmt_GSEA_measure@result$ID, suffix = RP)
      egmt_GSEA_measure@result$NES_rand_zscore[tochange] <- scale(egmt_GSEA_measure@result$NES[tochange])
    }

    if(dim(egmt_GSEA_measure@result)[1]>0){
      pdf(paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_",gset_tag,"_signed_GSEA_measure.pdf",sep = ""),height = 5,width = 5)
      for (i in which(!substring(egmt_GSEA_measure@result$Description,1,3)%in%c("MRf","FDf","Ran"))){
        if(abs(egmt_GSEA_measure@result$NES_rand_zscore[i])>zthr & egmt_GSEA_measure@result$p.adjust[i]<0.01)
          print(enrichplot::gseaplot2(egmt_GSEA_measure, geneSetID = i, title = paste(egmt_GSEA_measure$Description[i],"NES=",as.character(round(egmt_GSEA_measure$NES[i],2)),
                                                                        "; adjP=",as.character(round(egmt_GSEA_measure$p.adjust[i],4)),"; FDR=",as.character(round(egmt_GSEA_measure$qvalues[i],4)))))
      }
      dev.off()
      GSEA_result_df<-rbind(GSEA_result_df,
                            data.frame(measure="GSEA_measure", set=gset_tag,
                                       egmt_GSEA_measure@result[!substring(egmt_GSEA_measure@result$Description,1,3)%in%c("MRf","FDf","Ran"),
                                                                c("Description","NES","NES_rand_zscore","p.adjust","qvalues")]))
    }

    egmt_GSEA_measure_abs <- clusterProfiler::GSEA(geneList = GSEA_geneList_abs, TERM2GENE=gsea_sets_RP, verbose=TRUE, pvalueCutoff = 1)
    egmt_GSEA_measure_abs@result$NES_rand_zscore <- NA
    for (RP in colnames(org_RP_df)[c(-1,-2)]){
      tochange <- endsWith(x = egmt_GSEA_measure_abs@result$ID, suffix = RP)
      egmt_GSEA_measure_abs@result$NES_rand_zscore[tochange] <- scale(egmt_GSEA_measure_abs@result$NES[tochange])
    }

    if(dim(egmt_GSEA_measure_abs@result)[1]>0){
      pdf(paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_",gset_tag,"_abs_GSEA_measure.pdf",sep = ""),height = 5,width = 5)
      for (i in which(!substring(egmt_GSEA_measure_abs@result$Description,1,3)%in%c("MRf","FDf","Ran"))) {
        if(egmt_GSEA_measure_abs@result$NES_rand_zscore[i]>abs_zthr & egmt_GSEA_measure_abs@result$p.adjust[i]<0.01)
          print(enrichplot::gseaplot2(egmt_GSEA_measure_abs, geneSetID = i, title = paste(egmt_GSEA_measure_abs$Description[i],"NES=",as.character(round(egmt_GSEA_measure_abs$NES[i],2)),
                                                                            "; adjP=",as.character(round(egmt_GSEA_measure_abs$p.adjust[i],4)),"; FDR=",as.character(round(egmt_GSEA_measure_abs$qvalues[i],4)))))
      }
      dev.off()
      GSEA_result_df<-rbind(GSEA_result_df,
                            data.frame(measure="GSEA_measure_ABS", set=gset_tag,
                                       egmt_GSEA_measure_abs@result[!substring(egmt_GSEA_measure_abs@result$Description,1,3)%in%c("MRf","FDf","Ran"),
                                                                    c("Description","NES","NES_rand_zscore","p.adjust","qvalues")]))
    }

    write.csv(x =  GSEA_result_df, file = paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_results.csv",sep = ""), row.names = FALSE)
    all_GSEA_results <- rbind(all_GSEA_results, data.frame(GSEA_result_df,comp=paste(comp,collapse = "_vs_")))
  }

  to_highlight <- all_GSEA_results[(all_GSEA_results$measure=="GSEA_measure_ABS" & all_GSEA_results$NES_rand_zscore>abs_zthr) |
                                     (all_GSEA_results$measure=="GSEA_measure" & abs(all_GSEA_results$NES_rand_zscore)>zthr), ]
  to_highlight <- to_highlight[to_highlight$qvalues<.01, ]

  cols <- c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")
  g1 <- ggplot2::ggplot(all_GSEA_results, ggplot2::aes(x=NES, y=NES_rand_zscore))+
    ggplot2::geom_hline(yintercept = 0)+
    ggplot2::geom_hline(yintercept = c(1.25), linetype="dashed", col=cols[5])+
    ggplot2::geom_hline(yintercept = c(-1.5,1.5), linetype="dotted", col=cols[5])+
    ggplot2::geom_vline(xintercept = 0, col=cols[5])+
    ggplot2::facet_grid(comp~measure,scales = "free")+
    ggplot2::scale_y_reverse()+
    ggplot2::geom_point(ggplot2::aes(col=(p.adjust<.01)))+
    ggplot2::geom_point(data = to_highlight, inherit.aes = TRUE, col="black")+
    ggrepel::geom_text_repel(data = to_highlight,
                    inherit.aes = TRUE, col="black", ggplot2::aes(label=Description),max.overlaps = 100)+
    ggplot2::scale_color_manual(values = c(cols[3],cols[2],cols[5]))+
    ggplot2::theme_bw()
  print(g1)
  ggplot2::ggsave(plot = g1, filename = paste(targetDir,"/Pred_volcanos.pdf",sep=""), width = 8, height = length(comparisons)*3,limitsize = FALSE)

  return(all_GSEA_results)

}


#' RP4 wrapper
#' @description This function allows you to run the whole RP4 pipeline
#' @param samplesFile Samples.txt file that describes file locations and sample groupings
#' @keywords RP4 pipeline
#' @export
#' @examples
#' RP4("samples.txt")
RP4 <- function(samplesFile, organism="hs", QCplot=TRUE, targetDir=NA){

  # Check organism first
  if (!(organism %in% c("hs","mm"))) {
    print("Choose a valid organism. hs(human), mm(mouse), etc.")
    return(NA)
  }

  # Assign target directory
  if (is.na(targetDir)){
    targetDir=getwd()
  }

  samples_df <- read_RP4_samples_file(samplesFile)
  rRNA_counts_df <- RP4_read_rRNA_fragments(samples_df, organism = organism, QCplot = QCplot, targetDir = targetDir)

  return(RP4_predict_heterogenity(samples_df, rRNA_counts_df, compare="group", organism=organism, QCplot=QCplot, targetDir=targetDir))
}



