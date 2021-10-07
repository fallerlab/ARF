

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
  if (!(organism %in% c("hs","mm","rm","op","ch"))) {
    print("Choose a valid organism. hs (Homo sapiens, human), mm (Mus musculus, mouse), rm (rhesus macaque, Macaca mulatta),
          op (grey short-tailed opossum, Monodelphis domestica), ch (chicken, red junglefowl, Gallus gallus) etc.")
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
    rownames(df) <- c(paste("human_28S",0:5069,sep = "_"), paste("human_18S",0:1868,sep = "_"),
                    paste("human_5.8S",0:156,sep = "_"), paste("human_5S",0:120,sep = "_"))
  } else if (organism=="mm") {
    df <- setNames(data.frame(matrix(ncol = length(samples$sampleName), nrow = (4730+1870+157+121),data=0)), samples$sampleName)
    rownames(df) <- c(paste("NR_003279.1",0:4729,sep = "_"), paste("NR_003278.3",0:1869,sep = "_"),
                      paste("NR_003280.2",0:156,sep = "_"), paste("NR_030686.1",0:120,sep = "_"))
  } else if (organism=="rm") {
    df <- setNames(data.frame(matrix(ncol = length(samples$sampleName), nrow = (4810+1869+157+119),data=0)), samples$sampleName)
    rownames(df) <- c(paste("rhesus_28S",0:4809,sep = "_"), paste("rhesus_18S",0:1868,sep = "_"),
                      paste("rhesus_5.8S",0:156,sep = "_"), paste("rhesus_5S",0:118,sep = "_"))
  } else if (organism=="op") {
    df <- setNames(data.frame(matrix(ncol = length(samples$sampleName), nrow = (4955+1891+153+119),data=0)), samples$sampleName)
    rownames(df) <- c(paste("opossum_28S",0:4954,sep = "_"), paste("opossum_18S",0:1890,sep = "_"),
                      paste("opossum_5.8S",0:152,sep = "_"), paste("opossum_5S",0:118,sep = "_"))
  } else if (organism=="ch") {
    df <- setNames(data.frame(matrix(ncol = length(samples$sampleName), nrow = (4439+1823+153+120),data=0)), samples$sampleName)
    rownames(df) <- c(paste("chicken_28S",0:4438,sep = "_"), paste("chicken_18S",0:1822,sep = "_"),
                      paste("chicken_5.8S",0:152,sep = "_"), paste("chicken_5S",0:119,sep = "_"))
  } else if (organism=="sc") {
    df <- setNames(data.frame(matrix(ncol = length(samples$sampleName), nrow = (3396+1800+157+121),data=0)), samples$sampleName)
    rownames(df) <- c(paste("yeast_25S",0:3395,sep = "_"), paste("yeast_18S",0:1799,sep = "_"),
                      paste("yeast_5.8S",0:156,sep = "_"), paste("yeast_5S",0:120,sep = "_"))
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
    ggplot2::ggsave(filename = paste(targetDir, "/", "raw_rRNA_counts.png", sep = ""), plot = g1, limitsize = FALSE,
                    width = 8, height = 3+(floor(length(samples$sampleName)**0.5)*3))
  }

  return(df[rowSums(is.na(df))==0,])
}

#' Simplify result file
#' @description Simplify results based on thresholds
#' @param RP4results RP4 results data frame
#' @param randZscore_thr Default=1
#' @param ORA_adjP_thr Default=0.05
#' @param GSEA_adjP_thr Default=0.05
#' @keywords RP4 results candidates
#' @export
#' @examples
#' RP4_simplify_results(RP4results)
RP4_simplify_results <- function(RP4results, randZscore_thr=1, ORA_adjP_thr=0.05, GSEA_adjP_thr=0.05, ORA_sig_n=0) {
  temp <- NULL
  for (comp in unique(RP4results$comp)){
    ORA_pass <- RP4results$Description[RP4results$measure=="ORA" & RP4results$comp==comp &
                                         RP4results$p.adjust<ORA_adjP_thr & RP4results$NES>ORA_sig_n]

    randZscore_GSEA_pass <- RP4results$Description[RP4results$measure=="abs_GSEA_measure" & RP4results$comp==comp &
                                                     RP4results$p.adjust < GSEA_adjP_thr &
                                                     RP4results$NES_rand_zscore > randZscore_thr]

    temp <- rbind(temp, RP4results[(RP4results$comp==comp) &(RP4results$Description%in%intersect(ORA_pass, randZscore_GSEA_pass)) ,])
  }
  return(temp)
}


#' Draw heatmap for multiple comparisons.
#' @description Draw heatmap of NES and NES_randZscore based on different thresholds.
#' @param RP4_result_df RP4 result data.
#' @param title Title for the plot and files.
#' @param targetDir Directory to save the plots in.
#' @param addedRPs Add given RPs to the plot no matter if they are significant or not..
#' @keywords Differential Ribosome Heterogeneity Heatmap
#' @export
#' @examples
#' RP4_result_heatmap(RP4results,"Title","/Folder/to/save/in/")
RP4_result_heatmap <- function(RP4_result_df, title, targetDir, addedRPs=NULL,
                               randZscore_thr=c(1,1,1.5,0), ORA_adjP_thr=c(.05,2,.05,2), GSEA_adjP_thr=c(.05,.05,.05,.05),
                               ORA_sig_n=c(0,0,0,0)){

  thrlist <- list(z=randZscore_thr,p=GSEA_adjP_thr,o=ORA_adjP_thr,n=ORA_sig_n)

  palette1 = circlize::colorRamp2(c(0,1, 1.5, 2),
                                  c("white","#D9D0D3","#CCBA72","#0F0D0E"))

  for (i in 1:length(thrlist[[1]])) {
    print(i)

    filtered_tmp <- RP4_simplify_results(RP4_result_df,
                                         GSEA_adjP_thr=thrlist[["p"]][i],
                                         randZscore_thr = thrlist[["z"]][i],
                                         ORA_adjP_thr = thrlist[["o"]][i],
                                         ORA_sig_n=thrlist[["n"]][i])

    pdf(file = paste(targetDir,title,
                     "_randZ",as.character(thrlist[["z"]][i]*100),
                     "_GSEAp",as.character(thrlist[["p"]][i]*100),
                     "_ORAp",as.character(thrlist[["o"]][i]*100),
                     "_SigPosN",as.character(thrlist[["n"]][i]),
                     ".pdf",sep=""), width = 16, height = 16)
    set.seed(i)

    allstudies_abs_sig_RPS <- unique(filtered_tmp$Description[filtered_tmp$NES_rand_zscore > thrlist[["z"]][i] &
                                                                filtered_tmp$p.adjust < thrlist[["p"]][i] &
                                                                filtered_tmp$measure=="abs_GSEA_measure"])


    temp <- filtered_tmp[filtered_tmp$measure=="abs_GSEA_measure" &
                           filtered_tmp$Description%in%allstudies_abs_sig_RPS &
                           filtered_tmp$NES_rand_zscore > thrlist[["z"]][i] &
                           filtered_tmp$p.adjust < thrlist[["p"]][i], ]

    if (length(addedRPs)>0)
      temp <- unique(rbind(temp,RP4_result_df[RP4_result_df$Description %in% addedRPs & filtered_tmp$measure=="abs_GSEA_measure",]))

    comp_all_abs_matrix <- reshape2::acast(temp, formula=Description ~ comp, value.var = "NES")
    comp_all_abs_matrix[is.na(comp_all_abs_matrix)] <- 0
    ht <- ComplexHeatmap::Heatmap(comp_all_abs_matrix, name = "NES", na_col = "#0F0D0E",
                  column_title = paste(title,"absGSEA NES"), col = palette1)

    ComplexHeatmap::draw(ht, padding = grid::unit(c(30, 2, 2, 2), "mm"))

    ## NESrand
    comp_all_abs_matrix <- reshape2::acast(temp, formula=Description ~ comp, value.var = "NES_rand_zscore")
    comp_all_abs_matrix[is.na(comp_all_abs_matrix)] <- 0
    ht <- ComplexHeatmap::Heatmap(comp_all_abs_matrix, name = "NES_rand_zscore", na_col = "#0F0D0E",
                  column_title = paste(title,"absGSEA NES_rand_zscore"), col = palette1)

    ComplexHeatmap::draw(ht, padding = grid::unit(c(30, 2, 2, 2), "mm"))
    dev.off()
  }
}


#' Draw RP4 result scatterplot for multiple comparisons
#' @description Draw scatterplot for GSEA and ORA analyses
#' @param RP4_result_df Full RP4 result data.
#' @param targetDir Directory to save the plots in.
#' @param title Default is "DRH_prediction_volcanos"
#' @keywords Differential Ribosome Heterogeneity Heatmap
#' @export
#' @examples
#' RP4_result_scatterplot(RP4results,"/Folder/to/save/in/")
RP4_result_scatterplot <- function(RP4_result_df, targetDir, title="DRH_prediction_volcanos",
                                   randZscore_thr=1, ORA_adjP_thr=0.05, GSEA_adjP_thr=0.05, ORA_sig_n=-1){
  to_highlight <- RP4_simplify_results(RP4_result_df,
                                       randZscore_thr=randZscore_thr, ORA_adjP_thr=ORA_adjP_thr, GSEA_adjP_thr=GSEA_adjP_thr, ORA_sig_n=ORA_sig_n)

  cols <- c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")
  RP4_result_df$colorchange <- (RP4_result_df$p.adjust < min(ORA_adjP_thr,GSEA_adjP_thr))
  RP4_result_df$colorchange[RP4_result_df$measure!="ORA"] <-  (RP4_result_df$p.adjust[RP4_result_df$measure!="ORA"] < ORA_adjP_thr)
  RP4_result_df$colorchange[RP4_result_df$measure=="ORA"] <-  (RP4_result_df$p.adjust[RP4_result_df$measure=="ORA"] < GSEA_adjP_thr)

  g1 <- ggplot2::ggplot(RP4_result_df,
                        ggplot2::aes(x=NES_rand_zscore,y=NES))+
    ggplot2::geom_hline(yintercept = 0)+
    ggplot2::geom_hline(yintercept = c(1), linetype="dashed", col=cols[5])+
    ggplot2::geom_vline(xintercept = 0, col=cols[5])+
    ggplot2::facet_wrap(comp~measure,scales = "free")+
    ggplot2::geom_point(ggplot2::aes(col=colorchange))+
    ggplot2::geom_point(data = to_highlight, inherit.aes = TRUE, col="black")+
    ggrepel::geom_text_repel(data = to_highlight,
                             inherit.aes = TRUE, col="black", ggplot2::aes(label=Description),max.overlaps = 100)+
    ggplot2::scale_color_manual(values = c(cols[3],cols[2],cols[5]))+
    ggplot2::theme_bw()+
    ggplot2::labs(col=paste0("Significancy\nORA_padj<",as.character(ORA_adjP_thr),"\nGSEA_padj<",as.character(GSEA_adjP_thr)))+
    ggplot2::xlab("NES | rRNA proximity set size")+
    ggplot2::ylab("NES z-score | No of significant rRNA positions")
  print(g1)
  ggplot2::ggsave(plot = g1, filename = paste(targetDir, "/", title, ".pdf",sep=""),
                  width = 4+(ceiling(length(unique(RP4_result_df$comp))**0.5)*4),
                  height = ceiling(length(unique(RP4_result_df$comp))**0.5)*4,limitsize = FALSE)
}





#' Draw Differential rRNA abundance heatmap and its overlap with RP-rRNA proximity sets
#' @description Draw rRNA fragment change heatmaps to visualize position-specific differential rRNA fragment abundance
#' @param RP4_result_df Full RP4 result data.
#' @param RP4_DRF rRNA position specific differential rRNA fragment abundace results
#' @param organism Organism abbrevation. Pass "hs" for human and "mm" for mouse.
#' @param targetDir Directory to save the plots in.
#' @param title Default is "DRH_prediction_volcanos"
#' @keywords Differential Ribosome Heterogeneity Heatmap
#' @export
#' @examples
#' RP4_rRNApos_heatmaps(RP4results,"/Folder/to/save/in/")
RP4_rRNApos_heatmaps <- function(RP4_result_df, RP4_DRF, organism, RPs, targetDir, title="rRNA_pos_spec_heatmap",
                                   randZscore_thr=1, ORA_adjP_thr=0.05, GSEA_adjP_thr=0.05, ORA_sig_n=-1){
  org_RP_df <- NULL
  if (organism=="hs"){
    org_RP_df <- RP4:::RP_proximity_human_df
  } else if (organism=="mm") {
    org_RP_df <- RP4:::RP_proximity_mouse_df
  } else if (organism=="rm") {
    org_RP_df <- RP4:::RP_proximity_rhesus_df
  } else if (organism=="op") {
    org_RP_df <- RP4:::RP_proximity_opossum_df
  } else if (organism=="ch") {
    org_RP_df <- RP4:::RP_proximity_chicken_df
  } else if (organism=="sc") {
    org_RP_df <- RP4:::RP_proximity_yeast_df
  } else {
    print(paste(c("Organism", organism, "Not implemented yet!"), collapse = " "))
    return(NULL)
  }

  prox_col = circlize::colorRamp2(c(0,24.999,25.0,50,200,300),c("black","black","grey30","grey60","grey99","white"))

  profileplot <- RP4_DRF
  profileplot$sig <- (abs(profileplot$log2FoldChange)>0.5) & (profileplot$padj<0.05)
  if (organism =="mm") {
    profileplot$rrna <-vapply(strsplit(profileplot$pos,".",fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
  } else{
    profileplot$rrna <-vapply(strsplit(profileplot$pos,"_",fixed = TRUE), `[`, 2, FUN.VALUE=character(1))
  }
  profileplot$ipos <- as.numeric(vapply(strsplit(profileplot$pos,"_",fixed = TRUE), `[`, 3, FUN.VALUE=character(1)))
  profileplot$tomatrix <- abs(profileplot$log2FoldChange)*profileplot$sig
  profileplot$signedmatrix <- profileplot$log2FoldChange*profileplot$sig

  filtered_results <- RP4_simplify_results(RP4_result_df, randZscore_thr = randZscore_thr, ORA_adjP_thr=ORA_adjP_thr,
                                           GSEA_adjP_thr=GSEA_adjP_thr, ORA_sig_n=ORA_sig_n)
  #profileplot <- profileplot[profileplot$comp%in%filtered_results$comp[filtered_results$Description%in%RPs],]

  rrnas<-unique(org_RP_df$rRNA)
  focused_order <- c()
  for (rrna in unique(profileplot$rrna)){
    iorder <- order(profileplot$ipos)
    focused_order<-c(focused_order,unique(profileplot$pos[iorder[profileplot$rrna[iorder]==rrna]]))
  }

  temp <- org_RP_df[focused_order,]
  rids <- c()
  extended_rids <- c()
  for (RP in RPs){
    tmp_rids <- which(temp[,RP]<25)
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
    get_pos <- list(val_pos=val_pos,tags=tags)}

  whichPos <- get_pos[[1]]
  tags <- get_pos[[2]]

  foo <- NULL
  for (RP in RPs){
    foo<-cbind(foo, org_RP_df[focused_order[whichPos],RP])
  }
  colnames(foo)<-RPs

  ha <- ComplexHeatmap::HeatmapAnnotation(rRNA = org_RP_df[focused_order,"rRNA"][whichPos],
                         foo=foo, col=list(foo=prox_col,
                                           rRNA=setNames(wesanderson::wes_palette("Rushmore1")[c(1,3,4,5)],rrnas), na_col = "white"),
                         annotation_legend_param = list(foo=list(title="RP-rRNA\nproximity map\ndistance (Å)",at=c(0,25,100,200,300))))

  matrixplot <- reshape2::acast(profileplot, formula=comp ~ pos, value.var = "tomatrix")[,unique(profileplot$pos)]

  pdf(paste0(targetDir,"/",title,".pdf"),height = 6,width = 12)
  ht <- ComplexHeatmap::Heatmap(matrixplot[,focused_order[whichPos]], name = "abs(log2FC)", top_annotation = ha,
                                cluster_columns = FALSE,show_column_names = FALSE,
                col=circlize::colorRamp2(c(0,0.5,8),c("white","cornsilk","black")),na_col = "white",
                column_split = factor(tags,levels = unique(tags),labels = unique(tags)),
                column_title_rot = 90, heatmap_legend_param = list(at=c(0.5,2,4,6,8)))
  ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 10, 2), "mm"))
  dev.off()
  return(ht)
}

#' Predict Differential Ribosomal Heterogeneity candidates with rRNA count data
#' @description Overlap differential rRNA count data with 3d ribosome, rRNA-RP proximity data and predict heterogeneity candidates across groups.
#' @param samples Samples dataframe created by read_samples_file() function.
#' @param organism Organism abbrevation. Pass "hs" for human and "mm" for mouse.
#' @param QCplot TRUE or FALSE, whether to generate QC plots or not.
#' @param targetDir Direcory to save QC plots in.
#' @param comparisons List of comparisons to be included.
#' @keywords Diffferential Ribosome Heterogeneity rRNA ribosome RP
#' @export
#' @examples
#' RP4_predict_heterogenity(samples_df, rRNA_counts_df, organism="hs", QCplot=TRUE)
RP4_predict_heterogenity <- function(samples, rRNA_counts, compare="group", organism="hs", QCplot=FALSE, targetDir=NA, comparisons=NULL) {

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

  if(is.null(comparisons) || length(comparisons)==0) {
    comparisons <- list()
    for (i in 1:(s_l-1)) {
      for (j in (i+1):s_l) {
        print(paste("Comparing",s_n[i],"vs",s_n[j]))
        comparisons[[(length(comparisons) +1)]] <- c(s_n[i],s_n[j])
        res <- DESeq2::results(dds, contrast=c("DESEQcondition",s_n[i],s_n[j]), lfcThreshold = 0.5, alpha = 0.05, cooksCutoff = FALSE )
        print(DESeq2::summary(res))
      }
    }
  }

  # Read count transformations
  vsd <- DESeq2::vst(dds, blind=FALSE)

  if (QCplot) {
    g1 <- ggplot2::ggplot(reshape2::melt(SummarizedExperiment::assay(vsd)), ggplot2::aes(x=Var2, y=value)) +
      ggplot2::geom_violin()+ggplot2::coord_flip()+
      ggplot2::xlab("Positional rRNA fragment counts (vsd) after DESeq2 normalization")+ggplot2::ylab("Samples")
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "norm_rRNA_counts.png", sep = ""), limitsize = FALSE, plot = g1,
                    width = 8, height = 3+(floor(length(samples$sampleName)**0.5)*3))

    #PCA_analysis
    pca <- prcomp(SummarizedExperiment::assay(vsd)[order(matrixStats::rowVars(SummarizedExperiment::assay(vsd)),decreasing = TRUE),])
    pca_df <- data.frame(pca$rotation)
    pca_map <- paste(colnames(pca_df),' (Var.Exp.:%',round(100*((pca$sdev^2)/sum(pca$sdev^2)),1),')', sep='')
    pca_df$sample <- rownames(pca_df)
    pca_df$group <- sapply(pca_df$sample,FUN=function(x){return(samples$DESEQcondition[samples$sampleName==x])})

    g1 <- ggplot2::ggplot(pca_df, ggplot2::aes(x=PC1,y=PC2,color=group)) + ggplot2::geom_point() +
      ggplot2::labs(x=pca_map[1], y=pca_map[2]) + ggrepel::geom_text_repel(ggplot2::aes(label=sample), show.legend = FALSE)
    print(g1)
    ggplot2::ggsave(filename = paste(targetDir, "/", "vsd_based_PCA.png", sep = ""), plot = g1,limitsize = FALSE,
                    width = 2+(floor(length(samples$sampleName)**0.5)*2), height = 2+(floor(length(samples$sampleName)**0.5)*2))
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
  } else if (organism=="rm") {
    org_RP_df <- RP4:::RP_proximity_rhesus_df
    gsea_sets_RP <- RP4:::rhesus_gsea_sets_RP
  } else if (organism=="op") {
    org_RP_df <- RP4:::RP_proximity_opossum_df
    gsea_sets_RP <- RP4:::opossum_gsea_sets_RP
  } else if (organism=="ch") {
    org_RP_df <- RP4:::RP_proximity_chicken_df
    gsea_sets_RP <- RP4:::chicken_gsea_sets_RP
  } else if (organism=="sc") {
    org_RP_df <- RP4:::RP_proximity_yeast_df
    gsea_sets_RP <- RP4:::yeast_gsea_sets_RP
  } else {
    print(paste(c("Organism", organism, "Not implemented yet!"), collapse = " "))
    return(NULL)
  }

  ########## Overrepresentation Analysis ############
  RP_pathways <- sapply(colnames(org_RP_df)[c(-1,-2)],FUN=function(x){return(rownames(org_RP_df)[which(org_RP_df[,x]<25)])})

  ######### Gene set enrichment Analysis ############

  #library(clusterProfiler)
  #library(enrichplot)
  #library(gprofiler2)

  all_GSEA_results <- NULL
  # Separate for each DESEQcondition
  for (comp in comparisons) {
    print(paste("Running predictions for",comp[1],"vs",comp[2]))

    GSEA_result_df <- NULL

    res <- DESeq2::results(dds, contrast=c("DESEQcondition",comp[1],comp[2]), cooksCutoff = FALSE)

    ### Overrepresentation hook ####
    or_df <- fgsea::fora(pathways = RP_pathways,genes = rownames(res)[which(res$padj<.05 & abs(res$log2FoldChange)>0.5)], universe = rownames(res), minSize = 10)

    GSEA_result_df <- data.frame(measure="ORA", Description = or_df$pathway, NES=or_df$overlap, NES_rand_zscore=or_df$size,
                                 p.adjust=or_df$padj, qvalues=or_df$pval)

    temp_df <- res
    temp_df$weight <- scales::rescale(log10(temp_df$baseMean), to = c(0, 5))
    temp_df$padj[temp_df$padj<0.00001] = 0.00001

    for (measureID in c("abs_GSEA_measure")){ #},"GSEA_measure","w_GSEA_m", "abs_w_GSEA_m")){
      used_measure <- NULL
      if (measureID=="GSEA_measure")
        used_measure <- res$log2FoldChange*(-log10(temp_df$padj))
      else if (measureID=="abs_GSEA_measure")
        used_measure <- abs(res$log2FoldChange)*(-log10(temp_df$padj))
      else if (measureID=="w_GSEA_m")
        used_measure <- res$log2FoldChange*(-log10(temp_df$padj))*temp_df$weight
      else if (measureID=="abs_w_GSEA_m")
        used_measure <- abs(res$log2FoldChange)*(-log10(temp_df$padj))*temp_df$weight

      names(used_measure)<-rownames(res)
      used_measure <- used_measure[!sapply(used_measure, function(x) is.na(x))]
      used_geneList <- used_measure[order(used_measure, decreasing = TRUE)]
      used_geneList_abs <- abs(used_measure)[order(abs(used_measure), decreasing = TRUE)]

      egmt_used_measure <- clusterProfiler::GSEA(geneList = used_geneList, TERM2GENE=gsea_sets_RP, verbose=TRUE, pvalueCutoff = 1)
      egmt_used_measure@result$NES_rand_zscore <- NA
      for (RP in colnames(org_RP_df)[c(-1,-2)]){
        tochange <- endsWith(x = egmt_used_measure@result$ID, suffix = RP)
        egmt_used_measure@result$NES_rand_zscore[tochange] <- scale(egmt_used_measure@result$NES[tochange])
      }

      if(dim(egmt_used_measure@result)[1]>0){
        pdf(paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_", measureID,".pdf",sep = ""),height = 5,width = 5)
        for (i in which(!substring(egmt_used_measure@result$Description,1,3)%in%c("MRf","FDf","Ran"))){
          if(egmt_used_measure@result$NES_rand_zscore[i]>1 & egmt_used_measure@result$p.adjust[i]<0.01)
            print(enrichplot::gseaplot2(egmt_used_measure, geneSetID = i, title = paste(egmt_used_measure$Description[i],"NES=",as.character(round(egmt_used_measure$NES[i],2)),
                                                                          "; adjP=",as.character(round(egmt_used_measure$p.adjust[i],4)),"; FDR=",as.character(round(egmt_used_measure$qvalues[i],4)))))
        }
        dev.off()
        GSEA_result_df<-rbind(GSEA_result_df,
                              data.frame(measure=measureID,
                                         egmt_used_measure@result[!substring(egmt_used_measure@result$Description,1,3)%in%c("MRf","FDf","Ran"),
                                                                  c("Description","NES","NES_rand_zscore","p.adjust","qvalues")]))
      }
    }

    write.csv(x =  GSEA_result_df, file = paste(targetDir,"/",paste(comp,collapse = "_vs_"),"_results.csv",sep = ""), row.names = FALSE)
    all_GSEA_results <- rbind(all_GSEA_results, data.frame(GSEA_result_df,comp=paste(comp,collapse = "_vs_")))
  }

  return(all_GSEA_results)
}


#' Get lFC profiles of RP proximity sets
#' @description Fish out the Differential values for RP proximity sets
#' @param Samples.txt file that describes file locations and sample groupings
#' @param organism Organism abbrevation. Pass "hs" for human and "mm" for mouse.
#' @param compare If you want to compare samples based on other grouping, choose the columnname that is given in samplesFile (Default=group).
#' @param comparisons List of comparisons to be included.
#' @param exclude Sample ids to be excluded from analysis
#' @param save_sample_sim Create a sample similarity heatmap
#' @param savefolder where to save
#' @keywords lFCprofile RP rRNA proximity sets
#' @export
#' @examples
#' RP4("samples.txt", organism="mm")
RP4_report_RPspec_pos_results <- function(samplesFile, organism="hs", compare="group", comparisons=NULL, exclude=NULL, save_sample_sim=FALSE, savefolder=NA){

  # Check organism first
  if (!RP4_check_organism(organism))
    return(NA)

  samples <- read_RP4_samples_file(samplesFile)
  if(!is.null(exclude))
    samples <- samples[!samples$sampleName%in%exclude,]

  rRNA_counts <- RP4_read_rRNA_fragments(samples, organism = organism, QCplot = FALSE)

  s_n <- unique(samples[,compare])
  s_l <- length(s_n)

  samples$DESEQcondition <- samples[,compare]
  cts <- as.matrix(round(rRNA_counts, digits = 0))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = samples, design = ~DESEQcondition)

  keep <- rowSums(DESeq2::counts(dds)) >= 1000*dim(samples)[1]
  dds <- dds[keep,]
  dds <- DESeq2::DESeq(dds)

  if(is.null(comparisons) || length(comparisons)==0) {
    comparisons <- list()
    for (i in 1:(s_l-1)) {
      for (j in (i+1):s_l) {
        comparisons[[(length(comparisons) +1)]] <- c(s_n[i],s_n[j])
      }
    }
  }

  if(save_sample_sim){
    if (is.na(savefolder))
      savefolder=getwd()

    vsd <- DESeq2::vst(dds, blind=FALSE)
    sampleDists <- dist(t(SummarizedExperiment::assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- vsd$sampleName
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
    write.table(x = sampleDistMatrix, file = paste0(savefolder,'/Sample_similarity_distance_vsdbased.tsv'), sep = "\t",quote = FALSE,
              row.names = TRUE, col.names = rownames(sampleDistMatrix))
    pdf(paste0(savefolder,'/Sample_similarity_distance_vsdbased.pdf'),
        width = 3+(floor(length(samples$sampleName)**0.5)*3), height = 3+(floor(length(samples$sampleName)**0.5)*3))
    ht <- ComplexHeatmap::pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
    ComplexHeatmap::draw(ht)
    dev.off()
  }

  all_results <- NULL
  for (comp in comparisons){
    print(paste("Comparing",comp[1],"vs",comp[2]))
    res <- data.frame(DESeq2::results(dds, contrast=c("DESEQcondition",comp[1],comp[2]), cooksCutoff = FALSE ))
    res$pos <- rownames(res)
    res$comp <- paste0(comp[1],"_vs_",comp[2])
    all_results <- rbind(all_results, res)
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
  } else if (organism=="rm") {
    org_RP_df <- RP4:::RP_proximity_rhesus_df
    gsea_sets_RP <- RP4:::rhesus_gsea_sets_RP
  } else if (organism=="op") {
    org_RP_df <- RP4:::RP_proximity_opossum_df
    gsea_sets_RP <- RP4:::opossum_gsea_sets_RP
  } else if (organism=="ch") {
    org_RP_df <- RP4:::RP_proximity_chicken_df
    gsea_sets_RP <- RP4:::chicken_gsea_sets_RP
  } else if (organism=="sc") {
    org_RP_df <- RP4:::RP_proximity_yeast_df
    gsea_sets_RP <- RP4:::yeast_gsea_sets_RP
  } else {
    print(paste(c("Organism", organism, "Not implemented yet!"), collapse = " "))
    return(NULL)
  }
  rownames(all_results) <- 1:(dim(all_results)[1])
  return(all_results[all_results$pos%in%unique(gsea_sets_RP$gene),])
}


#' RP4 wrapper
#' @description This function allows you to run the whole RP4 pipeline
#' @param samplesFile Samples.txt file that describes file locations and sample groupings
#' @param organism Organism abbrevation. Pass "hs" for human and "mm" for mouse.
#' @param comparison If you want to compare samples based on other grouping, choose the columnname that is given in samplesFile (Default=group).
#' @param QCplot TRUE or FALSE, whether to generate QC plots or not.
#' @param targetDir Direcory to save QC plots in.
#' @param comparisons List of comparisons to be included.
#' @param exclude Sample ids to be excluded from analysis
#' @keywords RP4 pipeline
#' @export
#' @examples
#' RP4("samples.txt", organism="mm", targetDir="/target/directory/to/save/results")
RP4 <- function(samplesFile, organism="hs", comparison="group", QCplot=TRUE,  targetDir=NA, comparisons=NULL, exclude=NULL){

  # Check organism first
  if (!RP4_check_organism(organism))
    return(NA)

  # Assign target directory
  if (is.na(targetDir)){
    targetDir=getwd()
  }

  samples_df <- read_RP4_samples_file(samplesFile)

  if(!is.null(exclude))
    samples_df <- samples_df[!samples_df$sampleName%in%exclude,]

  rRNA_counts_df <- RP4_read_rRNA_fragments(samples_df, organism = organism, QCplot = QCplot, targetDir = targetDir)
  RP4_results <- RP4_predict_heterogenity(samples_df, rRNA_counts_df, compare="group", organism=organism, QCplot=QCplot, targetDir=targetDir, comparisons = comparisons)

  return(RP4_results)
}



