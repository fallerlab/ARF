
# Copyright (C) 2024  Ferhat Alkan
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


#' Draw dricARF result scatterplot for multiple comparisons
#' @description Draw scatterplot for GSEA and ORA analyses
#' @param dricARF_results Full dripARF result data.
#' @param targetDir Directory to save the plots in.
#' @param title Default is "dricARF (highlighted) & dripARF predictions"
#' @param addedRPs Add given RPs to the plot no matter if they are significant or not.
#' @param highlightRPs List of RPs to highlight instead of highlighting the top-predicted RPs.
#' @param randZscore_thr Default=1
#' @param ORA_adjP_thr Default=0.05
#' @param RPSEA_adjP_thr Default=0.05
#' @param ORA_sig_n Default=1
#' @keywords Differential Ribosome Heterogeneity Heatmap
#' @export
#' @examples
#' dripARF_result_scatterplot(dripARF_results, "/Folder/to/save/in/")
dricARF_result_scatterplot <- function(dricARF_results, targetDir,
                                       title="dricARF (highlighted) & dripARF predictions", addedRPs=NULL, highlightRPs=NULL,
                                       randZscore_thr=1, ORA_adjP_thr=0.05, RPSEA_adjP_thr=0.05, ORA_sig_n=1){
  `%>%` <- magrittr::`%>%`
  
  RP_results <- dricARF_results %>% dplyr::filter(!(Description %in% unique(ARF:::human_gsea_sets_Collision$ont)))
  final_coll <- dricARF_results %>% dplyr::filter(Description %in% unique(ARF:::human_gsea_sets_Collision$ont))
  #Exclude Col Int details
  final_coll <- final_coll %>% dplyr::filter(!Description%in%c( "sc_6I7O_Col.Int.", "sc_6T83_Col.Int.", "hs_7QVP_Col.Int."))
  cols <- c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")
  g1 <- (ggplot2::ggplot(NULL, ggplot2::aes(x=RPSEA.NES_randZ, y=RPSEA.NES))+
     ggplot2::geom_hline(yintercept = 0) + ggplot2::geom_vline(xintercept = c(randZscore_thr), linetype="dashed", col=cols[5]) + 
     ggplot2::geom_vline(xintercept = 0, col=cols[5])+ggplot2::facet_grid(comp~.)+
     ggplot2::geom_point(data=RP_results, ggplot2::aes(shape=RPSEA.padj<RPSEA_adjP_thr), col=cols[4], size=1)+
     ggplot2::geom_point(data=final_coll, ggplot2::aes(color=Description), alpha=.6,size=2.5,shape=3)+
     ggplot2::geom_point(data=final_coll%>%dplyr::filter(RPSEA.padj<RPSEA_adjP_thr), ggplot2::aes(color=Description), alpha=.6,size=2.5,shape=10)+
     ggplot2::geom_point(data=final_coll%>%dplyr::filter(RPSEA.padj>=RPSEA_adjP_thr), ggplot2::aes(color=Description), alpha=.6,size=2.5,shape=8)+
     ggrepel::geom_label_repel(data = final_coll%>%dplyr::filter(Description%in%c("Rib.Col.",addedRPs)), alpha=.8, size=2,
                               inherit.aes = TRUE, ggplot2::aes(label=Description),color="#000000", show.legend = F)+
     ggplot2::scale_color_manual(values = c("#e41a1c",  "#fb9a99",        "#d95f02", "#1a9641",           "#1f78b4", "#fb2ae9",          "#ab1be7", "#000000", "#fb9a99", rep("#110134",length(addedRPs))),
                                 breaks = c( "hs_7QVP_SAS", "hs_7QVP_Col.Int.","sc_6I7O_SAS", "sc_6I7O_Col.Int.", "sc_6T83_SAS", "sc_6T83_Col.Int.", "sc_6SV4_SAS", "Rib.Col.","Col.Int.", addedRPs))+
     ggplot2::scale_shape_manual(values = c(16,18), breaks = c(TRUE,FALSE))+
     ggplot2::theme_bw()+
     ggplot2::labs(col="Collision Prediction", shape=paste0("RPSEA padj<",as.character(RPSEA_adjP_thr)))+
     ggplot2::xlab("Enrichment Score 2\n(NES->Z-score within random sets)")+
     ggplot2::ylab("Enrichment Score 1 (RPSEA NES)"))
  
  print(g1)
  
  g2 <- (ggplot2::ggplot(NULL, ggplot2::aes(x=RPSEA.NES_randZ, y=-log10(RPSEA.padj)))+
           ggplot2::geom_hline(yintercept = c(-log10(RPSEA_adjP_thr)), linetype="dashed", col=cols[5]) + 
           #ggplot2::geom_vline(xintercept = 0, col=cols[5])+
           ggplot2::geom_vline(xintercept = c(randZscore_thr), linetype="dashed", col=cols[5]) + 
           ggplot2::facet_grid(comp~.)+
           ggplot2::geom_point(data=RP_results, ggplot2::aes(shape=RPSEA.padj<RPSEA_adjP_thr), col=cols[4], size=1)+
           ggplot2::geom_point(data=final_coll, ggplot2::aes(color=Description), alpha=.6,size=2.5,shape=3)+
           ggplot2::geom_point(data=final_coll%>%dplyr::filter(RPSEA.padj<RPSEA_adjP_thr), ggplot2::aes(color=Description), alpha=.6,size=2.5,shape=10)+
           ggplot2::geom_point(data=final_coll%>%dplyr::filter(RPSEA.padj>=RPSEA_adjP_thr), ggplot2::aes(color=Description), alpha=.6,size=2.5,shape=8)+
           ggrepel::geom_label_repel(data = final_coll%>%dplyr::filter(Description%in%c("Rib.Col.",addedRPs)), alpha=.8, size=2,
                                     inherit.aes = TRUE, ggplot2::aes(label=Description),color="#000000", show.legend = F)+
           ggplot2::scale_color_manual(values = c("#e41a1c",  "#fb9a99",        "#d95f02", "#1a9641",           "#1f78b4", "#fb2ae9",          "#ab1be7", "#000000", "#fb9a99", rep("#110134",length(addedRPs))),
                                       breaks = c( "hs_7QVP_SAS", "hs_7QVP_Col.Int.","sc_6I7O_SAS", "sc_6I7O_Col.Int.", "sc_6T83_SAS", "sc_6T83_Col.Int.", "sc_6SV4_SAS", "Rib.Col.","Col.Int.", addedRPs))+
           ggplot2::scale_shape_manual(values = c(16,18), breaks = c(TRUE,FALSE))+
           ggplot2::theme_bw()+
           ggplot2::labs(col="Collision Prediction", shape=paste0("RPSEA padj<",as.character(RPSEA_adjP_thr)))+
           ggplot2::xlab("Enrichment Score 2\n(NES->Z-score within random sets)")+
           ggplot2::ylab("-log10(RPSEA.padj)"))
  
  print(g2)
  
  gcombine <- cowplot::plot_grid(g1+ggplot2::theme(legend.position = "none"),g2,ncol=2,nrow=1,rel_widths = c(2,3))
  
  ggplot2::ggsave(plot = gcombine, filename = paste(targetDir, "/", title, ".pdf",sep=""),
                  width = 8, height = 3+(ceiling(length(unique(dricARF_results$comp)))*2),limitsize = FALSE)
  
  return(g1)
}

#' dricARF wrapper
#' @description This function allows you to run the whole dricARF pipeline
#' @param samplesFile File that describes file locations and sample groupings
#' @param rRNAs_fasta Fasta file for the 4 rRNAs of the organism.
#' @param samples Samples dataframe created by read_ARF_samples_file() function.
#' @param organism Organism abbrevation. Pass "hs" for human, "mm" for mouse, and "sc" for yeast.
#' @param compare If you want to compare samples based on other grouping, choose the columnname that is given in samplesFile (Default=group).
#' @param QCplot TRUE or FALSE, whether to generate QC plots or not.
#' @param targetDir Directory to save the QC plots in.
#' @param comparisons List of comparisons to be included.
#' @param exclude List of sample names to be excluded from the analysis.
#' @param GSEAplots Whether to produce and save the standard GSEA plots.
#' @param gsea_sets_RP RP-rRNA contact point sets to perform dripARF enrichments on. (preset for hs, mm, and sc predictions)
#' @param RP_proximity_df RP-rRNA proximity matrix that is calculated by ARF. (preset for hs, mm, and sc predictions)
#' @param gsea_sets_Collision Ribosome collision sets to perform dricARF enrichments on. (preset for hs, mm, and sc predictions)
#' @keywords dricARF pipeline
#' @export
#' @examples
#' dricARF("samples.txt", "rRNAs.fa", organism="mm", targetDir="/target/directory/to/save/results")
dricARF <- function(samplesFile, rRNAs_fasta, samples_df=NULL, organism=NULL, compare="group", QCplot=TRUE,  targetDir=NA,
                    comparisons=NULL, exclude=NULL, GSEAplots=FALSE, gsea_sets_RP=NULL, RP_proximity_df=NULL, gsea_sets_Collision=NULL){
  
  if(is.null(RP_proximity_df)){
    # Check organism first
    if (!ARF_check_organism(organism))
      return(NA)
    # preset organisms
    if (organism=="hs"){
      RP_proximity_df <- ARF:::RP_proximity_human_df
    } else if (organism=="mm") {
      RP_proximity_df <- ARF:::RP_proximity_mouse_df
    } else if (organism=="sc") {
      RP_proximity_df <- ARF:::RP_proximity_yeast_df
    } else {
      message(paste(c("Organism", organism, "is not preset/implemented yet! Please generate your own set of proximity matrix, RP-rRNA and ribosome collision sets using ARF.\n"), 
                    collapse = " "))
      return(NULL)
    }
  }
  
  if(is.null(gsea_sets_RP)){
    # Check organism first
    if (!ARF_check_organism(organism))
      return(NA)
    # preset organisms
    if (organism=="hs"){
      gsea_sets_RP <- ARF:::human_gsea_sets_RP
    } else if (organism=="mm") {
      gsea_sets_RP <- ARF:::mouse_gsea_sets_RP
    } else if (organism=="sc") {
      gsea_sets_RP <- ARF:::yeast_gsea_sets_RP
    } else {
      message(paste(c("Organism", organism, "is not preset/implemented yet! Please generate your own set of proximity matrix, RP-rRNA and ribosome collision sets using ARF.\n"), 
                    collapse = " "))
      return(NULL)
    }
  }
  
  if(is.null(gsea_sets_Collision)){
    # Check organism first
    if (!ARF_check_organism(organism))
      return(NA)
    # preset organisms
    if (organism=="hs"){
      gsea_sets_Collision <- ARF:::human_gsea_sets_Collision
    } else if (organism=="mm") {
      gsea_sets_Collision <- ARF:::mouse_gsea_sets_Collision
    } else if (organism=="sc") {
      gsea_sets_Collision <- ARF:::yeast_gsea_sets_Collision
    } else {
      message(paste(c("Organism", organism, "is not preset/implemented yet! Please generate your own set of proximity matrix, RP-rRNA and ribosome collision sets using ARF.\n"), 
                    collapse = " "))
      return(NULL)
    }
  } else{ # Add random sets for randZ the input gsea_sets_Collision
    
  }
  
  # Maybe check for coherence of rRNA ids ??
  ## TODOwork
  
  ## Lets integrate it all 
  added_sets <- unique(gsea_sets_Collision$ont)
  added_sets <- added_sets[!grepl("Rand",added_sets)]
  for (colset in added_sets){
    RP_proximity_df[,colset] <- 100
    RP_proximity_df[gsea_sets_Collision$gene[gsea_sets_Collision$ont==colset], colset] <- 1
  }
  gsea_sets_RP <- rbind(gsea_sets_RP, gsea_sets_Collision)

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
  
  results <- dripARF_predict_heterogenity(samples = samples_df, rRNAs_fasta=rRNAs_fasta, rRNA_counts = rRNA_counts_df,
                               compare="group", organism=organism, QCplot=QCplot, targetDir=targetDir,
                               comparisons = comparisons, GSEAplots=GSEAplots, 
                               gsea_sets_RP = gsea_sets_RP, RP_proximity_df = RP_proximity_df,
                               measureID = "abs_GSEA_measure_with_dynamic_p", runID="dricARF")
  
  dricARF_result_scatterplot(dricARF_results = results, targetDir = targetDir, title = "dricARF (highlighted) & dripARF predictions")
  # dripARF_result_heatmap(dripARF_results = results, targetDir = targetDir, title = "ALL dripARF predictions",
  #                        randZscore_thr = c(1), ORA_adjP_thr = c(0.05), RPSEA_adjP_thr = c(0.05), ORA_sig_n = 1)
  
  return(results)
}


#' Convert 3D Ribosome distance file to target organism through rRNA alignments!
#' @description Align rRNAs from hs/sc and target organism, then, convert the collision sets into target organism coordinates.
#' @param target_species ID for the target species, i.e. mm, sc, etc. 
#' @param target_rRNAs_fasta Fasta file for the rRNAs of the target organism. Same file used in rRNA fragment alignment.
#' @param rRNA_pairs List of rRNA ID pairs matching source and target rRNAs. (Use 28S, 18S, 5.8S, 5S as source ids) i.e. list(c("28S","species_28S"), c("18S","species_18S") etc.).
#' @keywords 3D ribosome analysis using PDB file
#' @export
#' @examples
#' ARF_convert_ribosome3D_rRNA_pos()
dricARF_liftover_collision_sets <- function(target_species, target_rRNAs_fasta, rRNA_pairs=list()) {
  # dplyr hack for %>%
  `%>%` <- magrittr::`%>%`
  
  yeast_rRNAs <- Biostrings::readBStringSet(file = system.file("extdata", "6T7I_yeast_rRNAs.fa", package = "ARF"), use.names = T)
  human_rRNAs <- Biostrings::readBStringSet(file = system.file("extdata", "4V6X_human_rRNAs.fa", package = "ARF"), use.names = T)
  names(yeast_rRNAs) <- sapply(sapply(names(yeast_rRNAs),strsplit,split=" ",fixed=T),"[",1)
  names(human_rRNAs) <- sapply(sapply(names(human_rRNAs),strsplit,split=" ",fixed=T),"[",1)
  
  target_rRNAs <- Biostrings::readBStringSet(file = target_rRNAs_fasta, use.names = T)
  names(target_rRNAs) <- sapply(sapply(names(target_rRNAs), strsplit, split=" ", fixed=T),"[",1)
  
  # Create rRNA pairs list
  yeast_rRNA_pairs=list()
  human_rRNA_pairs=list()
  for (rRNAs in rRNA_pairs){
    if (rRNAs[1]=="28S") {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_25S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_28S",rRNA)))
    } else if (rRNAs[1]=="18S") {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_18S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_18S",rRNA)))
    } else if (rRNAs[1]=="5.8S") {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_5.8S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_5.8S",rRNA)))
    } else if (rRNAs[1]=="5S") {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_5S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_5S",rRNA)))
    }
  }
  rRNA_yeast2t <- sapply(yeast_rRNA_pairs,"[",2)
  names(rRNA_yeast2t) <- sapply(yeast_rRNA_pairs,"[",1)
  rRNA_human2t <- sapply(human_rRNA_pairs,"[",2)
  names(rRNA_human2t) <- sapply(human_rRNA_pairs,"[",1)
  
  # list of yeast transformation vectors
  yeast_2_t <- list()
  for (pair in yeast_rRNA_pairs) {
    rRNAs <- yeast_rRNAs[pair[1]]
    alignment <- strsplit(as.character(msa::msaClustalW(Biostrings::RNAStringSet(append(rRNAs,target_rRNAs[pair[2]])))@unmasked),
                          split="")
    # Read the alignment into S2T vector 
    yeast_2_t[[pair[1]]] <- rep(NA,length(alignment[[1]])) 
    spos <- 0
    tpos <- 0
    for (i in 1:length(alignment[[1]])){
      if (alignment[[2]][i]!="-") {tpos=tpos+1}
      if (alignment[[1]][i]!="-") {
        spos = spos+1
        if (alignment[[2]][i]!="-") {yeast_2_t[[pair[1]]][spos] <- tpos }
      }
    }
  }
  
  # list of human transformation vectors
  human_2_t <- list()
  for (pair in human_rRNA_pairs) {
    rRNAs <- human_rRNAs[pair[1]]
    alignment <- strsplit(as.character(msa::msaClustalW(Biostrings::RNAStringSet(append(rRNAs,target_rRNAs[pair[2]])))@unmasked),
                          split="")
    # Read the alignment into S2T vector 
    human_2_t[[pair[1]]] <- rep(NA,length(alignment[[1]])) 
    spos <- 0
    tpos <- 0
    for (i in 1:length(alignment[[1]])){
      if (alignment[[2]][i]!="-") {tpos=tpos+1}
      if (alignment[[1]][i]!="-") {
        spos = spos+1
        if (alignment[[2]][i]!="-") {human_2_t[[pair[1]]][spos] <- tpos }
      }
    }
  }
  
  # This is the target ribosome collision set positions
  target_positions <- list()
  
  # These are the source collision sets
  yeast_source_positions <- sapply(c("sc_6I7O_Col.Int.", "sc_6I7O_SAS", "sc_6T83_Col.Int.", "sc_6T83_SAS", "sc_6SV4_SAS"),
                               function(x)(return(ARF:::yeast_gsea_sets_Collision$gene[ARF:::yeast_gsea_sets_Collision$ont==x]))) 
  for (setid in names(yeast_source_positions)){
    posset <- yeast_source_positions[[setid]]
    target_positions[[setid]] <- c()
    
    rRNAs_in_set <- sapply(strsplit(posset,split = "_[0-9]*$"),"[[",1)
    pos_in_set <- unname(sapply(posset, FUN=function(x){l <- strsplit(x, split = "_",fixed = T)[[1]]; return(as.numeric(l[length(l)]))}))
    
    for(rRNA in unique(rRNAs_in_set)){
      source_residues <- pos_in_set[rRNAs_in_set==rRNA]
      target_residues <- yeast_2_t[[rRNA]][source_residues]
      target_residues <- target_residues[!is.na(target_residues)]
      if (length(target_residues)>0){
        target_positions[[setid]] <- append(target_positions[[setid]], paste(rRNA_s2t[rRNA],target_residues,sep = "_"))
      }
    }
  }
  
  human_source_positions <- sapply(c("hs_7QVP_Col.Int.", "hs_7QVP_SAS"),
                                   function(x)(return(ARF:::human_gsea_sets_Collision$gene[ARF:::human_gsea_sets_Collision$ont==x]))) 
  for (setid in names(human_source_positions)){
    posset <- human_source_positions[[setid]]
    target_positions[[setid]] <- c()
    
    rRNAs_in_set <- sapply(strsplit(posset,split = "_[0-9]*$"),"[[",1)
    pos_in_set <- unname(sapply(posset, FUN=function(x){l <- strsplit(x, split = "_",fixed = T)[[1]]; return(as.numeric(l[length(l)]))}))
    
    for(rRNA in unique(rRNAs_in_set)){
      source_residues <- pos_in_set[rRNAs_in_set==rRNA]
      target_residues <- human_2_t[[rRNA]][source_residues]
      target_residues <- target_residues[!is.na(target_residues)]
      if (length(target_residues)>0){
        target_positions[[setid]] <- append(target_positions[[setid]], paste(rRNA_s2t[rRNA],target_residues,sep = "_"))
      }
    }
  }
  
  target_positions[["Col.Int."]] <- unique(c(target_positions[["sc_6I7O_Col.Int."]], target_positions[["sc_6T83_Col.Int."]], target_positions[["hs_7QVP_Col.Int."]]))
  target_positions[["Rib.Col."]] <- unique(c(target_positions[["sc_6I7O_SAS"]], target_positions[["sc_6T83_SAS"]], target_positions[["sc_6SV4_SAS"]],
                                             target_positions[["hs_7QVP_SAS"]]))
  
  # borrowed from dripARF_get_RP_proximity_sets for adding the random sets
  gsea_sets_RP <- do.call("rbind", lapply(names(target_positions), FUN = function(RP){
    tmp_df<-NULL
    proxpos <- target_positions[[RP]]
    
    tmp_df <- data.frame(ont=RP,gene=paste(RP_proximity_df$rRNA[proxpos], RP_proximity_df$resno[proxpos], sep = "_"))
    rands <- c((1:100)*(round(dim(RP_proximity_df)[1]/100,digits = 0)-1))
    
    tmp_df <- rbind(tmp_df, as.data.frame(do.call("rbind", lapply(1:99, FUN = function(i){
      randset <- ((proxpos+rands[i]) %% (dim(RP_proximity_df)[1]))+1
      return(data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
                        gene=paste(RP_proximity_df$rRNA[randset], RP_proximity_df$resno[randset], sep = "_")))
    }))))
    return(tmp_df)
  }))
  
  return(gsea_sets_RP) # all collision sets including randoms
}
