
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


#' Draw dricARF result scatterplot for multiple comparisons
#' @description Draws a two-panel scatterplot of RPSEA NES z-score vs NES (panel 1) and vs
#'   \code{-log10(RPSEA.padj)} (panel 2). Ribosome-collision structural sets (SAS, Col.Int.,
#'   Rib.Col.) are highlighted with distinct colours, while RP heterogeneity candidates are shown
#'   in the background.
#' @param dricARF_results Full dricARF result data.frame as returned by \code{dricARF()}.
#' @param targetDir Directory to save the PDF plot in.
#' @param title Title string and base name for the output PDF (default:
#'   \code{"dricARF (highlighted) & dripARF predictions"}).
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
#' @return A \code{cowplot} grid object combining two panels. A PDF is also saved to
#'   \code{targetDir}.
#' @keywords dricARF Differential Ribosome Collisions scatterplot
#' @export
#' @examples
#' \dontrun{
#' dricARF_result_scatterplot(dricARF_results, "/Folder/to/save/in/")
#' }
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
     ggplot2::labs(col="Collision Prediction", shape=paste0("RPSEA.padj < ",as.character(RPSEA_adjP_thr)))+
     ggplot2::xlab("Enrichment Score 2\n(NES->Z-score within random sets)")+
     ggplot2::ylab("Enrichment Score 1 (RPSEA NES)"))
  
  print(g1)
  
  g2 <- (ggplot2::ggplot(RP_results, ggplot2::aes(x=RPSEA.NES_randZ, y=-log10(RPSEA.padj)))+
           ggplot2::geom_rect(data = data.frame(comp=unique(RP_results$comp)),inherit.aes = F, ggplot2::aes(xmin=1,xmax=Inf,ymin=-log10(RPSEA_adjP_thr),ymax=Inf), alpha=0.1, fill="green")+
           ggplot2::geom_rect(data = data.frame(comp=unique(RP_results$comp)),inherit.aes = F, ggplot2::aes(xmin=-Inf,xmax=1,ymin=-log10(RPSEA_adjP_thr),ymax=Inf), alpha=0.1, fill="red")+
           ggplot2::geom_rect(data = data.frame(comp=unique(RP_results$comp)),inherit.aes = F, ggplot2::aes(xmin=-Inf,xmax=1,ymin=-Inf,ymax=-log10(RPSEA_adjP_thr)), alpha=0.3, fill="firebrick")+
           ggplot2::geom_rect(data = data.frame(comp=unique(RP_results$comp)),inherit.aes = F, ggplot2::aes(xmin=1,xmax=Inf,ymin=-Inf,ymax=-log10(RPSEA_adjP_thr)), alpha=0.3, fill="darkolivegreen3")+
           ggplot2::geom_hline(yintercept = c(-log10(RPSEA_adjP_thr)), linetype="dashed", col=cols[5]) + 
           ggplot2::geom_vline(xintercept = c(randZscore_thr), linetype="dashed", col=cols[5]) + 
           ggplot2::facet_grid(comp~.)+
           ggplot2::geom_point(col=cols[2], size=1.5, alpha=.5, shape=18)+
           ggplot2::geom_point(data=final_coll, ggplot2::aes(color=Description), alpha=.5,size=1.5, shape=20)+
           ggplot2::geom_point(data=final_coll, ggplot2::aes(color=Description), size=2.5, shape=10)+
           ggrepel::geom_label_repel(data = final_coll%>%dplyr::filter(Description%in%c("Rib.Col.",addedRPs)), alpha=.8, size=2,
                                     inherit.aes = TRUE, ggplot2::aes(label=Description), color="#000000", show.legend = F)+
           ggplot2::scale_color_manual(values = c("#e41a1c",  "#fb9a99",        "#d95f02", "#1a9641",           "#1f78b4", "#fb2ae9",          "#ab1be7", "#000000", "#fb9a99", rep("#110134",length(addedRPs))),
                                       breaks = c( "hs_7QVP_SAS", "hs_7QVP_Col.Int.","sc_6I7O_SAS", "sc_6I7O_Col.Int.", "sc_6T83_SAS", "sc_6T83_Col.Int.", "sc_6SV4_SAS", "Rib.Col.","Col.Int.", addedRPs))+
           ggplot2::theme_bw()+
           ggplot2::labs(col="Collision Prediction")+
           ggplot2::xlab("Enrichment Score 2\n(NES->Z-score within random sets)")+
           ggplot2::ylab("-log10(RPSEA.padj)"))
  
  print(g2)
  
  gcombine <- cowplot::plot_grid(g1+ggplot2::theme(legend.position = "none"),g2,ncol=2,nrow=1,rel_widths = c(2,3))
  
  ggplot2::ggsave(plot = gcombine, filename = paste(targetDir, "/", title, ".pdf",sep=""),
                  width = 8, height = 3+(ceiling(length(unique(dricARF_results$comp)))*2),limitsize = FALSE)
  
  return(gcombine)
}

#' dricARF wrapper
#' @description Runs the complete dricARF pipeline, which extends dripARF by additionally
#'   integrating ribosome collision structural sets derived from cryo-EM collision structures.
#'   Specifically, it incorporates stalled-adjacent-subunit (SAS) and collision-interface
#'   (Col.Int.) sets from available human and yeast collision structures into the RPSEA enrichment
#'   analysis, enabling simultaneous prediction of ribosomal heterogeneity candidates and ribosome
#'   collision enrichments in a single run.
#' @param samplesFile Path to the tab-separated samples file. Required columns: \code{sampleName},
#'   \code{bedGraphFile} (or BAM path), \code{group}. See \code{read_ARF_samples_file()}.
#' @param rRNAs_fasta FASTA file for the rRNAs of the organism.
#' @param samples_df Pre-loaded samples data.frame from \code{read_ARF_samples_file()} (optional).
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
#' @param gsea_sets_RP RP-rRNA contact point sets for dripARF enrichments (preset for \code{hs},
#'   \code{mm}, and \code{sc}; provide for custom organisms via
#'   \code{dripARF_get_RP_proximity_sets()}).
#' @param RP_proximity_df RP-rRNA proximity matrix (preset for \code{hs}, \code{mm}, and \code{sc};
#'   provide for custom organisms via \code{ARF_parse_PDB_ribosome()}).
#' @param gsea_sets_Collision Ribosome collision GSEA sets for dricARF enrichments (preset for
#'   \code{hs} and \code{sc}; provide for other organisms via
#'   \code{dricARF_liftover_collision_sets()}).
#' @return A data.frame in the same format as returned by \code{dripARF_predict_heterogenity()},
#'   with additional rows for the collision set enrichments (SAS, Col.Int., Rib.Col., etc.). A
#'   scatterplot PDF is also saved to \code{targetDir}.
#' @seealso \code{\link{dripARF}}, \code{\link{dricARF_liftover_collision_sets}},
#'   \code{\link{dricARF_result_scatterplot}}
#' @keywords dricARF pipeline
#' @export
#' @examples
#' \dontrun{
#' dricARF("samples.txt", "rRNAs.fa", organism="hs", targetDir="/target/directory/to/save/results")
#' }
dricARF <- function(samplesFile, rRNAs_fasta, samples_df=NULL, organism=NULL, compare="group", QCplot=TRUE,  targetDir=NA,
                    comparisons=NULL, exclude=NULL, GSEAplots=FALSE, ssRPSEAplots=FALSE, gsea_sets_RP=NULL, RP_proximity_df=NULL, gsea_sets_Collision=NULL){
  
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
                               compare=compare, organism=organism, QCplot=QCplot, targetDir=targetDir,
                               comparisons = comparisons, GSEAplots=GSEAplots, ssRPSEAplots=ssRPSEAplots, 
                               gsea_sets_RP = gsea_sets_RP, RP_proximity_df = RP_proximity_df,
                               measureID = "abs_GSEA_measure_with_dynamic_p", runID="dricARF")
  
  dricARF_result_scatterplot(dricARF_results = results, targetDir = targetDir, title = "dricARF (highlighted) & dripARF predictions")
  # dripARF_result_heatmap(dripARF_results = results, targetDir = targetDir, title = "ALL dripARF predictions",
  #                        randZscore_thr = c(1), ORA_adjP_thr = c(0.05), RPSEA_adjP_thr = c(0.05), ORA_sig_n = 1)
  
  return(results)
}


#' Lift over ribosome collision sets to a target organism via rRNA alignment
#' @description Aligns human (4V6X) and yeast (6T7I) rRNA sequences to the target organism's rRNAs
#'   via pairwise ClustalW alignment, then converts the built-in collision GSEA sets (SAS, Col.Int.,
#'   Rib.Col., etc.) to target-organism rRNA coordinates. Also generates 99 circular-shift
#'   randomised control sets per collision set for use in RPSEA.
#' @param target_species ID string for the target species, e.g. \code{"dm"} (Drosophila) or any
#'   custom identifier.
#' @param target_rRNAs_fasta FASTA file for the rRNAs of the target organism. Should be the same
#'   file used in rRNA fragment alignment.
#' @param rRNA_pairs List of 2-element character vectors matching source rRNA IDs to target rRNA
#'   IDs. Use \code{"28S"}, \code{"18S"}, \code{"5.8S"}, \code{"5S"} as source IDs. Example:
#'   \code{list(c("28S","dm_rRNA_28S"), c("18S","dm_rRNA_18S"))}.
#' @return A data.frame with two columns (\code{ont}, \code{gene}) containing lifted-over collision
#'   sets (e.g. \code{Rib.Col.}, \code{Col.Int.}, \code{sc_6I7O_SAS}, etc.) and their 99
#'   randomised control sets, ready for use as the \code{gsea_sets_Collision} argument in
#'   \code{dricARF()}.
#' @seealso \code{\link{dricARF}}
#' @keywords dricARF collision liftover rRNA
#' @export
#' @examples
#' \dontrun{
#' dricARF_liftover_collision_sets("dm", "drosophila_rRNAs.fa",
#'   rRNA_pairs=list(c("28S","dm_rRNA_28S"), c("18S","dm_rRNA_18S"),
#'   c("5.8S","dm_rRNA_5.8S"), c("5S","dm_rRNA_5S")))
#' }
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
    rRNA=rRNAs[2]
    if (rRNAs[1]%in%c("28S","25S")) {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_25S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_28S",rRNA)))
    } else if (rRNAs[1]%in%c("18S","16S")) {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_18S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_18S",rRNA)))
    } else if (rRNAs[1]=="5.8S") {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_5.8S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_5.8S",rRNA)))
    } else if (rRNAs[1]=="5S") {
      yeast_rRNA_pairs <- append(yeast_rRNA_pairs, list(c("rRNA_5S",rRNA)))
      human_rRNA_pairs <- append(human_rRNA_pairs, list(c("human_5S",rRNA)))
    } else{
      message('Problem with passed rRNA pairs list. Please use 28S, 18S, 5.8S, 5S as source ids, i.e. list(c("28S","species_25S"), c("18S","species_18S"), etc.)\n')
        
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
