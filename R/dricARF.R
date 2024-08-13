
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
  ggplot2::ggsave(plot = g1, filename = paste(targetDir, "/", title, ".pdf",sep=""),
                  width = 5, height = 3+(ceiling(length(unique(dricARF_results$comp)))*2),limitsize = FALSE)
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
#' @param gsea_sets_RP RP-rRNA contact point sets to perform enrichments on.
#' @param RP_proximity_df RP-rRNA proximity matrix that is calculated by ARF.
#' @keywords dricARF pipeline
#' @export
#' @examples
#' dricARF("samples.txt", "rRNAs.fa", organism="mm", targetDir="/target/directory/to/save/results")
dricARF <- function(samplesFile, rRNAs_fasta, samples_df=NULL, organism=NULL, compare="group", QCplot=TRUE,  targetDir=NA,
                    comparisons=NULL, exclude=NULL, GSEAplots=FALSE){
  
  # Check organism first
  if (!ARF_check_organism(organism))
    return(NA)
  
  # ADD the collision related stuff in
  if (organism=="hs"){
    RP_proximity_df <- ARF:::RP_proximity_human_df
    added_sets <- unique(ARF:::human_gsea_sets_Collision$ont)
    added_sets <- added_sets[!grepl("Rand",added_sets)]
    for (colset in added_sets){
      RP_proximity_df[,colset] <- 100
      RP_proximity_df[ARF:::human_gsea_sets_Collision$gene[ARF:::human_gsea_sets_Collision$ont==colset], colset] <- 1
    }
    gsea_sets_RP <- rbind(ARF:::human_gsea_sets_RP, ARF:::human_gsea_sets_Collision)
  } else if (organism=="mm") {
    RP_proximity_df <- ARF:::RP_proximity_mouse_df
    added_sets <- unique(ARF:::mouse_gsea_sets_Collision$ont)
    added_sets <- added_sets[!grepl("Rand",added_sets)]
    for (colset in added_sets){
      RP_proximity_df[,colset] <- 100
      RP_proximity_df[ARF:::mouse_gsea_sets_Collision$gene[ARF:::mouse_gsea_sets_Collision$ont==colset], colset] <- 1
    }
    gsea_sets_RP <- rbind(ARF:::mouse_gsea_sets_RP, ARF:::mouse_gsea_sets_Collision)
  } else if (organism=="sc") {
    RP_proximity_df <- ARF:::RP_proximity_yeast_df
    added_sets <- unique(ARF:::yeast_gsea_sets_Collision$ont)
    added_sets <- added_sets[!grepl("Rand",added_sets)]
    for (colset in added_sets){
      RP_proximity_df[,colset] <- 100
      RP_proximity_df[ARF:::yeast_gsea_sets_Collision$gene[ARF:::yeast_gsea_sets_Collision$ont==colset], colset] <- 1
    }
    gsea_sets_RP <- rbind(ARF:::yeast_gsea_sets_RP, ARF:::yeast_gsea_sets_Collision)
  } else {
    message(paste(c("Organism", organism, "Not implemented yet! Please generate your own set of proximity matrix and RP-rRNA sets using ARF.\n"), 
                  collapse = " "))
    return(NULL)
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
