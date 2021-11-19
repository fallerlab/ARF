

for (specie in c("human","mouse")){ #},"opossum","rhesus","chicken")){
  # HUMAN
  RP_proximity_df <- reshape2::dcast(
    read.csv(file = paste("/home/projects/ribosomal_heterogeneity/data/public/ribosome_structure/",specie,"_final_3D_distance_data.tsv",sep=""),
             header = TRUE,sep = ",",stringsAsFactors = FALSE), rRNA+resno~Species, value.var="mind",fun.aggregate = min)
  RP_proximity_df[RP_proximity_df == -Inf] <- NA
  RP_proximity_df[RP_proximity_df == Inf] <- NA
  rownames(RP_proximity_df) = paste(RP_proximity_df$rRNA, RP_proximity_df$resno, sep = "_")

  gsea_sets_RP <- NULL
  RPfocus <- colnames(RP_proximity_df)[c(-1,-2)]
  RPfocus <- RPfocus[-1*which(startsWith(x = RPfocus,"ES") | startsWith(x = RPfocus,"yeast"))]

  rrnas <- c("18S", "28S", "5.8S", "5S")

  poslists <- sapply(rrnas, FUN = function(x){return(which(sub(x = rownames(RP_proximity_df), pattern = "_[[:alnum:]]*$", replacement = "")==x))})
  rrna_lengths <- sapply(rrnas,function(x){return(length(poslists[[x]]))})

  RP_proxpos <- list()
  for (RP in RPfocus){
    proxpos <- which(RP_proximity_df[,RP]<25)

    # RP_proxpos[[RP]] <- proxpos
    # if (length(proxpos)>250 & !(startsWith(RP,"ES")))
    #   RP_proxpos[[paste0(RP,"_shrinked")]] <- sort(proxpos[order(RP_proximity_df[proxpos,RP])[1:250]])
    if (length(proxpos)>250 & !(startsWith(RP,"ES"))){
      RP_proxpos[[RP]] <- sort(proxpos[order(RP_proximity_df[proxpos,RP])[1:250]])
    } else if (length(proxpos)>0){
      RP_proxpos[[RP]] <- proxpos
    }

  }

  for (RP in names(RP_proxpos)){
    print(paste(RP,specie))
    proxpos <- RP_proxpos[[RP]]

    if(length(proxpos)>0) {
      gsea_sets_RP <- rbind(gsea_sets_RP,
                                  data.frame(ont=RP,gene=paste(RP_proximity_df$rRNA[proxpos], RP_proximity_df$resno[proxpos], sep = "_")))

      # nonnas <- sum(!is.na(RP_proximity_df[,RP]))
      # midrange <- which(RP_proximity_df[,RP] %in%
      #                     (sort(RP_proximity_df[,RP])[(round((nonnas/2))-round((length(proxpos)/2))):(round((nonnas/2))+round((length(proxpos)/2)))]))
      # gsea_sets_RP <- rbind(gsea_sets_RP, data.frame(ont=paste("MRf",RP,sep="_"),
      #                                                            gene=paste(RP_proximity_df$rRNA[midrange],
      #                                                                       RP_proximity_df$resno[midrange], sep = "_")))

      rands <- c((1:100)*(round(dim(RP_proximity_df)[1]/100,digits = 0)-1))
      for (i in 1:99) {
        randset <- ((proxpos+rands[i]) %% (dim(RP_proximity_df)[1]))+1
        gsea_sets_RP <- rbind(gsea_sets_RP, data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
                                                                   gene=paste(RP_proximity_df$rRNA[randset],
                                                                              RP_proximity_df$resno[randset], sep = "_")))
      }
    }
  }
  assign(paste("RP_proximity_",specie,"_df",sep = ""), RP_proximity_df)
  assign(paste(specie,"_gsea_sets_RP",sep = ""), gsea_sets_RP)
}


#
# # MOUSE
# RP_proximity_mouse_df <- read.csv(file = "/home/projects/ribosomal_heterogeneity/Rspace/mouse_RP.tsv",
#                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# rownames(RP_proximity_mouse_df) = paste(RP_proximity_mouse_df$rid, RP_proximity_mouse_df$resno, sep = "_")
#
# mouse_gsea_sets_RP <- NULL
# rrnas <- c("NR_003278.3","NR_003279.1","NR_030686.1","NR_003280.2")
#
# RPfocus <- colnames(RP_proximity_mouse_df)[3:84]
# poslists <- sapply(rrnas, FUN = function(x){return(which(sub(x = rownames(RP_proximity_mouse_df), pattern = "_[[:alnum:]]*$", replacement = "")==x))})
# rrna_lengths <- sapply(rrnas,function(x){return(length(poslists[[x]]))})
#
# #rrna shifts for "rand" controls
# rrna_rands <- list()
# for (rrna in rrnas) {
#   rrna_rands[[rrna]] <- c()
#   for (rand in rands){
#     if (rand < (rrna_lengths[rrna]/2))
#       rrna_rands[[rrna]] <- c(rrna_rands[[rrna]],rand,-rand)
#   }
# }
#
#
# for (RP in RPfocus){
#   print(RP)
#   proxpos <- which(RP_proximity_mouse_df[,RP]<25)
#   if(length(proxpos)>0) {
#     mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=RP,
#                                 gene=paste(RP_proximity_mouse_df$rid[proxpos], RP_proximity_mouse_df$resno[proxpos], sep = "_")))
#
#     # farpos <- which(RP_proximity_mouse_df[,RP]>= sort(RP_proximity_mouse_df[,RP],decreasing = TRUE)[length(proxpos)])
#     # mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=paste("FDf",RP,sep="_"),
#     #                             gene=paste(RP_proximity_mouse_df$rid[farpos],
#     #                                        RP_proximity_mouse_df$resno[farpos], sep = "_")))
#
#     nonnas <- sum(!is.na(RP_proximity_mouse_df[,RP]))
#     midrange <- which(RP_proximity_mouse_df[,RP] %in%
#                         (sort(RP_proximity_mouse_df[,RP])[(round((nonnas/2))-round((length(proxpos)/2))):(round((nonnas/2))+round((length(proxpos)/2)))]))
#     mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=paste("MRf",RP,sep="_"),
#                                 gene=paste(RP_proximity_mouse_df$rid[midrange],
#                                            RP_proximity_mouse_df$resno[midrange], sep = "_")))
#
#     rand_focus_rrnas <- unique(na.omit(sapply(rrnas, FUN=function(x){if (sum(proxpos%in%poslists[[x]])>0) return(x); return(NA)})))
#     for (i in 1:max(sapply(rand_focus_rrnas,FUN = function(x){return(length(rrna_rands[[x]]))}))){
#       randset <- unlist(sapply(rand_focus_rrnas, FUN = function(x){
#         ni = i
#         randl = length(rrna_rands[[x]])
#         if (i > randl)
#           ni = i%%randl
#
#         rs_pp <- proxpos[proxpos%in%poslists[[x]]] # rRNAspecific proxposs
#         return( min(rs_pp) + ( ( (rs_pp-min(rs_pp)) + rrna_rands[[x]][ni]) %% rrna_lengths[[x]] ) )
#         }))
#       mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
#                                                                  gene=paste(RP_proximity_mouse_df$rid[randset],
#                                                                             RP_proximity_mouse_df$resno[randset], sep = "_")))
#     }
#   }
# }

usethis::use_data(RP_proximity_mouse_df, RP_proximity_human_df, #RP_proximity_chicken_df, RP_proximity_opossum_df, RP_proximity_rhesus_df,
                  mouse_gsea_sets_RP, human_gsea_sets_RP, #chicken_gsea_sets_RP, opossum_gsea_sets_RP, rhesus_gsea_sets_RP,
                  internal = TRUE, overwrite = TRUE)

## TEST ##
load("R/sysdata.rda")
targetDir=paste0(getwd(),"/test/")
compare="group"
organism="mm"
QCplot=TRUE
comparisons=NULL

samples<-read_ARF_samples_file("/home/projects/ribosomal_heterogeneity/RP4_space/benchmark/SRP064202_Barna_Rpl10aRps25Rpl22_samples_mouse.tsv")
rRNA_counts <- DRIP_ARF_read_rRNA_fragments(samples, organism=organism, QCplot=QCplot, targetDir=targetDir)
org_RP_df <- RP_proximity_mouse_df
gsea_sets_RP <- mouse_gsea_sets_RP

DRIPARF_results <- DRIP_ARF_predict_heterogenity(samples,rRNA_counts,organism = "mm")
DRIP_ARF_result_heatmap(DRIPARF_results,"test",targetDir)
DRIPARF_rpspec <- DRIP_ARF_report_RPspec_pos_results(samplesFile = "/home/projects/ribosomal_heterogeneity/RP4_space/benchmark/SRP064202_Barna_Rpl10aRps25Rpl22_samples_mouse.tsv",
                                   organism = "mm",savefolder = targetDir)

DRIP_ARF_result_df <- DRIPARF_results
DRIP_ARF_DRF <- DRIPARF_rpspec
simpleresults <- DRIP_ARF_simplify_results(DRIPARF_results)
RPs<-unique(simpleresults$Description)
RPs<-RPs[startsWith(RPs,"RP")]
title="rRNA_pos_spec_heatmap"
randZscore_thr=1
ORA_adjP_thr=0.05
GSEA_adjP_thr=0.05
ORA_sig_n=-1

###
pca_df$figname <- factor(pca_df$sample,
  levels =  c("Rpl10a_rep1", "Rpl10a_rep2", "Rpl22_rep1", "Rpl22_rep2", "Rps25_rep1", "Rps25_rep2",
              "TotalRpl10a_rep1", "TotalRpl10a_rep2", "TotalRpl22_rep1", "TotalRpl22_rep2", "TotalRps25_rep1","TotalRps25_rep2"),
  labels = c("PD_Rpl10a_FLAG_1", "PD_Rpl10a_FLAG_2", "PD_Rpl22_HA_1", "PD_Rpl22_HA_2", "PD_Rps25_FLAG_1", "PD_Rps25_FLAG_2",
             "Total_Rpl10a_FLAG_1", "Total_Rpl10a_FLAG_2", "Total_Rpl22_HA_1", "Total_Rpl22_HA_2", "Total_Rps25_FLAG_1","Total_Rps25_FLAG_1"))
g1 <- ggplot2::ggplot(pca_df, ggplot2::aes(x=PC1,y=PC2,color=group)) + ggplot2::geom_point(show.legend = FALSE) +
  ggplot2::labs(x=pca_map[1], y=pca_map[2]) + ggrepel::geom_text_repel(ggplot2::aes(label=figname), show.legend = FALSE)
print(g1)
ggplot2::ggsave(filename = paste(targetDir, "/", "vsd_based_PCA.png", sep = ""), plot = g1,limitsize = FALSE,width = 4,height = 4)

tsne <- Rtsne::Rtsne(t(SummarizedExperiment::assay(vsd)[order(matrixStats::rowVars(SummarizedExperiment::assay(vsd)),decreasing = TRUE),]),check_duplicates = FALSE,perplexity = 1)
tsne_df <- cbind(data.frame(tsne$Y), colnames(vsd), sapply(colnames(vsd),FUN=function(x){return(pca_df$figname[pca_df$sample==x])}))
colnames(tsne_df)<-c("X","Y","sample","group")
g1<-ggplot2::ggplot(tsne_df,ggplot2::aes(x=X,y=Y,color=group))+ggplot2::geom_point(show.legend = FALSE)+
  ggrepel::geom_label_repel(ggplot2::aes(label=sample),show.legend = FALSE)+
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),axis.title.y = ggplot2::element_blank())
ggplot2::ggsave(filename = paste(targetDir, "/", "vsd_based_TSNE.png", sep = ""), plot = g1,limitsize = FALSE,width = 5,height = 5)

library(roxygen2); library(devtools); devtools::document(); devtools::install()

