

nomenclature <- read.csv2(file = "data-raw/nomenclature.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(RP_proximity_df)[!colnames(RP_proximity_df)%in%nomenclature$humanST]
nomenclature$humanST[!nomenclature$humanST%in%colnames(RP_proximity_df)]

for (specie in c("human","mouse")){ #},"opossum","rhesus","chicken")){
  # HUMAN
  RP_proximity_df <- reshape2::dcast(
    read.csv(file = paste("/home/projects/ribosomal_heterogeneity/data/public/ribosome_structure/",specie,"_final_3D_distance_data.tsv",sep=""),
             header = TRUE,sep = ",",stringsAsFactors = FALSE), rRNA+resno~Species, value.var="mind",fun.aggregate = min)
  RP_proximity_df[RP_proximity_df == -Inf] <- NA
  RP_proximity_df[RP_proximity_df == Inf] <- NA
  rownames(RP_proximity_df) = paste(RP_proximity_df$rRNA, RP_proximity_df$resno, sep = "_")

  rrnas <- c("18S", "28S", "5.8S", "5S")

  poslists <- sapply(rrnas, FUN = function(x){return(which(sub(x = rownames(RP_proximity_df), pattern = "_[[:alnum:]]*$", replacement = "")==x))})
  rrna_lengths <- sapply(rrnas,function(x){return(length(poslists[[x]]))})

  colnames(RP_proximity_df) <- as.character(factor(colnames(RP_proximity_df),
                                                   levels=colnames(RP_proximity_df),
                                                   labels=unlist(sapply(colnames(RP_proximity_df), function(x){
                                                     if(x %in% nomenclature$humanST)
                                                       return(nomenclature$new[nomenclature$humanST==x][1])
                                                     else
                                                       return(x)
                                                   }))))

  gsea_sets_RP <- NULL
  RPfocus <- colnames(RP_proximity_df)[c(-1,-2)]
  RPfocus <- RPfocus[-1*which(startsWith(x = RPfocus,"ES") | startsWith(x = RPfocus,"yeast"))]


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

library(roxygen2); library(devtools); devtools::document(); devtools::install()

