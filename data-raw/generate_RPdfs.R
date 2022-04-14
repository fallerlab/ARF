
# Copyright (C) 2021  Ferhat Alkan
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


nomenclature <- read.csv2(file = "data-raw/nomenclature.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(RP_proximity_df)[!colnames(RP_proximity_df)%in%nomenclature$humanST]
nomenclature$humanST[!nomenclature$humanST%in%colnames(RP_proximity_df)]

for (specie in c("human","mouse")){
  # HUMAN
  RP_proximity_df <- reshape2::dcast(
    read.csv(file = paste("data-raw/",specie,"_final_3D_distance_data.tsv",sep=""),
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

  # gsea_sets_RP <- NULL
  # RPfocus <- colnames(RP_proximity_df)[c(-1,-2)]
  # RPfocus <- RPfocus[-1*which(startsWith(x = RPfocus,"ES") | startsWith(x = RPfocus,"yeast"))]
  # 
  # 
  # RP_proxpos <- list()
  # for (RP in RPfocus){
  #   proxpos <- which(RP_proximity_df[,RP]<25)
  # 
  #   if (length(proxpos)>250 & !(startsWith(RP,"ES"))){
  #     RP_proxpos[[RP]] <- sort(proxpos[order(RP_proximity_df[proxpos,RP])[1:250]])
  #   } else if (length(proxpos)>0){
  #     RP_proxpos[[RP]] <- proxpos
  #   }
  # }
  # 
  # for (RP in names(RP_proxpos)){
  #   print(paste(RP,specie))
  #   proxpos <- RP_proxpos[[RP]]
  # 
  #   if(length(proxpos)>0) {
  #     gsea_sets_RP <- rbind(gsea_sets_RP,
  #                                 data.frame(ont=RP,gene=paste(RP_proximity_df$rRNA[proxpos], RP_proximity_df$resno[proxpos], sep = "_")))
  # 
  #     rands <- c((1:100)*(round(dim(RP_proximity_df)[1]/100,digits = 0)-1))
  #     for (i in 1:99) {
  #       randset <- ((proxpos+rands[i]) %% (dim(RP_proximity_df)[1]))+1
  #       gsea_sets_RP <- rbind(gsea_sets_RP, data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
  #                                                                  gene=paste(RP_proximity_df$rRNA[randset],
  #                                                                             RP_proximity_df$resno[randset], sep = "_")))
  #     }
  #   }
  # }
  
  gsea_sets_RP <- NULL
  RPfocus <- colnames(RP_proximity_df)[c(-1,-2)]
  RPfocus <- RPfocus[-1*which(startsWith(x = RPfocus,"ES") | startsWith(x = RPfocus,"yeast") | startsWith(x = RPfocus,"YCI"))]
  
  thr<-c(27.40514, 360)
  
  RP_proxpos <- list()
  for (RP in RPfocus){
    proxpos <- which(RP_proximity_df[,RP]<thr[1])
    
    if (length(proxpos)>thr[2] & !(startsWith(RP,"ES"))){
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

  assign(paste("RP_proximity_",specie,"_df",sep = ""), RP_proximity_df)
  assign(paste(specie,"_gsea_sets_RP",sep = ""), gsea_sets_RP)
}

#################################
df_4v6x <- read.csv("data-raw/4v6x_chain_data.edited.csv", sep=",")
get_chain_id <- function(x){
  temp <- human_gsea_sets_RP$gene[human_gsea_sets_RP$ont==x]
  
}
  

#################################

usethis::use_data(RP_proximity_mouse_df, RP_proximity_human_df, #RP_proximity_chicken_df, RP_proximity_opossum_df, RP_proximity_rhesus_df,
                  mouse_gsea_sets_RP, human_gsea_sets_RP, #chicken_gsea_sets_RP, opossum_gsea_sets_RP, rhesus_gsea_sets_RP,
                  internal = TRUE, overwrite = TRUE)

library(roxygen2); library(devtools); devtools::document(); devtools::install()



