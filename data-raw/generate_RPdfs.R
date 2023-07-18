
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
  print(specie)
  
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
    print(RP)
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

####################
ribosome_PDBs <- list(tomato=c("7QIZ"), wheat=c("4V7E"),
                      Staphylococcus_aureus=c("6S0X","6S13"),
                      Bacillus_subtilis=c("3J9W"),
                      ecoli=c("4V9D","7UG7","7OTC","7TOS","6XZA","6XZB","7K00"),
                      thermophilus=c("6QNQ","6QNR"),
                      TETRAHYMENA=c("4V8P"),
                      human=c("6QZP", "6Y0G", "6Y2L", "6Y57","5AJ0","7F5S","7TQL","4V6X"),
                      rabbit=c("7QWQ","7QWS","7QWR","7ZJX","7ZJW","7UCK","7UCJ","7TOR","7TOQ"),
                      yeast=c("5M1J","4V88","4V7R","7NRD","7NRC","6SNT","7TOP","7TOO","6T4Q","6T7I","6T7T"),
                      drosophila=c("4V6W","6XU6"))

ALL <- NULL
for (PDB_id in unlist(ribosome_PDBs)){
  PDB_file <- paste0(download_directory,"/",PDB_id,".cif")
  utils::download.file(paste0("https://files.rcsb.org/download/",PDB_id,".cif"),PDB_file)
  utils::download.file(paste0("https://www.rcsb.org/fasta/entry/",PDB_id),paste0(download_directory,"/",PDB_id,".fasta"))
  seqs <- data.frame(do.call(rbind, strsplit(names(seqinr::read.fasta(file = paste0(download_directory,"/",PDB_id,".fasta"),
                                                                      whole.header = T)), 
                                             split = "|",fixed=T)),PDB_id=PDB_id,stringsAsFactors = F)
  colnames(seqs) <- c("ID","chain_info","RP","organism","PDB_id")
  ALL <- rbind(ALL,seqs)
}
write.csv(ALL,"data-raw/ALL_PDB_RPS.csv",row.names = F)
# Do some editing on this doc then read back

RP_nomenclature_map <- read.csv(file = "data-raw/ALL_PDB_RPS_edited.csv",stringsAsFactors = F)
RP_nomenclature_map$chainID <- sapply(RP_nomenclature_map$chain_info,
            FUN = function(x){
              split_out <- c()
              splittedx <- strsplit(x,split = ",",fixed = T)[[1]]
              for (split in splittedx){
                if(grepl(pattern = "auth", x = split, fixed = T)){
                  splitted <- strsplit(gsub(pattern = "\\[|\\]",replacement = " ",split), split = " ", fixed=T)[[1]]
                  newx <- splitted[which(splitted=="auth")+1]
                }else{
                  newx <- gsub(pattern = "Chains*\\s*", replacement = '', split)
                }
                split_out <- c(split_out,newx)
              }
              return(paste(split_out,collapse = ","))
})

#################################

usethis::use_data(RP_proximity_mouse_df, RP_proximity_human_df,  #RP_proximity_yeast_df, #RP_proximity_chicken_df, RP_proximity_opossum_df, RP_proximity_rhesus_df,
                  mouse_gsea_sets_RP, human_gsea_sets_RP, #yeast_gsea_sets_RP, #chicken_gsea_sets_RP, opossum_gsea_sets_RP, rhesus_gsea_sets_RP,
                  RP_nomenclature_map,
                  internal = TRUE, overwrite = TRUE)

library(roxygen2); library(devtools); devtools::document(); devtools::install()


