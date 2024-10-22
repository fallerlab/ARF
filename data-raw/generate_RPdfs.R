
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

usethis::use_data(RP_proximity_mouse_df, RP_proximity_human_df, 
                  mouse_gsea_sets_RP, human_gsea_sets_RP, 
                  internal = TRUE, overwrite = TRUE)

### Above is old (before dricARF) ###
#####################################

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


############################################################
RP_proximity_human_df<-(ARF:::RP_proximity_human_df)
RP_proximity_mouse_df<-(ARF:::RP_proximity_mouse_df)
human_gsea_sets_RP<-(ARF:::human_gsea_sets_RP)
mouse_gsea_sets_RP<-(ARF:::mouse_gsea_sets_RP)

human_ARF_ribo <- read.table("data-raw/Ribosome.3D.4V6X.ARF.minDist_matrix.csv",
                             sep = ",",stringsAsFactors = F,header = 1)
mouse_ARF_ribo <- read.table("data-raw/Ribosome.3D.4V6X.ARF.minDist_matrix.mm.converted.csv",
                             sep = ",",stringsAsFactors = F,header = 1)
yeast_ARF_prox_df <- ARF_parse_PDB_ribosome(species = "sc",PDBid = "6T7I",download_directory="data-raw/",
                                            PDB_file = "data-raw/6T7I.cif",out_prefix = "data-raw/Ribosome.3D")
yeast_ARF_ribo <- read.table("data-raw/Ribosome.3D.6T7I.ARF.minDist_matrix.csv",
                             sep = ",",stringsAsFactors = F,header = 1)

row.names(human_ARF_ribo) <- paste(human_ARF_ribo$rRNA,human_ARF_ribo$resno,sep="_")
row.names(mouse_ARF_ribo) <- paste(mouse_ARF_ribo$rRNA,mouse_ARF_ribo$resno,sep="_")
row.names(yeast_ARF_ribo) <- paste(yeast_ARF_ribo$rRNA,yeast_ARF_ribo$resno,sep="_")

selected_col_sets <- c("human_7qvp_disome_noSAS_DistanceThr05_ext.1", "yeast_6i7o_disome_noSAS_DistanceThr05_ext.1", 
                       "yeast_6t83_disome_noSAS_DistanceThr05_ext.1", "both_4struct_isome_noSAS_DistanceThr05_ext.1",
                       "human_7qvp_disome_SAS_10_maxmin5_ext.1", "both_4struct_isome_SAS_10_maxmin5_ext.1",
                       "yeast_6i7o_disome_SAS_10_maxmin5_ext.1","yeast_6t83_disome_SAS_10_maxmin5_ext.1","yeast_6sv4_trisome_SAS_10_maxmin5_ext.1")
names(selected_col_sets) <- c("hs_7QVP_Col.Int.","sc_6I7O_Col.Int.", "sc_6T83_Col.Int.","Col.Int.", "hs_7QVP_SAS", "Rib.Col.", "sc_6I7O_SAS", "sc_6T83_SAS", "sc_6SV4_SAS")

human_collision_sets <- read.table("/home/ferro/GITspace/rRNA_distances/human.SASA.collision.sets.csv",
                                   sep="\t", header = 1,stringsAsFactors = F)
for (col_set in names(selected_col_sets)){
  human_ARF_ribo[,col_set] <- 100
  human_ARF_ribo[human_collision_sets$gene[human_collision_sets$ont==selected_col_sets[col_set]], col_set] <- 0
}

mouse_collision_sets <- read.table("/home/ferro/GITspace/rRNA_distances/mouse.SASA.collision.sets.csv",
                                   sep="\t", header = 1,stringsAsFactors = F)
for (col_set in names(selected_col_sets)){
  mouse_ARF_ribo[,col_set] <- 100
  mouse_ARF_ribo[mouse_collision_sets$gene[mouse_collision_sets$ont==selected_col_sets[col_set]], col_set] <- 0
}

yeast_collision_sets <- read.table("/home/ferro/GITspace/rRNA_distances/yeast.SASA.collision.sets.csv", 
                                   sep="\t", header = 1,stringsAsFactors = F)
for (col_set in names(selected_col_sets)){
  yeast_ARF_ribo[,col_set] <- 100
  yeast_ARF_ribo[yeast_collision_sets$gene[yeast_collision_sets$ont==selected_col_sets[col_set]], col_set] <- 0
}

human_RP_sets <- dripARF_get_RP_proximity_sets(RP_proximity_df = human_ARF_ribo, additional_RPcols = 85:(dim(human_ARF_ribo)[2]),
                                               rRNAs_fasta = "/home/ferro/GITspace/rRNA_distances/PDB/4V6X.rRNAs.fasta")
mouse_RP_sets <- dripARF_get_RP_proximity_sets(RP_proximity_df = mouse_ARF_ribo, additional_RPcols = 85:(dim(mouse_ARF_ribo)[2]),
                                               rRNAs_fasta = "/home/projects/ribosomal_heterogeneity/data/public/mouse_rRNAs.wU.fa")
yeast_RP_sets <- dripARF_get_RP_proximity_sets(RP_proximity_df = yeast_ARF_ribo, additional_RPcols = 55:(dim(yeast_ARF_ribo)[2]),
                                               rRNAs_fasta = "data-raw/6T7I.rRNAs.fasta")

RP_proximity_yeast_df <- yeast_ARF_ribo[,-76:-84]
RP_proximity_mouse_df <- mouse_ARF_ribo[,-85:-93]
RP_proximity_human_df <- human_ARF_ribo[,-85:-93]
rownames(RP_proximity_human_df) <- gsub(pattern = "rRNA",replacement = "human",rownames(RP_proximity_human_df))
human_RP_sets$gene <- gsub(pattern = "rRNA",replacement = "human",human_RP_sets$gene)

human_gsea_sets_Collision <- human_RP_sets[grepl("SAS|Rib.Col.|Col.Int.",human_RP_sets$ont),]
unique(gsub(pattern = "Rand[0-9]*_", replacement = "",unique(human_gsea_sets_Collision$ont)))
mouse_gsea_sets_Collision <- mouse_RP_sets[grepl("SAS|Rib.Col.|Col.Int.",mouse_RP_sets$ont),]
unique(gsub(pattern = "Rand[0-9]*_", replacement = "",unique(mouse_gsea_sets_Collision$ont)))
yeast_gsea_sets_Collision <- yeast_RP_sets[grepl("SAS|Rib.Col.|Col.Int.",yeast_RP_sets$ont),]
unique(gsub(pattern = "Rand[0-9]*_", replacement = "",unique(yeast_gsea_sets_Collision$ont)))

human_gsea_sets_RP <- human_RP_sets[!human_RP_sets$ont%in%human_gsea_sets_Collision$ont,]
unique(gsub(pattern = "Rand[0-9]*_", replacement = "",unique(human_gsea_sets_RP$ont)))
mouse_gsea_sets_RP <- mouse_RP_sets[!mouse_RP_sets$ont%in%mouse_gsea_sets_Collision$ont,]
unique(gsub(pattern = "Rand[0-9]*_", replacement = "",unique(mouse_gsea_sets_RP$ont)))
yeast_gsea_sets_RP <- yeast_RP_sets[!yeast_RP_sets$ont%in%yeast_gsea_sets_Collision$ont,]
unique(gsub(pattern = "Rand[0-9]*_", replacement = "",unique(yeast_gsea_sets_RP$ont)))

###########

usethis::use_data(RP_proximity_mouse_df, RP_proximity_human_df, RP_proximity_yeast_df, 
                  mouse_gsea_sets_RP, human_gsea_sets_RP, yeast_gsea_sets_RP,
                  human_gsea_sets_Collision, mouse_gsea_sets_Collision, yeast_gsea_sets_Collision,
                  RP_nomenclature_map,
                  internal = TRUE, overwrite = TRUE)

library(roxygen2); library(devtools); devtools::document(); devtools::install()


results <- ARF::dripARF(samplesFile = "test_data/Rpl10aRps25Rpl22_samples_mouse.tsv", rRNAs_fasta = "rRNAs/mouse_rRNAs.fa", organism = "mm", QCplot = F, targetDir = "test/",
             comparisons = list(c("uL1_FLAG","eL22_HA"),c("eS25_FLAG","eL22_HA"),c("uL1_FLAG","uL1_Total"),c("eS25_FLAG","eS25_Total")))
results2 <- ARF::dricARF(samplesFile = "test_data/Rpl10aRps25Rpl22_samples_mouse.tsv", rRNAs_fasta = "rRNAs/mouse_rRNAs.fa", organism = "mm", QCplot = F, targetDir = "test/",
                        comparisons = list(c("uL1_FLAG","eL22_HA"),c("eS25_FLAG","eL22_HA"),c("uL1_FLAG","uL1_Total"),c("eS25_FLAG","eS25_Total")))

results3 <- ARF::dricARF(samplesFile = "test_data/SRP043036_Lareau_Brown_2014.tsv", rRNAs_fasta = "rRNAs/6T7I_yeast_rRNAs.fa", organism = "sc", QCplot = F, targetDir = "test/")

