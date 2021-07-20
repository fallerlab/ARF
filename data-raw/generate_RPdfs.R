

rands <- c((1:100)*25)

# HUMAN
RP_proximity_human_df <- reshape2::dcast(read.csv(file = "/home/projects/ribosomal_heterogeneity/Rspace/human_ribosome_structure_summary.tsv",
                                                  header = TRUE,sep = ",",stringsAsFactors = FALSE),
                                         rRNA+resno~Species, value.var="mind",fun.aggregate = min)
RP_proximity_human_df[RP_proximity_human_df == -Inf] <- NA
rownames(RP_proximity_human_df) = paste(RP_proximity_human_df$rRNA, RP_proximity_human_df$resno, sep = "_")

human_gsea_sets_RP <- NULL
RPfocus <- colnames(RP_proximity_human_df)[3:84]
rrnas <- c("human_18S", "human_28S", "human_5.8S", "human_5S")
poslists <- sapply(rrnas, FUN = function(x){return(which(sub(x = rownames(RP_proximity_human_df), pattern = "_[[:alnum:]]*$", replacement = "")==x))})
rrna_lengths <- sapply(rrnas,function(x){return(length(poslists[[x]]))})

#rrna shifts for "rand" controls
rrna_rands <- list()
for (rrna in rrnas) {
  rrna_rands[[rrna]] <- c()
  for (rand in rands){
    if (rand < (rrna_lengths[rrna]/2))
      rrna_rands[[rrna]] <- c(rrna_rands[[rrna]],rand,-rand)
  }
}


for (RP in RPfocus){
  print(RP)
  proxpos <- which(RP_proximity_human_df[,RP]<25)
  if(length(proxpos)>0) {
    human_gsea_sets_RP <- rbind(human_gsea_sets_RP,
                                data.frame(ont=RP,gene=paste(RP_proximity_human_df$rRNA[proxpos], RP_proximity_human_df$resno[proxpos], sep = "_")))

    # farpos <- which(RP_proximity_human_df[,RP]>= sort(RP_proximity_human_df[,RP],decreasing = TRUE)[length(proxpos)])
    # human_gsea_sets_RP <- rbind(human_gsea_sets_RP, data.frame(ont=paste("FDf",RP,sep="_"),
    #                                                            gene=paste(RP_proximity_human_df$rRNA[farpos],
    #                                                                       RP_proximity_human_df$resno[farpos], sep = "_")))

    nonnas <- sum(!is.na(RP_proximity_human_df[,RP]))
    midrange <- which(RP_proximity_human_df[,RP] %in%
                        (sort(RP_proximity_human_df[,RP])[(round((nonnas/2))-round((length(proxpos)/2))):(round((nonnas/2))+round((length(proxpos)/2)))]))
    human_gsea_sets_RP <- rbind(human_gsea_sets_RP, data.frame(ont=paste("MRf",RP,sep="_"),
                                                               gene=paste(RP_proximity_human_df$rRNA[midrange],
                                                                          RP_proximity_human_df$resno[midrange], sep = "_")))


    rand_focus_rrnas <- unique(na.omit(sapply(rrnas, FUN=function(x){if (sum(proxpos%in%poslists[[x]])>0) return(x); return(NA)})))
    for (i in 1:max(sapply(rand_focus_rrnas,FUN = function(x){return(length(rrna_rands[[x]]))}))){
      randset <- unlist(sapply(rand_focus_rrnas, FUN = function(x){
        ni = i
        randl = length(rrna_rands[[x]])
        if (i > randl)
          ni = i%%randl

        rs_pp <- proxpos[proxpos%in%poslists[[x]]] # rRNAspecific proxposs
        return( min(rs_pp) + ( ( (rs_pp-min(rs_pp)) + rrna_rands[[x]][ni]) %% rrna_lengths[[x]] ) )
      }))
      human_gsea_sets_RP <- rbind(human_gsea_sets_RP, data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
                                                                 gene=paste(RP_proximity_human_df$rid[randset],
                                                                            RP_proximity_human_df$resno[randset], sep = "_")))
    }
  }
}

# MOUSE
RP_proximity_mouse_df <- read.csv(file = "/home/projects/ribosomal_heterogeneity/Rspace/mouse_RP.tsv",
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(RP_proximity_mouse_df) = paste(RP_proximity_mouse_df$rid, RP_proximity_mouse_df$resno, sep = "_")

mouse_gsea_sets_RP <- NULL
rrnas <- c("NR_003278.3","NR_003279.1","NR_030686.1","NR_003280.2")

RPfocus <- colnames(RP_proximity_mouse_df)[3:84]
poslists <- sapply(rrnas, FUN = function(x){return(which(sub(x = rownames(RP_proximity_mouse_df), pattern = "_[[:alnum:]]*$", replacement = "")==x))})
rrna_lengths <- sapply(rrnas,function(x){return(length(poslists[[x]]))})

#rrna shifts for "rand" controls
rrna_rands <- list()
for (rrna in rrnas) {
  rrna_rands[[rrna]] <- c()
  for (rand in rands){
    if (rand < (rrna_lengths[rrna]/2))
      rrna_rands[[rrna]] <- c(rrna_rands[[rrna]],rand,-rand)
  }
}


for (RP in RPfocus){
  print(RP)
  proxpos <- which(RP_proximity_mouse_df[,RP]<25)
  if(length(proxpos)>0) {
    mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=RP,
                                gene=paste(RP_proximity_mouse_df$rid[proxpos], RP_proximity_mouse_df$resno[proxpos], sep = "_")))

    # farpos <- which(RP_proximity_mouse_df[,RP]>= sort(RP_proximity_mouse_df[,RP],decreasing = TRUE)[length(proxpos)])
    # mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=paste("FDf",RP,sep="_"),
    #                             gene=paste(RP_proximity_mouse_df$rid[farpos],
    #                                        RP_proximity_mouse_df$resno[farpos], sep = "_")))

    nonnas <- sum(!is.na(RP_proximity_mouse_df[,RP]))
    midrange <- which(RP_proximity_mouse_df[,RP] %in%
                        (sort(RP_proximity_mouse_df[,RP])[(round((nonnas/2))-round((length(proxpos)/2))):(round((nonnas/2))+round((length(proxpos)/2)))]))
    mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=paste("MRf",RP,sep="_"),
                                gene=paste(RP_proximity_mouse_df$rid[midrange],
                                           RP_proximity_mouse_df$resno[midrange], sep = "_")))

    rand_focus_rrnas <- unique(na.omit(sapply(rrnas, FUN=function(x){if (sum(proxpos%in%poslists[[x]])>0) return(x); return(NA)})))
    for (i in 1:max(sapply(rand_focus_rrnas,FUN = function(x){return(length(rrna_rands[[x]]))}))){
      randset <- unlist(sapply(rand_focus_rrnas, FUN = function(x){
        ni = i
        randl = length(rrna_rands[[x]])
        if (i > randl)
          ni = i%%randl

        rs_pp <- proxpos[proxpos%in%poslists[[x]]] # rRNAspecific proxposs
        return( min(rs_pp) + ( ( (rs_pp-min(rs_pp)) + rrna_rands[[x]][ni]) %% rrna_lengths[[x]] ) )
        }))
      mouse_gsea_sets_RP <- rbind(mouse_gsea_sets_RP, data.frame(ont=paste(paste("Rand",as.character(i),sep = ""),RP,sep="_"),
                                                                 gene=paste(RP_proximity_mouse_df$rid[randset],
                                                                            RP_proximity_mouse_df$resno[randset], sep = "_")))
    }
  }
}



usethis::use_data(RP_proximity_mouse_df, RP_proximity_human_df, mouse_gsea_sets_RP, human_gsea_sets_RP, internal = TRUE, overwrite = TRUE)

## TEST ##
targetDir=paste(getwd(),"test",sep="/")
compare="group"
organism="mm"
QCplot=TRUE
samples<-read_RP4_samples_file("/home/projects/ribosomal_heterogeneity/RP4_space/Barna_Rpl10aRps25Rpl22_samples.tsv")
rRNA_counts <- RP4_read_rRNA_fragments(samples, organism=organism, QCplot=QCplot, targetDir=targetDir)
org_RP_df <- RP_proximity_mouse_df
gsea_sets_RP <- mouse_gsea_sets_RP
