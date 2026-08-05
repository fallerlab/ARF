---
output:
  pdf_document:
    pandoc_args: --listings
    includes:
      in_header: header.tex
---
# Analysis of Ribosomal rRNA Fragments (ARF)

This is a notebook that provides extensive documentation of the **ARF** pipeline. This documentation is divided in to:
-	[Installation instructions](#Installation instructions)
-	[Pipeline examples](#Pipeline examples)


# Installation instructions

In R (>= 4.0) environment

```
install.packages("devtools")
devtools::install_github("fallerlab/ARF@main")
```

### Required R libraries

Please make sure that you have the following packages installed as dripARF requires them:

-   bedr
-   DESeq2 (>= 1.30.1)
-   clusterProfiler
-   ComplexHeatmap
-   enrichplot
-   fgsea
-   grid
-   ggplot2
-   ggrepel
-   matrixStats
-   reshape2
-   scales
-   SummarizedExperiment
-   tidyverse
-   bio3d
-   Biostrings
-   msa

### Package installation in R
```R
install.packages('renv',
	dependencies = TRUE)

## initiate renv to manage R environment
renv.init()

## install packages in R environment
install.packages(
	'renv','remotes', 'targets', 'bedr', 'curl', 'ggrepel',
	'reshape2', 'tidyverse', 'bio3d', dependencies = TRUE)
	
remotes::install_bioc(
	c('DESeq2', 'matrixStats',
        'clusterProfiler', 'enrichplot', 'fgsea',
        'ComplexHeatmap',
        'msa', 'SummarizedExperiment'),
    dependencies = TRUE)

renv::install("fallerlab/ARF@main")
```



# Pipeline examples

### Running other organisms

#### Lifting over distances from PDB rRNA to study organism's  rRNA

To use ARF for any organism apart from those in the ARF structure database, the user must 

1. use structures from the ARF database of ribosome structures that is exactly for the organism of interest or closely related

2. have the structure of the organism of interest. 



**1. For an organism that is already in the structural database of ARF or an organism that is evolutionary close to the organism of interest, the PDB ID can be provided to ARF for downstream geneset generation.**

**Ribosome profiling reveals the fine-tuned response of *Escherichia coli* to mild and severe acid stress**

*The response to acidity is crucial for neutralophilic bacteria.  Escherichia coli has a well characterized regulatory network to induce  multiple defense mechanisms against excess of protons. Nevertheless,  systemic studies of the transcriptional and translational reprogramming  of E. coli to different acidic strengths have not yet been performed.  Here, we used ribosome profiling and mRNA sequencing to determine the  response of E.[ more...](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219022)*

*Organism: 	Escherichia coli str. K-12 substr. MG1655*



**NB:** Get the fasta file of the structure from PDB and use it for the alignment of the reads. 



Where the user wants to use ribosomes in the ARF ribosome database, the PDB ID **(6XZA)** is specified in the **ARF_parse_PDB_ribosome** function to generate the distances.

```R
############ Computing the distances between RP and rRNA: ARF_parse_PDB_ribosome
conversion_table_generator <- function(PDB_ID = NULL,
                                       ALL_PDB_RPS_edited_file = "./ARF/data-raw/ALL_PDB_RPS_edited.csv") {
  
  if(!is.null(PDB_ID)) {
    final_conversion_df <- all_pdb_rps_edited_file[all_pdb_rps_edited_file$PDB_id == PDB_ID, ] |>
      mutate(RPname = ifelse(grepl("rRNA",RPname), stringr::str_extract(RPname, "[0-9]+S"), RPname)) |>
      select(c("ID", RPname, "RP_new", "chain_info")) |>
      mutate(RP_name = "RPname", chainID = sub(".* (.*)\\[.*", "\\1", chain_info)) |>
      select(-c(chain_info, RPname))
    
  }   
  return(final_conversion_df)
}

final_conversion_df <- conversion_table_generator(PDB_ID = "6XZA") |> distinct()

> head(final_conversion_df)
       ID RP_new RP_name chainID
1 6XZA_29   uL10     L10      CA
2 6XZA_30   uL11     L11      DA
3 6XZA_31   uL13     L13      EA
4 6XZA_32   uL14     L14      FA
5 6XZA_33   uL15     L15      GA
6 6XZA_34   uL16     L16      HA
```



After generating the conversion table, the distances between the ribosomal proteins and the rRNA positions are computed to generate the **RP_proximity_df** table.

```R
############ ARF_parse_PDB_ribosome
RP_proximity_df <- ARF::ARF_parse_PDB_ribosome(species = "ec", PDBid = "6XZA",
                                               download_directory = "./Escherichia/ARF_results")

> head(RP_proximity_df)
               rRNA resno     bL17     bL19     bL20     bL21     bL25     bL27     bL28     bL32     bL33     bL34     bL35     bL36     bS16     bS18     bS20     bS21
rRNA_16S_1 rRNA_16S     1 124.6766 84.05763 167.6931 166.7970 157.9360 141.4937 154.0358 143.1039 155.1185 132.7517 166.9200 128.1418 19.13299 93.14598 60.17656 83.31601
rRNA_16S_2 rRNA_16S     2 121.3063 81.04034 164.6928 164.0205 156.0743 139.3862 151.4715 140.0738 153.1025 129.5194 164.7371 126.2618 19.42014 91.61643 57.53053 82.04972
rRNA_16S_3 rRNA_16S     3 118.5931 79.10021 162.2417 161.6628 154.9880 137.4868 148.5528 137.5285 150.9590 126.3271 162.4165 125.5205 20.86008 88.66694 54.95508 79.60290
rRNA_16S_4 rRNA_16S     4 115.8100 76.37127 154.7506 152.8312 146.7143 126.4935 136.4281 130.0562 138.6470 117.9102 150.2280 119.7526 30.05383 75.60108 57.87327 65.76929
rRNA_16S_5 rRNA_16S     5 123.7808 83.46115 161.1521 158.3869 149.4530 129.7835 141.2482 136.5964 141.5937 125.2058 153.7018 122.8201 33.53445 77.87900 65.70006 66.83129
rRNA_16S_6 rRNA_16S     6 123.8482 83.87247 159.3839 156.0171 146.2595 126.0625 137.9450 135.0728 137.4604 123.8948 149.7688 120.5839 36.55140 74.67179 67.84841 62.66192


```



#### Generating final rRNA position sets for RPs and collision sets

In obtaining the final genesets, the structural fasta and the rRNA sequence pairs (the same for both source and target) are provided to **ARF_convert_Ribo3D_pos**, **dripARF_get_RP_proximity_sets** and **dricARF_liftover_collision_sets** functions

**NB:** ARF automatically replaces non-IUPAC characters (e.g., X) in the rRNA FASTA with N prior to sequence alignment. No manual intervention is required.

```R
########### ARF_convert_Ribo3D_pos
LO.RP_proximity_df <- ARF::ARF_convert_Ribo3D_pos(
  source_distance_file = "./Ribosome.3D.6XZA.ARF.minimum_distances.csv",
  source_rRNAs_fasta = "./Escherichia/ARF_results/6XZA.rRNAs.fasta",
  target_species = "ec",
  target_rRNAs_fasta = "./Escherichia/ARF_results/6XZA.rRNAs.fasta",
  rRNA_pairs = list(c("rRNA_16S", "rRNA_16S"), c("rRNA_23S", "rRNA_23S"), c("rRNA_5S", "rRNA_5S")),
  source_positions = NULL,
  source_sets = NULL,
  type = "distances"
)

> head(LO.RP_proximity_df)
               rRNA resno     bL17     bL19     bL20     bL21     bL25     bL27     bL28     bL32     bL33     bL34     bL35     bL36     bS16     bS18     bS20     bS21
rRNA_16S_1 rRNA_16S     1 124.6766 84.05763 167.6931 166.7970 157.9360 141.4937 154.0358 143.1039 155.1185 132.7517 166.9200 128.1418 19.13299 93.14598 60.17656 83.31601
rRNA_16S_2 rRNA_16S     2 121.3063 81.04034 164.6928 164.0205 156.0743 139.3862 151.4715 140.0738 153.1025 129.5194 164.7371 126.2618 19.42014 91.61643 57.53053 82.04972
rRNA_16S_3 rRNA_16S     3 118.5931 79.10021 162.2417 161.6628 154.9880 137.4868 148.5528 137.5285 150.9590 126.3271 162.4165 125.5205 20.86008 88.66694 54.95508 79.60290
rRNA_16S_4 rRNA_16S     4 115.8100 76.37127 154.7506 152.8312 146.7143 126.4935 136.4281 130.0562 138.6470 117.9102 150.2280 119.7526 30.05383 75.60108 57.87327 65.76929
rRNA_16S_5 rRNA_16S     5 123.7808 83.46115 161.1521 158.3869 149.4530 129.7835 141.2482 136.5964 141.5937 125.2058 153.7018 122.8201 33.53445 77.87900 65.70006 66.83129
rRNA_16S_6 rRNA_16S     6 123.8482 83.87247 159.3839 156.0171 146.2595 126.0625 137.9450 135.0728 137.4604 123.8948 149.7688 120.5839 36.55140 74.67179 67.84841 62.66192

########### dripARF_get_RP_proximity_sets
gsea_sets_RP <- ARF::dripARF_get_RP_proximity_sets(
  RP_proximity_df = LO.RP_proximity_df,
  additional_RPcols = c(),
  rRNAs_fasta = "./Escherichia/ARF_results/6XZA.rRNAs.fasta",
  thresholds = NULL,
  cap_added_RPcols = F
)

> head(gsea_sets_RP)
   ont          gene
1 bL17  rRNA_23S_489
2 bL17  rRNA_23S_491
3 bL17  rRNA_23S_492
4 bL17 rRNA_23S_1266
5 bL17 rRNA_23S_1267
6 bL17 rRNA_23S_1268

###########  dricARF_liftover_collision_sets
gsea_sets_Collision <- ARF::dricARF_liftover_collision_sets(
  target_species = "ec",
  target_rRNAs_fasta = "./Escherichia/ARF_results/6XZA.rRNAs.fasta",
  rRNA_pairs = list(c("18S", "16S"), c("28S", "23S"), c("5S", "5S"))
)

> tail(gsea_sets_Collision)
                   ont          gene
184495 Rand99_Rib.Col. rRNA_23S_1083
184496 Rand99_Rib.Col. rRNA_23S_1084
184497 Rand99_Rib.Col. rRNA_23S_1085
184498 Rand99_Rib.Col. rRNA_23S_1196
184499 Rand99_Rib.Col. rRNA_5S_56
184500 Rand99_Rib.Col. rRNA_5S_57
```



#### Predicting rRNA changes

In determining rRNA positions changes and enrichment tests to predict likely changes in the populations of ribosomes, the samples file which contains the name of the sample, the path to the bedgraph files and the groups assigned to each sample that is used for differential analysis	

| sampleName   | bedGraphFile                                                 | group      |
| ---------------- | ------------------------------------------------------------ | ----- |
| RIBO_pH7.6_1 | PRJNA906596/riboseq/SRR22447291/tophat_align/accepted_hits.bedGraph | RIBO_pH7.6 |
| RIBO_pH5.8_1 | PRJNA906596/riboseq/SRR22447288/tophat_align/accepted_hits.bedGraph | RIBO_pH5.8 |
| RIBO_pH4.4_1 | PRJNA906596/riboseq/SRR22447286/tophat_align/accepted_hits.bedGraph | RIBO_pH4.4 |
| RIBO_pH7.6   | PRJNA906596/riboseq/SRR22447285/tophat_align/accepted_hits.bedGraph | RIBO_pH7.6 |
| RIBO_pH5.8_2 | PRJNA906596/riboseq/SRR22447284/tophat_align/accepted_hits.bedGraph | RIBO_pH5.8 |
| RIBO_pH4.4   | PRJNA906596/riboseq/SRR22447283/tophat_align/accepted_hits.bedGraph | RIBO_pH4.4 |
| RIBO_pH7.6_3 | PRJNA906596/riboseq/SRR22447282/tophat_align/accepted_hits.bedGraph | RIBO_pH7.6 |
| RIBO_pH5.8_3 | PRJNA906596/riboseq/SRR22447281/tophat_align/accepted_hits.bedGraph | RIBO_pH5.8 |
| RIBO_pH4.4_  | PRJNA906596/riboseq/SRR22447280/tophat_align/accepted_hits.bedGraph | RIBO_pH4.4 |



#### Run dripARF

```R
########## Run dricARF
dricARF_results <- ARF::dricARF(
  samplesFile = "./Escherichia/data/PRJNA906596/riboseq/samples.tsv",
  rRNAs_fasta = "./Escherichia/ARF_results/6XZA.rRNAs.fasta",
  samples_df = NULL,
  organism = NULL,
  compare = "group",
  QCplot = TRUE,
  targetDir = "./Escherichia/ARF_results/dricARF",
  comparisons = NULL,
  exclude = NULL,
  GSEAplots = TRUE,
  gsea_sets_RP = gsea_sets_RP,
  RP_proximity_df = LO.RP_proximity_df,
  gsea_sets_Collision = gsea_sets_Collision
)

> head(dricARF_results)
                      comp Description ORA.overlap ORA.setSize    ORA.padj        ORA.p RPSEA.NES RPSEA.NES_randZ   RPSEA.padj   RPSEA.pval      RPSEA.q C1.avg.read.c
1 RIBO_pH7.6_vs_RIBO_pH5.8        bS20          21         227 1.40943e-09 5.752778e-11  1.279929       1.1873237 5.980805e-07 1.220572e-08 1.039753e-07      14172.90
2 RIBO_pH7.6_vs_RIBO_pH5.8        uL24           0         158 1.00000e+00 1.000000e+00  1.279397       1.0531535 1.411524e-05 8.641985e-07 5.012537e-06      16624.91
3 RIBO_pH7.6_vs_RIBO_pH5.8        bS16          24         227 4.13438e-12 8.437509e-14  1.262678       1.2924111 2.427938e-06 9.909952e-08 7.009421e-07      12297.38
4 RIBO_pH7.6_vs_RIBO_pH5.8        bL33           2         119 1.00000e+00 6.020576e-01  1.219017       0.9005245 2.832980e-03 4.625274e-04 1.409363e-03      17501.75
5 RIBO_pH7.6_vs_RIBO_pH5.8        uL15           1         227 1.00000e+00 9.809963e-01  1.196271       1.0133013 1.084775e-03 8.855306e-05 3.249352e-04      17865.88
6 RIBO_pH7.6_vs_RIBO_pH5.8        bL35           2         227 1.00000e+00 9.030230e-01  1.187662       0.9754634 1.410927e-03 1.439722e-04 4.935334e-04      18550.52
  C2.avg.read.c
1      10466.22
2      20528.69
3      11700.20
4      14066.40
5      14744.64
6      16900.58
```





**2. For an organism whose structure is not in the database of ARF and you want to provide a PDB structure.**

**Arabidopsis HOT3/eIF5B1 constrains rRNA RNAi by facilitating 18S rRNA maturation during translation initiation**
*Ribosome biogenesis is essential for protein synthesis in gene expression. Yeast eIF5B has been shown biochemically to facilitate 18S rRNA 3' end maturation during late-40S ribosomal subunit assembly and gate the transition from translation initiation to elongation. But the effects of eIF5B have not been studied at the genome-wide level in any organism, and 18S rRNA 3' end maturation is poorly understood in plants. Arabidopsis HOT3/eIF5B1 was found to promote development and heat-stress acclimation by translational regulation, but its molecular function remained unknown. Here, we show that HOT3 is a late-stage ribosome biogenesis factor that facilitates 18S rRNA 3' end processing and is a translation initiation factor that globally impacts the transition from initiation to elongation. By developing and implementing 18S-ENDseq, we revealed previously unknown events in 18S rRNA 3' end maturation or metabolism. We quantitatively defined new processing hotspots and identified adenylation as the prevalent non-templated RNA modification at the 3' ends of pre-18S rRNAs. Aberrant 18S rRNA maturation in hot3 further activated RNAi to generate RDR1- and DCL2/4-dependent risiRNAs mainly from a 3' portion of 18S rRNA. We further showed that risiRNAs in hot3 were predominantly localized in ribosome-free fractions not responsible for the 18S rRNA maturation or translation initiation defects in hot3. Our study uncovered the molecular function of HOT3/eIF5B1 in 18S rRNA maturation at the late-40S assembly stage and revealed the regulatory crosstalk among ribosome biogenesis, mRNA translation initiation, and siRNA biogenesis in plants. Overall design: Comparative translation profiling analysis of Ribo-seq data for inflorescence of WT, hot3-2, hot3-3 and HOT-EYFP/hot3-2
Organism:	Arabidopsis thaliana*

###### 

**NB:** This is **not** automated in ARF so it has to be done manually due to inconsistencies in the PDB files. Chain names and ID tend to differ from structure to structure making it difficult to automate.

Using ARF with structures that are not inherent to it requires that the different chains in the ribosome structure (from PDB) are properly associated with their standard names. This can be done parsing the ribosome structure and the old names mapped to the new ones using a conversion of table which maps old names to the standard ones. The conversion table (**PDB_chains_2_RP_nomenclature**) should have columns **ID**, **RP_name**. **RP_new**, and **chainID**.

```R
############ Generating a conversion table for RPs with structural chain IDs from ARF structural database
conversion_table_generator <- function(rRNAs_file, organism, ALL_PDB_RPS_edited_file = "./ARF/data-raw/ALL_PDB_RPS_edited.csv") {
  
  ## Read rRNA fasta of PDB structure
  rRNA_seq_set <- Biostrings::readBStringSet(rRNAs_file)
  
  ## Clean PDB structure rRNA fasta names
  conversion_table <- as.data.frame(rRNA_seq_set@ranges) |>
    mutate(
      ID = sub("\\|.*", "", names),
      RP_name = sub(".*\\]\\|(.*)\\|.*", "\\1", names),
      chainID = paste0(sub(".*\\|(.*)\\]\\|.*", "\\1", names), "]")
    ) |>
    mutate(
      RP_name = ifelse(
        grepl("protein|factor|initiator", RP_name),
        sapply(strsplit(RP_name, " "), function(x) x[length(x)]),
        sapply(strsplit(RP_name, " "), function(x) x[1])
      ),
      chainID = gsub("^Chain.*auth |\\]", "", chainID)
    )
  
  ## Convert PDB structure rRNA fasta names to standard RP names
  if(!is.null(organism)) {
    all_pdb_rps_edited_file <- read.table(ALL_PDB_RPS_edited_file, sep = ",", header = TRUE)
    
    organism_selected_rp_conversion_df <- all_pdb_rps_edited_file |>
      filter(organism == organism) |>
      mutate(RPname = ifelse(grepl("rRNA",RPname), stringr::str_extract(RPname, "[0-9]+S"), RPname))
  } else {
    organism_selected_rp_conversion_df <- organism
  }
  
  ## Getting final dataframe and selecting columns needed by ARF
  final_conversion_df <- dplyr::left_join(organism_selected_rp_conversion_df[, c("RPname", "RP_new")], conversion_table,
                   by = join_by(RPname==RP_name)) |>
    mutate(RP_name = RPname) |>
    filter(!is.na(chainID)) |>
    select(ID, RP_name, RP_new, chainID)
  
  return(final_conversion_df)
}

## Run function to get standard RP names from PDB structure names
final_conversion_df <- conversion_table_generator(rRNAs_file = "./Arabidopsis/pdb/rcsb_pdb_8B2L.fasta",
                                                  organism = "Nicotiana tabacum") |>
                        distinct()
                            
> head(final_conversion_df)
       ID RP_name RP_new chainID
1 8B2L_38     L13   uL13      D3
2 8B2L_40     L15   uL15      F3
3 8B2L_42     L19   bL19      H3
4 8B2L_66     L23   uL23      e3
5 8B2L_49     L29   uL29      O3
6 8B2L_54     L34   bL34      T3
```

**NB:** If there is a mapping table between RP_names (**ID, RP_name, RP_new**) and chainIDs, downstream functions work fine especially for the **ARF::ARF_parse_PDB_ribosome** function.



After parsing the ribosome structure, the distances between the ribosomal proteins and the rRNA positions are computed to generate the **RP_proximity_df**.

```R          
############ Computing the distances between RP and rRNA: ARF_parse_PDB_ribosome
RP_proximity_df <- ARF::ARF_parse_PDB_ribosome(species = "AT", PDBid = "8B2L",
                                               download_directory = "./Arabidopsis/ARF_results/",
                                               PDB_chains_2_RP_nomenclature = filter(final_conversion_df, !grepl("RNA", RP_name))
                                            )

> head(RP_proximity_df)
               rRNA resno     bL19     bL34     bS16     bS18     bS21      bS6     eL13     eL15     eL19     eL29     eL34
rRNA_18S_1 rRNA_18S     1 77.98355 111.3755 80.66826 95.53819 27.88385 65.60532 153.5121 127.2228 77.98355 130.3601 111.3755
rRNA_18S_2 rRNA_18S     2 71.68672 104.5543 77.33179 94.14142 22.90794 72.84076 150.9881 123.0247 71.68672 130.3473 104.5543
rRNA_18S_3 rRNA_18S     3 82.59883 109.9158 68.69656 84.00509 28.98561 73.09128 147.7250 121.6217 82.59883 124.8351 109.9158 
rRNA_18S_4 rRNA_18S     4 80.18202 105.8229 62.25940 77.52245 32.41385 73.44756 140.6974 114.8993 80.18202 117.8812 105.8229 
rRNA_18S_5 rRNA_18S     5 82.71370 107.0446 57.66111 72.51181 33.76687 75.84136 138.3124 113.4135 82.71370 115.3651 107.0446 
rRNA_18S_6 rRNA_18S     6 79.25142 102.6017 53.45841 68.06296 38.13404 75.18194 131.2548 106.7126 79.25142 108.4252 102.6017 
```

The computed distances for the rRNA positions are lifted over the organism of interest using the **ARF_convert_Ribo3D_pos** function to produce the **LO.RP_proximity_df **.

```R
########### ARF_convert_Ribo3D_pos
LO.RP_proximity_df <- ARF::ARF_convert_Ribo3D_pos(
  source_distance_file = "./Ribosome.3D.8B2L.ARF.minimum_distances.csv",
  source_rRNAs_fasta = "./Arabidopsis/ARF_results/8B2L.rRNAs.fasta",
  target_species = "AT",
  target_rRNAs_fasta = "./Arabidopsis/organism/rRNA/Arabidopsis_thaliana.TAIR10.seq_20101214.mixed.ENA.23S.final.fa",
  rRNA_pairs = list(c("rRNA_18S", "rRNA_18S"), c("rRNA_25S", "rRNA_25S"), c("rRNA_5S", "rRNA_5S")),
  source_positions = NULL,
  source_sets = NULL,
  type = "distances"
)


> head(LO.RP_proximity_df)
               rRNA resno     bL19     bL34     bS16     bS18     bS21      bS6     eL13     eL15     eL19     eL29     eL34     
rRNA_18S_1 rRNA_18S     1 77.98355 111.3755 80.66826 95.53819 27.88385 65.60532 153.5121 127.2228 77.98355 130.3601 111.3755 
rRNA_18S_2 rRNA_18S     2 71.68672 104.5543 77.33179 94.14142 22.90794 72.84076 150.9881 123.0247 71.68672 130.3473 104.5543 
rRNA_18S_3 rRNA_18S     3 82.59883 109.9158 68.69656 84.00509 28.98561 73.09128 147.7250 121.6217 82.59883 124.8351 109.9158 
rRNA_18S_4 rRNA_18S     4 80.18202 105.8229 62.25940 77.52245 32.41385 73.44756 140.6974 114.8993 80.18202 117.8812 105.8229 
rRNA_18S_5 rRNA_18S     5 82.71370 107.0446 57.66111 72.51181 33.76687 75.84136 138.3124 113.4135 82.71370 115.3651 107.0446 
rRNA_18S_6 rRNA_18S     6 79.25142 102.6017 53.45841 68.06296 38.13404 75.18194 131.2548 106.7126 79.25142 108.4252 102.6017
```



#### Generating final rRNA position sets for RPs

rRNA position sets for RPs used in GSEA are finally generated with the **dripARF_get_RP_proximity_sets** function to generate the **gsea_sets_RP** dataframe.

**NB:** The sequences in the rRNA fasta file must have headers that correspond to the that in the gsea_sets_RP gene column. Therefore, in this case, the rRNA fasta should be ***>rRNA_23S*** , ***>rRNA_16S*** or ***>rRNA_5S***.

```R
########### dripARF_get_RP_proximity_sets
gsea_sets_RP <- ARF::dripARF_get_RP_proximity_sets(
  RP_proximity_df = LO.RP_proximity_df,
  additional_RPcols = c(),
  rRNAs_fasta = "./Arabidopsis/organism/rRNA/Arabidopsis_thaliana.TAIR10.seq_20101214.mixed.ENA.23S.final.fa",
  thresholds = NULL,
  cap_added_RPcols = F
)

> head(gsea_sets_RP)
   ont         gene
1 bL19 rRNA_18S_817
2 bL19 rRNA_18S_818
3 bL19 rRNA_18S_819
4 bL19 rRNA_18S_820
5 bL19 rRNA_18S_821
6 bL19 rRNA_18S_822
```



#### Generating final rRNA position sets for collisions

rRNA position sets for RPs used in GSEA are finally generated with the **dripARF_get_RP_proximity_sets** function to generate the **gsea_sets_RP** dataframe.

**NB:** Tweak rRNA sequence IDs in the fasta file (if they are different from the what is in the *rRNA pairs*) before mapping to so that IDs correspond with the rRNA pairs used here. **"18S", "28S", and "5S"** should not change in the pair list.

```R
###########  dricARF_liftover_collision_sets
gsea_sets_Collision <- ARF::dricARF_liftover_collision_sets(
  target_species = "AT",
  target_rRNAs_fasta = "./Arabidopsis/organism/rRNA/Arabidopsis_thaliana.TAIR10.seq_20101214.mixed.ENA.23S.final.fa",
  rRNA_pairs = list(c("18S", "rRNA_18S"), c("25S", "rRNA_25S"), c("5S", "rRNA_5S"))
)


> tail(gsea_sets_Collision)
                   ont               gene
221695 Rand99_Rib.Col. rRNA_rRNA_25S_1567
221696 Rand99_Rib.Col. rRNA_rRNA_25S_1568
221697 Rand99_Rib.Col. rRNA_rRNA_25S_1569
221698 Rand99_Rib.Col. rRNA_rRNA_25S_1688
221699 Rand99_Rib.Col.    rRNA_rRNA_5S_44
221700 Rand99_Rib.Col.    rRNA_rRNA_5S_45
```



#### Predicting rRNA changes

To predict changes in rRNA positions and position sets, the samples file which contains the name of the sample, the path to the bedgraph files and the groups used for differential analysis	

| sampleName        | bedGraphFile                                                 | group        |
| ----------------- | ------------------------------------------------------------ | ------------ |
| WT Ribo seq1      | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110909/tophat_align/accepted_hits.bedGraph | WT Ribo      |
| WT Ribo seq2      | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110908/tophat_align/accepted_hits.bedGraph | WT Ribo      |
| h32 Ribo seq1     | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110907/tophat_align/accepted_hits.bedGraph | h32 Ribo     |
| h32 Ribo seq2     | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110906/tophat_align/accepted_hits.bedGraph | h32 Ribo     |
| d2d4 Ribo seq1    | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110905/tophat_align/accepted_hits.bedGraph | d2d4 Ribo    |
| d2d4 Ribo seqq2   | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110904/tophat_align/accepted_hits.bedGraph | d2d4 Ribo    |
| h32d2d4 Ribo seq1 | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110903/tophat_align/accepted_hits.bedGraph | h32d2d4 Ribo |
| h32d2d4 Ribo seq2 | /home/edwin/test/Arabidopsis/data/PRJNA925168/riboseq/SRR23110902/tophat_align/accepted_hits.bedGraph | h32d2d4 Ribo |

#### Run dripARF

```R
########## Run dripARF
dripARF_results <- ARF::dripARF(
  samplesFile = "./Arabidopsis/data/PRJNA925168/riboseq/samples.tsv",
  rRNAs_fasta = "./Arabidopsis/organism/rRNA/Arabidopsis_thaliana.TAIR10.seq_20101214.mixed.ENA.23S.final.fa",
  samples_df = NULL,
  organism = "AT",
  compare = "group",
  QCplot = TRUE,
  targetDir = "./Arabidopsis/ARF_results/dripARF",
  comparisons = NULL,
  exclude = NULL,
  GSEAplots = TRUE,
  gsea_sets_RP = gsea_sets_RP,
  RP_proximity_df = LO.RP_proximity_df
)

> head(dripARF_results)
                 comp Description ORA.overlap ORA.setSize     ORA.padj        ORA.p RPSEA.NES RPSEA.NES_randZ   RPSEA.padj   RPSEA.pval      RPSEA.q C1.avg.read.c C2.avg.read.c
1 WT Ribo_vs_h32 Ribo        eL38          38          56 1.843017e-04 2.835410e-05  1.227182       0.9456723 5.703480e-03 1.316188e-03 3.941369e-03     13204.424      32137.14
2 WT Ribo_vs_h32 Ribo         eS8          46         238 1.000000e+00 1.000000e+00  1.214294       1.1423127 8.393819e-07 4.304523e-08 4.200291e-07    150467.415     109519.60
3 WT Ribo_vs_h32 Ribo         uS8          46         238 1.000000e+00 1.000000e+00  1.214294       1.1423127 8.393819e-07 4.304523e-08 4.200291e-07    150467.415     109519.60
4 WT Ribo_vs_h32 Ribo        bL19         136         177 1.799458e-22 9.227991e-24  1.196394       1.2830180 7.112028e-05 7.294388e-06 3.924466e-05      6817.121      13055.40
5 WT Ribo_vs_h32 Ribo        eL19         136         177 1.799458e-22 9.227991e-24  1.196394       1.2830180 7.112028e-05 7.294388e-06 3.924466e-05      6817.121      13055.40
6 WT Ribo_vs_h32 Ribo        bL34         120         178 1.085809e-12 1.113650e-13  1.173411       0.8794976 6.355478e-04 1.303688e-04 5.163078e-04     94137.795     101895.72
```

#### Run dricARF

```R
########## Run dricARF
dricARF_results <- ARF::dricARF(
  samplesFile = "./Arabidopsis/data/PRJNA925168/riboseq/samples.tsv",
  rRNAs_fasta = "./Arabidopsis/organism/rRNA/Arabidopsis_thaliana.TAIR10.seq_20101214.mixed.ENA.23S.final.fa",
  samples_df = NULL,
  organism = "AT",
  compare = "group",
  QCplot = TRUE,
  targetDir = "./Arabidopsis/ARF_results/dripARF",
  comparisons = NULL,
  exclude = NULL,
  GSEAplots = TRUE,
  gsea_sets_RP = gsea_sets_RP,
  RP_proximity_df = LO.RP_proximity_df,
  gsea_sets_Collision = gsea_sets_Collision
)

> head(dricARF_results)
                 comp Description ORA.overlap ORA.setSize     ORA.padj        ORA.p RPSEA.NES RPSEA.NES_randZ   RPSEA.padj   RPSEA.pval      RPSEA.q C1.avg.read.c C2.avg.read.c
1 WT Ribo_vs_h32 Ribo        eL38          38          56 1.843017e-04 2.835410e-05  1.220301       0.9421230 5.703480e-03 1.316188e-03 3.951069e-03     13204.424      32137.14
2 WT Ribo_vs_h32 Ribo         eS8          46         238 1.000000e+00 1.000000e+00  1.211576       1.1498347 8.393819e-07 4.304523e-08 4.210628e-07    150467.415     109519.60
3 WT Ribo_vs_h32 Ribo         uS8          46         238 1.000000e+00 1.000000e+00  1.211576       1.1498347 8.393819e-07 4.304523e-08 4.210628e-07    150467.415     109519.60
4 WT Ribo_vs_h32 Ribo        bL19         136         177 1.799458e-22 9.227991e-24  1.192329       1.2804162 7.112028e-05 7.294388e-06 3.934124e-05      6817.121      13055.40
5 WT Ribo_vs_h32 Ribo        eL19         136         177 1.799458e-22 9.227991e-24  1.192329       1.2804162 7.112028e-05 7.294388e-06 3.934124e-05      6817.121      13055.40
6 WT Ribo_vs_h32 Ribo        bL34         120         178 1.085809e-12 1.113650e-13  1.169130       0.8768647 6.355478e-04 1.303688e-04 5.175785e-04     94137.795     101895.72
```



```R
> sessionInfo()
R version 4.4.2 (2024-10-31)
Platform: x86_64-conda-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS/LAPACK: /home/<user name>/mambaforge/envs/<env name>/lib/libopenblasp-r0.3.28.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=nl_NL.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       

time zone: <region>/<city>
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ARF_2.1     tidyr_1.3.1 dplyr_1.1.4

loaded via a namespace (and not attached):
  [1] DBI_1.2.3                   gson_0.1.0                  rlang_1.1.4                 magrittr_2.0.3             
  [5] clue_0.3-66                 GetoptLong_1.0.5            DOSE_4.0.0                  matrixStats_1.4.1          
  [9] compiler_4.4.2              RSQLite_2.3.9               systemfonts_1.1.0           png_0.1-8                  
 [13] vctrs_0.6.5                 reshape2_1.4.4              stringr_1.5.1               pkgconfig_2.0.3            
 [17] shape_1.4.6.1               crayon_1.5.3                fastmap_1.2.0               XVector_0.46.0             
 [21] labeling_0.4.3              utf8_1.2.4                  tzdb_0.4.0                  enrichplot_1.26.3          
 [25] UCSC.utils_1.2.0            ragg_1.3.3                  purrr_1.0.2                 bit_4.5.0.1                
 [29] zlibbioc_1.52.0             cachem_1.1.0                aplot_0.2.3                 GenomeInfoDb_1.42.1        
 [33] jsonlite_1.8.9              blob_1.2.4                  DelayedArray_0.32.0         BiocParallel_1.40.0        
 [37] parallel_4.4.2              cluster_2.1.7               R6_2.5.1                    stringi_1.8.4              
 [41] RColorBrewer_1.1-3          GenomicRanges_1.58.0        GOSemSim_2.32.0             SummarizedExperiment_1.36.0
 [45] Rcpp_1.0.13-1               iterators_1.0.14            ggtangle_0.0.5              R.utils_2.12.3             
 [49] readr_2.1.5                 IRanges_2.40.1              Matrix_1.6-5                splines_4.4.2              
 [53] igraph_2.1.2                tidyselect_1.2.1            abind_1.4-8                 qvalue_2.38.0              
 [57] rstudioapi_0.17.1           doParallel_1.0.17           codetools_0.2-20            lattice_0.22-6             
 [61] tibble_3.2.1                plyr_1.8.9                  bio3d_2.4-5                 withr_3.0.2                
 [65] Biobase_2.66.0              treeio_1.30.0               KEGGREST_1.46.0             gridGraphics_0.5-1         
 [69] circlize_0.4.16             Biostrings_2.74.0           pillar_1.9.0                ggtree_3.14.0              
 [73] MatrixGenerics_1.18.0       renv_1.0.11                 foreach_1.5.2               stats4_4.4.2               
 [77] clusterProfiler_4.14.4      ggfun_0.1.8                 generics_0.1.3              vroom_1.6.5                
 [81] hms_1.1.3                   S4Vectors_0.44.0            ggplot2_3.5.1               munsell_0.5.1              
 [85] scales_1.3.0                tidytree_0.4.6              glue_1.8.0                  lazyeval_0.2.2             
 [89] tools_4.4.2                 data.table_1.15.4           locfit_1.5-9.10             fgsea_1.32.0               
 [93] fs_1.6.5                    fastmatch_1.1-4             cowplot_1.1.3               grid_4.4.2                 
 [97] ape_5.8                     AnnotationDbi_1.68.0        colorspace_2.1-1            nlme_3.1-165               
[101] GenomeInfoDbData_1.2.13     patchwork_1.3.0             msa_1.38.0                  cli_3.6.3                  
[105] textshaping_0.4.0           fansi_1.0.6                 S4Arrays_1.6.0              ComplexHeatmap_2.21.1      
[109] gtable_0.3.6                R.methodsS3_1.8.2           yulab.utils_0.1.8           DESeq2_1.46.0              
[113] digest_0.6.37               BiocGenerics_0.52.0         SparseArray_1.6.0           ggrepel_0.9.6              
[117] ggplotify_0.1.2             rjson_0.2.23                farver_2.1.2                memoise_2.0.1              
[121] R.oo_1.27.0                 lifecycle_1.0.4             httr_1.4.7                  GlobalOptions_0.1.2        
[125] GO.db_3.20.0                bit64_4.5.2                              
```