# Analysis of Ribosomal rRNA Fragments (ARF)

This is a notebook that provides extensive documentation of the **ARF** pipeline. This documentation is divided into:
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

-   DESeq2 (>= 1.30.1)
-   SummarizedExperiment
-   matrixStats
-   ggplot2
-   ggrepel
-   scales
-   reshape2
-   clusterProfiler
-   fgsea
-   ComplexHeatmap
-   grid
-   bio3d
-   Biostrings
-   msa
-   cowplot
-   dplyr
-   magrittr
-   readr
-   RColorBrewer
-   circlize
-   wesanderson
-   GSVA
-   limma
-   tidyr

### Package installation in R
```R
install.packages('renv',
	dependencies = TRUE)

## initiate renv to manage R environment
renv::init()

## install packages in R environment
install.packages(
	c('renv', 'remotes', 'curl', 'ggplot2', 'ggrepel', 'scales',
	  'reshape2', 'bio3d', 'cowplot', 'dplyr', 'magrittr', 'readr',
	  'RColorBrewer', 'circlize', 'wesanderson', 'tidyr'),
	dependencies = TRUE)

remotes::install_bioc(
	c('DESeq2', 'SummarizedExperiment', 'matrixStats',
	  'clusterProfiler', 'fgsea', 'ComplexHeatmap',
	  'msa', 'Biostrings', 'GSVA', 'limma'),
	dependencies = TRUE)

renv::install("fallerlab/ARF@main")
```



# Pipeline examples

### Running other organisms

#### Lifting over distances from PDB rRNA to study organism's  rRNA

To use ARF for any organism apart from those in the ARF structure database, the user must 

1. use structures from the ARF database of ribosome structures

2. either have the structure of or the closest structure to the organism of interest. 




1.**Ribosome profiling reveals the fine-tuned response of *Escherichia coli* to mild and severe acid stress**

*The response to acidity is crucial for neutralophilic bacteria.  Escherichia coli has a well characterized regulatory network to induce  multiple defense mechanisms against excess of protons. Nevertheless,  systemic studies of the transcriptional and translational reprogramming  of E. coli to different acidic strengths have not yet been performed.  Here, we used ribosome profiling and mRNA sequencing to determine the  response of E. [more...](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219022)*

*Organism: 	Escherichia coli str. K-12 substr. MG1655*



For an organism that is already in the structural database of dripARF, the PDB ID can be provided to ARF for downstream geneset generation.

> **Note:** Get the fasta file of the structure from PDB and use it for the alignment of the reads.



Where the user wants to use ribosomes in the ARF ribosome database, the PDB ID **(6XZA)** is specified in the **ARF_parse_PDB_ribosome** function to generate the distances.

```R
############ Computing the distances between RP and rRNA: ARF_parse_PDB_ribosome
conversion_table_generator <- function(PDB_ID = NULL,
                                       ALL_PDB_RPS_edited_file = "./ARF/data-raw/ALL_PDB_RPS_edited.csv") {
  # rRNAs_file <- "Pseudomonas/pdb/rcsb_pdb_7UNW.fasta"
  
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

> **Note:** Check for atypical (N, X, etc.) characters that may occur in the rRNA fasta file. These can generate an error in the ARF pipeline.

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



#### Run dricARF

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





2.**Hfq mediates transcriptome-wide RNA structurome reprogramming under virulence-inducing conditions in a phytopathogen**
*Although RNA structures play important roles in regulating gene expression, the mechanism and function of mRNA folding in plant bacterial pathogens remain elusive. Therefore, we perform dimethyl sulfate sequencing (DMS-seq) on the Pseudomonas syringae under nutrition-rich and deficient conditions, revealing that the mRNA structure changes substantially in the minimal medium (MM) that tunes global translation efficiency (TE), thereby inducing virulence.
Organism:	Pseudomonas savastanoi pv. phaseolicola 1448A*



For an organism whose structure is not in the database of ARF,

**a.** use a structure of an organism in the ARF database that is evolutionary close to the organism of interest.

**b**. provide a PDB structure 



###### a. Using ARF ribosome structure database

Where the user wants to use ribosomes in the ARF ribosome database, the PDB ID is specified in the **ARF_parse_PDB_ribosome** function to generate the distances.

```R
############ Computing the distances between RP and rRNA: ARF_parse_PDB_ribosome
RP_proximity_df <- ARF::ARF_parse_PDB_ribosome(species = "ps", PDBid = "3J9W",
                            download_directory = "./Pseudomonas/ARF_results/")

> head(RP_proximity_df)
                rRNA resno     bL17     bL19     bL20     bL21      bL27     bL28     bL31     bL32     bL33     bL34     bL35  
rRNA_16S_8  rRNA_16S     8 122.6946 83.35646 159.3836 157.4184 116.09391 137.6919 129.2263 134.1063 138.9464 123.1226 149.7557 
rRNA_16S_9  rRNA_16S     9 115.0248 78.22618 149.9908 147.5692 106.17896 125.7729 123.6571 124.8846 127.8942 112.8808 138.7025 
rRNA_16S_10 rRNA_16S    10 122.0681 82.35799 152.2817 148.5858 105.83935 129.4725 111.2955 128.2300 126.3189 120.0920 137.5317 
rRNA_16S_11 rRNA_16S    11 113.1197 75.09033 141.9361 138.1703  95.44860 118.9757 107.1514 118.0651 116.6108 109.8239 127.5232 
rRNA_16S_12 rRNA_16S    12 111.8091 74.61790 139.1231 135.0645  92.13528 115.3906 103.9626 115.5878 112.6200 107.5341 123.5958 
rRNA_16S_13 rRNA_16S    13 110.4894 74.40479 136.4506 132.1372  89.05201 111.6549 101.4642 113.2292 108.7133 105.1049 119.7709 
```

```R
LO.RP_proximity_df <- ARF::ARF_convert_Ribo3D_pos(
  source_distance_file = "./Ribosome.3D.3J9W.ARF.minimum_distances.csv",
  source_rRNAs_fasta = "./Pseudomonas/ARF_results/3J9W.rRNAs.fasta",
  target_species = "ps",
  target_rRNAs_fasta = "./Pseudomonas/organism/rRNA/Pseudomonas.savastanoi.fa",
  rRNA_pairs = list(c("rRNA_16S", "ps_rRNA_16S"), c("rRNA_23S", "ps_rRNA_23S"), c("rRNA_5S", "ps_rRNA_5S")),
  source_positions = NULL,
  source_sets = NULL,
  type = "distances"
)

> head(LO.RP_proximity_df)
                      rRNA resno     bL17     bL19     bL20     bL21      bL27     bL28     bL31     bL32     bL33     bL34     
ps_rRNA_16S_6  ps_rRNA_16S     6 122.6946 83.35646 159.3836 157.4184 116.09391 137.6919 129.2263 134.1063 138.9464 123.1226 
ps_rRNA_16S_7  ps_rRNA_16S     7 115.0248 78.22618 149.9908 147.5692 106.17896 125.7729 123.6571 124.8846 127.8942 112.8808 
ps_rRNA_16S_8  ps_rRNA_16S     8 122.0681 82.35799 152.2817 148.5858 105.83935 129.4725 111.2955 128.2300 126.3189 120.0920 
ps_rRNA_16S_9  ps_rRNA_16S     9 113.1197 75.09033 141.9361 138.1703  95.44860 118.9757 107.1514 118.0651 116.6108 109.8239 
ps_rRNA_16S_10 ps_rRNA_16S    10 111.8091 74.61790 139.1231 135.0645  92.13528 115.3906 103.9626 115.5878 112.6200 107.5341 
ps_rRNA_16S_11 ps_rRNA_16S    11 110.4894 74.40479 136.4506 132.1372  89.05201 111.6549 101.4642 113.2292 108.7133 105.1049 

```



###### 2. Providing a ribosome structure for your organism 

> **Note:** This is **not** automated in ARF and must be done manually due to inconsistencies in PDB files. Chain names and IDs tend to differ from structure to structure, making automation difficult.

Using ARF with structures not in its database requires that the different chains in the ribosome structure (from PDB) are properly associated with their standard names. This can be done by parsing the ribosome structure and mapping old names to the new ones using a conversion table that maps old names to standard ones. The conversion table (**PDB_chains_2_RP_nomenclature**) should have columns **ID**, **RP_name**, **RP_new**, and **chainID**.

```R
############ Generating a conversion table for RPs with structural chain IDs from ARF structural database
conversion_table_generator <- function(rRNAs_file, organism, ALL_PDB_RPS_edited_file = "./ARF/data-raw/ALL_PDB_RPS_edited.csv") {
  # rRNAs_file <- "Pseudomonas/pdb/rcsb_pdb_7UNW.fasta"
  
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
final_conversion_df <- conversion_table_generator(rRNAs_file = "./Pseudomonas/pdb/rcsb_pdb_7UNW.fasta",
                            organism = "Bacillus subtilis subsp. subtilis str. 168 (224308)")
                            
> head(final_conversion_df)
       ID RP_name RP_new chainID
1 7UNW_33     L10   uL10       I
2 7UNW_34     L11   uL11       J
3 7UNW_35     L13   uL13       L
4 7UNW_36     L14   uL14       M
5 7UNW_37     L15   uL15       N
6 7UNW_38     L16   uL16       O
```

> **Note:** If there is a mapping table between RP names (**ID, RP_name, RP_new**) and chainIDs, downstream functions work correctly, especially **ARF::ARF_parse_PDB_ribosome**.



After parsing the ribosome structure, the distances between the ribosomal proteins and the rRNA positions are computed to generate the **RP_proximity_df**.

```R          
############ Computing the distances between RP and rRNA: ARF_parse_PDB_ribosome
RP_proximity_df <- ARF::ARF_parse_PDB_ribosome(species = "ps", PDBid = "7UNW",
                            download_directory = "./Pseudomonas/ARF_results/",
                            PDB_chains_2_RP_nomenclature = filter(final_conversion_df, !grepl("RNA", RP_name))
)
```

The computed distances for the rRNA positions are lifted over the organism of interest using the **ARF_convert_Ribo3D_pos** function to produce the **LO.RP_proximity_df**.
```R
########### ARF_convert_Ribo3D_pos
LO.RP_proximity_df <- ARF::ARF_convert_Ribo3D_pos(
  source_distance_file = "./Ribosome.3D.7UNW.ARF.minimum_distances.csv",
  source_rRNAs_fasta = "./Pseudomonas/ARF_results/7UNW.rRNAs.fasta",
  target_species = "ps",
  target_rRNAs_fasta = "./Pseudomonas/organism/rRNA/Pseudomonas.savastanoi.fa",
  rRNA_pairs = list(c("rRNA_16S", "ps_rRNA_16S"), c("rRNA_23S", "ps_rRNA_23S"), c("rRNA_5S", "ps_rRNA_5S")),
  source_positions = NULL,
  source_sets = NULL,
  type = "distances"
)

```



#### Generating final rRNA position sets for RPs

rRNA position sets for RPs used in GSEA are finally generated with the **dripARF_get_RP_proximity_sets** function to generate the **gsea_sets_RP** dataframe.

> **Note:** The sequences in the rRNA fasta file must have headers that match the entries in the `gsea_sets_RP` gene column. In this case the rRNA fasta headers should be ***>ps_rRNA_23S***, ***>ps_rRNA_16S***, or ***>ps_rRNA_5S***.

```R
########### dripARF_get_RP_proximity_sets
gsea_sets_RP <- ARF::dripARF_get_RP_proximity_sets(
  RP_proximity_df = LO.RP_proximity_df,
  additional_RPcols = c(),
  rRNAs_fasta = "./Pseudomonas/organism/rRNA/Pseudomonas.savastanoi.fa",
  thresholds = NULL,
  cap_added_RPcols = F
)

> head(gsea_sets_RP)
   ont             gene
1 bL17  ps_rRNA_23S_480
2 bL17  ps_rRNA_23S_481
3 bL17  ps_rRNA_23S_482
4 bL17 ps_rRNA_23S_1255
5 bL17 ps_rRNA_23S_1256
6 bL17 ps_rRNA_23S_1257
```



#### Generating final rRNA position sets for collisions

rRNA position sets for RPs used in GSEA are finally generated with the **dripARF_get_RP_proximity_sets** function to generate the **gsea_sets_RP** dataframe.

> **Note:** Tweak rRNA sequence IDs in the fasta file (if they are different from what is in the *rRNA pairs*) before mapping, so that IDs correspond with the rRNA pairs used here. **"18S", "28S", and "5S"** should not change in the pair list.

```R
###########  dricARF_liftover_collision_sets
gsea_sets_Collision <- ARF::dricARF_liftover_collision_sets(
  target_species = "ps",
  target_rRNAs_fasta = "./Pseudomonas/organism/rRNA/Pseudomonas.savastanoi.fa",
  rRNA_pairs = list(c("18S", "ps_rRNA_16S"), c("28S", "ps_rRNA_23S"), c("5S", "ps_rRNA_5S"))
)

> head(gsea_sets_Collision)
               ont            gene
1 sc_6I7O_Col.Int. ps_rRNA_16S_296
2 sc_6I7O_Col.Int. ps_rRNA_16S_297
3 sc_6I7O_Col.Int. ps_rRNA_16S_298
4 sc_6I7O_Col.Int. ps_rRNA_16S_299
5 sc_6I7O_Col.Int. ps_rRNA_16S_300
6 sc_6I7O_Col.Int. ps_rRNA_16S_301
```



#### Predicting rRNA changes
To predict changes in rRNA positions and position sets, the samples file which contains the name of the sample, the path to the bedgraph files and the groups used for differential analysis	

|       sampleName         |bedGraphFile                          |group                         |
|----------------|-------------------------------|-----------------------------|
|Ribo-seq-WT-MM-1|PRJNA892464/riboseq/SRR21981107/tophat_align/accepted_hits.bedGraph|WT-MM            |
|Ribo-seq-WT-MM-2          |PRJNA892464/riboseq/SRR21981106/tophat_align/accepted_hits.bedGraph|WT-MM |
|Ribo-seq-WT-KB-1         |PRJNA892464/riboseq/SRR21981109/tophat_align/accepted_hits.bedGraph|WT-KB|
|Ribo-seq-WT-KB-2         |PRJNA892464/riboseq/SRR21981108/tophat_align/accepted_hits.bedGraph|WT-KB|

#### Run dripARF

```R
########## Run dripARF
dripARF_results <- ARF::dripARF(
  samplesFile = "./Pseudomonas/data/PRJNA892464/riboseq/samples.tsv",
  rRNAs_fasta = "./Pseudomonas/organism/rRNA/Pseudomonas.savastanoi.fa",
  samples_df = NULL,
  organism = NULL,
  compare = "group",
  QCplot = TRUE,
  targetDir = "./Pseudomonas/ARF_results/dripARF",
  comparisons = NULL,
  exclude = NULL,
  gsea_sets_RP = gsea_sets_RP,
  RP_proximity_df = LO.RP_proximity_df
)

> head(dripARF_results)
            comp Description ORA.overlap ORA.setSize  ORA.padj      ORA.p RPSEA.NES RPSEA.NES_randZ   RPSEA.padj   RPSEA.pval      RPSEA.q C1.avg.read.c C2.avg.read.c
1 WT-MM_vs_WT-KB        bS18           1         103 1.0000000 0.88123005  1.360593       1.0825019 6.699161e-05 2.734351e-06 1.241679e-05      46441.02      47565.55
2 WT-MM_vs_WT-KB        uL29           0          96 1.0000000 1.00000000  1.327863       1.1353856 3.808520e-04 2.331747e-05 9.162235e-05     348778.17    1074992.88
3 WT-MM_vs_WT-KB        bL17          10         227 0.7648441 0.01560906  1.300683       1.2126996 2.058084e-06 4.200172e-08 2.393615e-07      43475.23      77354.52
4 WT-MM_vs_WT-KB         uL5           0         135 1.0000000 1.00000000  1.248895       0.8168641 3.327019e-03 3.394917e-04 1.117176e-03      25554.95      29238.72
5 WT-MM_vs_WT-KB        uL23           4         206 1.0000000 0.60626532  1.226569       0.8454657 4.434675e-04 3.620143e-05 1.370837e-04     130670.95     416735.07
6 WT-MM_vs_WT-KB        uS11           1         178 1.0000000 0.97561035  1.179385       0.6180827 2.402107e-02 2.941356e-03 8.231266e-03      28808.16      29985.51
```

#### Run dricARF

```R
########## Run dricARF
dricARF_results <- ARF::dricARF(
  samplesFile = "./Pseudomonas/data/PRJNA892464/riboseq/samples.tsv",
  rRNAs_fasta = "./Pseudomonas/organism/rRNA/Pseudomonas.savastanoi.fa",
  samples_df = NULL,
  organism = "ps",
  compare = "group",
  QCplot = TRUE,
  targetDir = "./Pseudomonas/ARF_results/dricARF",
  comparisons = NULL,
  exclude = NULL,
  gsea_sets_RP = gsea_sets_RP,
  RP_proximity_df = LO.RP_proximity_df,
  gsea_sets_Collision = gsea_sets_Collision
)

> head(dricARF_results)
            comp Description ORA.overlap ORA.setSize  ORA.padj      ORA.p RPSEA.NES RPSEA.NES_randZ   RPSEA.padj   RPSEA.pval      RPSEA.q C1.avg.read.c C2.avg.read.c
1 WT-MM_vs_WT-KB        bS18           1         103 1.0000000 0.88123005  1.357511       1.0823762 3.399440e-05 1.172221e-06 5.845266e-06      46441.02      47565.55
2 WT-MM_vs_WT-KB        uL29           0          96 1.0000000 1.00000000  1.325801       1.1354463 5.588335e-04 2.890518e-05 1.153647e-04     348778.17    1074992.88
3 WT-MM_vs_WT-KB        bL17          10         227 0.9053257 0.01560906  1.298672       1.2126787 3.862281e-06 6.659106e-08 3.982719e-07      43475.23      77354.52
4 WT-MM_vs_WT-KB         uL5           0         135 1.0000000 1.00000000  1.244630       0.8168425 2.064576e-03 1.779807e-04 6.296391e-04      25554.95      29238.72
5 WT-MM_vs_WT-KB        uL23           4         206 1.0000000 0.60626532  1.225220       0.8454988 8.260331e-04 5.696780e-05 2.165598e-04     130670.95     416735.07
6 WT-MM_vs_WT-KB        uS11           1         178 1.0000000 0.97561035  1.177031       0.6181139 3.522836e-02 3.644313e-03 1.000334e-02      28808.16      29985.51
```

