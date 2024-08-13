## ARF (v 2.0)
R package for the <u>**A**</u>nalysis of <u>**R**</u>ibosomal RNA <u>**F**</u>ragments, generally generated by Ribosome Profiling (Ribo-Seq) experiments

Currently, it includes the following methods:

* **dripARF**, ribosome heterogeneity identification pipeline with <u>**d**</u>ifferential <u>**r**</u>ibosomal protein <u>**i**</u>ncorporation <u>**p**</u>redictions
* **dricARF**, <u>**d**</u>ifferential <u>**ri**</u>bosome <u>**c**</u>ollision prediction pipeline 

## Installation
In R (>= 4.0) environment

	install.packages("devtools")
    devtools::install_github("fallerlab/ARF@main")

### Required R libraries
Please make sure that you have the following packages installed as dripARF requires them:

* bedr
* DESeq2 (>= 1.30.1)
* clusterProfiler
* ComplexHeatmap
* enrichplot
* fgsea
* grid
* ggplot2
* ggrepel
* matrixStats
* reshape2
* scales
* SummarizedExperiment
* tidyverse
* bio3d
* Biostrings
* msa

## ARF package - How to use?
To be able to use the ARF package, you first need to align your <u>adapter-trimmed Ribo-seq reads</u> **to <u>given rRNA sequences</u>** which are provided under the **"rRNAs/"** folder for three organisms; human, mouse and yeast. Using these alignments, ARF quantifies the rRNA fragments at each rRNA position for further analysis.

Here is an example command on how to run your Ribo-seq rRNA-alignments.

    tophat -n 2 --no-novel-juncs --no-novel-indels --no-coverage-search --segment-length 25 -o test rRNAs/human_rRNAs example.fastq 

Currently, ARF platform contains the dripARF and dricARF methods. Both tools use these alignments as input. We therefore suggest to perform your positional rRNA abundance quantifications using the "bedtools" software, prior to your dripARF or dricARF run. Note that this is not mandatory since you can also use your indexed .bam files as input, however, it will significantly increase the computational speed of dripARF & dricARF. Following is the example on how you can create your .bedGraph files.

    bedtools genomecov -bg -ibam /your/bam/file.bam >  /your/bedGraph/file.bedGraph

<u>Please note that you have to specify the organism and align your reads strictly to given rRNA sequences when running the pipeline.</u> 

### The *dripARF* method (for ribosome heterogeneity predictions)
Based on position-specific rRNA abundances and the 3D structure of the ribosome, *dripARF* predicts which ribosomal proteins have a change in their ribosomal incorporation between given conditions. With such prediction, it identifies the prime candidates that contributes to the inter-sample heterogeneity of ribosomes across conditions.

The input for the dripARF pipeline is as follows:

* Read alignment files (bam or bedGraph files) that are generated by aligning adapter-trimmed Ribo-seq reads to <u>given</u> rRNA sequences
* <u>tab-seperated samples file</u>; this file contains the sample id as first column, file path to the rRNA read alignment file as second column, and a group/condition identifier as the third column (see test_data/samples.tsv as an example).

If your alignment files, whether in *.bam* or *.bedGraph* format, are ready and accessible through the paths given within *samples.tsv* file, and your target folder is accessible, you can run the dripARF wrapper function in one line within your R environment. 

    dripARF_results <- ARF::dripARF(samplesFile = "samples.tsv", rRNAs_fasta = "rRNAs/mouse_rRNAs.fa", organism="mm", QCplot=TRUE, targetDir="dripARF_results/")

This function will run the whole pipeline and will return a data frame with the ribosome heterogeneity prediction results for all possible comparisons. It will also create comparison-specific ".csv" files to store your results in addition to a few informative QC plots. Most importantly, the pipeline will save a scatterplot and heatmap visualization summarizing your results.

#### dripARF() function in detail

dripARF() wrapper function uses series of other dripARF functions to obtain your results. We summarize below what these functions are used for.

* read_ARF_samples_file()             : This function allows you to read the tab-seperated ARF-samples file, returning a dataframe.
* dripARF_read_rRNA_fragments()       : Read .bedgraph/.bam files to create the rRNA count data.
* dripARF_get_DESEQ_dds()             : The function that normalizes rRNA counts with DESEQ2.
* dripARF_report_RPset_group_counts() : Calculate the average read count of RP contact point sets.
* dripARF_predict_heterogenity()      : Overlap differential rRNA count data with 3d ribosome, rRNA-RP proximity data and predict heterogeneity candidates across groups.
* dripARF_simplify_results()          : Simplify dripARF results based on given thresholds.
* dripARF_result_heatmap()            : Draw heatmap of NES (ES1) and NES_randZscore (ES2) enrichments scores, where RPs are filtered based on given thresholds.
* dripARF_result_scatterplot()        : Draw comparison-specific scatter plots representing RPSEA and ORA analyses results.
* dripARF_report_RPspec_pos_results() : Obtain the positional differential abundance change values for RP proximity sets.
* dripARF_rRNApos_heatmaps()          : Draw rRNA fragment change heatmaps to visualize position-specific differential rRNA fragment abundances.
* dripARF_threshold_test()            : (Experimental function) This function allows you to run the whole dripARF pipeline with varying proximity thresholds.

### The *dricARF* method (for ribosome collision predictions)
Based on position-specific rRNA abundances and the 3D structure of the collided ribosomes, *dricARF* predicts if there is abundance change in ribosome collisions between Ribo-seq performed conditions.

The input for the dricARF pipeline is the same as dripARF function; **read alignment files** and **tab-seperated samples file**.

If your alignment files, whether in .bam or .bedGraph format, are ready and accessible through the paths given within *samples.tsv* file, and your target folder is accessible, you can run the dricARF function as follows. 

    dricARF_results <- ARF::dricARF(samplesFile = "samples.tsv", rRNAs_fasta = "rRNAs/6T7I_yeast_rRNAs.fa", organism="sc", QCplot=TRUE, targetDir="dricARF_results/")

This function will run the integrated dripARF & dricARF pipeline and will return a data frame with the ribosome collision and heterogeneity prediction results for all possible comparisons. It will also create comparison-specific ".csv" files to store your results in addition to a few informative QC plots. Most importantly, the pipeline will save a scatterplot summarizing your results.

#### dricARF() function in detail

dricARF() wrapper function uses series of dripARF functions to obtain your results. However, it also makes use of the following dricARF-specific functions.

* dricARF_result_scatterplot()        : Summarize dripARF & dricARF predictions visually, highlighting collision predictions in particular.


### TESTING ARF

If you have installed the ARF package successfully, you can test it with the data provided under the test_data/ folder.

## ARF package - How to use?
To be able to use the ARF package, you first need to align your <u>adapter-trimmed Ribo-seq reads</u> **to <u>given rRNA sequences</u>** which are provided under the "rRNAs/" folder for two organisms; human and mouse. Using these alignments, ARF quantifies the rRNA fragments at each rRNA position for further analysis.

Currently, ARF package contains only the *dripARF* method. Please stay tuned for future additions :)


### The *dripARF* method/pipeline
Based on position-specific rRNA abundances and the 3D structure of the ribosome, *dripARF* predicts which ribosomal proteins have a change in their ribosomal incorporation between given conditions. With such prediction, it identifies the prime candidates that contributes to the inter-sample heterogeneity of ribosomes across conditions.

The input for the dripARF pipeline is as follows:

* Read alignment files that are generated by aligning adapter-trimmed Ribo-seq reads to <u>given</u> rRNA sequences
* <u>tab-seperated samples file</u>; this file contains the sample id as first column, file path to the rRNA read alignment file as second column, and a group/condition identifier as the third column (see test_data/samples.tsv as an example).

Prior to using dripARF, assuming rRNA read alignments are ready, we also suggest to perform your positional rRNA abundance quantifications using the "bedtools" software. Note that this is not mandatory since you can also use your indexed *.bam* files as input. However, using pre-processed alignments (*.bedGraph*) will significantly increase the computational speed of dripARF. Following is the example on how you can create your *.bedGraph* files.

    bedtools genomecov -bg -ibam /your/bam/file.bam >  /your/bedGraph/file.bedGraph

If your alignment files, whether in *.bam* or *.bedGraph* format, are ready and accessible through the paths given within *samples.tsv* file, and your target folder is accessible, you can run the dripARF wrapper function in one line within your R environment. 

    dripARF_results <- ARF::dripARF(samplesFile = "samples.tsv", organism="mm", QCplot=TRUE, targetDir="dripARF_results/")

This function will run the whole pipeline and will return a data frame with the ribosome heterogeneity prediction results for all possible comparisons. It will also create comparison-specific ".csv" files to store your results in addition to a few informative QC plots. Most importantly, the pipeline will save a scatterplot and heatmap visualization summarizing your results.

Within this function, there exists a few other options as well. For example, you can specify which comparisons you want to focus on, or which samples to exclude from the analysis. 

<u>Please note that you have to specify the organism and align your reads strictly to given rRNA sequences when running the pipeline.</u> 

#### dripARF() function in detail

dripARF() wrapper function uses series of other dripARF functions to obtain your results. We summarize below what these functions are used for.

* ARF_check_organism()                : Checks if given organism is supported within ARF.
* read_ARF_samples_file()             : This function allows you to read the tab-seperated ARF-samples file, returning a dataframe.
* dripARF_read_rRNA_fragments()       : Read .bedgraph/.bam files to create the rRNA count data.
* dripARF_get_DESEQ_dds()             : The function that normalizes rRNA counts with DESEQ2.
* dripARF_report_RPset_group_counts() : Calculate the average read count of RP contact point sets.
* dripARF_predict_heterogenity()      : Overlap differential rRNA count data with 3d ribosome, rRNA-RP proximity data and predict heterogeneity candidates across groups.
* dripARF_simplify_results()          : Simplify dripARF results based on given thresholds.
* dripARF_result_heatmap()            : Draw heatmap of NES (ES1) and NES_randZscore (ES2) enrichments scores, where RPs are filtered based on given thresholds.
* dripARF_result_scatterplot()        : Draw comparison-specific scatter plots representing RPSEA and ORA analyses results.
* dripARF_report_RPspec_pos_results() : Obtain the positional differential abundance change values for RP proximity sets.
* dripARF_rRNApos_heatmaps()          : Draw rRNA fragment change heatmaps to visualize position-specific differential rRNA fragment abundances.
* dripARF_threshold_test()            : (Experimental function) This function allows you to run the whole dripARF pipeline with varying proximity thresholds.

### TESTING ARF

If you have installed the ARF package successfully, you can test it with the data provided under the test_data/ folder.

## Copyright

Copyright 2024 by Ferhat Alkan

ARF is released under the GNU General Public License.

GNU GENERAL PUBLIC LICENSE

This is a free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License, either version 3 of the License, or
(at your option) any later version. You should have received a copy of the GNU General Public License
along with ARF, see file LICENSE. If not, see <http://www.gnu.org/licenses/>.

This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citations

If you actively use the ARF platform, please cite the following publications:

Ferhat Alkan, Oscar Wilkins, Santiago Hernandez-Perez, Sofia Ramalho, Joana Silva, Jernej Ule,  William J Faller; Identifying ribosome heterogeneity using ribosome profiling data, *Nucleic Acids Research*, 2022; gkac484, [https://doi.org/10.1093/nar/gkac484](https://doi.org/10.1093/nar/gkac484)

Edwin Sakyi Kyei-Baffour, Jitske Bak, Joana Silva, William J Faller, Ferhat Alkan; Detecting ribosome collisions with differential rRNA fragment analysis in ribosome profiling data, 2024, (Under Review).

## Contact

For questions, queries, bug reports; please use the **<u>github issues tab</u>** or contact: <fallerlab@gmail.com> and/or <feralkan@gmail.com>.

