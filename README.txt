This directory contains all the scripts and data necessary to produce generate the SCNA_Proyect analysis.

To run the analysis you must:

1) Make sure that all the raw data is contained in the appropriate directory inside the SCNA_Proyect repository. You only need to manually add the files marked with ***.

2) Pre install all the necessary R libraries.

3) Run the scripts in the specified order.

-----------------------------------------------------------------------------------------------------------------------

1) The necessary raw data which is needed to run the analyses, as well as their location in this repository is indicated below. The README.txt file in each of these directories contains the information of the nature and origin of the raw data. 

data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12.txt
data/COSMIC_data/tissue_cnv Tue Nov 26 23_27_19 2019.csv
data/GISTIC_raw_data/basal/all_data_by_genes.txt
data/GISTIC_raw_data/her2/all_data_by_genes.txt
data/GISTIC_raw_data/luma/all_data_by_genes.txt
data/GISTIC_raw_data/lumb/all_data_by_genes.txt
data/GISTIC_raw_data/basal/scores.gistic
data/GISTIC_raw_data/her2/scores.gistic
data/GISTIC_raw_data/luma/scores.gistic
data/GISTIC_raw_data/lumb/scores.gistic
data/GISTIC_raw_data/basal/amp_genes.conf_99.txt
data/GISTIC_raw_data/her2/amp_genes.conf_99.txt
data/GISTIC_raw_data/luma/amp_genes.conf_99.txt
data/GISTIC_raw_data/lumb/amp_genes.conf_99.txt
data/GISTIC_raw_data/basal/del_genes.conf_99.txt
data/GISTIC_raw_data/her2/del_genes.conf_99.txt
data/GISTIC_raw_data/luma/del_genes.conf_99.txt
data/GISTIC_raw_data/lumb/del_genes.conf_99.txt
data/network_data/basal/deg-basal.tsv
data/network_data/her2/deg-her2.tsv
data/network_data/luma/deg-luma.tsv
data/network_data/lumb/deg-lumb.tsv
data/network_data/basal/basal.sif.sorted  ***
data/network_data/her2/her2.sif.sorted  ***
data/network_data/luma/luma.sif.sorted  ***
data/network_data/lumb/lumb_norm.sif.sorted  ***

 ***  WARNING ******************************************************************************************************

The marked files (***) (the input files for the sif_processing_pipeline) data files are not contained within these directories beacuse of their size, they need to be added localy for the scripts to run. They can be found in the following directory of INMEGENS's remote server:

username@castillo.cluster.inemen.gob.mx:/home/labs/csbig/subtipos_mama_2018/networks/

	- her2.sif.sorted
	- basal.sif.sorted
	- luma.sif.sorted
	- lumb.sif.sorted
	- healthy.sif.sorted

********************************************************************************************************************

2) The scripts were designed in RStudio Version 1.1.463, and will make use of the following packages:

install.packages(“tidyverse”)
install.packages(“vroom”)
install.packages(“janitor”)
install.packages(“stringr”)
install.packages(“bc3net”)
install.packages(“ggsci”)
install.packages(“scales”)
install.packages(“cowplot”)

3) Once all of the raw data is in place, you can run the scripts in the specified order (to generate all the intermediate data that will be used for the analysis). Also before running any of the scripts the setwd(dir = “”) command on each of the individual scripts must be modified to match the local location of the SCNA_Proyect directory. 
 
data/BioMart_data/Biomart_processing_pipeline.R
data/COSMIC_data/COSMIC_SCNA_processing_pipeline.R
data/GISTIC_processed_data/SCNAs_processing_pipeline.R
data/GISTIC_processed_data/Gscore_precessing_pipeline.R
data/GISTIC_processed_data/GISTIC_filtering_pipeline.R
data/GISTIC_processed_data/SCNAs_compilation_pipeline.R
data/GISTIC_processed_data/SCNAs_compilation_pipeline_for_default_SCNAs.R
data/network_data/sif_processing_pipeline.R

Finally to run the key analysis scripts that have been developed so far run the following scripts in the specified order (an explanation of what the script dose and its inputs and outputs is available in the README.txt that is in each of the specified directories of the respective directory):

results/part_1_presence_significance/enrichment_analysis/enrichment_analysis.R
results/part_1_presence_significance/enrichment_analysis/enrichment_analysis_for_default_SCNAs.R
results/GISTIC_SCNAs_discovery_congruence/GISTIC_SCNAs_discovery_congruence_script.R
results/GISTIC_SCNAs_discovery_congruence/lower_q_value_congruence_analysis.R
results/GISTIC_SCNAs_discovery_congruence/manual_default_SCNAs_comparison_report.Rmd
