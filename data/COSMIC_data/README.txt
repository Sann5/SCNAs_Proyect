This repository contains the raw data downloaded from the COSMIC database for breast cancer SCNA's. It aslo contains the script that is used to clean and adapt the data. The tissue_cnv Tue Nov 26 23_27_19 2019.csv (the file with the raw data) contaisn the summerised infomation of the number of patients with messured deletions and amlifications for every gene in the genome and the amount of samples tested. The results shown in this table are derived from high value (numeric) copy number data. Data where the copy number and minor allele information is unknown are excluded. Gain and loss (amplifications and deletions) are messured as follows: 

Gain:
average genome ploidy <= 2.7 AND total copy number >= 5
OR average genome ploidy > 2.7 AND total copy number >= 9

Loss:
average genome ploidy <= 2.7 AND total copy number = 0
OR average genome ploidy > 2.7 AND total copy number < ( average genome ploidy - 2.7 )

This file can be downloaded from this link: https://cancer.sanger.ac.uk/cosmic/browse/tissue?wgs=off&sn=breast&ss=all&hn=all&sh=all&in=t&src=tissue&all_data=n

This file contains all breast cancer samples, independently of the subtype or histological categorization. It is worth mentioning that there is no availabe SCNA for the specific subtypes which are stuided in this proyect.

Pipeline description:
	- COSMIC_SCNA_processing_pipeline.R: This script cleans the data from tissue_cnv Tue Nov 26 23_27_19 2019.csv and it joins it with the Biomart dicctionary

		Input:
		1) tissue_cnv Tue Nov 26 23_27_19 2019.csv
		2) Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt (in BioMart_data direcotry)

		Outputs:
		1) COSMIC_porcessed_SCNA.tsv
