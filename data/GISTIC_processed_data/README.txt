This file contains the output files of the script also contained in this file. Each subtype directory has four output files:

	- genes_w_Gscore.tsv: List of protein coding genes with an associated G score (aMplification and deletions). It also contains all the genes that did not fall into any Gscore interval.
	- proccesed_SCNAs_summerised.tsv: List of genes with statistical summary of all the samples SCNA messures for each specifiC gene.
	- processed_SCNAs_by_gene_by_sample.tsv: Large format data of every SCNA messure for every sample for every protein coding gene.
	- SCNAs_processing_summary.tsv: summary of the SCNAs_processing_pipeline.R pipeline for each subtype.

Pipeline description:

- SCNAs_processing_pipeline.R: This pipeline summarizes the all_data_by_genes.txt file in the GISTIC_raw_data directory for each subtype, generating three files for each cancer subtype:

		Input:
		1) all_data_by_genes.txt

		Outputs:
		1) proccesed_SCNAs_summerised.tsv
		2) processed_SCNAs_by_gene_by_sample.tsv
		3) SCNAs_processing_summary.tsv

- Gscore_precessing_pipeline.R: This pipeline takes the list of human protein coding genes and the G score data for each subtype,  and it joins them, producing as output a list of protein coding genes with an associated G score. The pipeline produces one output for each cancer subtype:

		Input:
		1) scores.gistic (in GISTIC_raw_data directory)
		2) data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt

		Outputs:
		1) genes_w_Gscore.tsv

- GISTIC_filtering_pipeline.R: This pipeline filters out the non significan SCNA's found through the GISTIC algorithm for each subtype. The filtering criteria was that the product of the q-value times the number of genes with a q-value equal or lower to any given q-value was not higher than 1. This translates into the fact that we would expect by chance alone to have less than one false positive among the filtered SCNA's. Then out of the individual threshold q values calculated for each subtype, the smallest was chosen as a the threshold value for all four subtypes.

		Input:
		1) genes_w_Gscore.tsv

		Outputs:
		1) genes_w_Gscore_filtered.tsv

- SCNAs_compilation_pipeline.R: Makes a list of all unique SCNAs detected through the GISTIC2.0 algorithm in each of the subtypes and saves list of unique elements into healthy directory for its use in enrichment_analysis.R

		Input:
		1) *subtypes*/genes_w_Gscore_filtered.tsv

		Outputs:
		1) healthy/genes_w_Gscore_filtered.tsv

- SCNAs_compilation_pipeline_for_default_SCNAs.R: Same as above but for default SCNAs: Makes a list of all unique SCNAs detected through the GISTIC2.0 algorithm in each of the subtypes and saves list of unique elements into healthy directory for its use in enrichment_analysis_for_default_SCNAs.R

		Input:
		1) *subtypes*/default_SCNAs_w_q_values_filtered.tsv

		Outputs:
		1) healthy/default_SCNAs_w_q_values_filtered.tsv
