Analysis description: 

- enrichment_analysis.R: Computes the fdr value in an enrichment analysis (one-sided Fisher's exact test) for significantly amplified/deleted genes (GISTIC2.0 algorithm, q_value = 2.8E-04)' for different network sizes of the corresponding subtype MI correlation networks. Outputs two .rds figures plotting -log10(p value) vs. network size. The analysis also uses the complete set of detected SCNAs and probes weather these genes are enriched in the healthy network. This was included because the is a possibility that these SCNAs are enriched in the network just because they are important regulatory elements independently of their SCNA status. 

		Input:
		1) data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt
		2) data/GISTIC_processed_data/*subtype*/genes_w_Gscore_filtered.tsv
		3) data/network_data/*subtype*/*subtype*_norm_133k_interactions.sif

		Outputs:
		1) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh.rds
		2) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh.png
		3) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh_splited.rds

enrichment_analysis_for_default_SCNAs.R: Same as above but uses default SCNAs GISTIC set with a q value threshold of 3.7E-04

		Input:
		1) data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt
		2) data/GISTIC_processed_data/*subtype*/default_SCNAs_w_q_values_filtered.tsv
		3) data/network_data/*subtype*/*subtype*_norm_133k_interactions.sif

		Outputs:
		1) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh_default.rds
		2) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh_default.png
		3) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh_splited_default.rds