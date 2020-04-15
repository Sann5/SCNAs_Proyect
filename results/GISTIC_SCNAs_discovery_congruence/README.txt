Analysis description: 

- GISTIC_SCNAs_discovery_congruence_script.R: Extracts the SCNAs reported in the del_genes.conf_99.txt and amp_genes.conf_99.txt files for each subtype, extracts the q_value and residual q_values and compares the items in this list to the list manually obtained from other GISTIC output files. It generates a summary of both the q values that are used by default by GISTIC and a set analysis of the genes that a are shared in both pipelines.

		Input:
		1) data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt
		2) data/GISTIC_raw_data/*subtype*/amp_genes.conf_99.txt
		3) data/GISTIC_raw_data/*subtype*/del_genes.conf_99.txt
		4) data/GISTIC_processed_data/*subtype*/genes_w_Gscore.tsv

		Outputs:
		1) results/GISTIC_SCNAs_discovery_congruence/default_q_values.tsv
		2) results/GISTIC_SCNAs_discovery_congruence/gene_sets_comparison_2.5E-01.tsv
		3) data/GISTIC_processed_data/*subtype*/default_SCNAs_w_q_values.tsv

- lower_q_value_congruence_analysis.R: Performs two primary tasks: 1) creates filtered versions of the default SCNAs gene sets for each subtype and 2) performed a comparison of the manual curated GISTIC SCNAs set and the default GISTIC SCANs set (identical to the on that is outputted in GISTIC_SCNAs_discovery_congruence_script.R, but here the q value threshold is 3.7E-04, instead of 0.25.

		Input:
		1) data/GISTIC_processed_data/*subtype*/default_SCNAs_w_q_values.tsv
		2) data/GISTIC_processed_data/*subtype*/genes_w_Gscore.tsv

		Outputs:
		1) results/GISTIC_SCNAs_discovery_congruence/gene_sets_comparison_3.7E-04.tsv
		2) data/GISTIC_processed_data/*subtype*/default_SCNAs_w_q_values_filtered.tsv

- manual_default_SCNAs_comparison_report.Rmd: Reportes the findings of the previous analysis and of the enrichment analysis. (Report in spanish).


		Input:
		1) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh.rds
		2) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh_splited.rds
		3) results/GISTIC_SCNAs_discovery_congruence/gene_sets_comparison_2.5E-01.tsv
		4) results/GISTIC_SCNAs_discovery_congruence/gene_sets_comparison_3.7E-04.tsv
		5) results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh_default.rds

		Outputs:
		1) manual_default_SCNAs_comparison_report.html