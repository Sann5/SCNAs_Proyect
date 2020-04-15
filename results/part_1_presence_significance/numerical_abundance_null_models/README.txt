Analysis description: 

- COSMIC_GISTIC_congruence_analysis.R: This script analysises the congruence between the SCNA's anotated in the COSMIC database and the significanly aplified/deleted genes detected by the GISTIC2.0 algorithm. It generates a bar graph witht the number of genes found by GISTIC and the number of genes that are present both in the most significan aplifications/deletions in COSMIC nad in GISTIC (for each subtype).
		
		Input:
		1) data/COSMIC_data/COSMIC_porcessed_SCNA.tsv
		2) data/GISTIC_processed_data/subtype/genes_w_Gscore.tsv

		Outputs:
		1) GISTIC_COSMIC_fig.rds

- GISTIC_networks_enrichment_analysis.R: This script analysises whether the GISTIC SCNA's are significantly enriched in our networks, independently of network size. It generates a bar grapgh for each network size (in terms of edges) where each graph has the information of the amount of genes in the network, the amount of genes in the network with significant SCNA's and three random sample set fo protein codeing genes which we use to graphically compare weather the presence of these genes is above random ocurrences (for each subtype). 
		
		Input:
		1) data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt
		2) data/GISTIC_processed_data/subtype/genes_w_Gscore_filtered.tsv
		3) data/network_data/subtype/subtype_norm_133k_interactions.sif"

		Outputs:
		1) GISTIC_network_ 1332 _fig.rds
		2) GISTIC_network_ 13317 _fig.rds
		3) GISTIC_network_ 133170 _fig.rds

- GISTIC_networks_enrichment_analysis_part_2.R: This script is similar to its predecesor but instead of generating a bargrapgh for each network size it generates line graphs for each subtype where we plot the ratio of the intersection between the network and the GISTIC set for each network size, along with the intersection of random samles of genes of the same size as the GISTIC set.  
		
		Input:
		1) data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt
		2) data/GISTIC_processed_data/subtype/genes_w_Gscore_filtered.tsv
		3) data/network_data/subtype/subtype_norm_133k_interactions.sif"

		Outputs:
		1) basal_network_enrichment_ratio.rds
		2) her2_network_enrichment_ratio.rds
		3) luma_network_enrichment_ratio.rds
		4) lumb_network_enrichment_ratio.rds

- healthy_network_contrast_analysis.R: This analysis probes weather the significant SCNA found through the GISTIC algorithm are equaly present in the healthy network as compared to their presence in the cancerous networks of the same edge size. If we observe that this is not so, it could mean that the fact that they are geting more presence in the cancerus networks could be due to the fact that they are geting copy number mutated. It counts the number of SCNA's found through GISTIC in all the cancerous networks and in the healthy network.

*This script was used to develope healthy_network_contrast_analysis_part_2.R but it no longer holds analyisis value. It is kept in the directory in case of later use but its output files have been deleted from the directory. 
		
		Input:
		1) data/network_data/healthy/healthy_norm_133k_interactions.sif
		2) data/GISTIC_processed_data/subtype/genes_w_Gscore_filtered.tsv
		3) data/network_data/subtype/subtype_norm_133k_interactions.sif"

		Outputs:
		1) GISTIC_healthy_network_ 1332 _fig.png
		2) GISTIC_healthy_network_ 13317 _fig.png
		3) GISTIC_healthy_network_ 133170 _fig.png

- healthy_network_contrast_analysis_part_2.R: This analyisis makes a line graph of the ratio of intersecting SCNAs and the nodes of each network for each subtype (same as GISTIC_networks_enrichment_analysis_part_2.R), plus a a line which represents the intersection of the GSITIC set with the healthy network, and one more line whch represents the intersection of random sets of genes with the healthy network.
		
		Input:
		1) data/network_data/healthy/healthy_norm_133k_interactions.sif
		2) data/GISTIC_processed_data/subtype/genes_w_Gscore_filtered.tsv
		3) data/network_data/subtype/subtype_norm_133k_interactions.sif"

		Outputs:
		1) basal_healthy_ratio.rds
		2) her2_healthy_ratio.rds
		3) luma_healthy_ratio.rds
		4) lumb_healthy_ratio.rds

- numerical_abundance_analysis_report.Rmd: This is the R Markdown report of all the previously mentiones analysis. 

		Input:
		1) GISTIC_COSMIC_fig.rds
		2) GISTIC_network_ 1332 _fig.rds
		3) GISTIC_network_ 13317 _fig.rds
		4) GISTIC_network_ 133170 _fig.rds
		5) basal_network_enrichment_ratio.rds
		6) her2_network_enrichment_ratio.rds
		7) luma_network_enrichment_ratio.rds
		8) lumb_network_enrichment_ratio.rds
		9) data/GISTIC_processed_data/(basal, her2, luma, lumb)/genes_w_Gscore_filtered.tsv
		10) data/network_data/(basal, her2, luma, lumb)/(basal, her2, luma, lumb)_norm_133k_interactions.sif
		11) basal_healthy_ratio.rds
		12) her2_healthy_ratio.rds
		13) luma_healthy_ratio.rds
		14) lumb_healthy_ratio.rds
              

		Outputs:
		1) merical_abundance_analysis_report.html


