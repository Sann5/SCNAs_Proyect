################################################################################
# title: 'Enrichment Analysis'
# author: "CASTRO DAU Santiago Manuel"
# date: '19 MAR 2020'
# last modified: '19 MAR 2020'
# description: 
# 'Computes the fdr value in an enrichment analysis (one-sided Fisher's exact 
# test) for significantly amplified/deleted genes (GISTIC2.0 algorithm, 
# q_value = 2.8E-04)' for different network sizes of the corresponding subtype 
# MI correlation networks. Outputs two .rds figures plotting -log10(p value) vs.
# network size.
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(bc3net)
library(ggsci)

# Set working directory
setwd(dir = "/Users/Santiago/SCNA_Proyect")

# Loading proteing coding genes list to sample from
Biomart_path <- paste("data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt")

Biomart <- vroom(file = Biomart_path,
                 delim = "\t",
                 col_types = c(col_character(), col_double())) %>%
  pull(gene_stable_id) 

# Crating a list where we specify the subtypes and the network sizes. 
# The list will be used to iterate the pipeline for each of the subtypes and
# network sizes.
subtypes_list <- c("basal", "her2", "luma", "lumb", "healthy")
network_sizes <- seq(1332, 133170, by = 10000)

# Creating working object outside the loop
working_tibble <- tibble()

# 2nd Loop start
for(subtype in subtypes_list) {
  
  # Loading SCNAs
  path_SCNAs <- paste("data/GISTIC_processed_data/",
                      subtype,
                      "/genes_w_Gscore_filtered.tsv",
                      sep = "")
  
  SCNAs <- vroom(file = path_SCNAs,
                 delim = "\t",
                 col_types = c(col_character(), col_double())) %>% 
    select(gene_stable_id)
  
  # Loading network
  cancer_path <- paste("data/network_data/",
                       subtype,
                       "/",
                       subtype,
                       "_norm_133k_interactions.sif",
                       sep = "")
  
  cancer_net <- vroom(file = cancer_path,
                      delim = "\t",
                      col_types = c(col_character(), col_double()))
  
  # Network size loop start
  for(size in network_sizes) {
  
    # Filtering networks and creating list of unique elements in the network
    cancer_net_filtered <- cancer_net %>% 
      top_n(size, MI)
    
    cancer_nodes <- bind_rows(cancer_net_filtered %>% 
                                select(gene_stable_id = target),
                              cancer_net_filtered %>% 
                                select(gene_stable_id = source)) %>%
      distinct() %>% 
      pull(gene_stable_id)
    
    # Changing the name of the column in the current SCNA_in_loop object so
    # it identifies the iteration in the wroking_tibble object
    SCNAs_in_loop <- SCNAs
    names(SCNAs_in_loop)[names(SCNAs_in_loop) == "gene_stable_id"] <- paste(subtype,
                                                                            size,
                                                                            sep = "_")
    
    # Calculating and adding enrichment analysis to the working_tibble
    working_tibble <- bind_rows(working_tibble,
                                enrichment(genes = cancer_nodes,
                                 reference = Biomart,
                                 genesets = SCNAs_in_loop,
                                 adj = "fdr",
                                 verbose = FALSE)) 
    } 
}

# Making plot
  working_figure <- working_tibble %>% 
    separate(TermID, c("subtype", "net_size"), "_") %>% 
    mutate(net_size = as.double(net_size)) %>% 
    mutate(log_10_padj = -log10(padj)) %>% 
    ggplot(mapping = aes(y = log_10_padj,
                         x = net_size,
                         colour = subtype)) +
    geom_line(aes(color = subtype, linetype = subtype)) +
    geom_point(aes(color = subtype)) +
    labs(title = paste("-log10(p) vs network size for all four subtypes and the healthy network"),
         y = "-log10(p)",
         x = "Network Size",
         fill = "") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_minimal() +
    ggsci::scale_color_lancet() 
  
# Saving figures as R objects and png onject
saveRDS(working_figure,
        file = "results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh.rds")
  
ggsave(filename = "enrichment_line_grapgh.png",
       plot = working_figure,
       path = "results/part_1_presence_significance/enrichment_analysis/",
       scale = 1,
       device = "png",
       width = 7.08,
       height = 5.08,
       units = c("in"))

# Spliting graph
working_figure <- working_figure +
    facet_grid(subtype ~ ., scales = "free") + 
    theme(legend.position="none")

saveRDS(working_figure,
        file = "results/part_1_presence_significance/enrichment_analysis/enrichment_line_grapgh_splited.rds")
    
