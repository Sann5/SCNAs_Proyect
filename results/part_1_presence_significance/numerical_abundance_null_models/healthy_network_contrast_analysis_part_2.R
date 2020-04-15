################################################################################
# title: 'GISTIC healthy network contrast analysis second part'
# author: "CASTRO DAU Santiago Manuel"
# date: '04 MAR 2020'
# last modified: '10 MAR 2020'
# description: 
# 'Adds line graphs fo the ratio of intersecting SCNAs and the nodes of each
# network' 
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)
library(ggsci)

# Set working directory
setwd(dir = "/Users/Santiago/SCNA_Proyect")

# Crating working tibble outside loop
working_tibble <- tibble()

# Loading healthy network
healthy_path <- paste("data/network_data/healthy/healthy_norm_133k_interactions.sif")

healthy <- vroom(file = healthy_path,
                 delim = "\t",
                 col_types = c(col_character(), col_double())) %>%
  mutate(is_ = "is_network")

# Loading proteing coding genes list to sample from
Biomart_path <- paste("data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt")

Biomart <- vroom(file = Biomart_path,
                 delim = "\t",
                 col_types = c(col_character(), col_double())) %>%
  select(gene_stable_id) %>%
  mutate(is_ = "is_Biomart")

# Creating list of network sizes to filter
network_sizes <- c(1332, 13317, 133170)
size <- 1332
# 1st Loop start
for(size in network_sizes) {
  
  # Crating a list where we specify the subtypes. The list will be used to iterate
  # the pipeline for each of the subtypes.
  subtypes_list <- c("basal", "her2", "luma", "lumb")
  subtype <- "basal"
  # 2nd Loop start
  for(subtype in subtypes_list) {
    
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
    
    # Filtering networks and creating list of unique elements in the network
    cancer_net <- cancer_net %>% 
      top_n(size, MI)
    
    healthy_net <- healthy %>%
      top_n(size, MI)
    
    cancer_nodes <- bind_rows(cancer_net %>% 
                                select(gene_stable_id = target),
                              cancer_net %>% 
                                select(gene_stable_id = source)) %>%
      distinct() %>%
      mutate(is_ = "is_network")
    
    n_cancer_nodes <- cancer_nodes %>% 
      group_by(is_) %>% 
      summarise(count = n()) %>% 
      pull(count)
    
    healthy_nodes <- bind_rows(healthy_net %>% 
                                 select(gene_stable_id = target, is_),
                               healthy_net %>% 
                                 select(gene_stable_id = source, is_)) %>%
      distinct()
    
    n_healthy_nodes <- healthy_nodes %>% 
      group_by(is_) %>% 
      summarise(count = n()) %>% 
      pull(count)
    
    # Loading SCNAs
    path_SCNAs <- paste("data/GISTIC_processed_data/",
                        subtype,
                        "/genes_w_Gscore_filtered.tsv",
                        sep = "")
    
    SCNAs <- vroom(file = path_SCNAs,
                   delim = "\t",
                   col_types = c(col_character(), col_double()))
    
    n_SCNAs <- SCNAs %>% 
      group_by(is_) %>% 
      summarise(count = n()) %>% 
      pull(count)
    
    # Joining SCNAs with nodes
    cancer_intersec <- inner_join(x = cancer_nodes %>% 
                                    select(-is_),
                                  y = SCNAs %>% 
                                    select(gene_stable_id) %>% 
                                    mutate(is_ = "is_GISTIC"),
                                  by = c("gene_stable_id" = "gene_stable_id"))
    
    healthy_intersec <- inner_join(x = healthy_nodes %>% 
                                     select(-is_),
                                   y = SCNAs %>% 
                                     select(gene_stable_id) %>% 
                                     mutate(is_ = "is_GISTIC"),
                                   by = c("gene_stable_id" = "gene_stable_id"))
    
    working_tibble <- bind_rows(working_tibble,
                                cancer_intersec %>% 
                                  group_by(is_) %>%
                                  summarise(count = n()/n_cancer_nodes) %>%
                                  mutate(net_type = "cancer") %>% 
                                  mutate(subtype = subtype) %>%
                                  mutate(net_size = size),
                                healthy_intersec %>% 
                                  group_by(is_) %>%
                                  summarise(count = n()/n_healthy_nodes) %>%
                                  mutate(net_type = "healthy") %>% 
                                  mutate(subtype = subtype) %>%
                                  mutate(net_size = size))
    
    # Sampleing proteing coding genes from Biomart and joining them with nodes
    # 3rd loop start
    for (i in 1:100){
      intersection_2 <- inner_join(x = healthy_nodes %>% 
                                     select(-is_),
                                   y = Biomart %>% 
                                     select(gene_stable_id, is_) %>% 
                                     sample_n(n_SCNAs, replace = TRUE),
                                   by = c("gene_stable_id" = "gene_stable_id"))
      
      # Writing sampling found genes to outside object
      working_tibble <- bind_rows(working_tibble,
                                  intersection_2 %>% 
                                    group_by(is_) %>%
                                    summarise(count = n()/n_healthy_nodes) %>%
                                    mutate(subtype = subtype) %>%
                                    mutate(net_size = size) %>% 
                                    mutate(net_type = "healthy_random") %>% 
                                    mutate(is_ = "is_GISTIC"))
    }
  } 
}

# Loop start to figures
for (subtype in subtypes_list) {
  
  working_figure <- working_tibble %>%
    rename(subtype_ = subtype) %>% 
    filter(subtype_ == subtype) %>% 
    mutate(ratio = count*100) %>% 
    group_by(subtype_, net_size, net_type) %>% 
    summarise(median = median(ratio),
              sd = sd(ratio)) %>% 
    ggplot(mapping = aes(y = median,
                         x = net_size,
                         color = net_type)) +
    geom_errorbar(aes(ymin=median-sd, ymax=median+sd),
                  width=.01) +
      geom_line(aes(color = net_type, linetype = net_type)) +
      geom_point(aes(color = net_type)) +
    labs(title = paste("Ratio of ",
                       subtype,
                       "SCNAs that intersect with networks",
                       sep = ""),
         y = "Percentage %",
         x = "Network Size") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_startrek() +
    theme_minimal()

  # Saving R objects
  saveRDS(working_figure,
          file = paste("results/part_1_presence_significance/numerical_abundance_null_models/",
                       subtype,
                       "_healthy_ratio.rds",
                       sep = ""))
  
}
