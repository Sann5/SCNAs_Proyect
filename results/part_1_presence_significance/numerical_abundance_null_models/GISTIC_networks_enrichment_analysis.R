################################################################################
# title: 'GISTIC networks enrichment analysis'
# author: "CASTRO DAU Santiago Manuel"
# date: '04 MAR 2020'
# last modified: '04 MAR 2020'
# description: 
# 'This script analysises whether the GISTIC SCNA's are significantly enriched 
# in out networks, independently of network size         ' 
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

# Loading proteing coding genes list to sample from
Biomart_path <- paste("data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt")

Biomart <- vroom(file = Biomart_path,
                 delim = "\t",
                 col_types = c(col_character(), col_double())) %>%
  select(gene_stable_id) %>%
  mutate(is_ = "is_Biomart")

# Creating list of network sizes to filter
network_sizes <- c(1332, 13317, 133170)

# 1st Loop start
for(size in network_sizes) {
  
  # Crating a list where we specify the subtypes. The list will be used to iterate
  # the pipeline for each of the subtypes.
  subtypes_list <- c("basal", "her2", "luma", "lumb")
  
  # 2nd Loop start
  for(subtype in subtypes_list) {
    
    # Loading network
    net_path <- paste("data/network_data/",
                      subtype,
                      "/",
                      subtype,
                      "_norm_133k_interactions.sif",
                      sep = "")
    
    net <- vroom(file = net_path,
                 delim = "\t",
                 col_types = c(col_character(), col_double()))
    
    # Filtering network and creating list of unique elements in the network
    net <- net %>% 
      top_n(size, MI)
    
    nodes <- bind_rows(net %>% select(gene_stable_id = target),
                       net %>% select(gene_stable_id = source)) %>%
      distinct() %>%
      mutate(is_ = "is_node")
    
    # Loading SCNAs an counting number of SCNA's for each subtype
    path_SCNAs <- paste("data/GISTIC_processed_data/",
                        subtype,
                        "/genes_w_Gscore_filtered.tsv",
                        sep = "")

    SCNAs <- vroom(file = path_SCNAs,
                   delim = "\t",
                   col_types = c(col_character(), col_double()))
    
    n_nodes <- SCNAs %>% 
      group_by(is_) %>% 
      summarise(count = n()) %>% 
      pull(count)
    
    # Joining SCNAs with nodes
    intersection <- inner_join(x = nodes %>% select(-is_),
                               y = SCNAs %>% select(gene_stable_id, is_),
                               by = c("gene_stable_id" = "gene_stable_id"))
    
    working_tibble <- bind_rows(working_tibble,
                                intersection %>% 
                                  group_by(is_) %>%
                                  summarise(count = n()) %>%
                                  mutate(subtype = subtype) %>%
                                  mutate(net_size = size),
                                nodes %>% 
                                  group_by(is_) %>%
                                  summarise(count = n()) %>%
                                  mutate(subtype = subtype) %>%
                                  mutate(net_size = size))
    
    # Sampleing proteing coding genes from Biomart and joining them with nodes
    # 3rd loop start
    for (i in 1:3){
      intersection_2 <- inner_join(x = nodes %>% 
                                     select(-is_),
                                   y = Biomart %>% 
                                     select(gene_stable_id, is_) %>% 
                                     sample_n(n_nodes, replace = TRUE),
                                   by = c("gene_stable_id" = "gene_stable_id"))
      
      # Writing sampling found genes to outside object
      working_tibble <- bind_rows(working_tibble,
                                  intersection_2 %>% 
                                    group_by(is_) %>%
                                    summarise(count = n()) %>%
                                    mutate(subtype = subtype) %>%
                                    mutate(net_size = size) %>% 
                                    mutate(is_ = paste("is_Biomart",
                                                       "_",
                                                       i,
                                                       sep = "")))
    }
  } 
}

# Crating and saving figures
for(size in network_sizes) {
working_figure <- working_tibble %>% 
  filter(net_size == size) %>% 
  mutate(is_ = as_factor(is_)) %>%
  mutate(is_ = fct_relevel(is_, 
                           "is_node",
                           "is_GISTIC",
                           "is_Biomart_1",
                           "is_Biomart_2",
                           "is_Biomart_3")) %>%
  ggplot(mapping = aes(fill = is_,
                       y = count,
                       x = as.factor(subtype))) + 
  geom_bar(position = "dodge",
           stat = "identity") +
  labs(title = paste("Number of SCNA's present in network of size ",
                     size,
                     sep = ""),
       y = "Number of genes",
       x = "Subtypes",
       fill = "") + 
  scale_fill_discrete(labels = c("In netowrk", 
                                 "In network & in GISTIC",
                                 "Random sample No.1",
                                 "Random sample No.2",
                                 "Random sample No.3")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() 

# ggsave(filename = paste("GISTIC_network_",
#                         size,
#                         "_fig.png"),
#        plot = working_figure,
#        path = "results/part_1_presence_significance/numerical_abundance_null_models/",
#        scale = 1,
#        device = "png",
#        width = 7.08,
#        height = 5.08,
#        units = c("in"))

saveRDS(working_figure,
        file = paste("results/part_1_presence_significance/numerical_abundance_null_models/",
                     "GISTIC_network_",
                     size,
                     ".rds",
                     sep = ""))
}
