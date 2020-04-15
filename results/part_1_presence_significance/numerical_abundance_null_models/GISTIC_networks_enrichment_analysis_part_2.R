################################################################################
# title: 'GISTIC networks enrichment analysis part 2'
# author: "CASTRO DAU Santiago Manuel"
# date: '06 MAR 2020'
# last modified: '06 MAR 2020'
# description: 
# 'Generates a line grapgh that quantifies the ratio of SCNA compared to the 
# ratio of ramndomly smpled genes' 
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)

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
size <- 1332
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
    
    net_nodes <- nodes %>% 
      group_by(is_) %>% 
      summarise(count = n()) %>% 
      pull(count)
    
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
                                  mutate(net_size = size) %>% 
                                  mutate(net_nodes = net_nodes))
    
    # Sampleing proteing coding genes from Biomart and joining them with nodes
    # 3rd loop start
    for (i in 1:100){
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
                                    mutate(is_ = "is_random_sample") %>% 
                                    mutate(net_nodes = net_nodes))
    }
  } 
}

# Crating and saving figures
for(subtype in subtypes_list) {
  
working_figure <- working_tibble %>%
  rename(subty = subtype) %>%
  filter(subty == subtype) %>% 
  mutate(ratio = count*100/net_nodes) %>% 
  group_by(is_, subty, net_size) %>% 
  summarise(median = median(ratio),
            sd = sd(ratio)) %>% 
  ggplot(mapping = aes(y = median,
                       x = net_size,
                       colour = is_)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), 
                width=.01) +
  geom_line(aes(color = is_, linetype = is_)) +
  geom_point(aes(color = is_)) +
  labs(title = paste("Ratio of",
                     subtype,
                     "SCNA's and randomly sampled genes in networks"),
       y = "Percentage %",
       x = "Network Size",
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  theme_minimal()

  # # Saving figure
  # ggsave(filename = paste(subtype,
  #                         "_network_enrichment_ratio_fig.png",
  #                         sep = ""),
  #        plot = working_figure,
  #        path = "results/part_1_presence_significance/numerical_abundance_null_models/",
  #        scale = 1,
  #        device = "png",
  #        width = 8.08,
  #        height = 5.08,
  #        units = c("in"))
  
  # Saving figures as R objects
  saveRDS(working_figure,
          file = paste("results/part_1_presence_significance/numerical_abundance_null_models/",
                       subtype,
                       "_network_enrichment_ratio.rds",
                       sep = ""))
  
}
