################################################################################
# title: 'GISTIC healthy network contrast analysis'
# author: "CASTRO DAU Santiago Manuel"
# date: '04 MAR 2020'
# last modified: '04 MAR 2020'
# description: 
# 'This analysis probes weather the significant SCNA found through the GISTIC  
# algorithm are equaly present in the healthy network. If we observe that this
# is not so it could mean that the fact that they are geting more presence in 
# the cancerus network could be due to the fact that they are geting copy 
# number mutated.' 
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)

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
                                  summarise(count = n()) %>%
                                  mutate(subtype = subtype) %>%
                                  mutate(net_size = size),
                                cancer_nodes %>% 
                                  group_by(is_) %>%
                                  summarise(count = n() - cancer_intersec %>% 
                                              group_by(is_) %>%
                                              summarise(count = n()) %>% 
                                              pull(count)) %>%
                                  mutate(subtype = subtype) %>%
                                  mutate(net_size = size))

  } 
  
  # Writing healthy data for the network outsid the loop
  working_tibble <- bind_rows(working_tibble,
                              healthy_intersec %>% 
                                group_by(is_) %>%
                                summarise(count = n()) %>%
                                mutate(subtype = "healthy") %>%
                                mutate(net_size = size),
                              healthy_nodes %>% 
                                group_by(is_) %>%
                                summarise(count = n() - healthy_intersec %>% 
                                            group_by(is_) %>%
                                            summarise(count = n()) %>% 
                                            pull(count)) %>%
                                mutate(subtype = "healthy") %>%
                                mutate(net_size = size))
  
}

# Crating figures
for(size in network_sizes) {
  working_figure <- working_tibble %>% 
      filter(net_size == size) %>% 
      mutate(is_ = as_factor(is_)) %>%
      mutate(is_ = fct_relevel(is_, 
                               "is_network",
                               "is_GISTIC")) %>%
    mutate(subtype = fct_relevel(subtype, 
                             "healthy",
                             "basal",
                             "her2",
                             "luma",
                             "lumb")) %>%
      ggplot(mapping = aes(fill = is_,
                           y = count,
                           x = subtype)) + 
      geom_bar(position = "stack",
               stat = "identity") +
      labs(title = paste("Number of SCNA's in Healthy and Neoplasic Networks of Size ",
                         size,
                         sep = ""),
           y = "Number of genes",
           x = "Subtypes",
           fill = "") + 
      scale_fill_discrete(labels = c("Genes in network", 
                                     "SCNA's in network")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_minimal()
  
  # Saving figure
  ggsave(filename = paste("GISTIC_healthy_network_",
                          size,
                          "_fig.png"),
         plot = working_figure,
         path = "results/part_1_presence_significance/numerical_abundance_null_models/",
         scale = 1,
         device = "png",
         width = 7.08,
         height = 5.08,
         units = c("in"))

  
}
