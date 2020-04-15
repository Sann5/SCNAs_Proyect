################################################################################
# title: 'GISTIC filtering pipeline'
# author: "CASTRO DAU Santiago Manuel"
# date: '04 MAR 2020'
# last modified: '04 MAR 2020'
# description: 
# 'Filters out the non significant SCNA's found throught the GISTIC algorithm' 
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)

# Set working directory
setwd(dir = "/Users/Santiago/SCNA_Proyect")

# Crating a list where we specify the subtypes. The list will be used to iterate
# the pipeline for each of the subtypes.
subtypes_list <- c("basal", "her2", "luma", "lumb")

# Creating outside loop working object
working_tibble <- c()

# Previous loop to determine threshold q vlaue
for(subtype in subtypes_list) {
  
  # Set file path
  path_SCNAs <- paste("data/GISTIC_processed_data/",
                      subtype,
                      "/genes_w_Gscore.tsv",
                      sep = "")
  
  # Read file and clean names
  SCNAs <- vroom(file = path_SCNAs,
                 delim = "\t",
                 col_types = c(col_character(), col_double()))
  
  # Filtering out SCNAs so that the cumulative frequency times its respective 
  # q value dose not surpasre 1 (less that one fals epostivie per GISTIC set)
  SCNAs <- SCNAs %>% 
    filter(!is.na(g_score)) %>% 
    mutate(q_value = 1/(10**log10_q_value)) %>%
    arrange(q_value) %>%
    rownames_to_column(var = "cum_frequency") %>%
    mutate_at("cum_frequency", ~as.double(.)) %>%
    mutate(q_x_cum_freq = q_value*cum_frequency) %>%
    filter(q_x_cum_freq < 1)
  
  working_tibble <- c(working_tibble,
                      SCNAs %>% 
                        pull(q_value) %>% 
                        max(., na.rm = FALSE))
}

# Asigning variable for threshold q value
q_val_th <- min(working_tibble)

# Loop start
for(subtype in subtypes_list) {
  
  # Set file path
  path_SCNAs <- paste("data/GISTIC_processed_data/",
                      subtype,
                      "/genes_w_Gscore.tsv",
                      sep = "")
  
  # Read file and clean names
  SCNAs <- vroom(file = path_SCNAs,
                 delim = "\t",
                 col_types = c(col_character(), col_double()))
  
  # Filtering out non significan aplifications/deletions
  # and counting number of aplifications/deletions
  SCNAs <- SCNAs %>% 
    filter(!is.na(g_score)) %>% 
    mutate(q_value = 1/(10**log10_q_value)) %>%
    filter(q_value <= q_val_th) %>% 
    distinct(gene_stable_id, type, .keep_all = TRUE) %>% 
    mutate(is_ = "is_GISTIC")
  
  # Writing filtered genes with GISTIC SCNA's
  filtered_path <- paste("data/GISTIC_processed_data/",
                        subtype,
                        "/genes_w_Gscore_filtered.tsv",
                        sep = "")
  
  vroom_write(SCNAs,
              delim = "\t",
              path = filtered_path)
}  
