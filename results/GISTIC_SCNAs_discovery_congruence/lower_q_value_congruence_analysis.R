################################################################################
# title: 'Lower q value congruence script'
# author: "CASTRO DAU Santiago Manuel"
# date: '24 MAR 2020'
# last modified: '24 MAR 2020'
# description: 
# 'Performs two primary tasks: 1) creates filtered versions of the default SCNAs
# gene sets for each subtype and 2) performed a comparison of the manual curated
# GISTIC SCNAs set and the default GISTIC SCANs set (identical to the on that is
# outputted in GISTIC_SCNAs_discovery_congruence_script.R, but here the q value
# threshold is 3.7E-04, instead of 0.25.'
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)

# Set working directory
setwd(dir = "/Users/Santiago/SCNA_Proyect")

# Creating working object outside the loop
working_tibble <- tibble()

# Crating a list where we specify the subtypes.
subtypes_list <- c("basal", "her2", "luma", "lumb")
subtype <- "basal"

# Creating outside loop working object
working_list <- c()

# Previous loop to determine threshold q vlaue
for(subtype in subtypes_list) {
  
  # Set file path
  path_def_SCNAs <- paste("data/GISTIC_processed_data/",
                          subtype,
                          "/default_SCNAs_w_q_values.tsv",
                          sep = "")
  
  # Read file and clean names
  def_SCNAs <- vroom(file = path_def_SCNAs,
                     delim = "\t",
                     col_types = c(col_character(), col_double()))
  
  # Filtering out SCNAs so that the cumulative frequency times its respective 
  # q value dose not surpasre 1 (less that one fals epostivie per GISTIC set)
  def_SCNAs <- def_SCNAs %>% 
    arrange(q_value) %>%
    rownames_to_column(var = "cum_frequency") %>%
    mutate_at("cum_frequency", ~as.double(.)) %>%
    mutate(q_x_cum_freq = q_value*cum_frequency) %>%
    filter(q_x_cum_freq < 1)
  
  working_list <- c(working_list,
                    def_SCNAs %>% 
                      pull(q_value) %>% 
                      max(., na.rm = FALSE))
}

# Asigning variable for threshold q value
q_val_th <- min(working_list)

# Loop start
for(subtype in subtypes_list) {
  
  # Set file path for manually curated SCNAs
  path_man_SCNAs <- paste("data/GISTIC_processed_data/",
                          subtype,
                          "/genes_w_Gscore.tsv",
                          sep = "")
  
  # Read file and clean names for manually curated SCNAs
  man_SCNAs <- vroom(file = path_man_SCNAs,
                     delim = "\t",
                     col_types = c(col_character(), col_double()))
  
  # Filtering out non significan aplifications/deletions
  # and counting number of aplifications/deletions 
  # for manually curated SCNAs
  man_SCNAs <- man_SCNAs %>%
    filter(!is.na(g_score)) %>% 
    mutate(q_value = 1/(10**log10_q_value)) %>%
    distinct(gene_stable_id, type, .keep_all = TRUE) %>% 
    filter(q_value <= q_val_th)
  
  # Set file path for default SCNAs
  path_def_SCNAs <- paste("data/GISTIC_processed_data/",
                          subtype,
                          "/default_SCNAs_w_q_values.tsv",
                          sep = "")
  
  # Read file and clean names
  def_SCNAs <- vroom(file = path_def_SCNAs,
                     delim = "\t",
                     col_types = c(col_character(), col_double()))
  
  # Filtering out non significan aplifications/deletions
  # and counting number of aplifications/deletions 
  # for manually curated SCNAs
  def_SCNAs <- def_SCNAs %>%
    filter(q_value <= q_val_th)
  
  # Writtig default FILTERED SNACs to file in data/GISTIC_processed_data/subtype
  path_def_filt_SCNAs <- paste("data/GISTIC_processed_data/",
                               subtype,
                               "/default_SCNAs_w_q_values_filtered.tsv",
                               sep = "")
  
  vroom_write(def_SCNAs,
              delim = "\t",
              path = path_def_filt_SCNAs)
  
  # Compare
  intersection <- inner_join(x = def_SCNAs %>% 
                               select(gene_stable_id),
                             y = man_SCNAs %>% 
                               select(gene_stable_id),
                             by = c("gene_stable_id" = "gene_stable_id")) %>% 
    mutate(is_ = "is_both") %>% 
    group_by(is_) %>% 
    summarise(is_both = n()) %>% 
    select(is_both)
  
  left <- left_join(x = def_SCNAs %>% 
                      select(gene_stable_id),
                    y = man_SCNAs %>% 
                      select(gene_stable_id) %>% 
                      mutate(is_ = "is_manual"),
                    by = c("gene_stable_id" = "gene_stable_id")) %>% 
    filter(is.na(is_)) %>% 
    group_by(is_) %>% 
    summarise(is_only_GISTIC_default = n()) %>% 
    select(is_only_GISTIC_default)
  
  right <- right_join(x = def_SCNAs %>% 
                        select(gene_stable_id) %>% 
                        mutate(is_ = "is_default"),
                      y = man_SCNAs %>% 
                        select(gene_stable_id),
                      by = c("gene_stable_id" = "gene_stable_id")) %>% 
    filter(is.na(is_)) %>% 
    group_by(is_) %>% 
    summarise(is_only_GISTIC_manual = n()) %>% 
    select(is_only_GISTIC_manual)
  
  # Bind results
  temp_object <- bind_cols(left, 
                           intersection,
                           right) %>% 
    mutate(subtype = subtype) %>% 
    select(subtype, everything())
  
  working_tibble <- bind_rows(working_tibble,
                                temp_object)
  
}  

# Writing output files
comparison_path <- "results/GISTIC_SCNAs_discovery_congruence/gene_sets_comparison_3.7E-04.tsv"
vroom_write(working_tibble,
            delim = "\t",
            path = comparison_path)
