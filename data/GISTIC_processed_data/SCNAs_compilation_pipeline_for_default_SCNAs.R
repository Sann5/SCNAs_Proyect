################################################################################
# title: 'All subtypes SCNAs compilation script for default SCNAs'
# author: "CASTRO DAU Santiago Manuel"
# date: '24 MAR 2020'
# last modified: '24 MAR 2020'
# description: 
# 'Makes a list of all unique SCNAs detected through the GISTIC2.0 algorithm
# in each of the subtypes and saves list of unique elements into healthy 
# directory for its use in enrichment_analysis.R'
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)

# Set working directory
setwd(dir = "/Users/Santiago/SCNA_Proyect")

# Crating a list where we specify the subtypes and the network sizes. 
# The list will be used to iterate the pipeline for each of the subtypes and
# network sizes.
subtypes_list <- c("basal", "her2", "luma", "lumb")

# Creating working object outside the loop
working_tibble <- tibble()

# 2nd Loop start
for(subtype in subtypes_list) {
  
  # Loading SCNAs
  path_SCNAs <- paste("data/GISTIC_processed_data/",
                      subtype,
                      "/default_SCNAs_w_q_values_filtered.tsv",
                      sep = "")
  
  SCNAs <- vroom(file = path_SCNAs,
                 delim = "\t",
                 col_types = c(col_character(), col_double())) %>% 
    select(gene_stable_id)
  
  # Adding SCNAs to outside loop object
  working_tibble <- bind_rows(working_tibble,
                              SCNAs)
}

# Creating and saving list of unique elements into healthy directory for its 
# use in enrichment_analysis.R

working_tibble <- working_tibble %>% 
  distinct(gene_stable_id)

vroom_write(working_tibble,
            path = "data/GISTIC_processed_data/healthy/default_SCNAs_w_q_values_filtered.tsv",
            delim = "\t")