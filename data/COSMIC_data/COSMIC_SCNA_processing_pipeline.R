################################################################################
# title: 'COSMIC SCNA processing PIPELINE'
# author: "CASTRO DAU Santiago Manuel"
# date: '20 FEB 2109'
# last modified: '20 FEB 2020'
# description: 
# 'This script cleans the data from tissue_cnv Tue Nov 26 23_27_19 2019.csv and
# it joins it with the Biomart dicctionary'
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)

# Set working directory
setwd(dir = "/home/chewbacca/SCNA_Proyect")

# Reading COSMIC data and creating columns with gain/loss ratios and celaning
# gene name column
ASNC_cosmic <- vroom(file = "data/COSMIC_data/tissue_cnv Tue Nov 26 23_27_19 2019.csv",
                     delim = ",",
                     col_names = c("gene", 
                                   "gain",
                                   "loss",
                                   "tested",
                                   "over",
                                   "under",
                                   "tested2"),
                     col_types = c(col_character(),col_double())) %>%
  mutate(ratio_gain = gain/tested) %>%
  mutate(ratio_loss = loss/tested)  %>% 
  mutate(gene = paste(gene, "_", sep = "")) %>%
  mutate(gene = str_match(gene, pattern = "[:alnum:]+")) %>%
  distinct(gene, .keep_all = TRUE) %>%
  select(-c(over, under, tested2))

# This script can be usful in the future to filter in significant
# amplifications and deletions
# ASNC_gain <- ASNC_cosmic %>%
#   filter(gain != 1) %>%
#   filter(ratio_loss < ratio_gain) %>%
#   top_n(650, ratio_gain)
# ASNC_loss <- ASNC_cosmic %>%
#   filter(loss != 1) %>%
#   filter(ratio_loss > ratio_gain) %>%
#   top_n(650, ratio_loss)
# ASNC_top <- bind_rows(ASNC_gain, ASNC_loss)

# Reading BioMart dictionary to match gene_id with gene name
Ensbl_IDs <- vroom(file = "data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt",
                   col_names = TRUE,
                   col_types = c(col_double(), col_character())) %>% 
  janitor::clean_names(case = "snake")

# Jouining data sets
ASNC_to_write <- inner_join(x = ASNC_cosmic, 
                            y = Ensbl_IDs,
                            by = c("gene" = "hgnc_symbol")) %>% 
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  mutate(is_ = "is_COSMIC")

# Writting 
vroom_write(ASNC_to_write, 
            path = "data/COSMIC_data/COSMIC_porcessed_SCNA.tsv",
            delim = "\t")
