################################################################################
# title: 'BioMart processing PIPELINE'
# author: "CASTRO DAU Santiago Manuel"
# date: '21 FEB 2109'
# last modified: '21 FEB 2020'
# description: 
# 'Takes the Biomart_EnsemblG94_GRCh38_p12 file and filters out non protein 
# codeing genes and genes not anottated in the nuclear 
# genome filtered out.'
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)

# Set working directory
setwd(dir = "/home/chewbacca/SCNA_Proyect")

# Reading the data of the BioMart dicctionary and filtering out non 
# coding genes and genens not annotated in nuclear DNA or 
# known genomic loci
path_Ensbl <- "data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12.txt"
Ensbl_IDs <- vroom::vroom(file = path_Ensbl,
                          col_names = TRUE,
                          col_types = c(col_double(), col_character())) %>% 
  janitor::clean_names(case = "snake") %>% 
  filter(gene_type == "protein_coding") %>% 
  mutate(chromosome_scaffold_name = str_replace(chromosome_scaffold_name,
                                                "X",
                                                "23")) %>%
  filter(chromosome_scaffold_name %in% c(1:23)) %>% 
  mutate(chromosome_scaffold_name = as.double(chromosome_scaffold_name)) %>%
  select(gene_stable_id,
         hgnc_symbol,
         everything()) %>% 
  distinct(gene_stable_id, .keep_all = TRUE)

# Writing new file
vroom_write(Ensbl_IDs,
            path = "data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt",
            delim = "\t")
