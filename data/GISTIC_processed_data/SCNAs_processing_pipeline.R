################################################################################
# title: 'SCNA's processing PIPELINE'
# author: "CASTRO DAU Santiago Manuel"
# date: '15 NOV 2109'
# last modified: '18 FEB 2020'
# description: 
# 'This pipeline summerises the all_data_by_genes.txt file in the 
# GISTIC_raw_data directory for each subtype, generating three files:
# 1) proccesed_SCNAs_summerised.tsv: List of genes with statistical summary of 
#    all the samples SCNA messures for each specifi gene.
# 2) processed_SCNAs_by_gene_by_sample.tsv: Large format data of every SCNA 
#    messure for every sample for every protein coding gene.
# 3) SCNAs_processing_summary.tsv: summary of the SCNAs_processing_pipeline.R 
#    pipeline for each subtype.'
################################################################################
# Calling libreries
library(dplyr)
library(tidyr)
library(stringr)
library(vroom)
library(janitor)

# Set working directory
setwd(dir = "/home/chewbacca/SCNA_Proyect")

# Reading the data of the BioMart dicctionary to join with the SCNA data
path_Ensbl <- "data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12.txt"
Ensbl_IDs <- vroom::vroom(file = path_Ensbl,
                          col_names = TRUE,
                          col_types = c(col_double(), col_character())) %>% 
  janitor::clean_names(case = "snake")

# Crating a list where we specify the subtypes. The list will be used to iterate
# the pipeline for each of the subtypes.
subtypes_list <- c("basal", "her2", "luma", "lumb")

# Loop start
for(subtype in subtypes_list) {
  # Set file path
  path_SCNAs <- paste("data/GISTIC_raw_data/",
                      subtype,
                      "/all_data_by_genes.txt",
                      sep = "")
  
  # Read file and clean names
  SCNAs <- vroom::vroom(file = path_SCNAs,
                       col_names = TRUE,
                       col_types = c(col_double(), col_character())) %>% 
    janitor::clean_names(case = "snake")
  
  # Count number of genes in the original SCNAs file
  original_genes <- dim(SCNAs)[1]
  num_samples <- dim(SCNAs)[2] - 3
  orignial_duplicates <- SCNAs %>% 
    select(gene_symbol) %>%
    duplicated() %>% 
    which() %>% 
    length()
  
  # Replace gene symbols with Ensembl IDs 
  SCNAs <- SCNAs %>% dplyr::left_join(x = SCNAs,
                                    y = Ensbl_IDs,
                                    by = c("gene_symbol" = "hgnc_symbol"))
  
  # Count number of genes in SCNAs that were not found in Ensbl_IDs, 
  # (multriplications and contractions)
  not_found <- SCNAs %>% 
    select(gene_stable_id) %>%
    is.na %>% 
    which() %>% 
    length()
  multriplications <- SCNAs %>% 
    select(gene_symbol) %>%
    duplicated() %>% 
    which() %>% 
    length()
  contractions <- SCNAs %>% 
    select(gene_stable_id) %>%
    duplicated() %>% 
    which() %>% 
    length()
  
  # Removal of not found values, duplicates and all genes that are not 
  # protein coding or not in nuclear DNA. Then we count remaining genes.
  chromosome_list <- c(1:22, "X")
  SCNAs <- SCNAs %>% 
    distinct(gene_stable_id, .keep_all = TRUE) %>%
    drop_na(gene_stable_id) %>% 
    filter(gene_type == "protein_coding") %>% 
    filter(chromosome_scaffold_name %in% chromosome_list)
  
  processed_genes <- dim(SCNAs)[1]
  
  # Change X for 23 in chromosome_scafold_name and change to number format
  SCNAs <- SCNAs %>% 
    mutate(chromosome_scaffold_name = str_replace(chromosome_scaffold_name,
                                                  "X",
                                                  "23")) %>%
    mutate(chromosome_scaffold_name = as.double(chromosome_scaffold_name))
  
  # Pivot data and discard gene_id column
  SCNAs <- SCNAs %>% 
    select(-gene_id) %>%
    pivot_longer(names_to = "sample", 
                 values_to = "copy_number_change",
                 cols = -c(gene_stable_id,
                           gene_symbol,
                           chromosome_scaffold_name,
                           cytoband,
                           gene_start_bp,
                           gene_end_bp,
                           gene_percent_gc_content,
                           gene_type)) %>% 
    select(gene_stable_id,
           gene_symbol,
           chromosome_scaffold_name,
           everything()) %>% 
    mutate(type = ifelse(copy_number_change > 0,
                         "amp", 
                         "del"))
  
  # We write the SCNAs tibble as an intermediate data product in the pipeline 
  # because it will be useful for other analysis
  vroom_write(SCNAs,
              path = paste("data/GISTIC_processed_data/",
                           subtype,
                           "/processed_SCNAs_by_gene_by_sample.tsv",
                           sep = ""),
              delim = "\t")
  
  # We procede to summerise all of the samples SCNA infromation by gene
  SCNAs <- SCNAs %>%
    group_by(gene_stable_id,
             gene_symbol,
             chromosome_scaffold_name,
             cytoband,
             gene_start_bp,
             gene_end_bp,
             gene_percent_gc_content,
             gene_type) %>% 
    summarise(mean = mean(copy_number_change),
              abs_mean = abs(mean(copy_number_change)),
              sd = sd(copy_number_change), 
              median = median(copy_number_change),
              max = max(copy_number_change),
              min = min(copy_number_change),
              num_amp = sum(type == "amp"),
              num_del = sum(type == "del"))
  # This section wil be left commented because the Amp/Del tagg will be attached
  # in the subsequent pipeline: Gscore_pipeline
    # mutate(ratio = ifelse(num_amp > num_del,
    #                        num_amp/num_samples,
    #                        num_del/num_samples),
    #        type = ifelse(num_amp > num_del,
    #                      "Amp",
    #                      "Del"))
  
  # Write SCNAs file
  vroom_write(SCNAs,
              path = paste("data/GISTIC_processed_data/",
                           subtype, 
                           "/proccesed_SCNAs_summerised.tsv",
                           sep = ""),
              delim = "\t")
  
  #Write summary of the process
  col_1 <- c(subtype,
             num_samples,
             original_genes,
             orignial_duplicates,
             contractions,
             multriplications,
             not_found,
             processed_genes)
  col_2 <- c("Subtype",
             "Samples",
             "Genes before process",
             "Duplicates before process",
             "Contractions",
             "Multriplications",
             "Gene symbols with no ID",
             "Genes after process")
  summary <- bind_cols(as.data.frame(col_2),
                       as.data.frame(col_1)) %>% 
    as_tibble()
  
  path_SCNAs_summary <- paste("data/GISTIC_processed_data/",
                             subtype, 
                             "/SCNAs_processing_summary.tsv",
                             sep = "")
  vroom_write(summary,
              path = path_SCNAs_summary,
              delim = "\t")
}

