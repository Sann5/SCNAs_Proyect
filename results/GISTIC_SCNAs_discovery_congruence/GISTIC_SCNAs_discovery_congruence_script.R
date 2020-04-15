################################################################################
# title: 'GISTIC SCNAs discovery congruence script'
# author: "CASTRO DAU Santiago Manuel"
# date: '24 MAR 2020'
# last modified: '24 MAR 2020'
# description: 
# 'Extracts the SCNAs reported in the del_genes.conf_99.txt and 
# amp_genes.conf_99.txt files for each subtype, extracts the q_value and 
# residual q_values and compares the items in this list to the list manually
# obtained from other GISTIC output files.'
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)

# Set working directory
setwd(dir = "/Users/Santiago/SCNA_Proyect")

# Loading proteing coding genes
biomart_path <- paste("data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt")

biomart <- vroom(file = biomart_path,
                 delim = "\t",
                 col_types = c(col_character(), col_double()))

# Creating working object outside the loop
working_tibble <- tibble()
working_tibble_2 <- tibble()

# Crating a list where we specify the subtypes.
subtypes_list <- c("basal", "her2", "luma", "lumb")
subtype <- "basal"

# For each subtype load the genes.conf_99.txt files, extract q values and get
# SCNAs list
for (subtype in subtypes_list) {
  
  # Defining path to genes.conf_99.txt files
  amps_path <- paste("data/GISTIC_raw_data/",
                subtype,
                "/amp_genes.conf_99.txt",
                sep = "")
  
  dels_path <- paste("data/GISTIC_raw_data/",
                subtype,
                "/del_genes.conf_99.txt",
                sep = "")
  
  # Reading files
  amps <- read_tsv(file = amps_path, col_names = FALSE)
  dels <- read_tsv(file = dels_path, col_names = FALSE)
  
  # Celaning data
  amps <- amps %>% 
    t() %>%                                                # Transpose data
    as_tibble(.name_repair = "minimal") %>%                # Transform to tibble
    row_to_names(row_number = 1) %>%                       # First row to col names
    clean_names(case = "snake") %>%                        # Clean col names
    pivot_longer(cols = -c("cytoband",                     # Pivot longer genes and removing na values
                           "q_value", 
                           "residual_q_value", 
                           "wide_peak_boundaries"),
      names_to = "trash_col",
      names_prefix = "na",
      values_to = "genes_in_wide_peak",
      values_drop_na = TRUE) %>% 
    select(-trash_col) %>%                                 # Get rid of trash col
    separate(wide_peak_boundaries, c("chromosome",         # Separate wide_peak_boundaries col
                                     "temp_col"), ":") %>% 
    separate(temp_col, c("wide_peak_start_bp", 
                         "wide_peak_end_bp"), "-") %>% 
    mutate(chromosome = str_remove(chromosome, "chr")) %>% # Remove chr from chromosome col
    mutate(chromosome = str_replace(chromosome,            # Change X for 23
                                    "X",
                                    "23")) %>% 
    mutate_at(c("q_value",                                 # Cols as double
                "residual_q_value", 
                "chromosome", 
                "wide_peak_start_bp", 
                "wide_peak_end_bp"), ~as.double(.)) %>% 
    select(genes_in_wide_peak,                             # Reorder cols
           q_value, 
           residual_q_value, 
           chromosome, 
           cytoband, 
           wide_peak_start_bp, 
           wide_peak_end_bp) %>%
    mutate(type = "amp")                                   # Identifie SCNA as amp or del
  
  dels <- dels %>% 
    t() %>%                                                # Transpose data
    as_tibble(.name_repair = "minimal") %>%                # Transform to tibble
    row_to_names(row_number = 1) %>%                       # First row to col names
    clean_names(case = "snake") %>%                        # Clean col names
    pivot_longer(cols = -c("cytoband",                     # Pivot longer genes and removing na values
                           "q_value", 
                           "residual_q_value", 
                           "wide_peak_boundaries"),
                 names_to = "trash_col",
                 names_prefix = "na",
                 values_to = "genes_in_wide_peak",
                 values_drop_na = TRUE) %>% 
    select(-trash_col) %>%                                 # Get rid of trash col
    separate(wide_peak_boundaries, c("chromosome",         # Separate wide_peak_boundaries col
                                     "temp_col"), ":") %>% 
    separate(temp_col, c("wide_peak_start_bp", 
                         "wide_peak_end_bp"), "-") %>% 
    mutate(chromosome = str_remove(chromosome, "chr")) %>% # Remove chr from chromosome col
    mutate(chromosome = str_replace(chromosome,            # Change X for 23
                                    "X",
                                    "23")) %>% 
    mutate_at(c("q_value",                                 # Cols as double
                "residual_q_value", 
                "chromosome", 
                "wide_peak_start_bp", 
                "wide_peak_end_bp"), ~as.double(.)) %>% 
    select(genes_in_wide_peak,                             # Reorder cols
           q_value, 
           residual_q_value, 
           chromosome, 
           cytoband, 
           wide_peak_start_bp, 
           wide_peak_end_bp) %>%
    mutate(type = "del")                                   # Identifie SCNA as amp or del
  
  # Binding amps and dels
  SCNAs <- bind_rows(amps, dels)
  
  # Filtering out of SCNAs non protein coding genes and 
  # adding gene_stable_id, keeping only unique elements
  SCNAs <- left_join(x = SCNAs,
                     y = biomart %>% select(gene_stable_id, hgnc_symbol),
                     by = c("genes_in_wide_peak" = "hgnc_symbol")) %>% 
    filter(!is.na(gene_stable_id)) %>% 
    select(gene_stable_id, everything()) %>% 
    distinct(gene_stable_id, type, .keep_all = TRUE)
  
  # Creating summary
  working_tibble <- bind_rows(working_tibble,
                              SCNAs %>% 
                                summarise(q_value_max = max(q_value), 
                                          residual_q_value_max = max(residual_q_value)) %>% 
                                mutate(subtype = subtype))
  
  # Load manually curated SCNAs to compare both gene set composition
  # Set file path
  manual_SCNAs_path <- paste("data/GISTIC_processed_data/",
                      subtype,
                      "/genes_w_Gscore.tsv",
                      sep = "")
  
  # Read file and clean names
  manual_SCNAs <- vroom(file = manual_SCNAs_path,
                 delim = "\t",
                 col_types = c(col_character(), col_double())) %>% 
    filter(!is.na(g_score)) %>% 
    distinct(gene_stable_id, type, .keep_all = TRUE) %>% 
    mutate(q_value = 1/(10**log10_q_value)) %>% 
    filter(q_value <= SCNAs %>% pull(q_value) %>% max(.)) 
  
  # Compare
  intersection <- inner_join(x = SCNAs %>% select(gene_stable_id),
                             y = manual_SCNAs %>% 
                               select(gene_stable_id) %>% 
                               mutate(is_ = "is_"),
                             by = c("gene_stable_id" = "gene_stable_id")) %>% 
    group_by(is_) %>% 
    summarise(is_both = n()) %>% 
    select(is_both)
  
  left <- left_join(x = SCNAs %>% select(gene_stable_id),
                    y = manual_SCNAs %>% 
                      select(gene_stable_id) %>% 
                      mutate(is_ = "is_"),
                    by = c("gene_stable_id" = "gene_stable_id")) %>% 
    filter(is.na(is_)) %>% 
    group_by(is_) %>% 
    summarise(is_only_GISTIC_default = n()) %>% 
    select(is_only_GISTIC_default)
  
  right <- right_join(x = SCNAs %>% 
                        select(gene_stable_id) %>% 
                        mutate(is_ = "is_"),
                      y = manual_SCNAs %>% 
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
  
  working_tibble_2 <- bind_rows(working_tibble_2,
                                temp_object)
  
  # Writtig default SNACs to file in data/GISTIC_processed_data/subtype
  default_SCNAs_path <- paste("data/GISTIC_processed_data/",
                              subtype,
                              "/default_SCNAs_w_q_values.tsv",
                              sep = "")
  
  vroom_write(SCNAs,
              delim = "\t",
              path = default_SCNAs_path)
  
} 

# Writing output files

q_values_path <- "results/GISTIC_SCNAs_discovery_congruence/default_q_values.tsv"
vroom_write(working_tibble,
            delim = "\t",
            path = q_values_path)

comparison_path <- "results/GISTIC_SCNAs_discovery_congruence/gene_sets_comparison_2.5E-01.tsv"
vroom_write(working_tibble_2,
            delim = "\t",
            path = comparison_path)
