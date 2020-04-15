################################################################################
# title: 'Gscore processing PIPELINE'
# author: "CASTRO DAU Santiago Manuel"
# date: '17 NOV 2109'
# last modified: '19 FEB 2020'
# description: 
# 'This pipeline takes the list of human protein coding genes and 
# the Gscore data for each subtype, AS INPUT and it joins them, producing as  
# output a list of protein coding genes with an associated Gscore. 
# The pipeline produces one ouput for each cancer subtype:
# -The list of genes with associated Gscores (aplifications and deletions) and
#  genes that did not fall into any Gscore interval
################################################################################
# Calling libreries
library(tidyverse)
library(stringr)
library(vroom)
library(janitor)

# Set working directory
setwd(dir = "/home/chewbacca/SCNA_Proyect")

# Reading the data of the BioMart dicctionary to join with the Gscore data and
# filtering out non coding genes and genens not annotated in nuclear DNA or 
# known genomic loci
path_Ensbl <- "data/BioMart_data/Biomart_EnsemblG94_GRCh38_p12_only_protein_coding_genes.txt"
Ensbl_IDs <- vroom::vroom(file = path_Ensbl,
                          col_names = TRUE,
                          col_types = c(col_double(), col_character())) 

# Crating a list where we specify the subtypes. The list will be used to iterate
# the pipeline for each of the subtypes.
subtypes_list <- c("basal", "her2", "luma", "lumb")
subtype <- "basal"

# 1st Loop start
for(subtype in subtypes_list) {
  # Set file path
  path_Gscore <- paste("data/GISTIC_raw_data/",
                       subtype, 
                       "/scores.gistic",
                       sep = "")
  
  # Read file and clean names
  Gscores <- vroom::vroom(file = path_Gscore,
                          col_names = TRUE,
                          col_types = c(col_double(), col_character())) %>% 
    janitor::clean_names(case = "snake") 
  
  # Generate new columns names to add to the genes dictionary
  column_names <- c(colnames(Gscores)[c(1, 5:8)], "g_score_sd") 
  
  # Creating object outside of the loop to which we will asign the products
  # of the in-loop objects to row_bind them and write them at the end
  genes_Gscore <- as_tibble()
  
  # Generating a list to do the join two times, one for amplifications and one
  # for deletions
  types <- c("Amp", "Del")
  type_ <- "Del"
  
  # 2nd loop start
  for (type_ in types) {
    
    # Asign Gscores to an object side the loop and fileting by SCNA type
    Gscores_in_loop <- Gscores %>%
      filter(type == type_)
    
    # Creating an in-loop object witht the protein coding genes plus the
    # new columns to be filld in the next loop
    genes_in_loop <- Ensbl_IDs
    genes_in_loop[, column_names] <- NA
      
    i <- 1
      # 3rd loop start: for every gene in the dicctionary  
      for (i in 1:dim(genes_in_loop)[1]) {
        
        # Obtaining filtering data form the gene in question
        chr_gene <- genes_in_loop$chromosome_scaffold_name[i]
        start_gene <- genes_in_loop$gene_start_bp[i]
        end_gene <- genes_in_loop$gene_end_bp[i]
     
        # Generating new row witht the median information of all the intervals
        # that the gene feel into
        new_row <- Gscores_in_loop %>% 
          filter(chromosome == chr_gene &
                   (start_gene <= start | start_gene <= end) &
                   (end_gene >= end | end_gene >= start)) %>% 
          summarise(g_score_median = median(g_score),
                    g_score_sd = sd(g_score),
                    log10_q_value_median = median(log10_q_value),
                    average_amplitude_median = median(average_amplitude),
                    frequency_median = median(frequency),
                    type = type_)
        
        # Adding each of the components of the new_row to the empty columns of
        # the in-loop genes_in_loop data
        genes_in_loop$g_score[i] <- new_row$g_score_median
        genes_in_loop$g_score_sd[i] <- new_row$g_score_sd
        genes_in_loop$log10_q_value[i] <- new_row$log10_q_value_median
        genes_in_loop$average_amplitude[i] <- new_row$average_amplitude_median
        genes_in_loop$frequency[i] <- new_row$frequency_median
        genes_in_loop$type[i] <- new_row$type
      }
      
    # Ordering the columns in the desired order
    genes_in_loop <- genes_in_loop %>% 
      select("gene_stable_id",
             "hgnc_symbol",
             "chromosome_scaffold_name",
             "gene_start_bp",
             "gene_end_bp",
             "gene_percent_gc_content",             
             "gene_type",
             "g_score",
             "g_score_sd", 
             "average_amplitude",
             "frequency",
             "log10_q_value",
             "type")
    
    # Binding amp and del data into the same object wich will be written 
    # before starting the iteration with the next subtype
    genes_Gscore <- bind_rows(genes_in_loop,
                              genes_Gscore)
    
  }
  # Write data in the corresponding subtype directorty with the corresponding
  # type_ name
  vroom_write(genes_Gscore,
              path = paste("data/GISTIC_processed_data/",
                           subtype, 
                           "/genes_w_Gscore.tsv",
                           sep = ""),
              delim = "\t")
}
  