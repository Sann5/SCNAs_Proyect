################################################################################
# title: 'sif files processing PIPELINE'
# author: "CASTRO DAU Santiago Manuel"
# date: '24 FEB 2109'
# last modified: '24 FEB 2020'
# description: 
# 'Filter and normalizeing the top 133,170 interactions for each subtype'
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)

# Set working directory
setwd(dir = "/home/chewbacca/SCNA_Proyect")

# Crating a list where we specify the subtypes. The list will be used to iterate
# the pipeline for each of the subtypes.
subtypes_list <- c("basal", "her2", "luma", "lumb", "healthy")
subtype <- "healthy"

# Loop start
for(subtype in subtypes_list) {
  
  # Set file path
  sif_path <- paste("data/network_data/",
                    subtype,
                    "/",
                    subtype,
                    ".sif.sorted",
                    sep = "")
  
  # Read file
  working_sif <- vroom(file = sif_path,
                       delim = "\t",
                       col_types = c(col_character(), col_double()),
                       col_names = c("source",
                                     "MI",
                                     "target"))
  
  # Geting the max MI value and defining the fun_max function
  max_MI <- max(working_sif$MI)
  fun_max <- function(x) (x/max_MI)
  
  # Filtering top 133k interactions and normalizeing
  working_sif <- working_sif %>% 
    top_n(133170, MI) %>%
    mutate_at(c("MI"), fun_max)
  
  # Setting path and writting data
  new_file_path <- paste("data/network_data/",
                         subtype,
                         "/",
                         subtype,
                         "_norm_133k_interactions.sif",
                         sep = "")
  
  vroom_write(working_sif,
              path = new_file_path,
              delim = "\t")
  
}

