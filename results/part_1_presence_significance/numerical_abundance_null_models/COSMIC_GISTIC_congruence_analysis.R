################################################################################
# title: 'COSMIC, GISTIC congruence analysis'
# author: "CASTRO DAU Santiago Manuel"
# date: '03 MAR 2020'
# last modified: '09 MAR 2020'
# description: 
# 'This script analysises the congruence between the SCNA's anotated in the
# COSMIC database and the significanly aplifies/deleted genes detected by the
# GISTIC2.0 algorithm'
################################################################################
# Calling libreries
library(tidyverse)
library(vroom)
library(janitor)
library(scales)

# Set working directory
setwd(dir = "/Users/Santiago/SCNA_Proyect")

# Crating a list where we specify the subtypes. The list will be used to iterate
# the pipeline for each of the subtypes.
subtypes_list <- c("basal", "her2", "luma", "lumb")

# Crating working tibble outside loop
working_tibble <- tibble()

# Loop start
for(subtype in subtypes_list) {

  # Set file path
  path_SCNAs <- paste("data/GISTIC_processed_data/",
                      subtype,
                      "/genes_w_Gscore_filtered.tsv",
                      sep = "")
  
  # Read file and clean names
  SCNAs <- vroom(file = path_SCNAs,
                   delim = "\t",
                   col_types = c(col_character(), col_double()))
  
  # Counting number of aplifications/deletions
  amps <- SCNAs %>% 
    group_by(type) %>%
    summarise(count = n()) %>%
    filter(type == "Amp") %>%
    pull(count)
  
  dels <- SCNAs %>% 
    group_by(type) %>%
    summarise(count = n()) %>%
    filter(type == "Del") %>%
    pull(count)
  
  # Jioning most significant GISTIC SCNA's in the COSMIC database
    # Seting file path
    path_COSMIC <- "data/COSMIC_data/COSMIC_porcessed_SCNA.tsv"
    
    # Reading data
    COSMIC <- vroom(file = path_COSMIC,
                    delim = "\t",
                    col_types = c(col_character(), col_double()))
    
    # Filtering COSMIC SCNA's
    COSMIC_top <- bind_rows(COSMIC %>%
                              filter(gain > loss) %>%
                              filter(gain != 1) %>%
                              top_n(amps, ratio_gain) %>%
                              arrange(-ratio_gain),
                            COSMIC %>%
                              filter(gain < loss) %>%
                              filter(loss != 1) %>%
                              top_n(dels, ratio_loss) %>%
                              arrange(-ratio_loss))
    
    intersection <- inner_join(SCNAs %>% select(-is_),
                               COSMIC_top %>% select(gene_stable_id, is_),
                               by = c("gene_stable_id" = "gene_stable_id"))
    
    # Writing data into outside object
    working_tibble <- bind_rows(working_tibble,
                                intersection %>%
                                  group_by(is_) %>%
                                  summarise(count = n()) %>%
                                  mutate(subtype = subtype),
                                SCNAs %>%
                                  group_by(is_) %>%
                                  summarise(count = n()) %>%
                                  mutate(subtype = subtype))
}  

# Procesing data and making graph
GISTIC_COSMIC_fig <- working_tibble %>% 
  pivot_wider(names_from = is_,
              values_from = count) %>% 
  mutate(ratio = is_COSMIC/is_GISTIC) %>% 
  pivot_longer(-c(subtype, ratio), 
               names_to = "is_",
               values_to = "count") %>% 
  mutate(ratio = replace(ratio, is_=="is_GISTIC", NA)) %>% 
  mutate(is_ = as_factor(is_)) %>%
  mutate(is_ = fct_relevel(is_, "is_GISTIC", "is_COSMIC")) %>%
  ggplot(mapping = aes(fill = is_,
                       y = count,
                       x = as.factor(subtype),
                       label = scales::percent(ratio))) + 
  geom_bar(position = "dodge",
           stat = "identity") +
  geom_text(position = position_dodge(width = .9),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 3) +
  labs(title = "SCNA's found both throught the GISTIC algorithm and in COSMIC",
       y = "Number of genes",
       x = "Subtypes",
       fill = "") + 
  scale_fill_discrete(labels = c("In GISTIC", "In GISTIC & in COSMIC")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal()

# Saving figure
saveRDS(GISTIC_COSMIC_fig,
        file = "results/part_1_presence_significance/numerical_abundance_null_models/GISTIC_COSMIC_fig.rds")
