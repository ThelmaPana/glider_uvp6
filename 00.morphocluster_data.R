#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Extract data that needs to be imported to morphocluster
# Date: 22/09/2021
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

library(tidyverse)

# Initiate empty df to store samples for each project
morpho_samples <- tibble()
morpho_samples_back <- tibble()

# List csv files to read
files <- list.files(path = "data/morphocluster", full.names = TRUE)

# Loop over csv files and process them
for (file in files){

  ## Before renaming samples
  #df <- read_csv(file, col_types = cols()) %>% # read file
  #  filter(part != "transit") %>%  # discard transit part
  #  summarise(samples = paste(samples, collapse = ",")) %>% # paste samples from out and back
  #  mutate(project = basename(file) %>% str_remove(".csv") %>% str_remove("morpho_"), .before = samples) # generate project name

  ## After renaming samples
  # Rename samples
  df <- read_csv(file, col_types = cols()) %>% # read file
    filter(part != "transit") %>%  # discard transit part
    separate_rows(samples, convert = TRUE, sep = ",") %>% 
    mutate(samples = paste(samples, str_split_fixed(basename(file) %>% str_remove(".csv"), "_", n = 5)[5], sep = "_")) # rename samples
  
  # Extract back transects for scientific analysis
  df_back <- df %>% 
    filter(part == "back") %>% 
    summarise(samples = paste(samples, collapse = ",")) %>%  # paste samples from out and back
    mutate(project = basename(file) %>% str_remove(".csv") %>% str_remove("morpho_"), .before = samples) # generate project name
    
  # Out and back transects for morphocluster subset
  df <- df %>% 
    summarise(samples = paste(samples, collapse = ",")) %>%  # paste samples from out and back
    mutate(project = basename(file) %>% str_remove(".csv") %>% str_remove("morpho_"), .before = samples) # generate project name
  
  # bind to other projects
  morpho_samples <- bind_rows(morpho_samples, df)
  morpho_samples_back <- bind_rows(morpho_samples_back, df_back)
}

# Save as a csv file
write_csv(morpho_samples, "data/00.morphocluster_samples_subset.csv")

# List back samples for science
morpho_samples_back <- morpho_samples_back %>% 
  summarise(samples = paste(samples, collapse = ","))
write_lines(morpho_samples_back$samples, "data/00.morphocluster_samples_science.txt")


## Subsample to assess prediction quality ----
#--------------------------------------------------------------------------#
prop <- 0.1 # proportion of samples to validate

samples <- tibble(samples = str_split(morpho_samples_back$samples, ",")[[1]]) 

test_set <- samples %>%
  filter((row_number() - 1) %% (1/prop) == 0) 

train_set <- anti_join(samples, test_set)

test_set <- test_set %>% summarise(samples = paste(samples, collapse = ","))
train_set <- train_set %>% summarise(samples = paste(samples, collapse = ","))

write_lines(test_set$samples, "data/00.morphocluster_samples_test.txt")
write_lines(train_set$samples, "data/00.morphocluster_samples_train.txt")
