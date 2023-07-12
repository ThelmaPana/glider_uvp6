#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Set taxonomy to higher levels for plankton data
# Date: 06/09/2022
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(arrow)
library(googlesheets4)
library(data.tree)

gs4_auth(use_oob = TRUE)
2

# Link to spreadsheet
ss <- "https://docs.google.com/spreadsheets/d/1uzxfHgg56bpfK2u12-DC-K0hUuZxjsu8rnYIHD8MJik/edit?usp=sharing"
# Taxonomy level to be used
level_use <- "level2"


## Read data ----
#--------------------------------------------------------------------------#
# Read objects 
df <- read_parquet("data/08.uvp6_objects_science.parquet")

# Read score thresholds
thr <- read_csv("data/06.taxa_threshold.csv") %>% select(taxon, score_threshold)

# List of non planktonic objects
taxa <- sort(unique(thr$taxon))
non_plankton <- c("artefact", "detritus", "other_living")


## Regroup objects ----
#--------------------------------------------------------------------------#
tcd <- read_sheet(ss, sheet = "learn")

# Get match between level0 (EcoTaxa taxonomy), level1 (taxonomy to use for classif) and level2 (ecological group)
taxo_match <- tcd %>% 
  select(level0, all_of(level_use), plankton = plankton2) %>% 
  mutate(plankton = as.logical(plankton)) %>% 
  drop_na(level0) %>% 
  drop_na(all_of(level_use))

names(taxo_match) <- c("level0", "new_taxon", "plankton")


df <- df %>% 
  rename(level0=taxon) %>% 
  left_join(taxo_match, by = "level0") %>% 
  select(-level0) %>% 
  rename(taxon=new_taxon) %>% 
  drop_na(taxon) %>% 
  select(objid:longitude, classif_qual, rf_score, taxon, plankton, area)


df %>% count(taxon) %>% filter(!taxon %in% c("detritus", "artefact")) %>% arrange(desc(n)) %>% pull(n) %>% sum()

# Number of predicted objects
df %>% 
  filter(classif_qual == "P") %>% 
  count(taxon)


## Threshold predictions ----
#--------------------------------------------------------------------------#
# Keep objects whether validated, whether with prediction score above threshold

df <- df %>% 
  left_join(thr) %>% 
  mutate(keep = (classif_qual == "V") | (score_threshold > rf_score)) %>% 
  filter(keep) %>% 
  select(-keep)



## Proportion of living in 1-2 mm ----
#--------------------------------------------------------------------------#
df_living <- df %>% 
  mutate(
    esd = 2 * sqrt(area/pi), # compute esd from area
    esd = esd * 0.073 # 1 px = 0.073 mm
    ) %>% 
  filter(between(esd, 1, 2)) %>% 
  select(sample_id, depth, datetime, latitude, longitude, taxon, esd, plankton)

sum(df_living$plankton) / nrow(df_living)


## Dataset composition ----
#--------------------------------------------------------------------------#
df_counts <- df %>% 
  select(objid, taxon, plankton) %>% 
  count(plankton, taxon) %>% 
  arrange(desc(plankton))

df_counts %>% pull(n) %>% sum()
df_counts %>% filter(plankton) %>%  pull(n) %>% sum()

df_counts %>% 
  filter(plankton) %>% 
  select(-plankton) %>% 
  arrange(desc(n)) %>% 
  mutate(taxon = fct_inorder(taxon)) %>% 
  ggplot() +
  geom_col(aes(x = taxon, y = n)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Taxonomic group", y = "Total number of images") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## Discard non relevant objects ----
#--------------------------------------------------------------------------#
df <- df %>% filter(!(taxon %in% c("detritus", "other_living", "artefact")))


# Save data
write_csv(df, file = "data/09.uvp6_grouped_objects.csv")


## Size of houses ----
#--------------------------------------------------------------------------#
sort(unique(df$taxon))
df %>% 
  filter(taxon == "Appendicularia") %>% 
  select(objid:taxon, area) %>% 
  mutate(
    esd = 2*sqrt(area/pi),
    esd = esd *0.073
  ) %>% 
  summary()


## Size of salps ----
#--------------------------------------------------------------------------#
sort(unique(df$taxon))
df %>% 
  filter(taxon == "Salpida") %>% 
  select(objid:taxon, area) %>% 
  mutate(
    esd = 2*sqrt(area/pi),
    esd = esd *0.073
  ) %>% 
  summary()
# median esd = 2.7 mm