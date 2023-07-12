#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Prepare learning and test set for RF training
# Date: 06/01/2023
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(arrow)
library(googlesheets4)
library(data.tree)
library(ecotaxar)

gs4_auth(use_oob = TRUE)
2

# Link to spreadsheet
ss <- "https://docs.google.com/spreadsheets/d/1uzxfHgg56bpfK2u12-DC-K0hUuZxjsu8rnYIHD8MJik/edit?usp=sharing"
# Taxonomy level to be used
level_use <- "level2"

# Proportion of training data, remaining is used for validation
train_prop <- 0.8
# Number of detritus for training
n_detritus <- 100000


## Read data ----
#--------------------------------------------------------------------------#
# Read data and drop useless columns
df_learn <- read_parquet("data/dataset/learn/01.uvp6_objects.parquet") %>% select(-c(txo.id, x, y, xm, ym, angle, morphocluster_label))
df_test <- read_parquet("data/dataset/test/01.uvp6_objects.parquet") %>% select(-c(txo.id, x, y, xm, ym, angle, morphocluster_label))
# All data is the combination of learning and test set
df_new <- bind_rows(df_learn, df_test)

# Keep only validated objects for learning data
df_learn <- df_learn %>% filter(classif_qual == "V") 


## Build the taxonomy tree ----
#--------------------------------------------------------------------------#
tc <- df_learn %>% 
  count(taxon, lineage) %>% 
  # convert it into a tree
  rename(pathString=lineage) %>%
  arrange(pathString) %>%
  as.Node()

print(tc, "taxon","n", limit = 50)
# Convert to dataframe
tcd <- ToDataFrameTree(tc, "taxon", "n")%>% 
  as_tibble() %>% 
  rename(level0=taxon, nb_level0=n)


## Write tree into GSS ----
#--------------------------------------------------------------------------#
# Start by erasing previous data (3 first columns) in spreadsheet
range_flood(ss, sheet = "learn", range = "learn!A:C", reformat = FALSE)
# Write new tree
range_write(ss, sheet = "learn", data = tcd) 
# Open it in browser tab to make sure everything is ok
gs4_browse(ss)


## Read tree count from Google Spread Sheet (GSS) and create table for taxonomy match ----
#--------------------------------------------------------------------------#
tcd <- read_sheet(ss, sheet = "learn")

# Get match between level0 (EcoTaxa taxonomy), level1 (taxonomy to use for classif) and level2 (ecological group)
taxo_match <- tcd %>% 
  select(level0, all_of(level_use), plankton = plankton2) %>% 
  mutate(plankton = as.logical(plankton)) %>% 
  drop_na(level0) %>% 
  drop_na(all_of(level_use))

names(taxo_match) <- c("level0", "new_taxon", "plankton")


## Match taxonomy between EcoTaxa export and taxonomy to use ----
#--------------------------------------------------------------------------#
df_learn <- df_learn %>% 
  rename(level0=taxon) %>% 
  left_join(taxo_match, by = "level0") %>% 
  select(-level0) %>% 
  rename(taxon=new_taxon) %>% 
  drop_na(taxon) %>% 
  select(objid:longitude, plankton, taxon, area:skeleton_area)

df_test <- df_test %>% 
  rename(level0=taxon) %>% 
  left_join(taxo_match, by = "level0") %>% 
  select(-level0) %>% 
  rename(taxon=new_taxon) %>% 
  drop_na(taxon) %>% 
  select(objid:longitude, plankton, taxon, area:skeleton_area)

df_new <- df_new %>% 
  rename(level0=taxon) %>% 
  left_join(taxo_match, by = "level0") %>% 
  select(-level0) %>% 
  rename(taxon=new_taxon) %>% 
  drop_na(taxon) %>% 
  select(objid:longitude, area:skeleton_area)


## Reduce number of detritus in learning set ----
#--------------------------------------------------------------------------#
to_drop <- df_learn %>% 
  filter(taxon == "detritus") %>% 
  slice_sample(n = nrow(.) - n_detritus) %>% 
  pull(objid)

df_learn <- df_learn %>% 
  filter(!(objid %in% to_drop))



## Split training and validation in the learning set ----
#--------------------------------------------------------------------------#
df_learn <- df_learn %>% 
  group_by(taxon) %>% 
  # shuffle rows
  sample_frac(1) %>%
  mutate(
    # percent rank (which is random)
    r=(1:n()) / n(),
    # assign in set based in this
    set=case_when(
      r >= train_prop ~ "val",
      TRUE ~ "train"
    )
  ) %>% 
  select(-r) %>%
  ungroup() %>% 
  relocate(set, .after = plankton)

df_test <- df_test %>% mutate(set = "test", .after = plankton)


## Get deep features from ecotaxa ----
#--------------------------------------------------------------------------#
# Connect to database
db <- db_connect_ecotaxa()

# Extract deep features
d_feat <- tbl(db, "objects") |> 
  filter(projid %in% c(7544L, 7545L)) |> 
  select(objid) |> 
  left_join(
    tbl(db, "obj_cnn_features"),
    by=c("objid"="objcnnid")
  ) %>% 
  collect()

# Disconnect
db_disconnect_ecotaxa(db)


# Add deep features to datasets
df_learn <- df_learn %>% left_join(d_feat)
df_test <- df_test %>% left_join(d_feat)
df_new <- df_new %>% left_join(d_feat)


# Save data
write_csv(df_learn, file = "data/dataset/learn/04.learn_dataset.csv")
write_csv(df_test, file = "data/dataset/test/04.test_dataset.csv")
write_csv(df_new, file = "data/dataset/test/04.new_data.csv")

