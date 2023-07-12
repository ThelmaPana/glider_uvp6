#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Process uvp6 data
# Date: 06/09/2022
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#



library(tidyverse)
library(arrow)
library(lubridate)
library(oce)
library(parallel)
library(ecotaxar)
library(scales)
library(castr)
library(multidplyr)
library(cmocean)
library(akima)
library(broom)
library(suncalc)  
library(pastecs)

source("lib/lib.R")

n_cores <- 10

coast <- read_csv("data/coast.csv", col_types = cols())


## Read data ----
#--------------------------------------------------------------------------#
bins <- read_parquet("data/10.bins.parquet")
parts <- read_parquet("data/10.parts.parquet")
obj <- read_parquet("data/10.obj.parquet")


## Compute plankton abundance per depth bin ----
#--------------------------------------------------------------------------#
obj_bins <- obj %>% 
  mutate(
    mission = str_split_fixed(sample, "_", n = 3)[,3] %>% tolower(),
    sample = str_split_fixed(sample, "_sea", n = 2)[,1],
    depth = roundp(depth, precision = 5, f = floor) + 2.5
  ) %>% 
  select(-c(mission, mission_part)) %>% 
  group_by(mission_det, sample, depth, taxon) %>% 
  summarise(abund = n()) %>% 
  ungroup()

obj_counts <- obj_bins %>% 
  group_by(taxon) %>% 
  summarise(n = sum(abund)) %>% 
  ungroup() %>% 
  arrange(desc(n))
sum(obj_counts$n)
# total number of objects

# Plot number of objects per taxon
obj_counts %>% 
  mutate(taxon = fct_inorder(taxon)) %>% 
  ggplot() +
  geom_col(aes(x =  taxon, y = n)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8000)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file = "figures/11.sup_dataset_comp.png", width = 164, height = 80, unit = "mm", dpi = 300, bg = "white")


## Compute concentration per bin ----
#--------------------------------------------------------------------------#
# Generate all combination of bins and taxa
taxa <- sort(unique((obj$taxon)))
obj_conc <- bins %>% 
  select(-mission) %>% 
  crossing(taxon = taxa) %>%
  # Join abundances per bin
  left_join(obj_bins) %>% 
  # Compute concentration
  mutate(
    conc = abund / watervolume,
    conc = ifelse(is.na(conc), 0, conc),
    conc = conc * 1000 # from ind per L to ind per m3
    ) %>% 
  select(mission_det, sample:depth, taxon, conc) %>% 
  # One column per taxon
  pivot_wider(names_from = taxon, values_from = conc)


## Compute particles concentration from abundances ----
#--------------------------------------------------------------------------#
parts_conc <- parts %>% 
  pivot_longer(class21:class33) %>% 
  mutate(value = value / watervolume) %>% 
  pivot_wider()


## Save ----
#--------------------------------------------------------------------------#
write_parquet(obj_conc, sink = "data/11.obj_conc.parquet")
write_parquet(parts_conc, sink = "data/11.parts_conc.parquet")


