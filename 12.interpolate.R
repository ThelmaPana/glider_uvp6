#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Interpolate and plot data
# Date: 14/09/2022  
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
library(patchwork)
library(ggforce)
library(sp)
library(ggrepel)

source("lib/lib.R")

coast <- read_csv("data/coast.csv", col_types = cols())



## Read data ----
#--------------------------------------------------------------------------#
obj_conc <- read_parquet("data/11.obj_conc.parquet")
ctd <- read_parquet("data/10.ctd.parquet")
bins <- read_parquet("data/10.bins.parquet")
parts_conc <- read_parquet("data/11.parts_conc.parquet")

# List all missions
all_mission <- unique(ctd$mission_det)

# Beginning date of each transect
dates <- obj_conc %>% 
  group_by(mission_det) %>% 
  summarise(date = date(min(datetime))) %>% 
  ungroup()

# List taxa
taxa <- obj_conc %>% 
  select(Annelida:Salpida) %>% 
  colnames()
# Ignore unwanted taxa
det_taxa <- c("artefact", "badfocus", "bubble", "detritus" , "feces", "fiber", "filament", "othertocheck", "reflection", "tentacle", "turbid")
temp_taxa <- c("t003", "t004", "t005", "t008", "t011")

taxa <- taxa %>% base::setdiff(det_taxa) %>% base::setdiff(temp_taxa)


## Interpolate CTD data ----
#--------------------------------------------------------------------------#
# Group data to ease interpolation
ctd_g <- ctd %>% 
  select(mission_det, dist, depth, temp:cdom) %>% 
  mutate(dist = round(dist)) %>% 
  group_by(mission_det, depth, dist) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup() 

ctd_g %>% 
  ggplot() +
  geom_point(aes(x = dist, y = -depth, color = sal)) +
  scale_color_cmocean(name = "haline") +
  facet_wrap(~mission_det, ncol = 4)


## Coarse interpolation
step_x <- 1 # in km
step_y <- 5 # in m

# Variables to interpolates
vars <- c("temp", "sal", "dens", "chla", "oxy", "bb700", "cdom") 


ctd_int <- tibble()

for (id in all_mission){
  #print(id)
  transect_int <- tibble()
  
  # Keep only CTD data for m493 part 1
  df_t <- ctd_g %>% filter(mission_det == id)
  
  # Generate output grid for transect
  xo <- seq(floor(min(df_t$dist)),  ceiling(max(df_t$dist)), by = step_x)
  yo <- seq(0, ceiling(max(df_t$depth)), by = step_y)
  
  # Loop over variables to interpolate
  for (variable in vars){
    # Compute variable interpolation
    if (all(is.na(df_t %>% pull(all_of(variable))))){# in everything is NA, then create an output table of NA
      var_int <- crossing(dist = xo, depth = yo) %>% 
        mutate(mission_det = id) %>% 
        mutate(name = variable, value = NA) %>% 
        pivot_wider(names_from = name, values_from = value)
      
    } else { # else proceed to interpolation
      var_int <- coarse_interp_1v(df_t, variable, xo=xo, yo=yo)  
    }
    
    
    
    # Join to table with other variables
    if (length(transect_int) == 0) { # If first interpolated variable on this transect
      transect_int <- var_int # Replace transect table by newly computed interpolation
    } else { # Else perform a left join with previously interpolated variables
      transect_int <- left_join(transect_int, var_int, by = c("mission_det", "dist", "depth"))
    } 
    
  }
  ctd_int <- bind_rows(ctd_int, transect_int)
}


ctd_int %>% 
  left_join(dates) %>% 
  ggplot() +
  geom_tile(aes(x = dist, y = -depth, fill = temp)) +
  facet_wrap(~mission_det, ncol = 4) +
  scale_fill_cmocean(na.value = NA) +
  theme_minimal()

write_parquet(ctd_int, sink = "data/12.ctd_int.parquet")


## Fine interpolation
step_x <- 0.2 # fine interpolation step in X axis: 0.2 km or 200 m
step_y <- 0.5 # fine interpolation step in Y axis: 0.5 m
theta <- 0.5 # bandwidth or scale parameter


# Initiate empty tibble for this transect
ctd_int_fine <- tibble()

for (id in all_mission){
  #print(mission_ctd)
  
  transect_int <- tibble()
  
  # Keep only CTD data for m493 part 1
  df_t <- ctd_int %>% filter(mission_det == id)
  
  # Generate output grid for transect
  xo <- seq(floor(min(df_t$dist)),  ceiling(max(df_t$dist)), by = step_x)
  yo <- seq(0, ceiling(max(df_t$depth)), by = step_y)
  
  # Loop over variables to interpolate
  for (variable in vars){
    # Compute variable interpolation
    if (all(is.na(df_t %>% pull(all_of(variable))))){# in everything is NA, then create an output table of NA
      var_int <- crossing(dist = xo, depth = yo) %>% 
        mutate(mission_det = id) %>% 
        mutate(name = variable, value = NA) %>% 
        pivot_wider(names_from = name, values_from = value)
      
    } else { # else proceed to interpolation
      var_int <- fine_interp_1v(df_t, variable, plankton = c(), step_x_fine = step_x, step_y_fine = step_y, theta = theta)
    }
    
    
    
    # Join to table with other variables
    if (length(transect_int) == 0) { # If first interpolated variable on this transect
      transect_int <- var_int # Replace transect table by newly computed interpolation
    } else { # Else perform a left join with previously interpolated variables
      transect_int <- left_join(transect_int, var_int, by = c("mission_det", "dist", "depth"))
    } 
    
  }
  ctd_int_fine <- bind_rows(ctd_int_fine, transect_int)
}

ctd_int_fine <- ctd_int_fine %>% left_join(dates)

write_parquet(ctd_int_fine, sink = "data/12.ctd_int_fine.parquet")



## Interpolate plankton data ----
#--------------------------------------------------------------------------#
## Coarse interpolation
# Variables to interpolates
vars <- taxa

step_x <- 5 # in km
step_y <- 30 # in m

conc_summary <- obj_conc %>% 
  select(mission_det, dist, depth, all_of(taxa)) %>% 
  mutate(
    dist = roundp(dist, precision = step_x, f = floor) + step_x/2,
    depth = roundp(depth, precision = step_y, f = floor) + step_y/2
  ) %>% 
  group_by(mission_det, dist, depth) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup()


conc_int <- tibble()

for (id in all_mission){
  #print(id)
  
  transect_int <- tibble()
  
  # Keep only CTD data for m493 part 1
  df_t <- conc_summary %>% filter(mission_det == id)
  
  # Generate output grid for transect
  xo <- seq(floor(min(df_t$dist)),  ceiling(max(df_t$dist)), by = step_x)
  yo <- seq(0, ceiling(max(df_t$depth)), by = step_y)
  
  # Loop over variables to interpolate
  for (variable in taxa){
    # Compute variable interpolation
    var_int <- coarse_interp_1v(df_t, variable, xo=xo, yo=yo)
    
    # Join to table with other variables
    if (length(transect_int) == 0) { # If first interpolated variable on this transect
      transect_int <- var_int # Replace transect table by newly computed interpolation
    } else { # Else perform a left join with previously interpolated variables
      transect_int <- left_join(transect_int, var_int, by = c("mission_det", "dist", "depth"))
    } 
    
  }
  conc_int <- bind_rows(conc_int, transect_int)
}

plankton_int <- conc_int %>% left_join(dates)


write_parquet(plankton_int, sink = "data/12.plankton_int.parquet")

plankton_int %>% 
  ggplot() +
  geom_tile(aes(x = dist, y = -depth, fill = Copepoda)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "white", size = 0.5, data = ctd_int_fine) +
  facet_wrap(~date, ncol = 4) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(-300, 0)) +
  #facet_wrap(~mission_det) +
  #theme(text = element_text(size = 18), legend.key.height = unit(1.5, "cm")) +
  theme(legend.key.height = unit(1.5, "cm")) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Interpolated Copepoda", fill = "#.m-3")




## Interpolate particle data ----
#--------------------------------------------------------------------------#
parts_summary <- parts_conc %>% 
  filter(mission_det %in% all_mission) %>% 
  mutate(dist = round(dist)) %>% 
  group_by(mission_det, dist, depth) %>% 
  summarise(across(class21:class33, mean, na.rm = TRUE)) %>% 
  ungroup()


## Coarse interpolation
step_x <- 1 # in km
step_y <- 5 # in m

# Variables to interpolates
vars <- parts_summary %>% select(contains("class")) %>% colnames()

part_int <- tibble()

for (id in all_mission){
  #print(id)  
  transect_int <- tibble()
  
  df_t <- parts_summary %>% filter(mission_det == id)
  
  # Generate output grid for transect
  xo <- seq(floor(min(df_t$dist)),  ceiling(max(df_t$dist)), by = step_x)
  yo <- seq(0, ceiling(max(df_t$depth)), by = step_y)
  
  # Loop over variables to interpolate
  for (variable in vars){
    # Compute variable interpolation
    var_int <- coarse_interp_1v(df_t, variable, xo=xo, yo=yo)
    
    # Join to table with other variables
    if (length(transect_int) == 0) { # If first interpolated variable on this transect
      transect_int <- var_int # Replace transect table by newly computed interpolation
    } else { # Else perform a left join with previously interpolated variables
      transect_int <- left_join(transect_int, var_int, by = c("mission_det", "dist", "depth"))
    } 
    
  }
  part_int <- bind_rows(part_int, transect_int)
}

part_int <- part_int %>% left_join(dates)

write_parquet(part_int, sink = "data/12.part_int.parquet")

## Fine interpolation
step_x <- 0.2 # fine interpolation step in X axis: 0.2 km or 200 m
step_y <- 0.5 # fine interpolation step in Y axis: 0.5 m
theta <- 0.3 # bandwidth or scale parameter


# Initiate empty tibble for this transect
part_int_fine <- tibble()

for (id in all_mission){
  #print(id)  
  transect_int <- tibble()
  
  # Keep only CTD data for m493 part 1
  df_t <- part_int %>% filter(mission_det == id)
  
  # Generate output grid for transect
  xo <- seq(floor(min(df_t$dist)),  ceiling(max(df_t$dist)), by = step_x)
  yo <- seq(0, ceiling(max(df_t$depth)), by = step_y)
  
  # Loop over variables to interpolate
  for (variable in vars){
    # Compute variable interpolation
    var_int <- fine_interp_1v(df_t, variable, plankton = c(), step_x_fine = step_x, step_y_fine = step_y, theta = theta)
    
    # Join to table with other variables
    if (length(transect_int) == 0) { # If first interpolated variable on this transect
      transect_int <- var_int # Replace transect table by newly computed interpolation
    } else { # Else perform a left join with previously interpolated variables
      transect_int <- left_join(transect_int, var_int, by = c("mission_det", "dist", "depth"))
    } 
    
  }
  part_int_fine <- bind_rows(part_int_fine, transect_int)
}

part_int_fine <- part_int_fine %>% left_join(dates)

write_parquet(part_int_fine, sink = "data/12.part_int_fine.parquet")






