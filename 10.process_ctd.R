#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Download and process CTD and particle data
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
library(pastecs)
library(cmocean)
library(measurements)
library(gsw)
library(patchwork)

n_cores <- 24

coast <- read_csv("data/coast.csv", col_types = cols())


## Get samples and CTD data from ecotaxa and ecopart ----
#--------------------------------------------------------------------------#

# Connect to ecotaxa
db <- db_connect_ecotaxa()

# Get list of projects
projects <- tbl(db, "projects") %>% 
  select(projid, title) %>% 
  filter(str_detect(title, 'uvp6.*sea002')) %>% 
  arrange(projid) %>% 
  collect()

# Get list of samples
samples <- tbl(db, "samples") %>% 
  filter(projid %in% !!projects$projid) %>% 
  select(projid, sampleid, orig_id, lat=latitude, lon=longitude, datetime = t18) %>% 
  collect()

# Disconnect from ecotaxa
db_disconnect_ecotaxa(db)

# Connect to ecopart
db <- db_connect_ecopart()

# Get projects in ecopart
part_projects <- tbl(db, "part_projects") %>% 
  filter(projid %in% !!projects$projid) %>% 
  select(pprojid, ptitle, projid) %>% 
  collect()

# Get samples in ecopart
part_samples <- tbl(db, "part_samples") %>% 
  filter(pprojid %in% !!part_projects$pprojid) %>% 
  select(pprojid, psampleid, sampleid, ctd_desc) %>% 
  collect()

# Prepare CTD names
desc <- parse_mapping(iconv(part_samples$ctd_desc[1], from="latin1", to="utf8"))
# remove special characters in the variable names
clean_names <- names(desc) %>% 
  #str_replace("µ", "u") %>%
  make.names() %>%
  tolower()
  #str_replace_all("\\.{2,}", ".") %>% # remove double dots
  #str_replace_all("\\.$" , "")        # remove dots at end f name
# prepare the variable description for select()
desc <- str_c("extrames", desc)
names(desc) <- clean_names

# Get watervolumes from particles histo in ecopart
histopart <-  tbl(db, "part_histopart_det") %>% 
  filter(psampleid %in% !!part_samples$psampleid) %>% 
  select(psampleid, depth, watervolume, class21:class33) %>% 
  collect() %>% 
  unique()

# Get CTD data in relevant psamples
raw_ctd_data <- tbl(db, "part_ctd") %>% 
  filter(psampleid %in% !!part_samples$psampleid) %>% 
  select(psampleid, depth, datetime, temperature, chloro_fluo, conductivity, practical_salinity, oxygen_mass, oxygen_vol, cpar, fcdom_factory, desc) %>% 
  #glimpse()
  collect()

# Disconnect from ecopart
db_disconnect_ecopart(db)



## Process CTD data to generate clean 5 meter bins ----
#--------------------------------------------------------------------------#

not_all_na <- function(x) any(!is.na(x))
na_prop <- function(x) sum(is.na(x))/length(x)

# Clean raw data
ctd_data <- raw_ctd_data %>% 
  # Ignore variables which are all na
  select(where(not_all_na)) %>% 
  left_join(part_samples) %>% 
  left_join(part_projects) %>% 
  left_join(samples %>% select(-datetime)) %>% 
  arrange(datetime) %>% 
  select(title = ptitle, sample = orig_id, depth, datetime, temperature:flbbcd_cdom_scaled) %>% 
  # and we have to separate downcast from upcast
  # keep only downcasts: down and upcast are duplicated, d+u is already in d
  filter(str_detect(sample, "d")) %>% 
  # ignore values with no depth
  drop_na(depth) %>% 
  filter(depth > 0) %>% 
  # ignore data with weird date
  filter(year(datetime) == 2021) %>% 
  # round datetime to second
  mutate(datetime = round_date(datetime)) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  # reorganise columns
  select(title:datetime, lon = nav_longitude, lat = nav_latitude, everything())

unique(ctd_data$title)
rm(raw_ctd_data)


# Convert coordinates in deg/min/sec to decimal degree
# First we need to convert coordinates to character and add a space between degrees and minutes
ctd_data <- ctd_data %>% 
  mutate(
    lon = as.character(lon),
    lat = as.character(lat),
    lon = sub("(.{1})(.*)", "\\1 \\2", lon),
    lat = sub("(.{2})(.*)", "\\1 \\2", lat)
  )

# Convert coordinates to decimal degree
# NB not in mutate as `conv_unit` does not work with mutate
ctd_data$lon = conv_unit(ctd_data$lon, "deg_dec_min", "dec_deg") %>% as.numeric()
ctd_data$lat = conv_unit(ctd_data$lat, "deg_dec_min", "dec_deg") %>% as.numeric()


## Read and join SMRU data ----
#--------------------------------------------------------------------------#

# Empty tibble to store SMRU data
smru_data <- tibble()

# List SMRU files
smru_files <- list.files(path = "data/ctd_smru", full.names = TRUE)
# Loop over SMRU files and read them
for(file in smru_files){
  
  # Extract mission number
  mission_nb <- str_split_fixed(file, "[.]", n = 3)[,2]
  
  # Read file
  df <- read_tsv(file, col_types = cols_only( # Read only a few columns
    "PLD_REALTIMECLOCK" = "c", 
    "SMRU_TEMPERATURE" = "d", 
    "SMRU_CONDUCTIVITY" = "d", 
    "SMRU_SALINITY" = "d" 
  )) %>% 
    rename_with(tolower)
  
  # Clean data 
  df <- df %>% 
    mutate(
      mission = mission_nb, 
      datetime = round_date(dmy_hms(pld_realtimeclock)), 
      .before = pld_realtimeclock
      ) %>% # convert datetime to a datetime object and round to second
    arrange(datetime) %>% 
    select(-pld_realtimeclock) %>% 
    distinct() # keep distinct rows (duplicates induced by rounding datetime at second)
  
  # Append to previous files
  smru_data <- bind_rows(smru_data, df)
}

# Join CTD and SMRU data where SMRU data is required
smru_needed <- ctd_data %>% 
  group_by(title) %>% 
  summarise_all(na_prop) %>% 
  filter(temperature == 1) %>% 
  pull(title)

ctd_data <- ctd_data %>% 
  left_join(smru_data %>% select(-mission)) %>% 
  mutate(
    temperature = ifelse(title %in% smru_needed, smru_temperature, temperature),
    conductivity = ifelse(title %in% smru_needed, smru_conductivity, conductivity),
    practical_salinity = ifelse(title %in% smru_needed, smru_salinity, practical_salinity)
  ) %>% select(-contains("smru"))
rm(smru_data)


## Interpolate depth ----
#--------------------------------------------------------------------------#
# Depth is recorded at discrete intervals while time is continuous: interpolate depth from time
# Datetime at which depth is known
k_depth <- ctd_data %>% 
  group_by(title) %>% 
  select(title, sample, depth, datetime) %>% 
  mutate(
    lag_depth = lag(depth),
    change = depth != lag_depth,
    change = ifelse(is.na(change), TRUE, change)
  ) %>% 
  ungroup() %>% 
  filter(change) %>% 
  select(-c(lag_depth, change))

ctd_data <- ctd_data %>% 
  group_by(title) %>% 
  mutate(depth = castr::interpolate( x = k_depth$datetime, y = k_depth$depth, xout = datetime)) %>% 
  ungroup()


## Separate up and downcasts ----
#--------------------------------------------------------------------------#
# Separate down and upcasts using turning points
cluster <- new_cluster(n_cores)
cluster_library(cluster, "castr")
ctd_data <- ctd_data %>% 
  # compute very smoothed depth by sample to separate down and upcasts
  group_by(title, sample) %>% 
  partition(cluster) %>% 
  mutate(depth_sm = smooth(depth, k = 50, n = 5), .after = depth) %>% 
  collect() %>% 
  ungroup()

# Compute turning points on smoothed depth
tp <- turnpoints(ctd_data$depth_sm)

# Compute up and downcast from turning points
ctd_data_ad <- ctd_data %>% 
  #select(title:datetime) %>% 
  mutate(
    pits = tp$pits,
    peaks = tp$peaks
  ) %>% 
  #filter(title == "uvp6_sn000003lp_2021_sea002_m493" & sample == "Yo_208d") %>% 
  group_by(title, sample) %>% 
  mutate(
    rank = row_number(), # compute rank within each sample
    last_peak = max(which(peaks)), # find lask peak = end of dive
    last_pit = ifelse(length(which(pits)) > 0, max(which(pits)), 1), # find last pit = beginning of dive
    up = ifelse(rank >= last_peak, TRUE, FALSE),
    down = ifelse(rank >= last_pit, TRUE, FALSE),
    ud = ifelse(up & down, "a", ifelse(down, "d", NA))
  ) %>% 
  ungroup() %>% 
  drop_na(ud) %>% 
  mutate(
    yo_nb = str_extract_all(sample, "[0-9]+") %>% as.numeric(),
    sample = str_c("Yo_", yo_nb, ud),
    .after = sample
  ) %>% 
  select(-c(depth_sm, pits:down))

unique(ctd_data_ad$title)


# Plot to make sure that up and downcast were properly computed
ctd_data_ad %>% 
  filter(title == "uvp6_sn000003lp_2021_sea002_m489") %>% 
  mutate(batch = floor(yo_nb / 10)) %>% 
  #filter(between(yo_nb, 0, 20)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = ud)) +
  facet_wrap(~batch, scales = "free")


## Select columns ----
#--------------------------------------------------------------------------#
# Keep only columns of interest
ctd_data_ad <- ctd_data_ad %>% 
  mutate(mission = str_split_fixed(title, "_", n = 4)[,4] %>% tolower(), .after = title) %>% 
  select(
    mission:lat, 
    temp = temperature, sal = practical_salinity, chla = chloro_fluo, oxy = arod_ft_do, bb700 = flbbcd_bb_700_scaled, cdom = flbbcd_cdom_scaled
  )


## Remove aberrant data ----
#--------------------------------------------------------------------------#
# Remove aberrant data
ctd_data_ad <- ctd_data_ad %>% 
  mutate(
    sal = ifelse(sal < 37, NA, sal),
    chla = ifelse(chla < 0, NA, chla),
    bb700 = ifelse(bb700 < 0, NA, bb700),
    bb700 = ifelse(bb700 > 0.001, NA, bb700),
    cdom = ifelse(cdom < 0.7, NA, cdom),
    cdom = ifelse(cdom > 0.9, NA, cdom),
  )


## Compute density ----
#--------------------------------------------------------------------------#
ctd_data_ad <- ctd_data_ad %>% 
  mutate(
    pressure = gsw_p_from_z(z = -depth, latitude = lat),
    SA = gsw_SA_from_SP(SP = sal, p = pressure, longitude = lon, latitude = lat),
    CT = gsw_CT_from_t(SA = SA, t = temp, p = pressure),
    dens = gsw_sigma0(SA = SA, CT = CT),
    .after = sal
  ) %>% 
  select(-c(pressure, SA, CT))


## Convert oxygen from µmol/L to µmol/kg ----
#--------------------------------------------------------------------------#
ctd_data_ad <- ctd_data_ad %>% mutate(oxy = oxy *1000/(dens + 1000))


## Correct oxygen values from bottle data ----
#--------------------------------------------------------------------------#
# Need to remove 52 µmol/L
ctd_data_ad <- ctd_data_ad %>% mutate(oxy = oxy - 52)


## Bin data at 1 meter ----
#--------------------------------------------------------------------------#
# To reduce number of observations

cluster_library(cluster, "tidyverse")
ctd_1m <- ctd_data_ad %>% 
  mutate(depth = roundp(depth, precision = 1, f = round)) %>% 
  group_by(mission, sample, yo_nb, depth) %>% 
  partition(cluster) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  collect() %>% 
  ungroup() %>% 
  arrange(datetime) %>% 
  select(mission, sample, yo_nb, lon, lat, datetime, depth, everything())


## Despike CTD data ----
#--------------------------------------------------------------------------#
# Group only by mission to despike 
ctd_1m <- ctd_1m %>% 
  group_by(mission, sample, yo_nb) %>% 
  partition(cluster) %>% 
  mutate(
    temp  = despike(temp),
    sal   = despike(sal),
    dens  = despike(dens), 
    chla  = despike(chla),
    oxy   = despike(oxy),
    bb700 = despike(bb700),
    cdom  = despike(cdom),
  ) %>% 
  collect() %>% 
  ungroup()


## Bin data at 5 meter ----
#--------------------------------------------------------------------------#
# Bin data at 5 meters to match with Ecopart bins for particles
ctd_5m <- ctd_1m %>% 
  mutate(depth = roundp(depth, precision = 5, f = floor) + 2.5) %>% 
  group_by(mission, sample, yo_nb, depth) %>% 
  partition(cluster) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  collect() %>% 
  ungroup() %>% 
  arrange(datetime) %>% 
  select(mission, sample, yo_nb, lon, lat, datetime, depth, everything())


## Despike CTD data at 5 meters ----
#--------------------------------------------------------------------------#
ctd_5m <- ctd_5m %>% 
  group_by(mission) %>% 
  partition(cluster) %>% 
  mutate(
    temp  = despike(temp),
    sal   = despike(sal),
    dens  = despike(dens), 
    chla  = despike(chla),
    oxy   = despike(oxy),
    bb700 = despike(bb700),
    cdom  = despike(cdom),
  ) %>% 
  collect() %>% 
  ungroup()
summary(ctd_5m)

# Plot first mission
mission_plot <- "sea002_m486"
ctd_5m %>% 
  filter(mission == mission_plot) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = temp)) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "thermal")

ctd_5m %>% 
  filter(mission == mission_plot) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = sal)) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "haline")

ctd_5m %>% 
  filter(mission == mission_plot) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = dens)) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "dense")

ctd_5m %>% 
  filter(mission == mission_plot) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = chla)) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "algae")

ctd_5m %>% 
  filter(mission == mission_plot) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = oxy)) +
  facet_wrap(~mission, scales = "free") +
  scale_color_distiller(palette = "Blues", direction = 1)

ctd_5m %>% 
  filter(mission == mission_plot) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = bb700)) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "turbid")

ctd_5m %>% 
  #filter(mission == mission_plot) %>% 
  #filter(between(yo_nb, 10, 30)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = cdom)) +
  geom_hline(yintercept = -50) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "matter")

ctd_5m %>% 
  filter(mission == mission_plot) %>% 
  #filter(between(yo_nb, 10, 30)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = cdom)) +
  geom_hline(yintercept = -50) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "matter")


## Interpolate missing CTD data ----
#--------------------------------------------------------------------------#
ctd_na_prop <- ctd_5m %>% 
  group_by(mission, sample, yo_nb) %>% 
  summarise(across(temp:cdom, na_prop)) %>% 
  ungroup()

# Max proportion of missing value to interpolate a profile
prop_thr <- 0.25

# Initiate empty tibble
ctd_int <- tibble()

# Variables to interpolate
vars <- c("temp", "sal", "dens", "chla", "oxy", "bb700", "cdom")

for (my_var in vars){
  # Interpolation
  int <- ctd_5m %>% 
    select(mission:depth, all_of(my_var)) %>% 
    left_join(ctd_na_prop %>% select(mission, sample, yo_nb, na_prop = all_of(my_var))) %>% 
    pivot_longer(all_of(my_var), names_to = "variable") %>%
    group_by(mission, sample, yo_nb) %>% 
    mutate(
      value = ifelse(na_prop > 0 & na_prop <= prop_thr, interpolate(x = depth, y = value, xout = depth), value)
    ) %>% 
    ungroup() %>% 
    select(-na_prop) %>% 
    pivot_wider(names_from = variable, values_from = value)
  
  # Join with other variables
  if (nrow(ctd_int) == 0){
    ctd_int <- int
  } else {
    ctd_int <- ctd_int %>% left_join(int)
  }
  
}


## Smooth CTD data ----
#--------------------------------------------------------------------------#
ctd_int <- ctd_int %>% 
  group_by(mission) %>% 
  partition(cluster) %>% 
  mutate(
    temp  = smooth(temp,  k = 2, n = 2),
    sal   = smooth(sal,   k = 2, n = 2),
    dens  = smooth(dens,  k = 2, n = 2),
    chla  = smooth(chla,  k = 2, n = 2),
    oxy   = smooth(oxy,   k = 2, n = 2),
    bb700 = smooth(bb700, k = 2, n = 2),
    cdom  = smooth(cdom,  k = 2, n = 2)
  ) %>% 
  collect() %>% 
  ungroup() %>% 
  arrange(datetime)


mission_plot <- "sea002_m486"
ctd_int %>% 
  filter(mission == mission_plot) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, color = temp)) +
  facet_wrap(~mission, scales = "free") +
  scale_color_cmocean(name = "thermal")


## Compute distance from Villefranche ----
#--------------------------------------------------------------------------#
# Coordinates of Nice Cap
lon_vlfr <- 7.304707
lat_vlfr <- 43.685812


# Compute distance along path and from villefranche
ctd_int <-  ctd_int %>% 
 mutate(dist = geodDist( # then distance is computed from villefranche
   longitude1 = lon, 
   latitude1  = lat, 
   longitude2 = lon_vlfr, 
   latitude2  = lat_vlfr, 
 ), .after = datetime)

# Plot first mission with distance
ctd_int %>% 
 filter(mission == "sea002_m486") %>% 
 #filter(depth < 10) %>% 
 filter(yo_nb < 50) %>% 
 ggplot() +
 geom_point(aes(x = dist, y = -depth, color = yo_nb))

# Plot map
ctd_int %>% 
 filter(mission == "sea002_m486") %>% 
 filter(depth < 10) %>% 
 ggplot() +
 geom_point(aes(x = lon, y = lat, color = yo_nb)) +
 coord_quickmap()


## Add sampled watervolumes to CTD bins ----
#--------------------------------------------------------------------------#
# Get water volumes with samples and projects
watervolumes <- histopart %>% 
  select(psampleid, depth, watervolume) %>% 
  left_join(part_samples) %>% 
  left_join(samples) %>% 
  left_join(projects) %>% 
  select(title, sample = orig_id, depth, watervolume) %>% 
  mutate(
    mission = str_split_fixed(title, "_", n = 4)[,4] %>% tolower()
  )  

# Join with CTD data
sampled_bins <- watervolumes %>% 
  left_join(ctd_int %>% select(mission:depth)) %>% 
  select(mission, sample, yo_nb, lon, lat, datetime, dist, depth, watervolume) %>% 
  drop_na(datetime) %>% 
  arrange(datetime)

  
## Particles ----
#--------------------------------------------------------------------------#
# Join particles on bins
parts <- histopart %>% 
  left_join(part_samples) %>% 
  left_join(part_projects) %>% 
  left_join(samples %>% select(-datetime)) %>% 
  select(title = ptitle, sample = orig_id, depth, watervolume, class21:class33) %>% 
  rename(mission = title) %>% 
  mutate(mission = str_split_fixed(mission, "_", n = 4)[,4] %>% tolower())

parts <- sampled_bins %>% left_join(parts)


## Keep only backs in ctd, bins and parts tables ----
#--------------------------------------------------------------------------#
# First we need to read objects data
obj <- read_csv("data/09.uvp6_grouped_objects.csv")
samp_rem <- c("Yo_61a_sea002_m493", "Yo_61d_sea002_m493") # One yo deeper than 300m in m493 has to be removed
obj <- obj %>% filter(!sample_id %in% samp_rem)

# Clean mission and sample name
obj <- obj %>% 
  rename(sample = sample_id) %>% 
  mutate(
    mission = str_split_fixed(sample, "_", n = 3)[,3] %>% tolower(),
    sample = str_split_fixed(sample, "_sea", n = 2)[,1],
    yo_nb = str_extract(sample, "[0-9]+") %>% as.numeric(),
    .after = sample
  ) 

# Compute back samples
back_samples <- obj %>% 
  arrange(datetime) %>% 
  select(mission, sample, yo_nb) %>% 
  unique() %>% 
  group_by(mission) %>% 
  mutate(
    yo_diff = yo_nb - lag(yo_nb),
    rank = row_number(),
    mission_part = ifelse(rank < which.max(yo_diff), 1, 2),
    mission_part = as.factor(mission_part),
    mission_det = str_c(mission, "part", mission_part, sep = "_"),
    .after = mission
  ) %>% 
  ungroup() %>% 
  select(-(yo_diff:rank)) %>% 
  mutate(back = TRUE)

ctd <- ctd_int %>% left_join(back_samples) %>% filter(back) %>% select(-c(mission_part, back)) %>% select(mission, mission_det, everything())
bins <- sampled_bins %>% left_join(back_samples) %>% filter(back) %>% select(-c(mission_part, back)) %>% select(mission, mission_det, everything())
parts <- parts %>% left_join(back_samples) %>% filter(back) %>% select(-c(mission_part, back)) %>% select(mission, mission_det, everything())
obj <- obj %>% left_join(back_samples) %>% select(-c(mission_part, back)) %>% select(mission, mission_det, everything())


## Save data ----
#--------------------------------------------------------------------------#
write_parquet(ctd_int, sink = "data/10.ctd_all.parquet")
write_parquet(ctd, sink = "data/10.ctd.parquet")
write_parquet(bins, sink = "data/10.bins.parquet")
write_parquet(parts, sink = "data/10.parts.parquet")
write_parquet(obj, sink = "data/10.obj.parquet")



