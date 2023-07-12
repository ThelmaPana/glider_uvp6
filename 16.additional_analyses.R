#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Perform supplementary analyses: wind, surface chla…
# Date: 07/12/2022
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(lubridate)
library(arrow)
library(tidync)
library(cmocean)
source("lib/lib.R")


## Read data ----
#--------------------------------------------------------------------------#
ctd <- read_parquet("data/10.ctd.parquet")
obj_conc <- read_parquet("data/11.obj_conc.parquet")
coast <- read_csv("data/coast.csv", col_types = cols())


## Wind speed ----
#--------------------------------------------------------------------------#
# Read weather data
weather <- tibble()
files <- list.files("data/weather", full.names = TRUE)

for (file in files){
  d <- read_delim(file, 
             delim = ";", escape_double = FALSE, na = c("", "NA", "mq"), col_types = cols(date = col_datetime(format = "%Y%m%d%H%M%S")), 
             trim_ws = TRUE) %>% 
    filter(numer_sta == "07690") %>% 
    select(date, avg_wind = ff)
  
  weather <- bind_rows(weather, d)

}

weather <- weather %>% mutate(month = month(date))

# Get mission dates
missions <- ctd %>% 
  group_by(mission_det) %>% 
  summarise(
    beg = min(datetime),
    end = max(datetime)
  ) %>% 
  ungroup() %>% 
  mutate(month = month(beg))

# Plot wind speed by month and highlight missions
weather %>% 
  ggplot() +
  coord_cartesian(clip = "off", ylim = c(0, 13)) +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = Inf), fill = "#4b7bac", alpha = 0.8, data = missions) +
  geom_path(aes(x = date, y = avg_wind)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~month, scales = "free_x") +
  labs(x = "Date", y = expression(paste("Average wind 10 min (", m^{-1}, ")"))) +
  theme_minimal() +
  theme(
    strip.text.x = element_blank(), text = element_text(size = 10), panel.spacing = unit(1, "lines"), 
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    plot.margin = unit(c(20, 0,0,0), "points")
    )
ggsave(file = "figures/15.sup_wind.png", width = 164, height = 80, unit = "mm", dpi = 300, bg = "white")


## Surface chla ----
#--------------------------------------------------------------------------#
# Read surface chla data
chla_surf <- tidync("data/satellite/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_1674220699768.nc") %>% 
  hyper_tibble() %>% 
  mutate(date = as_date(time, origin = "1900-01-01")) %>% 
  select(-time)

# Plot one chla map to see the data
chla_surf %>% 
  mutate(month = month(date)) %>% 
  filter(month == 1) %>% 
  ggplot() +
  geom_polygon(aes(x = lon, y = lat), fill = "gray", data = coast) +
  geom_raster(aes(x = lon, y = lat, fill = CHL)) +
  scale_fill_viridis_c() +
  coord_quickmap() +
  facet_wrap(~date, ncol = 7)

# Plot a time series of surface chla values on the transect
chla_surf %>% 
  ggplot(aes(x = date, y = CHL), alpha = 0.01) +
  geom_point(size = 0.5) +
  geom_smooth(se = FALSE) +
  geom_vline(xintercept = as_date("2021-02-24"), color = "red", linewidth = 0.5, alpha = 0.5) +
  scale_x_date(date_breaks = "2 week") +
  scale_y_continuous(trans = "log1p", expand = c(0, 0.1)) +
  labs(x = "Date", y = bquote('Surface Chl a (mg m'~m^-3~')')) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 10))
ggsave(file = "figures/15.sup_surface_chla.png", width = 164, height = 80, unit = "mm", dpi = 300, bg = "white")



## Lenses (SCV) ----
#--------------------------------------------------------------------------#
ctd <- read_parquet("data/10.ctd_all.parquet") %>% 
  filter(mission == "sea002_m491") %>% 
  filter(datetime <= ymd_hms("2021-04-24 00:00:00"))

## Interpolate data
# Coarse interpolation
step_x <- 3600 # in s
step_y <- 5 # in m
# Convert datetime to numeric for interpolation
ctd_g <- ctd %>% mutate(datetime = as.numeric(datetime))

ctd_g %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = -depth, colour = oxy))


# Variables to interpolates
vars <- c("temp", "sal", "dens", "chla", "oxy", "bb700", "cdom") 

# Initiate empty tibble
transect_int <- tibble()

# Generate output grid for transect
xo <- seq(min(ctd_g$datetime), max(ctd_g$datetime), by = step_x)
yo <- seq(0, ceiling(max(ctd_g$depth)), by = step_y)

# Loop over variables to interpolate
for (variable in vars){
  print(variable)
  # Compute variable interpolation
  if (all(is.na(ctd_g %>% pull(all_of(variable))))){# in everything is NA, then create an output table of NA
    var_int <- crossing(datetime = xo, depth = yo) %>% 
      mutate(mission = id) %>% 
      mutate(name = variable, value = NA) %>% 
      pivot_wider(names_from = name, values_from = value)
    
  } else { # else proceed to interpolation
    var_int <- coarse_interp_1v_dt(ctd_g, variable, xo=xo, yo=yo)  
  }
  
  # Join to table with other variables
  if (length(transect_int) == 0) { # If first interpolated variable on this transect
    transect_int <- var_int # Replace transect table by newly computed interpolation
  } else { # Else perform a left join with previously interpolated variables
    transect_int <- left_join(transect_int, var_int, by = c("mission", "datetime", "depth"))
  } 
}

# Datetime to datetime format
transect_int <- transect_int %>% mutate(datetime = as_datetime(datetime))

summary(transect_int)

transect_int %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = oxy)) +
  scale_fill_viridis_c()



# Fine interpolation
step_x <- 300 # fine interpolation step in X axis: 300 seconds
step_y <- 0.5 # fine interpolation step in Y axis: 0.5 m
theta <- 5 # bandwidth or scale parameter


transect_int_fine <- tibble()

# Generate output grid for transect
xo <- seq(min(as.numeric(transect_int$datetime)), max(as.numeric(transect_int$datetime)), by = step_x)
yo <- seq(0, ceiling(max(transect_int$depth)), by = step_y)



# Loop over variables to interpolate
for (variable in vars){
  print(variable)
  # Compute variable interpolation
  if (all(is.na(transect_int %>% pull(all_of(variable))))){# in everything is NA, then create an output table of NA
    var_int <- crossing(datetime = xo, depth = yo) %>% 
      mutate(mission = id) %>% 
      mutate(name = variable, value = NA) %>% 
      pivot_wider(names_from = name, values_from = value)
    
  } else { # else proceed to interpolation
    var_int <- fine_interp_1v_dt(transect_int %>% mutate(datetime = as.numeric(datetime)), variable, plankton = c(), xo=xo, yo=yo, theta = theta)
  }
  
  
  
  # Join to table with other variables
  if (length(transect_int_fine) == 0) { # If first interpolated variable on this transect
    transect_int_fine <- var_int # Replace transect table by newly computed interpolation
  } else { # Else perform a left join with previously interpolated variables
    transect_int_fine <- left_join(transect_int_fine, var_int, by = c("mission", "datetime", "depth"))
  } 
  
}

## Limit maximum interpolation depth
# Get back transects
backs <- obj_conc %>% 
  filter(str_detect(mission_det, "sea002_m491")) %>% 
  group_by(mission_det) %>% 
  summarise(
    beg = min(datetime),
    end = max(datetime)
  )
# extract beginning and end of each back transect
begs <- backs %>% pull(beg)
ends <- backs %>% pull(end)

transect_int_fine <- transect_int_fine %>% 
  pivot_longer(temp:cdom, names_to = "variable", values_to = "value") %>% 
  mutate(datetime = as_datetime(datetime)) %>% 
  mutate(
    value = ifelse(datetime > begs[1] & datetime < ends[1] & depth > 300, NA, value),
    value = ifelse(datetime > begs[2] & datetime < ends[2] & depth > 300, NA, value)
  ) %>% 
  pivot_wider(names_from = "variable", values_from = "value")


## Investigate size of SCV, and time between two samplings
# L1
ctd %>% 
  filter(datetime > as.POSIXct("2021-04-20 00:59:00") & datetime < as.POSIXct("2021-04-20 01:01:00")) 
# dist = 21.1 km from shore for the first lense

ctd %>% 
  mutate(dist = roundp(dist, precision = 0.1)) %>% 
  filter(dist == 21.1) %>% 
  filter(yo_nb == 85) %>% 
  pull(datetime) %>% 
  mean()
as.numeric(as.POSIXct("2021-04-18 23:27:43") %--% as.POSIXct("2021-04-20 01:00:00"), "hours")
# area of L1 was sampled 25.5 h before
ctd %>% filter(datetime > as.POSIXct("2021-04-19 21:59:00") & datetime < as.POSIXct("2021-04-19 22:01:00")) # L1 starts at 19 km
ctd %>% filter(datetime > as.POSIXct("2021-04-20 05:59:00") & datetime < as.POSIXct("2021-04-20 06:01:00")) # L1 ends at 25 km

# L2 and L3
transect_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = sal)) +
  scale_fill_viridis_c() +
  annotate("text", x = as.POSIXct("2021-04-21 15:00:00"), y = -220, label = "L2", size = 2.5, color = "white") +
  annotate("text", x = as.POSIXct("2021-04-22 08:00:00"), y = -220, label = "L3", size = 2.5, color = "white") +
  geom_vline(xintercept = ymd_hms("2021-04-21 13:00:00")) +
  geom_vline(xintercept = ymd_hms("2021-04-21 19:00:00")) +
  geom_vline(xintercept = ymd_hms("2021-04-22 06:00:00"), color = "red") +
  geom_vline(xintercept = ymd_hms("2021-04-22 13:00:00"), color = "red")

as.numeric(as.POSIXct("2021-04-21 15:00:00") %--% as.POSIXct("2021-04-22 08:00:00"), "hours")
# 17 hours between the two sampling

ctd %>% filter(datetime > as.POSIXct("2021-04-21 12:59:00") & datetime < as.POSIXct("2021-04-21 13:01:00")) # L2 starts at 52 km
ctd %>% filter(datetime > as.POSIXct("2021-04-21 18:59:00") & datetime < as.POSIXct("2021-04-21 19:01:00")) # L2 ends at 55 km
ctd %>% filter(datetime > as.POSIXct("2021-04-22 05:59:00") & datetime < as.POSIXct("2021-04-22 06:01:00")) # L3 ends at 55 km
ctd %>% filter(datetime > as.POSIXct("2021-04-22 12:59:00") & datetime < as.POSIXct("2021-04-22 13:01:00")) # L3 starts at 46 km

## And generate plot for supmat
p1 <- transect_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = temp)) +
  annotate("text", x = as.POSIXct("2021-04-14 12:00:00"), y = 100, label = "▿", size = 5) +
  annotate("text", x = as.POSIXct("2021-04-17 18:00:00"), y = 100, label = "▾", size = 5) +
  annotate("text", x = as.POSIXct("2021-04-19 09:00:00"), y = 100, label = "▿", size = 5) +
  annotate("text", x = as.POSIXct("2021-04-22 01:00:00"), y = 100, label = "▾", size = 5) +
  annotate("text", x = as.POSIXct("2021-04-24 00:00:00"), y = 100, label = "▿", size = 5) +
  coord_cartesian(clip = "off", ylim = c(-615, 0)) +
  scale_fill_cmocean(name = "thermal", na.value = NA) + 
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", x = as.POSIXct("2021-04-20 01:00:00"), y = -220, label = "L1", size = 2.5, color = "white") +
  annotate("text", x = as.POSIXct("2021-04-21 15:00:00"), y = -220, label = "L2", size = 2.5, color = "white") +
  annotate("text", x = as.POSIXct("2021-04-22 08:00:00"), y = -220, label = "L3", size = 2.5, color = "white") +
  labs(x = "", y = "Depth (m)", fill = "Temperature (°C)") +
  theme_minimal() +
  theme(legend.key.height = unit(0.4, "cm"), plot.margin = unit(c(15, 5.5, 5.5, 5.5), "points"), axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10), legend.title = element_text(size = 8))

p2 <- transect_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = sal)) +
  coord_cartesian(clip = "off", ylim = c(-615, 0)) +
  scale_fill_cmocean(name = "haline", na.value = NA) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Depth (m)", fill = "Salinity") +
  theme_minimal() +
  theme(legend.key.height = unit(0.4, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10), legend.title = element_text(size = 8))

p3 <- transect_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = oxy)) +
  coord_cartesian(clip = "off", ylim = c(-615, 0)) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Oxygen (µmol ", kg^{-1}, ")"))) +
  theme_minimal() +
  theme(legend.key.height = unit(0.4, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10), legend.title = element_text(size = 8))

p4 <- transect_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = cdom)) +
  coord_cartesian(clip = "off", ylim = c(-615, 0)) +
  scale_fill_cmocean(name = "matter", na.value = NA) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("CDOM (µg ", L^{-1}, ")"))) +
  theme_minimal() +
  theme(legend.key.height = unit(0.4, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10), legend.title = element_text(size = 8))

p5 <- transect_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = dens)) +
  coord_cartesian(clip = "off", ylim = c(-615, 0)) +
  scale_fill_cmocean(name = "dense", na.value = NA) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Density (kg ", m^{-3}, ")"))) +
  theme_minimal() +
  theme(legend.key.height = unit(0.4, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10), legend.title = element_text(size = 8))

p6 <- transect_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = datetime, y = -depth, fill = bb700)) +
  coord_cartesian(clip = "off", ylim = c(-615, 0)) +
  scale_fill_cmocean(name = "turbid", na.value = NA) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Datetime", y = "Depth (m)", fill = expression(paste("BB700 (", m^{-1}, ")"))) +
  theme_minimal() +
  theme(legend.key.height = unit(0.4, "cm"), text = element_text(size = 10), legend.title = element_text(size = 8))

p <- p1 / p3 / p2 / p4 / p5 / p6
p

ggsave(plot = p, file = "figures/15.sup_lenses.pdf", width = 164, height = 180, unit = "mm", dpi = 300)



## Small and large particles ----
#--------------------------------------------------------------------------#
# Investigate if the decrease of ratio of small to large particles in April is caused by:
# - a decrease in small particles
# - an increase in large particles

part_int <- read_parquet("data/12.part_int_fine.parquet")

# Definition of size classes
#class~21: 102-128~µm, 
#class~22: 128-161~µm, 
#class~23: 161-203~µm, 
#class~24: 203-256~µm, 
#class~25: 256-323~µm, 
#class~26: 323-406~µm, 
#class~27: 406-512~µm, 
#class~28: 512-645~µm, 
#class~29: 645-813~µm, 
#class~30: 813-1020~µm, 
#class~31: 1020-1290~µm, 
#class~32: 1290-1630~µm, 
#class~33: 1630-2050~µm.

# Regroup to large classes (same as in Ecopart)
class_names <- tribble(
  ~class,  ~label,
  "class08", "128 - 256 µm",
  "class09", "256 - 512 µm",
  "class10", "512 µm - 1.02 mm",
  "class11", "1.02 mm - 2.05 mm"
) %>% 
  mutate(label = fct_inorder(label))


part_int %>% 
  filter(depth < 200 & dist > 20) %>% ## Look offshore and not too deep
  mutate(
    class08 = class22 + class23 + class24,
    class09 = class25 + class26 + class27,
    class10 = class28 + class29 + class30,
    class11 = class31 + class32 + class33
  ) %>% 
  select(mission_det, dist, date, class08:class11) %>% 
  pivot_longer(cols = class08:class11, names_to = "class", values_to = "conc") %>% 
  # Average by day and by size class
  group_by(date, class) %>% 
  summarise(conc = mean(conc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # Add nice names for classes
  left_join(class_names) %>% 
  ggplot() +
  geom_path(aes(x = date, y = conc, color = class, group = class), show.legend = F) +
  facet_wrap(~label, scales = "free") +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(x = "Date", y = expression(paste("Averaged particle concentration (", m^{-3}, ")")))
ggsave(file = "figures/15.sup_small_big_parts.png", width = 164, height = 80, unit = "mm", dpi = 300, bg = "white")



## Distance for 1 m3 ----
#--------------------------------------------------------------------------#
# Compute the distance covered by the glider for the UVP6 to image 1m3 of water
# 1m3 takes about 20 min
# look at distance covered in 20 min
bins <- read_parquet("data/10.bins.parquet")

bins %>% 
  mutate(dt_round = round_date(datetime, unit = "20 min")) %>% 
  group_by(dt_round) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank == 1 | rank == max(rank)) %>% 
  ungroup()


## Tests of differences in concentrations at day VS night ----
#--------------------------------------------------------------------------#
## Wilcoxon test to compare mean yos conc at day VS night
dn_conc_wc <- plankton %>% 
  filter(depth < 315) %>% 
  group_by(mission_det, sample) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # Compute day and night
  mutate(date = date(datetime)) %>% 
  # compute datetime for nauticalDawn (beginning), sunrise (beginning), sunset (end) and nauticalDusk (end) from UTC datetime
  left_join(getSunlightTimes(data = select(., date, lat, lon), keep = c("nauticalDawn", "sunrise", "sunset", "nauticalDusk")),  by = c("lon", "lat", "date")) %>% 
  mutate(
    period = ifelse(datetime < nauticalDawn, # if UTC datetime is before nautical dawn, it is night
                    "night",
                    ifelse(datetime < sunrise, # if UTC datetime is before sunrise, it is dawn
                           "dawn",
                           ifelse(datetime < sunset, # if UTC datetime is before sunset, it is day
                                  "day",
                                  ifelse(datetime < nauticalDusk, # if UTC datetime is before nautical dusk, it is dusk
                                         "dusk",
                                         "night")))) # else if UTC datetime is after nautical dusk, it is night
  ) %>% 
  select(-c("date", "nauticalDawn", "sunrise", "sunset", "nauticalDusk")) %>% 
  filter(period %in% c("day", "night")) %>% 
  pivot_longer(Annelida:Salpida, names_to = "taxon", values_to = "conc")


# Empty list to store p-values
p_values <- c()
# Loop over taxa
for (ta in taxa){
  day_conc <- dn_conc_wc %>% filter(taxon == ta) %>% filter(period == "day") # Day concentrations
  night_conc <- dn_conc_wc %>% filter(taxon == ta) %>% filter(period == "night") # Night concentrations
  
  # Perform test
  res <- wilcox.test(day_conc$conc, night_conc$conc)
  # Store p value
  p_values <- c(p_values, res$p.value)
}

dn_wilcox <- tibble(
  taxon = taxa,
  p_value = p_values
) %>% 
  mutate(signif = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", ""))))


dn_conc_wc %>% 
  ggplot() +
  geom_boxplot(aes(x = taxon, y = conc, color = period), size = 0.5, outlier.size = 0.5) +
  geom_text(aes(x = taxon, y = 300, label = signif), data = dn_wilcox, size = 5) +
  scale_y_continuous(trans = "log1p", expand = c(0,0)) +
  scale_color_manual(values = c("gray", "black")) +
  coord_cartesian(clip = "off") +
  labs(x = "Taxon", y = expression(paste("Concentration (ind.", m^{-3}, ")")), color = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 10), plot.margin = unit(c(15, 5.5, 5.5, 5.5), "points"))
ggsave(file = "figures/15.sup_dn_conc_wc.pdf", width = 164, height = 80, unit = "mm", dpi = 300)


## Abundances on 10 vertical bins with Kolmogorov 
step_y <- 30
step_x <- 5
dn_conc_ks <- plankton %>% 
  filter(depth < 315) %>% 
  select(-sample) %>% 
  mutate(
    dist = roundp(dist, precision = step_x, f = floor) + step_x/2,
    depth = roundp(depth, precision = step_y, f = floor) + step_y/2
  ) %>% 
  group_by(mission_det, dist, depth) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # Compute day and night
  mutate(date = date(datetime)) %>% 
  # compute datetime for nauticalDawn (beginning), sunrise (beginning), sunset (end) and nauticalDusk (end) from UTC datetime
  left_join(getSunlightTimes(data = select(., date, lat, lon), keep = c("nauticalDawn", "sunrise", "sunset", "nauticalDusk")),  by = c("lon", "lat", "date")) %>% 
  mutate(
    period = ifelse(datetime < nauticalDawn, # if UTC datetime is before nautical dawn, it is night
                    "night",
                    ifelse(datetime < sunrise, # if UTC datetime is before sunrise, it is dawn
                           "dawn",
                           ifelse(datetime < sunset, # if UTC datetime is before sunset, it is day
                                  "day",
                                  ifelse(datetime < nauticalDusk, # if UTC datetime is before nautical dusk, it is dusk
                                         "dusk",
                                         "night")))) # else if UTC datetime is after nautical dusk, it is night
  ) %>% 
  select(-c("date", "nauticalDawn", "sunrise", "sunset", "nauticalDusk")) %>% 
  filter(period %in% c("day", "night")) %>% 
  pivot_longer(Annelida:Salpida, names_to = "taxon", values_to = "conc") %>% 
  group_by(depth, period, taxon) %>% 
  summarise(conc = mean(conc)) %>% 
  ungroup()


# Empty list to store p values
p_values <- c()
# Loop over taxa
for (ta in taxa){
  day_conc <- dn_conc_ks %>% filter(taxon == ta) %>% filter(period == "day") # Day concentrations
  night_conc <- dn_conc_ks %>% filter(taxon == ta) %>% filter(period == "night") # Night concentrations
  
  # Perform test
  res <- ks.test(day_conc$conc, night_conc$conc)
  # Store p value
  p_values <- c(p_values, res$p.value)
}

dn_ks <- tibble(
  taxon = taxa,
  p_value = p_values
) %>% 
  mutate(signif = ifelse(p_value < 0.001, " (***)", ifelse(p_value < 0.01, " (**)", ifelse(p_value < 0.05, " (*)", ""))))

dn_conc_ks %>% 
  left_join(dn_ks) %>% 
  mutate(taxon = str_c(taxon, signif)) %>% 
  mutate(
    conc = ifelse(period == "day", -conc, conc)
  ) %>% 
  ggplot() +
  geom_col(aes(x = -depth, y = conc, color = period), fill = "white") +
  scale_color_manual(values = c("gray", "black")) +
  coord_flip() +
  theme_minimal() +
  labs(x = "Depth (m)", y = expression(paste("Concentration (ind.", m^{-3}, ")")), color = "") +
  facet_wrap(~taxon, scales = "free") +
  theme(text = element_text(size = 10))
ggsave(file = "figures/15.sup_dn_conc_ks.pdf", width = 164, height = 100, unit = "mm", dpi = 300)


## Distance between surfacing events ----
#--------------------------------------------------------------------------#
surf_dist <- ctd %>% 
  select(mission_det, sample, yo_nb, depth, dist, lat, lon) %>% 
  filter(str_detect(sample, "d")) %>% # keep only downcasts
  group_by(mission_det, sample, yo_nb) %>% 
  filter(depth == min(depth)) %>% # keep only minimum depth of each downcast
  ungroup() %>% 
  group_by(mission_det) %>% 
  mutate(
    lon_end = lead(lon),
    lat_end = lead(lat),
    surf_dist = abs(dist - lead(dist))
  ) %>% 
  ungroup()


surf_dist %>% 
  ggplot() +
  geom_density(aes(x = surf_dist, color = mission_det)) +
  labs(x = "Distance between two surfacing events")

