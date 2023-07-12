#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Make all the plots
# Date: 11/07/2023
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(patchwork)
source("lib/lib.R")

## Read data ----
#--------------------------------------------------------------------------#
coast <- read_csv("data/coast.csv")

ctd <- read_parquet("data/10.ctd.parquet")
ctd_int_fine <- read_parquet("data/12.ctd_int_fine.parquet")
plankton_int <- read_parquet("data/12.plankton_int.parquet")

# Transects for the paper
paper_transects <- c("sea002_m487_part_2", "sea002_m489_part_1", "sea002_m491_part_2", "sea002_m496_part_2")


## Compute nights ----
#--------------------------------------------------------------------------#
# To generate grey rectangles in plot backgrounds to show nighttime
dist_dt <- ctd %>% 
  mutate(
    depth = roundp(depth, precision = 10),
    dist = roundp(dist, precision = 1)
  ) %>% 
  group_by(mission_det, dist) %>% 
  summarise(
    datetime = mean(datetime),
    lon = mean(lon),
    lat = mean(lat),
  ) %>% 
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
  select(mission_det, dist, period) %>%  
  filter(period == "night") %>% 
  group_by(mission_det) %>% 
  mutate(
    diff = dist - lag(dist),
    diff = ifelse(is.na(diff), 1, diff),
    night = cumsum(diff > 2) + 1
  ) %>% 
  ungroup()

dates <- ctd %>% 
  select(mission_det, datetime) %>% 
  group_by(mission_det) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  mutate(date = date(datetime)) %>% 
  select(-datetime)

nights <- dist_dt %>% 
  left_join(dates) %>% 
  select(-diff) %>% 
  group_by(mission_det, date, night) %>% 
  summarise(
    beg = min(dist),
    end = max(dist)
  ) %>% 
  ungroup()


## Plot env data for JO ----
#--------------------------------------------------------------------------#

ctd_plot_jo <- ctd %>% 
  filter(mission_det == "sea002_m487_part_2") %>% 
  select(dist, depth, temp, sal, chla, oxy) 

p1 <- ctd_plot_jo %>% 
  ggplot() +
  geom_point(aes(x = dist, y = -depth, color = temp), size = 0.5) +
  scale_color_cmocean(name = "thermal") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", color = "Temperature (°C)") +
  theme_minimal()

p2 <- ctd_plot_jo %>% 
  ggplot() +
  geom_point(aes(x = dist, y = -depth, color = sal), size = 0.5) +
  scale_color_cmocean(name = "haline") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", color = "Salinity") +
  theme_minimal()

p3 <- ctd_plot_jo %>% 
  ggplot() +
  geom_point(aes(x = dist, y = -depth, color = chla), size = 0.5) +
  scale_color_cmocean(name = "algae") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", color = expression(paste("Chlorophyll (mg ", m^{-3}, ")"))) +
  theme_minimal()

p4 <- ctd_plot_jo %>% 
  ggplot() +
  geom_point(aes(x = dist, y = -depth, color = oxy), size = 0.5) +
  scale_color_distiller(palette = "Blues", direction = 1, na.value = NA) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", color = expression(paste("Oxygen (µmol ", kg^{-1}, ")"))) +
  theme_minimal()


g <- p1 + p2 + p3 + p4 +  plot_layout(ncol = 2)
g
ggsave(plot = g, file = "figures/13.env_no_int_jo.png", width = 350, height = 200, unit = "mm", dpi = 300)
#ggsave(plot = g, file = "figures/13.env_no_int_jo.pdf", width = 350, height = 200, unit = "mm", dpi = 300)

ctd_int_plot_jo <- ctd_int_fine %>% 
  filter(mission_det == "sea002_m487_part_2") %>% 
  select(dist, depth, temp, sal, chla, oxy) 

p1 <- ctd_int_plot_jo %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = temp)) +
  scale_fill_cmocean(name = "thermal", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", fill = "Temperature (°C)") +
  theme_minimal()

p2 <- ctd_int_plot_jo %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = sal)) +
  scale_fill_cmocean(name = "haline", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", fill = "Salinity") +
  theme_minimal()

p3 <- ctd_int_plot_jo %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = chla)) +
  scale_fill_cmocean(name = "algae", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", fill = expression(paste("Chlorophyll (mg ", m^{-3}, ")"))) +
  theme_minimal()

p4 <- ctd_int_plot_jo %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth", fill = expression(paste("Oxygen (µmol ", kg^{-1}, ")"))) +
  theme_minimal()


g <- p1 + p2 + p3 + p4 +  plot_layout(ncol = 2)
ggsave(plot = g, file = "figures/env_int_jo.png", width = 350, height = 200, unit = "mm", dpi = 300)
ggsave(plot = g, file = "figures/env_int_jo.pdf", width = 350, height = 200, unit = "mm", dpi = 300)



## Plot env interpolated (paper) ----
#--------------------------------------------------------------------------#
p1 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = temp)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
  scale_fill_cmocean(name = "thermal", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = "Temperature (°C)", tag = "A") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_text(hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p2 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = sal)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
  scale_fill_cmocean(name = "haline", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = "Salinity", tag = "B") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p3 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = dens)) +
  geom_contour(aes(x = dist, y = -depth, z = dens), breaks = c(28.6, 28.8, 29), colour = "hotpink", linewidth = 0.2) +
  scale_fill_cmocean(name = "dense", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Density (kg ", m^{-3}, ")")), tag = "C") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p4 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = chla)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
  scale_fill_cmocean(name = "algae", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Chlorophyll (mg ", m^{-3}, ")")), tag = "D") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p5 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Oxygen (µmol ", kg^{-1}, ")")), tag = "E") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p6 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = cdom)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
  scale_fill_cmocean(name = "matter", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("CDOM (µg ", L^{-1}, ")")), tag = "F") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )


p7 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = bb700)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
  scale_fill_cmocean(name = "turbid", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(y = "Depth (m)", fill = expression(paste("BB700 (", m^{-1}, ")")), tag = "G") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal()+
  theme(
    legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )


p8 <- ctd_int_fine %>% 
  filter(mission_det %in% paper_transects) %>% 
  mutate(chla_bbp = chla/bb700) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = chla_bbp)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
  #scale_fill_cmocean(name = "matter", na.value = NA) +
  scale_fill_viridis_c(option = "mako", na.value = NA) +
  #scale_fill_distiller(palette = "YlGnBu") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Chl / BB700 (mg ", m^{-2}, ")")), tag = "H") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    legend.key.height = unit(0.25, "cm"), 
    text = element_text(size = 10), strip.text.x = element_blank(),
    axis.title.y = element_text(size = 8), legend.title = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )


p <- p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8
ggsave(plot = p, file = "figures/13.transects_env.pdf", width = 164, height = 190, unit = "mm", dpi = 300)


## Plot env interpolated (supp) ----
#--------------------------------------------------------------------------#
ctd_int_fine <- ctd_int_fine %>% mutate(chla_bbp = chla/bb700) %>% filter(depth < 300)

all_dates <- ctd_int_fine %>% pull(date) %>% unique()

lims = function(x){
  # Define limits so that colour scales are shared among plots
  return(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  }

for (batch in c(1:5)){
  print(batch)
  tag <- LETTERS[batch]
  dates_plot <- all_dates[((batch-1)*4+1):(batch*4)]
  nights_plot <- nights %>% filter(date %in% dates_plot)
  
  
  p1 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = temp)) +
    geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
    geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    scale_fill_cmocean(name = "thermal", na.value = NA, limits = lims(ctd_int_fine$temp)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(-300, 0)) +
    labs(x = "", y = "Depth (m)", fill = "Temperature (°C)", tag = tag) +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
      axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_text(hjust = 1), 
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  p2 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = sal)) +
    geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
    geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    scale_fill_cmocean(name = "haline", na.value = NA, limits = lims(ctd_int_fine$sal)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(-300, 0)) +
    labs(x = "", y = "Depth (m)", fill = "Salinity", tag = " ") +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
      axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  p3 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = dens)) +
    geom_contour(aes(x = dist, y = -depth, z = dens), breaks = c(28.6, 28.8, 29), colour = "hotpink", linewidth = 0.2) +
    geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    scale_fill_cmocean(name = "dense", na.value = NA, limits = lims(ctd_int_fine$dens)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(-300, 0)) +
    labs(x = "", y = "Depth (m)", fill = expression(paste("Density (kg ", m^{-3}, ")")), tag = " ") +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
      axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  p4 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = chla)) +
    geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
    geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    scale_fill_cmocean(name = "algae", na.value = NA, limits = lims(ctd_int_fine$chla)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    labs(x = "", y = "Depth (m)", fill = expression(paste("Chlorophyll (mg ", m^{-3}, ")")), tag = " ") +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
      axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  p5 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
    geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
    geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA, limits = lims(ctd_int_fine$oxy)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(-300, 0)) +
    labs(x = "", y = "Depth (m)", fill = expression(paste("Oxygen (µmol ", kg^{-1}, ")")), tag = " ") +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
      axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  p6 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = cdom)) +
    geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
    geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    scale_fill_cmocean(name = "matter", na.value = NA, limits = lims(ctd_int_fine$cdom)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(-300, 0)) +
    labs(x = "", y = "Depth (m)", fill = expression(paste("CDOM (µg ", L^{-1}, ")")), tag = " ") +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
      axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  p7 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = bb700)) +
    geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
    geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    scale_fill_cmocean(name = "turbid", na.value = NA, limits = lims(ctd_int_fine$bb700)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(-300, 0)) +
    labs(x = "", y = "Depth (m)", fill = expression(paste("BB700 (", m^{-1}, ")")), tag = "") +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal()+
    theme(
      legend.key.height = unit(0.25, "cm"), text = element_text(size = 10),
      axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), strip.text.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  
  p8 <- ctd_int_fine %>% 
    filter(date %in% dates_plot) %>% 
    ggplot() +
    geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights_plot) +
    geom_raster(aes(x = dist, y = -depth, fill = chla_bbp)) +
    geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2) +
    #geom_contour(aes(x = dist, y = -depth, z = chla), breaks = c(0.1), colour = "black", linewidth = 0.2) +
    #scale_fill_cmocean(name = "matter", na.value = NA) +
    scale_fill_viridis_c(option = "mako", na.value = NA, limits = lims(ctd_int_fine$chla_bbp)) +
    #scale_fill_distiller(palette = "YlGnBu") +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off", ylim = c(-300, 0)) +
    labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Chl / BB700 (mg ", m^{-2}, ")")), tag = "") +
    facet_wrap(~date, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.25, "cm"), 
      text = element_text(size = 10), strip.text.x = element_blank(),
      axis.title.y = element_text(size = 8), legend.title = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  
  p <- p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8
  p
  ggsave(plot = p, file = paste0("figures/13.transects_env_supp", batch, ".pdf"), width = 164, height = 190, unit = "mm", dpi = 300)
}



## Plot plankton interpolated (paper) ----
#--------------------------------------------------------------------------#

p1 <- plankton_int %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = Copepoda)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, breaks = c(0, 10, 100)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Copepoda (", m^{-3}, ")"))) +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  guides(fill = guide_colourbar(barwidth = 0.5, barheigh = 2.5)) +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p2 <- plankton_int %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = Appendicularia)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, breaks = c(0, 10, 30)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Appendicularia (", m^{-3}, ")"))) +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  guides(fill = guide_colourbar(barwidth = 0.5, barheigh = 2.5)) +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p3 <- plankton_int %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = Salpida)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "", y = "Depth (m)", fill = expression(paste("Salpida (", m^{-3}, ")"))) +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  guides(fill = guide_colourbar(barwidth = 0.5, barheigh = 2.5)) +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    axis.text.x = element_blank(), axis.title.x = element_blank(), 
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p4 <- plankton_int %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = Rhizaria)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Rhizaria (", m^{-3}, ")"))) +
  facet_wrap(~date, ncol = 4) +
  theme_minimal()+
  guides(fill = guide_colourbar(barwidth = 0.5, barheigh = 2.5)) +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    #axis.text.x = element_blank(), axis.title.x = element_blank(), 
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

p <- p1 / p2 / p3 / p4 #+ plot_layout(guides = 'collect')
p
ggsave(plot = p, file = "figures/13.transects_plankton.pdf", width = 160, height = 100, unit = "mm", dpi = 300)


## Plot plankton interpolated (supp) ----
#--------------------------------------------------------------------------#

# All transects
# Copepoda
p <- plankton_int %>% 
  filter(depth <= 300) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = Copepoda)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, breaks = c(0, 10, 100)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Copepoda (", m^{-3}, ")")), tag = "A") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
ggsave(plot = p, file = "figures/13.transects_pl_cop_supp.png", width = 160, height = 120, unit = "mm", dpi = 300)

# Appendicularia
p <- plankton_int %>% 
  filter(depth <= 300) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = Appendicularia)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, breaks = c(0, 10, 30)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Append. (", m^{-3}, ")")), tag = "B") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
ggsave(plot = p, file = "figures/13.transects_pl_app_supp.png", width = 160, height = 120, unit = "mm", dpi = 300)

# Salpida
p <- plankton_int %>% 
  filter(depth <= 300) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = Salpida)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, breaks = c(0, 10, 100)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Salpida (", m^{-3}, ")")), tag = "C") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
ggsave(plot = p, file = "figures/13.transects_pl_salp_supp.png", width = 160, height = 120, unit = "mm", dpi = 300)

# Rhizaria
p <- plankton_int %>% 
  filter(depth <= 300) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -320, ymax = 15), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = Rhizaria)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.2, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, breaks = c(0, 10, 20)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", ylim = c(-300, 0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Rhizaria (", m^{-3}, ")")), tag = "D") +
  facet_wrap(~date, ncol = 4) +
  theme_minimal() +
  theme(
    text = element_text(size = 10), legend.title = element_text(size = 10),
    strip.text.x = element_text(size = 10, hjust = 1),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
ggsave(plot = p, file = "figures/13.transects_pl_rhiz_supp.png", width = 160, height = 120, unit = "mm", dpi = 300)


## Current maps ----
#--------------------------------------------------------------------------#
# Formula to compute current projects (thanks Laurent!)
#phi = GLI_current2['heading'].values
#rho = GLI_current2['speed'].values
#
#Unorth = rho * np.cos(np.deg2rad(phi))
#Ueast = rho * np.sin(np.deg2rad(phi))


# List yos to keep for each mission plot of a mission
yo_back <- tibble(mission_det = paper_transects)  %>% 
  mutate(
    mission = str_split_fixed(mission_det, "_", n = 3)[, 2]
  ) %>% 
  left_join(
    ctd %>% 
      group_by(mission_det, yo_nb) %>% 
      arrange(datetime) %>% 
      slice(1) %>% 
      ungroup() %>% 
      select(mission_det, yo_nb, datetime) 
  ) %>% 
  left_join(dates)


# Read files containing trajectory data
files <- list.files("data/trajectory", full.names = TRUE)
df <- lapply(files, function(file){
  read_csv(file) %>% 
    mutate(mission = str_split_fixed(str_split_fixed(file, "/", n = 3)[3], "_", n = 3)[2] %>% tolower())
}) %>% 
  bind_rows() %>% 
  rename(yo_nb = YO_NUMBER, lat = Lat, lon = Lon, phi = heading, rho = speed)


# Compute current projections
traj <- yo_back %>% 
  left_join(df) %>% 
  mutate(
    u_north = rho * cos(deg2rad(phi)),
    u_east = rho * sin(deg2rad(phi))
  )
# Set a projection scaling factor
factor <- 250

# Define two reference points
sp_points <- tibble(
  name = c("Nice", "Dyf"),
  lat = c(43.695833, 43.418),
  lon = c(7.271389,  	7.869333)
)

traj %>% 
  group_by(mission_det, mission) %>% 
  mutate(time_since_beg = as.period(datetime - min(datetime))) %>% 
  ungroup()

# Plot it!
traj %>% 
  drop_na(lat) %>% 
  group_by(mission_det, mission) %>% 
  mutate(time_since_beg = as.numeric(datetime - min(datetime), "hours")) %>% 
  ungroup() %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_polygon(data = coast, fill = "grey80") +
  geom_segment(aes(xend = lon + u_east/factor, yend = lat + u_north/factor), arrow = arrow(length = unit(0.1,"cm"))) +
  geom_point(aes(colour = time_since_beg), size = 0.8) +
  geom_point(data = sp_points) +
  geom_text(aes(x = lon + 0.03, y = lat + 0.04, label = name), data = sp_points, size = 3) +
  #geom_path() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) + scale_y_continuous(expand = expansion(mult = c(0.05, 0))) + 
  scale_color_viridis_c() +
  labs(x = "Longitude", y = "Latitude", alpha = "Current \nspeed \n(cm s⁻¹)", colour = "Time since \nbeginning of \ntransect (h)") +
  coord_quickmap() +
  facet_wrap(~date) +
  theme_minimal() +
  theme(
    text = element_text(size = 10), strip.text.x = element_text(hjust = 1),
    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
ggsave(file = "figures/13.map_currents.pdf", width = 164, height = 108, unit = "mm", dpi = 300, bg = "white")


