#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Exploratory analysis of glider + uvp6 data
# Date: 16/11/2022
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(arrow)
library(vegan)
library(ggrepel)
library(castr)
library(lubridate)
library(factoextra)
library(NbClust)
library(patchwork)
library(suncalc)
library(gbm)

source("lib/lib.R")



## Read data ----
#--------------------------------------------------------------------------#
coast <- read_csv("data/coast.csv")

ctd <- read_parquet("data/10.ctd.parquet")
plankton <- read_parquet("data/11.obj_conc.parquet")
parts <- read_parquet("data/11.parts_conc.parquet")

ctd_int_fine <- read_parquet(file = "data/12.ctd_int_fine.parquet")

ctd <- ctd %>% 
  select(-mission) %>% 
  select(mission_det, everything())

plankton <- plankton %>% 
  select(-yo_nb) %>% 
  select(mission_det, sample, lon, lat, datetime, dist, depth, everything())

parts <- parts %>% 
  select(mission_det, sample, lon, lat, datetime, dist, depth, class21:class33)


all_data <- ctd %>% 
  left_join(plankton) %>% 
  left_join(parts) %>% 
  # convert datetime to julian day
  mutate(datetime = yday(datetime)) %>% 
  filter(depth < 315)

taxa <- plankton %>% select(Annelida:Salpida) %>% colnames()

dates <- ctd %>% 
  select(mission_det, datetime) %>% 
  group_by(mission_det) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  mutate(date = date(datetime)) %>% 
  select(-datetime)

# Missions to plot
paper_transects <- c("sea002_m487_part_2", "sea002_m489_part_1", "sea002_m491_part_2", "sea002_m496_part_2")


## Compute nights ----
#--------------------------------------------------------------------------#
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



nights <- dist_dt %>% 
  left_join(dates) %>% 
  select(-diff) %>% 
  group_by(mission_det, date, night) %>% 
  summarise(
    beg = min(dist),
    end = max(dist)
  ) %>% 
  ungroup()


## PCA on log transformed particles data ----
#--------------------------------------------------------------------------#
# Group env, particles and plankton on same bins
step_x <- 1 # in km
step_y <- 10 # in m
all_data_g <- all_data %>% 
  select(mission_det, datetime, dist, depth, lat, lon, temp:cdom, class21:class33, all_of(taxa)) %>% 
  mutate(
    dist = roundp(dist, precision = step_x, f = floor) + step_x/2,
    depth = roundp(depth, precision = step_y, f = floor) + step_y/2
  ) %>% 
  group_by(mission_det, dist, depth) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  drop_na(Annelida) # delete bins where plankton is missing

env_g <- all_data_g %>% select(lat, lon, dist, depth, datetime, temp:cdom)
meta_g <- all_data_g %>% select(mission_det:datetime)
parts_g <- all_data_g %>% select(class21:class33) %>% drop_na(class21)
parts_g <- log1p(parts_g)


# Perform pca
pca_par <- rda(parts_g, scale = TRUE) 

# Extract eigenvalues
eig_par <- as.data.frame(t(summary(eigenvals(pca_par)))) %>% 
  rownames_to_column(var = "componant") %>% 
  as_tibble() %>% 
  rename(eigenvalue = Eigenvalue, prop_exp = `Proportion Explained`, cum_prop = `Cumulative Proportion`) %>% 
  mutate(componant = factor(componant, levels = componant))

# Plot eigenvalues
eig_par %>% 
  ggplot() +
  geom_col(aes(x = componant, y = prop_exp)) +
  theme_classic()

# Fit standardized env data on PCA axes
sup_fit_env <- envfit(pca_par ~ ., data = env_g, perm = 999, na.rm = T, choices=c(1:4))

# Extract coordinates of env data projection
sup_proj_env <- as.data.frame(sup_fit_env$vectors$arrows*sqrt(sup_fit_env$vectors$r)) %>% 
  rownames_to_column(var = "variable") %>% 
  as_tibble() %>% 
  mutate(orig = 0) %>% 
  mutate(across(PC1:PC4, ~opposite(.)))

# Extract PCA scores of profiles for plots (scaling 2)
sta_par_plot <- vegan::scores(pca_par, display="sites", choices=c(1:4), scaling=2) %>% 
  as_tibble() %>% 
  bind_cols(meta_g, .) %>% 
  mutate(across(PC1:PC4, ~opposite(.)))

# Extract PCA scores of profiles for computations (scaling 1)
sta_par_comp <- vegan::scores(pca_par, display="sites", choices=c(1:4), scaling=1) %>% 
  as_tibble() %>% 
  bind_cols(meta_g, .) %>% 
  mutate(across(PC1:PC4, ~opposite(.)))

# Extract variables scores in scaling 2
var_par_plot <- vegan::scores(pca_par, display="species", choices=c(1:4), scaling=2) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "variable") %>% 
  as_tibble() %>% 
  mutate(orig = 0) %>% 
  mutate(across(PC1:PC4, ~opposite(.)))

# PCA biplot
k <- 15
ggplot() +
  # Vertical and horizontal lines
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  theme_classic() +
  # sites
  geom_point(data = sta_par_plot, aes(x = k*PC1, y = k*PC2), alpha = 0.2, color = "gray") +
  # variables
  #geom_point(data = var_par_plot, aes(x = PC1, y = PC2)) +
  #geom_label_repel(data = var_par_plot, aes(x = PC1, y = PC2, label = variable), fill="white") +
  geom_segment(data = var_par_plot, aes(x = orig, y = orig, xend = PC1, yend = PC2),  colour = "black", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text_repel(data = var_par_plot, aes(x = PC1, y = PC2, label = variable), colour = "black") +
  # Segments for projected environmental data
  geom_segment(data = sup_proj_env, aes(x = orig, y = orig, xend = k*PC1, yend = k*PC2),  colour = "red", arrow = arrow(length = unit(0.2, "cm"))) +
  # Labels for projected environmental data
  geom_text(data = sup_proj_env, aes(x = k*PC1, y = k*PC2, label = variable), colour = "red", hjust=ifelse(sup_proj_env$PC1>0, -0.05, 1.05)) +
  # Labels
  ggtitle("PCA of log transformed particles concentrations") +
  xlab(paste0("PC 1", " (", format(round(100*eig_par$prop_exp[1], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC 2", " (", format(round(100*eig_par$prop_exp[2], 1), nsmall = 1), "%)")) +
  # Fixed ratio between axes
  coord_fixed() 
  

# PCA biplot for paper
k <- 5
ggplot() +
  # Vertical and horizontal lines
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  theme_classic() +
  # sites
  #geom_point(data = sta_par_plot, aes(x = k*PC1, y = k*PC2), alpha = 0.1, color = "gray80", size = 0.1) +
  # variables
  #geom_point(data = var_par_plot, aes(x = PC1, y = PC2)) +
  #geom_label_repel(data = var_par_plot, aes(x = PC1, y = PC2, label = variable), fill="white") +
  geom_segment(data = var_par_plot, aes(x = orig, y = orig, xend = PC1, yend = PC2),  colour = "black", arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text_repel(data = var_par_plot, aes(x = PC1, y = PC2, label = variable), colour = "black", size = 3) +
  # Segments for projected environmental data
  geom_segment(data = sup_proj_env, aes(x = orig, y = orig, xend = k*2*PC1, yend = k*2*PC2),  colour = "#1f78b4", arrow = arrow(length = unit(0.2, "cm"))) +
  # Labels for projected environmental data
  geom_text(data = sup_proj_env, aes(x = k*2*PC1, y = k*2*PC2, label = variable), colour = "#1f78b4", hjust=ifelse(sup_proj_env$PC1>0, -0.05, 1.05), size = 3) +
  # Labels
  xlab(paste0("PC 1", " (", format(round(100*eig_par$prop_exp[1], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC 2", " (", format(round(100*eig_par$prop_exp[2], 1), nsmall = 1), "%)")) +
  # Fixed ratio between axes
  coord_fixed(xlim = c(-5, 8), ylim = c(-4, 5)) +
  theme(text = element_text(size = 10))
ggsave(file = "figures/15.pca_particles_log.pdf", width = 164, height = 123, unit = "mm", dpi = 300)

ggplot() +
  # Vertical and horizontal lines
  #geom_vline(xintercept = 0, color = "gray") +
  #geom_hline(yintercept = 0, color = "gray") +
  theme_classic() +
  # sites
  geom_point(data = sta_par_plot, aes(x = k*PC1, y = k*PC2), alpha = 0.1, color = "gray80", size = 0.1) +
  # variables
  #geom_point(data = var_par_plot, aes(x = PC1, y = PC2)) +
  #geom_label_repel(data = var_par_plot, aes(x = PC1, y = PC2, label = variable), fill="white") +
  #geom_segment(data = var_par_plot, aes(x = orig, y = orig, xend = PC1, yend = PC2),  colour = "black", arrow = arrow(length = unit(0.2, "cm"))) +
  #geom_text_repel(data = var_par_plot, aes(x = PC1, y = PC2, label = variable), colour = "black") +
  # Segments for projected environmental data
  #geom_segment(data = sup_proj_env, aes(x = orig, y = orig, xend = k*PC1, yend = k*PC2),  colour = "#1f78b4", arrow = arrow(length = unit(0.2, "cm"))) +
  # Labels for projected environmental data
  #geom_text(data = sup_proj_env, aes(x = k*PC1, y = k*PC2, label = variable), colour = "#1f78b4", hjust=ifelse(sup_proj_env$PC1>0, -0.05, 1.05)) +
  # Labels
  xlab(paste0("PC 1", " (", format(round(100*eig_par$prop_exp[1], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC 2", " (", format(round(100*eig_par$prop_exp[2], 1), nsmall = 1), "%)")) +
  # Fixed ratio between axes
  coord_fixed(xlim = c(-5, 8), ylim = c(-4, 5)) +
  theme(text = element_text(size = 10))
ggsave(file = "figures/15.pca_particles_log_points.png", width = 164, height = 123, unit = "mm", dpi = 300)

# Bin stations at 10 m and 1 km
sta_par_plot_bin <- sta_par_plot %>% 
  left_join(dates)

# Plot PCA axes on transects
# Plot plankton PC1 on transects

lim_PC1 <- quantile(sta_par_plot_bin$PC1, c(0.001, 0.999)) %>% abs() %>% max()
lim_PC2 <- quantile(sta_par_plot_bin$PC2, c(0.001, 0.999)) %>% abs() %>% max()

sta_par_plot_bin %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -Inf, ymax = Inf), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = PC1)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "black", linewidth = 0.5, data = ctd_int_fine) +
  geom_contour(aes(x = dist, y = -depth, z = dens), color = "black", breaks = c(28.6, 28.8, 29), linewidth = 0.5, linetype = "dotted", data = ctd_int_fine) + #breaks = c(38.2, 38.3), 
  geom_contour(aes(x = dist, y = -depth, z = chla), color = "darkgreen", breaks = c(0.05), linewidth = 0.5, linetype = "dashed", data = ctd_int_fine %>% filter(date %in% c(as.Date("2021-03-20")))) +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "#f7f7f7", limits = c(-lim_PC1, lim_PC1), na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(-300, 0)) +
  facet_wrap(~date, ncol = 4) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", tag = "A") +
  #ggtitle("PC1 distribution across transects") + 
  theme_minimal() +
  theme(text = element_text(size = 8), strip.text.x = element_text(size = 10, hjust = 1))
ggsave(file = "figures/15.sup_pc1_particles_log.pdf", width = 164, height = 160, unit = "mm", dpi = 300)


# Plot plankton PC2 on transects
sta_par_plot_bin %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -Inf, ymax = Inf), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = PC2)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "black", linewidth = 0.5, data = ctd_int_fine) +
  geom_contour(aes(x = dist, y = -depth, z = dens), color = "black", breaks = c(28.6, 28.8, 29), linewidth = 0.5, linetype = "dotted", data = ctd_int_fine) + #breaks = c(38.2, 38.3), 
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "#f7f7f7", limits = c(-lim_PC2, lim_PC2), na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(-300, 0)) +
  facet_wrap(~date, ncol = 4) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", tag = "B") +
  #ggtitle("PC2 distribution across transects") + 
  theme_minimal() +
  theme(text = element_text(size = 8), strip.text.x = element_text(size = 10, hjust = 1))
ggsave(file = "figures/15.sup_pc2_particles_log.pdf", width = 164, height = 160, unit = "mm", dpi = 300)


# Plots for paper
p1 <- sta_par_plot_bin %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -Inf, ymax = Inf), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = PC1)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "black", linewidth = 0.2, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) +
  geom_contour(aes(x = dist, y = -depth, z = dens), breaks = c(28.6, 28.8, 29), linetype = "dotted", color = "black", linewidth = 0.5, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) + 
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "#f7f7f7", limits = c(-lim_PC1, lim_PC1)) +
  #scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~date, ncol = 4) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "") +
  #ggtitle("Particles richness") + 
  theme_minimal() +
  guides(fill = guide_colourbar(ticks = FALSE, label = FALSE, barwidth = 0.5, barheigh = 3)) +
  theme(legend.key.height = unit(0.5, "cm"), axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10), strip.text.x = element_text(size = 10, hjust = 1))

p2 <- sta_par_plot_bin %>% 
  filter(mission_det %in% paper_transects) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -Inf, ymax = Inf), fill = "gray80", data = nights %>% filter(mission_det %in% paper_transects)) +
  geom_raster(aes(x = dist, y = -depth, fill = PC2)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "black", linewidth = 0.2, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) +
  geom_contour(aes(x = dist, y = -depth, z = dens), breaks = c(28.6, 28.8, 29), linetype = "dotted", color = "black", linewidth = 0.5, data = ctd_int_fine %>% filter(mission_det %in% paper_transects)) + 
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "#f7f7f7", limits = c(-lim_PC2, lim_PC2)) +
  #scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~date, ncol = 4) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "") +
  #ggtitle("Particles size") + 
  theme_minimal() +
  guides(fill = guide_colourbar(ticks = FALSE, label = FALSE, barwidth = 0.5, barheigh = 3)) +
  theme(strip.text.x = element_blank(), text = element_text(size = 10))


p <- p1 / p2

ggsave(p, file = "figures/15.transects_particles_log.pdf", width = 164, height = 70, unit = "mm", dpi = 300)



## PCA on log transformed plankton data ----
#--------------------------------------------------------------------------#
# Generate larger bins for plankton
step_x <- 5 # in km
step_y <- 30 # in m
all_data_g <- all_data %>% 
  select(mission_det, datetime, dist, depth, lat, lon, temp:cdom, class21:class33, all_of(taxa)) %>% 
  mutate(
    dist = roundp(dist, precision = step_x, f = floor) + step_x/2,
    depth = roundp(depth, precision = step_y, f = floor) + step_y/2
  ) %>% 
  group_by(mission_det, dist, depth) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  drop_na(Annelida) # delete bins where plankton is missing

plankton_g <- all_data_g %>% select(mission_det:lon, all_of(taxa))
env_g <- all_data_g %>% select(lat, lon, dist, depth, datetime, temp:class33)
meta_g <- all_data_g %>% select(mission_det:datetime)

plankton_g <- plankton_g %>% mutate(across(all_of(taxa), ~log1p(.))) 

# Replace empty bins by averaged concentration so they do not contribute to PCA space
plankton_g <- plankton_g %>% 
  pivot_longer(cols = all_of(taxa), names_to = "taxon", values_to = "conc") %>% 
  mutate(conc = ifelse(conc == 0, NA, conc)) %>% 
  group_by(taxon) %>% 
  mutate(conc = ifelse(is.na(conc), mean(conc, na.rm = TRUE), conc)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "taxon", values_from = "conc") %>% 
  select(all_of(taxa))
  

# Perform pca
pca_plh <- rda(plankton_g, scale = TRUE) # do not scale after log transform


# Extract eigenvalues
eig_plh <- as.data.frame(t(summary(eigenvals(pca_plh)))) %>% 
  rownames_to_column(var = "componant") %>% 
  as_tibble() %>% 
  rename(eigenvalue = Eigenvalue, prop_exp = `Proportion Explained`, cum_prop = `Cumulative Proportion`) %>% 
  mutate(componant = factor(componant, levels = componant))

# Plot eigenvalues
eig_plh %>% 
  ggplot() +
  geom_col(aes(x = componant, y = prop_exp)) +
  ggtitle("PCA on log transformed plankton concentrations") +
  theme_classic()

# Fit standardized env data on PCA axes
sup_fit_env <- envfit(pca_plh ~ ., data = env_g, perm = 999, na.rm = T, choices=c(1:5))


# Extract coordinates of env data projection
sup_proj_env <- as.data.frame(sup_fit_env$vectors$arrows*sqrt(sup_fit_env$vectors$r)) %>% 
  rownames_to_column(var = "variable") %>% 
  as_tibble() %>% 
  mutate(orig = 0)

# Extract PCA scores of profiles for plots (scaling 2)
sta_plh_plot <- vegan::scores(pca_plh, display="sites", choices=c(1:5), scaling=2) %>% 
  as_tibble() %>% 
  bind_cols(meta_g, .)

# Extract PCA scores of profiles for computations (scaling 1)
sta_plh_comp <- vegan::scores(pca_plh, display="sites", choices=c(1:5), scaling=1) %>% 
  as_tibble() %>% 
  bind_cols(meta_g, .)

# Extract variables scores in scaling 2
var_plh_plot <- vegan::scores(pca_plh, display="species", choices=c(1:5), scaling=2) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "variable") %>% 
  as_tibble() %>% 
  mutate(orig = 0) 

var_plh_plot <- var_plh_plot %>% 
  mutate(
    contrib = sqrt(PC1^2 + PC2^2),
    contrib = contrib - mean(contrib)
  ) #%>% 
  #filter(contrib > 0)



# PCA biplot
k <- 10
ggplot() +
  # Vertical and horizontal lines
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  theme_classic() +
  geom_point(data = sta_plh_plot, aes(x = PC1, y = PC2), alpha = 0.5, color = "gray", size = 0.2) +
  geom_segment(data = var_plh_plot, aes(x = orig, y = orig, xend = PC1, yend = PC2),  colour = "black", arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.2) +
  geom_text_repel(data = var_plh_plot, aes(x = PC1, y = PC2, label = variable), colour = "black", size = 2) +
  #geom_point(data = var_plh_plot, aes(x = PC1, y = PC2)) +
  #geom_label_repel(data = var_plh_plot, aes(x = PC1, y = PC2, label = variable), fill="white") +
  # Segments for projected environmental data
  geom_segment(data = sup_proj_env, aes(x = orig, y = orig, xend = 10*PC1, yend = 10*PC2),  colour = "#1f78b4", arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.2) +
  # Labels for projected environmental data
  geom_text(data = sup_proj_env, aes(x = k*PC1, y = k*PC2, label = variable), colour = "#1f78b4", hjust=ifelse(sup_proj_env$PC1>0, -0.05, 1.05), size = 2) +
  # Labels
  xlab(paste0("PC 1", " (", format(round(100*eig_plh$prop_exp[1], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC 2", " (", format(round(100*eig_plh$prop_exp[2], 1), nsmall = 1), "%)")) +
  # Fixed ratio between axes
  coord_fixed() +
  #coord_fixed(xlim = c(-15, 12), ylim = c(-10, 5)) +
  theme(text = element_text(size = 8))
ggsave(file = "figures/15.sup_pca_plankton_log_12.pdf", width = 90, height = 90, unit = "mm", dpi = 300)

# PCA biplot 2-3
k <- 10
ggplot() +
  # Vertical and horizontal lines
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  theme_classic() +
  geom_point(data = sta_plh_plot, aes(x = PC2, y = PC3), alpha = 0.5, color = "gray", size = 0.2) +
  geom_segment(data = var_plh_plot, aes(x = orig, y = orig, xend = PC2, yend = PC3),  colour = "black", arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.2) +
  geom_text_repel(data = var_plh_plot, aes(x = PC2, y = PC3, label = variable), colour = "black", size = 2) +
  #geom_point(data = var_plh_plot, aes(x = PC2, y = PC3)) +
  #geom_label_repel(data = var_plh_plot, aes(x = PC2, y = PC3, label = variable), fill="white") +
  # Segments for projected environmental data
  geom_segment(data = sup_proj_env, aes(x = orig, y = orig, xend = k*PC2, yend = k*PC3),  colour = "#1f78b4", arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.2) +
  # Labels for projected environmental data
  geom_text(data = sup_proj_env, aes(x = k*PC2, y = k*PC3, label = variable), colour = "#1f78b4", hjust=ifelse(sup_proj_env$PC2>0, -0.05, 1.05), size = 2) +
  # Labels
  xlab(paste0("PC 2", " (", format(round(100*eig_plh$prop_exp[2], 1), nsmall = 1), "%)")) +
  ylab(paste0("PC 3", " (", format(round(100*eig_plh$prop_exp[3], 1), nsmall = 1), "%)")) +
  # Fixed ratio between axes
  coord_fixed() +
  #coord_fixed(xlim = c(-15, 12), ylim = c(-10, 5)) +
  theme(text = element_text(size = 8))
ggsave(file = "figures/15.sup_pca_plankton_log_23.pdf", width = 90, height = 90, unit = "mm", dpi = 300)

lim_PC1 <- max(abs(sta_plh_plot$PC1))
lim_PC2 <- max(abs(sta_plh_plot$PC2))
lim_PC3 <- max(abs(sta_plh_plot$PC3))

# Plot plankton PC1 on transects
sta_plh_plot %>% 
  left_join(dates) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -Inf, ymax = Inf), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = PC1)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "black", linewidth = 0.5, data = ctd_int_fine) +
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low = "#018571", mid = "#f5f5f5", high = "#a6611a", limits = c(-lim_PC1, lim_PC1)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(-300, 0)) +
  facet_wrap(~date, ncol = 4) +
  labs(x = "Distance from shore (km)", y = "Depth (m)") +
  #ggtitle("PC1 distribution across transects") + 
  theme_minimal() +
  theme(text = element_text(size = 10))
ggsave(file = "figures/15.sup_pc1_plankton_log.pdf", width = 164, height = 160, unit = "mm", dpi = 300)


# Plot plankton PC2 on transects
sta_plh_plot %>% 
  left_join(dates) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -Inf, ymax = Inf), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = PC2)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "black", linewidth = 0.5, data = ctd_int_fine) +
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low = "#018571", mid = "#f5f5f5", high = "#a6611a", limits = c(-lim_PC2, lim_PC2)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~date, ncol = 4) +
  labs(x = "Distance from shore (km)", y = "Depth (m)") +
  ggtitle("PC2 distribution across transects") + 
  theme_minimal()

# Plot plankton PC3 on transects
sta_plh_plot %>% 
  left_join(dates) %>% 
  ggplot() +
  geom_rect(aes(xmin = beg, xmax = end, ymin = -Inf, ymax = Inf), fill = "gray80", data = nights) +
  geom_raster(aes(x = dist, y = -depth, fill = PC3)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), color = "black", linewidth = 0.5, data = ctd_int_fine) +
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low = "#018571", mid = "#f5f5f5", high = "#a6611a", limits = c(-lim_PC3, lim_PC3)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~date, ncol = 4) +
  labs(x = "Distance from shore (km)", y = "Depth (m)") +
  ggtitle("PC3 distribution across transects") + 
  theme_minimal()


## Part of variance explained by day/night variations ----
#--------------------------------------------------------------------------#
# Prepara data, separate day and night
plankton_dn <- all_data_g %>% 
  select(-datetime) %>% 
  left_join(dates) %>% 
  mutate(dist = roundp(dist, precision = 1)) %>% 
  left_join(dist_dt %>% select(mission_det:period)) %>% 
  mutate(period = ifelse(is.na(period), "day", period)) %>% 
  select(date, period, dist, depth, temp:oxy, Copepoda, Salpida, Appendicularia, Rhizaria) %>% 
  mutate(
    period = factor(period)
  ) %>% 
  drop_na() %>% 
  filter(depth == 45) # In the 30-60 m depth range

# Initate empty tibble
inf_dn <- tibble()
# Loop over taxa and fit model
for (taxon in c("Copepoda", "Salpida", "Appendicularia", "Rhizaria")){
  m <- gbm(as.formula(paste(taxon, "~ temp + sal + chla + oxy + period")), data = plankton_dn, distribution="gaussian")
  rel_inf <- relative.influence(m, n.trees=100)
  inf_dn <- inf_dn %>% 
    bind_rows(
      as_tibble_row(rel_inf / sum(rel_inf)) %>% mutate(taxon = taxon))
}

inf_dn
