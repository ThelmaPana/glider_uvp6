#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Generate exploratory plots
# Date: 11/07/2023
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


## Read data ----
#--------------------------------------------------------------------------#
ctd <- read_parquet("data/10.ctd.parquet")
ctd_int_fine <- read_parquet("data/12.ctd_int_fine.parquet")


## Plot interpolated CTD data ----
#--------------------------------------------------------------------------#

## Clean output directory
output_dir <- "plots/env/interpolated_ctd"
dir.create(output_dir)
files <- list.files(output_dir, full.names = TRUE)
file.remove(files)

p <- ctd_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = temp)) +
  scale_fill_cmocean(name = "thermal", na.value = NA) +
  theme_minimal() +
  facet_wrap(~date, ncol = 4) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Temperature", fill = "°C")
ggsave(p, file = file.path(output_dir, "01_temperature.pdf"), width=15, height=10)


p <- ctd_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = sal)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.5) +
  scale_fill_cmocean(name = "haline", na.value = NA) +
  theme_minimal() +
  facet_wrap(~date, ncol = 4) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Salinity", fill = "")
ggsave(p, file = file.path(output_dir, "02_salinity.pdf"), width=15, height=10)

p <- ctd_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = dens)) +
  geom_contour(aes(x = dist, y = -depth, z = dens), breaks = c(28.6, 28.8, 29), colour = "hotpink", linewidth = 0.5) +
  scale_fill_cmocean(name = "dense", na.value = NA) +
  theme_minimal() +
  facet_wrap(~date, ncol = 4) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Density", fill = "kg.m-3")
ggsave(p, file = file.path(output_dir, "03_density.pdf"), width=15, height=10)

p <- ctd_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = chla)) +
  scale_fill_cmocean(name = "algae", na.value = NA) +
  theme_minimal() +
  facet_wrap(~date, ncol = 4) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Chl a", fill = "mg.m-3")
ggsave(p, file = file.path(output_dir, "04_chla.pdf"), width=15, height=10)

p <- ctd_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
  scale_fill_distiller(palette = "Blues", na.value = NA, direction = 1) +
  theme_minimal() +
  facet_wrap(~date, ncol = 4) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Oxygen", fill = "µmol.kg-1")
ggsave(p, file = file.path(output_dir, "05_oxygen.pdf"), width=15, height=10)

p <- ctd_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = cdom)) +
  scale_fill_cmocean(name = "matter", na.value = NA, trans = "log1p") +
  theme_minimal() +
  facet_wrap(~date, ncol = 4) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "CDOM scaled", fill = "µg.L-1")
ggsave(p, file = file.path(output_dir, "06_cdom.pdf"), width=15, height=10)


p <- ctd_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = bb700)) +
  scale_fill_cmocean(name = "turbid", na.value = NA) +
  theme_minimal() +
  facet_wrap(~date, ncol = 4) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "BB700 scaled", fill = "m-1")
ggsave(p, file = file.path(output_dir, "07_bb700.pdf"), width=15, height=10)


# List of newly created plots 
my_plots <- list.files(output_dir, full.names = TRUE)

# Combine all pages in one pdf
qpdf::pdf_combine(input = my_plots, output = paste0(output_dir, ".pdf"))

# Delete single page plots
file.remove(my_plots)

# Delete temporary dir
unlink(output_dir, recursive = TRUE)


## Plot particles ----
#--------------------------------------------------------------------------#
# Clean output directory
output_dir <- "plots/particles/interpolated_parts"
dir.create(output_dir)
files <- list.files(output_dir, full.names = TRUE)
file.remove(files)

p <- part_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = class07)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", size = 0.5, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, labels = label_scientific())+#, breaks = c(1000, 10000, 50000)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Particles 64-128 µm", fill = "#.L-1") +
  facet_wrap(~date, ncol = 4) +
  #theme(text = element_text(size = 18), legend.key.height = unit(1, "cm"))
  theme(legend.key.height = unit(1, "cm"))
ggsave(p, file = file.path(output_dir, "class_07.pdf"), width=15, height=10)


p <- part_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = class08)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", size = 0.5, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Particles 128-256 µm", fill = "#.L-1") +
  facet_wrap(~date, ncol = 4) +
  #theme(text = element_text(size = 18), legend.key.height = unit(1.5, "cm"))
  theme(legend.key.height = unit(1.5, "cm"))
ggsave(p, file = file.path(output_dir, "class_08.pdf"), width=15, height=10)


p <- part_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = class09)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", size = 0.5, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, labels = label_scientific(), breaks = c(10, 100, 1000, 2000)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Particles 256-512 µm", fill = "#.L-1") +
  facet_wrap(~date, ncol = 4) +
  #theme(text = element_text(size = 18), legend.key.height = unit(1, "cm"))
  theme(legend.key.height = unit(1, "cm"))
ggsave(p, file = file.path(output_dir, "class_09.pdf"), width=15, height=10)


p <- part_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = class10)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", size = 0.5, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Particles 0.512 - 1.02 mm", fill = "#.L-1") +
  facet_wrap(~date, ncol = 4) +
  #theme(text = element_text(size = 18), legend.key.height = unit(1.5, "cm"))
  theme(legend.key.height = unit(1.5, "cm"))
ggsave(p, file = file.path(output_dir, "class_10.pdf"), width=15, height=10)


p <- part_int_fine %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = class11)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", size = 0.5, data = ctd_int_fine) +
  scale_fill_viridis_c(trans = "log1p", na.value = NA, labels = label_scientific()) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Particles 1.02 - 2.05 mm", fill = "#.L-1") +
  facet_wrap(~date, ncol = 4) +
  #theme(text = element_text(size = 18), legend.key.height = unit(1, "cm"))
  theme(legend.key.height = unit(1, "cm"))
ggsave(p, file = file.path(output_dir, "class_11.pdf"), width=15, height=10)


# List of newly created plots 
my_plots <- list.files(output_dir, full.names = TRUE)

# Combine all pages in one pdf
qpdf::pdf_combine(input = my_plots, output = paste0(output_dir, ".pdf"))

# Delete single page plots
file.remove(my_plots)

# Delete temporary dir
unlink(output_dir, recursive = TRUE)



## Plot watervolumes ----
#--------------------------------------------------------------------------#
bins %>% 
  #filter(mission_det == "sea002_m495_part_1") %>% 
  #mutate(dist = roundp(dist, precision = 0.5)) %>% 
  #group_by(mission_det, dist, depth) %>% 
  #summarise(watervolume = mean(watervolume, na.rm = TRUE)) %>% 
  #ungroup() %>% 
  ggplot() +
  geom_point(aes(x = dist, y = -depth, color = watervolume)) +
  scale_color_viridis_c() +
  theme_minimal() +
  facet_wrap(~mission_det) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  #theme(text = element_text(size = 18)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", title = "Sampled volume per 5-m bin", color = "L")

bins %>% filter(mission_det == "sea002_m495_part_1") %>% filter(depth < 200) %>% summary()

# Volume of larger bins
step_x <- 5 # in km
step_y <- 30 # in m

large_bins <- bins %>% 
  mutate(
    dist = roundp(dist, precision = step_x, f = floor) + step_x/2,
    depth = roundp(depth, precision = step_y, f = floor) + step_y/2
  ) %>% 
  group_by(mission_det, dist, depth) %>% 
  summarise(
    tot_vol = sum(watervolume),
    datetime_avg = mean(datetime)
  ) %>% 
  ungroup()
summary(large_bins)

large_bins %>% 
  ggplot() +
  geom_histogram(aes(x = tot_vol)) +
  labs(x = "Volume of larger bins (L)", title = "Plankton bins volume") +
  theme_classic()


