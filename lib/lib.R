## Dir to save data on complex ----
#--------------------------------------------------------------------------#
save_dir <- "/remote/complex/home/tpanaiotis/datasets/glider_uvp6/data"

## Useful functions ----
#--------------------------------------------------------------------------#


scale2 <- function(x, na.rm = T){(x - mean(x, na.rm = na.rm)) / sd(x, na.rm)}

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

opposite <- function(x){return(-x)}


## Interpolate (coarse) a variable across a mission ----
#--------------------------------------------------------------------------#
coarse_interp_1v = function(df, variable, xo, yo) {
  #' Interpolate a variable across a mission.
  #' 
  #' Compute the linear interpolation of a variable on a given mission 
  #' with onto a specified output grid. 
  #' @param df Dataframe containing ISIIS data to interpolate
  #' @param variable name of variable to interpolate
  #' @param xo vector of x-coordinate of output grid (dist)
  #' @param yo vector of y-coordinate of output grid (depth)
  
  library(akima)
  library(broom)
  library(tidyverse)
  
  
  # Get mission name
  mission_det <- df %>% pull(mission_det) %>% unique()
  
  # Keep relevant columns and drop NA values
  df <- df %>% 
    select(dist, depth, all_of(variable)) %>% 
    drop_na(all_of(variable))
  
  # Perform the interpolation 
  df_int <- akima::interp(
    x=df$dist, 
    y=df$depth, 
    z=df[variable] %>% pull(), 
    linear=TRUE, 
    duplicate = "mean",
    xo = xo,
    yo = yo,
  ) %>% 
    tidy() %>% # convert list to dataframe
    rename(dist = x, depth = y, value = z) %>% # rename columns
    # change name of column containing newly interpolated variable with mutate and spread
    mutate(
      variable = variable,
      mission_det = mission_det,
    ) %>% 
    spread(variable, value)
  
  return(df_int)
}


## Interpolate (coarse) a variable across a mission with datetime ----
#--------------------------------------------------------------------------#
coarse_interp_1v_dt = function(df, variable, xo, yo) {
  #' Interpolate a variable across a mission.
  #' 
  #' Compute the linear interpolation of a variable on a given mission 
  #' with onto a specified output grid. 
  #' @param df Dataframe containing ISIIS data to interpolate
  #' @param variable name of variable to interpolate
  #' @param xo vector of x-coordinate of output grid (datetime)
  #' @param yo vector of y-coordinate of output grid (depth)
  
  library(akima)
  library(broom)
  library(tidyverse)
  
  
  # Get mission name
  mission <- df %>% pull(mission) %>% unique()
  
  # Keep relevant columns and drop NA values
  df <- df %>% 
    select(datetime, depth, all_of(variable)) %>% 
    drop_na(all_of(variable))
  
  # Perform the interpolation 
  df_int <- akima::interp(
    x=df$datetime, 
    y=df$depth, 
    z=df[variable] %>% pull(), 
    linear=TRUE, 
    duplicate = "mean",
    xo = xo,
    yo = yo,
  ) %>% 
    tidy() %>% # convert list to dataframe
    rename(datetime = x, depth = y, value = z) %>% # rename columns
    # change name of column containing newly interpolated variable with mutate and spread
    mutate(
      variable = variable,
      mission = mission,
    ) %>% 
    spread(variable, value)
  
  return(df_int)
}


## Interpolate (fine) a variable across a mission_det ----
#--------------------------------------------------------------------------#
fine_interp_1v = function(df_int, variable, plankton, step_x_fine, step_y_fine, theta = 0.5) {
  #' Interpolate a variable across a mission_det.
  #'
  #' Compute the linear interpolation of a variable on a given mission_det
  #' with onto a specified output grid.
  #' @param df_int dataframe of coarse interpolated data
  #' @param variable name of variale to interpolate 
  #' @param plankton list of plankton taxa (to ignore data around the thermocline)
  #' @param step_x_fine dimension of fine output grid in x direction (dist)
  #' @param step_y_fine dimension of fine output grid in y direction (depth)
  #' @param theta bandwidth or scale parameter
  
  
  library(akima)
  library(broom)
  library(tidyverse)
  
  # Get mission_det name
  mission_det <- df_int %>% pull(mission_det) %>% unique()
  
  ## Keep relevant columns and drop NA values
  df_int <- df_int %>%
    select(dist, depth, all_of(variable)) %>% 
    arrange(dist, depth)
  
  # Reformat data as akima::interp output
  x <- unique(df_int$dist)
  y <- unique(df_int$depth)
  z <- matrix(df_int[[variable]], nrow = length(x), byrow = TRUE)
  list_int <- list(x=x, y=y, z=z)
  
  # Perform fine interpolation
  df_int_fine <- list_int %>% 
    fields::interp.surface.grid(grid.list = list(
      x = seq(floor(min(x)), ceiling(max(x)), by = step_x_fine),
      y = seq(floor(min(y)), ceiling(max(y)), by = step_y_fine)
    )) %>%
    fields::image.smooth(theta = theta) %>%
    broom::tidy() %>%
    dplyr::rename(dist = x, depth = y, value = z) %>% # rename columns
    # change name of column containing newly interpolated variable with mutate and spread
    mutate(
      variable = variable,
      mission_det = mission_det
    ) %>%
    spread(variable, value)
  
  # Delete data in holes for plankton
  if (variable %in% plankton){
    
    df_int <- df_int %>% 
      rename(dist_round = dist, depth_round = depth, present = variable)
    
    # computed rounded distance and depth for fine interpolation grid
    df_int_fine <- df_int_fine %>% 
      mutate(
        dist_round = round(dist),
        depth_round = round(depth)
      ) %>% 
      # match with coarse grid
      left_join(df_int, by = c("dist_round", "depth_round")) %>% 
      mutate_at(variable, ~ ifelse(is.na(present), NA, .)) %>% 
      select(mission_det, dist, depth, all_of(variable))
    
  }
  
  
  return(df_int_fine)
}


## Interpolate (fine) a variable across a mission_det ----
#--------------------------------------------------------------------------#
fine_interp_1v_dt = function(df_int, variable, plankton, xo, yo, theta = 0.5) {
  #' Interpolate a variable across a mission_det.
  #'
  #' Compute the linear interpolation of a variable on a given mission_det
  #' with onto a specified output grid.
  #' @param df_int dataframe of coarse interpolated data
  #' @param variable name of variale to interpolate 
  #' @param plankton list of plankton taxa (to ignore data around the thermocline)
  #' @param xo output grid in x direction (datetime)
  #' @param yo output grid in y direction (depth)
  #' @param theta bandwidth or scale parameter
  
  
  library(akima)
  library(broom)
  library(tidyverse)
  
  # Get mission_det name
  mission <- df_int %>% pull(mission) %>% unique()
  
  ## Keep relevant columns and drop NA values
  df_int <- df_int %>%
    select(datetime, depth, all_of(variable)) %>% 
    arrange(datetime, depth)
  
  # Reformat data as akima::interp output
  x <- unique(df_int$datetime)
  y <- unique(df_int$depth)
  z <- matrix(df_int[[variable]], nrow = length(x), byrow = TRUE)
  list_int <- list(x=x, y=y, z=z)
  
  # Perform fine interpolation
  df_int_fine <- list_int %>% 
    fields::interp.surface.grid(grid.list = list(
      x = xo,
      y = yo
    )) %>%
    fields::image.smooth(theta = theta) %>%
    broom::tidy() %>%
    dplyr::rename(datetime = x, depth = y, value = z) %>% # rename columns
    # change name of column containing newly interpolated variable with mutate and spread
    mutate(
      variable = variable,
      mission = mission
    ) %>%
    spread(variable, value)
  
  # Delete data in holes for plankton
  if (variable %in% plankton){
    
    df_int <- df_int %>% 
      rename(datetime_round = datetime, depth_round = depth, present = variable)
    
    # computed rounded datetimeance and depth for fine interpolation grid
    df_int_fine <- df_int_fine %>% 
      mutate(
        datetime_round = round(datetime),
        depth_round = round(depth)
      ) %>% 
      # match with coarse grid
      left_join(df_int, by = c("datetime_round", "depth_round")) %>% 
      mutate_at(variable, ~ ifelse(is.na(present), NA, .)) %>% 
      select(mission_det, datetime, depth, all_of(variable))
    
  }
  
  
  return(df_int_fine)
}


## Plot interpolated values of concentration  per taxon and save plot ----
#--------------------------------------------------------------------------#
plot_interp_conc <- function(df, df_env, taxa, output_file){
  #' Plot interpolated values of concentration per taxon along a transect and save plot.
  #'
  #' Values are plotted as -depth = f(dist) with concentration in color. Plots are saved in a multipage pdf. 
  #' @param df dataframe with interpolated concentration data to plot
  #' @param df_env dataframe with interpolated environmental data
  #' @param taxa list of taxa to plot 
  #' @param output_file name of saved file 
  
  library(ggtext)
  library(qpdf)
  
  # Create output_dir
  output_dir <- tools::file_path_sans_ext(output_file)
  dir.create(output_dir)
  
  # Empty output_dir in case it already exists
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  for (my_taxon in taxa){
    
    # Plot 
    p <- df %>% 
      select(beg, dist, depth, all_of(my_taxon)) %>% 
      pivot_longer(all_of(my_taxon), names_to = "taxon", values_to = "conc") %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = conc)) +
      geom_contour(aes(x = dist, y = -depth, z = practical_salinity), breaks = c(38.2, 38.3), color = "white", size = 0.5, data = df_env) +
      scale_fill_viridis_c(trans = "log1p", na.value = NA) +
      facet_wrap(~beg, ncol = 4) +
      labs(x = "Distance (km)", y = "Depth (m)", fill = "Concentration<br>(ind.m<sup>-3</sup> )", title = my_taxon) +
      theme_minimal() +
      theme(legend.title = element_markdown(), text = element_text(size = 12)) 
    
    # Save
    ggsave(p, file=str_c(output_dir, "/", my_taxon, ".pdf"), width=15, height=10)
  }
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_taxa.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
  
  # Move multi transect file out of temporary dir
  file.rename(from = str_c(output_dir, "/all_taxa.pdf"), to = output_file)
  
  # Delete temporary dir
  unlink(output_dir, recursive = TRUE)
  
}


## Plot interpolated temperature, salinity, density, fluo and oxy for each transect and save as pdf ----
#--------------------------------------------------------------------------#
plot_interp_env <- function(df, env_vars, output_file){
  #' Make a combined pdf of plots of interpolated temp, sal, dens, fluo and oxy variables for each transect.
  #'
  #' The pdf has one page per transect, 5 plots per page, using the cmocean color map. 
  #' @param df Dataframe containing interpolated ISIIS data to plot
  #' @param output_file name of output file
  
  
  # Create output_dir
  output_dir <- tools::file_path_sans_ext(output_file)
  dir.create(output_dir)
  
  # Empty output_dir in case it already exists
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  for (my_var in env_vars){
    
    # Plot 
    p <- df %>% 
      select(beg, dist, depth, all_of(my_taxon)) %>% 
      pivot_longer(all_of(my_taxon), names_to = "taxon", values_to = "conc") %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = conc)) +
      geom_contour(aes(x = dist, y = -depth, z = practical_salinity), breaks = c(38.2, 38.3), color = "white", size = 0.5, data = df_env) +
      scale_fill_viridis_c(trans = "log1p", na.value = NA) +
      facet_wrap(~beg, ncol = 4) +
      labs(x = "Distance (km)", y = "Depth (m)", fill = "Concentration<br>(ind.m<sup>-3</sup> )", title = my_taxon) +
      theme_minimal() +
      theme(legend.title = element_markdown(), text = element_text(size = 12)) 
    
    # Save
    ggsave(p, file=str_c(output_dir, "/", my_taxon, ".pdf"), width=15, height=10)
  }  
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_transects.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
  
  # Move multi transect file out of temporary dir
  file.rename(from = str_c(output_dir, "/all_transects.pdf"), to = output_file)
  
  # Delete temporary dir
  unlink(output_dir, recursive = TRUE)
}


## Plot examples of objects in detritus clusters ----
#--------------------------------------------------------------------------#

plot_clust = function(objs, clust_plot, n_row, seed=NA) {
  #' Plot examples of objects for a detritus cluster
  #' @param objs df of all detritus objects, with "path_to_img" column
  #' @param clust_plot cluster to be plotted
  #' @param n_row number of rows in the plot, the number of plotted objects will be n_row*n_row
  #' @param seed optional seed for repeatability
  
  library(png)
  library(grid)
  library(gridExtra)
  
  n_img = n_row * n_row
  
  # Optionally set seed
  if (!is.na(seed)) {set.seed(seed)}
  
  # Select objects to plot
  objs_plot <- objs %>% 
    mutate(hac_clust = ind_scores$hac_clust) %>% 
    filter(hac_clust == clust_plot) %>% 
    slice_sample(n = n_img)
  # Extract file names
  filenames <- objs_plot$path_to_img
  # Read and store images
  imgs <- list()
  for (j in 1:n_img) {
    imgs[[j]] <- readPNG(filenames[j])
  }
  # Set plot layout
  par(mar = rep(0.5, 4))
  layout(matrix(1:n_img, nr=n_row, byr=T))
  # Plot images
  for (j in 1:n_img) {image(imgs[[j]], col=grey(seq(0, 1, length = 256)), axes=F, asp=1)}
  # Add title
  mtext(str_c("cluster ", clust_plot), side = 3, line = -2, outer = TRUE)
  # Get plot
  my_plot <- recordPlot()   
  
  # and return it
  return(my_plot)
}
