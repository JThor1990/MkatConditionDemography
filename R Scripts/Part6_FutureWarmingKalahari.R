###########################################################

# LINKING CLIMATE VARIABILITY TO DEMOGRAPHY IN ARID ENVIRONMENTS: A MECHANISTIC STUDY IN KALAHARI MEERKATS 

# AUTHORS: Thorley, Duncan, ... and Clutton-Brock (2024) 

# PART 6: PROJECTING THE FUTURE NUMBER OF HOT DAYS ACROSS THE DISTRIBUTION OF MEERKATS

###########################################################

# Load required packages
lapply(c("ggh4x", "ggplot2", "ggpattern", "ggstar", "patchwork", "RColorBrewer", "rnaturalearth", "rnaturalearthdata", "sf", "tidyverse", "terra"), FUN = library, character.only = TRUE)

# Set working directory
setwd("INSERT FILE PATH HERE")

# Load in the meerkat shapefiles via terra
mkat_shp <- vect("Data//data_0.shp")
plot(mkat_shp, col = "darkorange", alpha = 0.5) # basic plot using terra

# Obtain the Country Boundaries
# Get the country boundaries for Africa
africa <- ne_countries(continent = "Africa", returnclass = "sf")
# Filter for Southern Africa (adjust countries as necessary)
southern_africa_countries <- c("Angola", "Zambia", "South Africa", "Namibia", "Botswana", "Zimbabwe", "Mozambique", "Lesotho", "eSwatini", "Uganda", "Malawi")
southern_africa <- africa[africa$name %in% southern_africa_countries, ]
# crop Southern Africa to below = 14.5 degrees  
bbox <- st_bbox(c(xmin = 10, xmax = 41, ymin = -90,  ymax = -14.5))  
southern_africa <- st_crop(southern_africa, bbox)
# Convert to terra SpatVector object
southern_africa <- vect(southern_africa)
# Clip the species distribution shapefile to the boundaries of Southern Africa to smooth off the edges where the distribution falls into the sea
mkat_shp <- intersect(mkat_shp, southern_africa)
# and make sure it is a single polygon again (union outside boundary)
mkat_shp <- aggregate(mkat_shp, dissolve = TRUE)
# also aggregate southern Africa for later used 
southern_africa_agg <- aggregate(southern_africa, dissolve = TRUE)
  
# Convert the two to sf objects
mkat_sf <- st_as_sf(mkat_shp)
southern_africa <- st_as_sf(southern_africa)

# Get coordinates of the Kalahari
krr_lat <- -26.979
krr_lon <- 21.832

# set up a base map to adjust later
base_map <- ggplot() +
    geom_sf(data = southern_africa, fill = "white", color = "black", lwd = 0.6) +  # Plot country
  geom_sf_pattern(
  data = mkat_sf, aes(), pattern = "stripe",
  fill = 'white',
  alpha = 0,
  pattern_colour = "black",
  pattern_fill = "black",
  pattern_size = 0.1,
  pattern_density = 0.01,  
  pattern_angle = 45) +
    geom_star(aes(x = krr_lon, y = krr_lat), fill = "red", color = "black", size = 4) + 
    theme_minimal() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.title = element_blank()) + 
    annotate ("rect", xmin= 12, xmax = 41, ymin = -14, ymax= -15.1,
              color="white", fill="white")
base_map

#-------------------------------------
  
# Incorporate the CMIP6 historic data and forecasts and plot them out

# Load in the data 
cmip6 <- read.csv("Data//MeerkatCMIP6Tmax.csv", header = TRUE)
cmip6_con <- read.csv("Data//MeerkatCMIP6ConsecutiveTmax.csv", header = TRUE)

# Rework the data before aggregating
 # For each period, work out the mean number of days per year / pixel over the threshold temperatures 
 # and the number of consecutive days over the given thresholds / pixel 
 # (i.e. longest continuous period during which daily max temp never decreases below the threshold)
 # Note that the latter was calculated for breeding seasons (Jul - Jul) so as not to divide the austral summer 
cmip6_agg <- cmip6 %>% 
  mutate(period = case_when(year %in% 2000:2009 ~ "2000-2010", 
                            year %in% 2050:2059 ~ "2050-2060",
                            year %in% 2090:2099 ~ "2090-2100"), 
         period_experiment = paste0(period, "_", experiment)) %>% 
  group_by(period, experiment, period_experiment, Longitude, Latitude) %>% 
  summarise(tmax35 = mean(tmax35, na.rm = T), 
            tmax37 = mean(tmax37, na.rm = T)) %>% 
  data.frame()

cmip6_agg2 <- cmip6_con %>% 
  # remove the breeding seaons without full runs of data through Dec/Jan
  filter(!(breeding_season %in% c("1999/2000", "2009/2010", 
                                  "2049/2050", "2059/2060", 
                                  "2080/2090", "2099/2100"))) %>% 
  mutate(year = substring(breeding_season, 1, 4),
         period = case_when(year %in% 2000:2009 ~ "2000-2010", 
                            year %in% 2050:2059 ~ "2050-2060",
                            year %in% 2090:2099 ~ "2090-2100")) %>% 
  filter(!is.na(period)) %>% 
  mutate(period_experiment = paste0(period, "_", experiment)) %>% 
  group_by(period, experiment, period_experiment, Longitude, Latitude) %>% 
  summarise(consecutive_tmax35 = mean(max_streak_35, na.rm = T), 
            consecutive_tmax37 = mean(max_streak_37, na.rm = T)) %>% 
  data.frame()


#-------------------------------------

 # Restrict the spatial data to southern Africa and interpolate to a higher resolution
   # First we deal with the total number of days over a given threshold
 cmip6_highres_ndays <- list() 
 exp <- unique(cmip6_agg$period_experiment)
 
 for (k in exp) {

  # first run for days > 35 deg
 sub_tmax35 <- filter(cmip6_agg, period_experiment == k) %>% 
   dplyr::select(Longitude, Latitude, tmax35)
 points <- vect(sub_tmax35, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
 extent_df <- ext(points)
 r_temp <- rast(extent_df, resolution = 1) 
 r <- rasterize(points, r_temp, field = "tmax35", fun = mean)

 factor <- 10 #  Define the resolution increase factor
 res <- res(r) / factor # Calculate new resolution
 higher_res_r <- rast(ext(r), res = res, crs = crs(r))  
 higher_res_r <- resample(r, higher_res_r, method = "bilinear")
 masked_r <- mask(higher_res_r, southern_africa_agg)
 masked_r_df_tmax35 <- as.data.frame(masked_r, xy = TRUE, na.rm = TRUE) %>% 
   mutate(period_experiment = k) %>% 
   rename(Longitude = x, Latitude = y, tmax35 = mean)
 
 # now run for days > 37 deg
 sub_tmax37 <- filter(cmip6_agg, period_experiment == k) %>% 
   dplyr::select(Longitude, Latitude, tmax37)
 
 points <- vect(sub_tmax37, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
 extent_df <- ext(points)
 r_temp <- rast(extent_df, resolution = 1)  # 1 degree resolution
 r <- rasterize(points, r_temp, field = "tmax37", fun = mean)
 res <- res(r) / factor # Calculate new resolution
 higher_res_r <- rast(ext(r), res = res, crs = crs(r))  
 higher_res_r <- resample(r, higher_res_r, method = "bilinear")
 masked_r <- mask(higher_res_r, southern_africa_agg)
 masked_r_df_tmax37 <- as.data.frame(masked_r, xy = TRUE, na.rm = TRUE) %>% 
   mutate(period_experiment = k) %>% 
   rename(Longitude = x, Latitude = y, tmax37 = mean)
 
 # join the two temp thresholds together 
 masked_r_df <- left_join(masked_r_df_tmax35, masked_r_df_tmax37)
 
 # store in the list
 cmip6_highres_ndays[[which(k == exp)]] <- masked_r_df
 
}
 
 # bind the list into a single data frame and reconstitute the period and experiment info
 cmip6_highres_ndays <- do.call(rbind, cmip6_highres_ndays) %>% 
  separate(period_experiment, into = c("period", "experiment"), sep = "_") 
 
 cmip6_highres_ndays <- cmip6_highres_ndays %>% 
  mutate(experiment = case_when(experiment == "historical" ~ "historical",
                                experiment == "ssp245" ~ "SSP2-4.5", 
                                experiment == "ssp585" ~ "SSP5-8.5"))

 # Now repeat the procedure for the number of consecutive days over the thresholds. 
 cmip6_highres_consecutivedays <- list() 
 exp <- unique(cmip6_agg2$period_experiment) # should be the same as above 
 
 for (k in exp) {
   
   # first run for days > 35 deg
   sub_tmax35 <- filter(cmip6_agg2, period_experiment == k) %>% 
     dplyr::select(Longitude, Latitude, consecutive_tmax35)
   points <- vect(sub_tmax35, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
   extent_df <- ext(points)
   r_temp <- rast(extent_df, resolution = 1) 
   r <- rasterize(points, r_temp, field = "consecutive_tmax35", fun = mean)
   
   factor <- 10 #  Define the resolution increase factor
   res <- res(r) / factor # Calculate new resolution
   higher_res_r <- rast(ext(r), res = res, crs = crs(r))  
   higher_res_r <- resample(r, higher_res_r, method = "bilinear")
   masked_r <- mask(higher_res_r, southern_africa_agg)
   masked_r_df_tmax35 <- as.data.frame(masked_r, xy = TRUE, na.rm = TRUE) %>% 
     mutate(period_experiment = k) %>% 
     rename(Longitude = x, Latitude = y, consecutive_tmax35 = mean)
   
   # now run for days > 37 deg
   sub_tmax37 <- filter(cmip6_agg2, period_experiment == k) %>% 
     dplyr::select(Longitude, Latitude, consecutive_tmax37)
   
   points <- vect(sub_tmax37, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
   extent_df <- ext(points)
   r_temp <- rast(extent_df, resolution = 1)  # 1 degree resolution
   r <- rasterize(points, r_temp, field = "consecutive_tmax37", fun = mean)
   res <- res(r) / factor # Calculate new resolution
   higher_res_r <- rast(ext(r), res = res, crs = crs(r))  
   higher_res_r <- resample(r, higher_res_r, method = "bilinear")
   masked_r <- mask(higher_res_r, southern_africa_agg)
   masked_r_df_tmax37 <- as.data.frame(masked_r, xy = TRUE, na.rm = TRUE) %>% 
     mutate(period_experiment = k) %>% 
     rename(Longitude = x, Latitude = y, consecutive_tmax37 = mean)
   
   # join the two temp thresholds together 
   masked_r_df <- left_join(masked_r_df_tmax35, masked_r_df_tmax37)
   
   # store in the list
   cmip6_highres_consecutivedays[[which(k == exp)]] <- masked_r_df
   
 }
 
 # bind the list into a single data frame and reconstitute the period and experiment info
 cmip6_highres_consecutivedays <- do.call(rbind, cmip6_highres_consecutivedays) %>% 
   separate(period_experiment, into = c("period", "experiment"), sep = "_") 
 
 cmip6_highres_consecutivedays <- cmip6_highres_consecutivedays %>% 
   mutate(experiment = case_when(experiment == "historical" ~ "historical",
                                 experiment == "ssp245" ~ "SSP2-4.5", 
                                 experiment == "ssp585" ~ "SSP5-8.5"))
 
 
#----------------------
 
# Plot the data. 
 # (Need to ensure a big plotting window or R/RStudio may crash)
 
#-----------------------

 # First I will plot all of the future scenarios to place in the supplementary, before then focusing on a medium level SSP to plot for the main text 
 
# set up the colour palette 
brewer_colors <- c("white", brewer.pal(9, "YlOrRd"))  # Reds palette with 9 colors

# Set the colour or the label text for the different SSP "pathways". 
ssp_cols <-  c("black", "goldenrod", "darkred")

# Put the coordinate for the reserve into a data frame for easier handling
krr_coord <- cmip6_highres_ndays %>% 
  dplyr::select(experiment, period) %>% 
 distinct() %>% 
 mutate(Longitude = krr_lon, Latitude = krr_lat)

# n total days over 35 degrees
p1 <- ggplot() +
  geom_tile(data = cmip6_highres_ndays, aes(x = Longitude, y = Latitude, fill = tmax35), colour = NA) + 
  geom_contour(data = cmip6_highres_ndays, aes(x = Longitude, y = Latitude, z = tmax35),
               breaks = seq(0, 150, 25), 
               colour = "white", alpha = 0.5, linewidth = 0.1) +
  geom_sf(data = southern_africa, fill = "white", color = "black",
          alpha = 0, lwd = 0.3) + 
  geom_sf_pattern(
    data = mkat_sf,
    aes(),
    pattern = "stripe",
    fill = 'white',
    alpha = 0,
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_size = 0.1,
    pattern_density = 0.01,  
    pattern_angle = 45) +
  geom_star(data = krr_coord, aes(x = Longitude, y = Latitude), fill = "black", size = 2) +
  annotate("rect", xmin= 12, xmax = 41, ymin = -14, ymax= -15.1,
           color="white", fill="white") + 
  facet_grid2(experiment ~ period, 
              strip = strip_themed(text_y = elem_list_text(colour = ssp_cols))) + 
  scale_fill_gradientn(colors = brewer_colors,  breaks = seq(0, 200, 25)) +
  scale_x_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), abs(x), "°")) +
  scale_y_continuous(labels = function(y) paste0(ifelse(y < 0, "-", ""), abs(y), "°")) + 
  labs(fill = "Days ≥ 35°C") +  
  theme_minimal() +
  theme(panel.grid = element_line(linewidth = 0.3),
        axis.title = element_blank(), 
        axis.text = element_text(colour = "black", size = 7.5),
        strip.text.x.top =  element_text(size = 11.5),
        strip.text.y.right = element_text(size = 10.5),
        legend.title = element_text(size = 10.5, hjust = 0.5), 
        legend.text = element_text(size = 9.5)) + 
  coord_sf(expand = F)
p1   
 
#-----------------

# Repeat for n total days over 37 degrees

p2 <- ggplot() +
  geom_tile(data = cmip6_highres_ndays, aes(x = Longitude, y = Latitude, fill = tmax37), colour = NA) + 
  geom_contour(data = cmip6_highres_ndays, aes(x = Longitude, y = Latitude, z = tmax37),
               breaks = seq(0, 150, 25), 
               colour = "white", alpha = 0.5, linewidth = 0.1) +
  geom_sf(data = southern_africa, fill = "white", color = "black",
          alpha = 0, lwd = 0.3) + 
  geom_sf_pattern(
    data = mkat_sf,
    aes(),
    pattern = "stripe",
    fill = 'white',
    alpha = 0,
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_size = 0.1,
    pattern_density = 0.01,  
    pattern_angle = 45) +
  geom_star(data = krr_coord, aes(x = Longitude, y = Latitude), fill = "black", size = 2) +
  annotate("rect", xmin= 12, xmax = 41, ymin = -14, ymax= -15.1,
           color="white", fill="white") + 
  facet_grid2(experiment ~ period, 
              strip = strip_themed(text_y = elem_list_text(colour = ssp_cols))) + 
  scale_fill_gradientn(colors = brewer_colors, breaks = seq(0, 160, 20)) + 
  scale_x_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), abs(x), "°")) +
  scale_y_continuous(labels = function(y) paste0(ifelse(y < 0, "-", ""), abs(y), "°")) + 
  labs(fill = "Days ≥ 37°C") +  
  theme_minimal() +
  theme(panel.grid = element_line(linewidth = 0.3),
        axis.title = element_blank(), 
        axis.text = element_text(colour = "black", size = 7.5),
        strip.text.x.top =  element_text(size = 11.5),
        strip.text.y.right = element_text(size = 10.5),
        legend.title = element_text(size = 10.5, hjust = 0.5), 
        legend.text = element_text(size = 9.5)) + 
  coord_sf(expand = F)
p2


#-----------------
 
# n consecutive days over 35 degrees

p3 <- ggplot() +
  geom_tile(data = cmip6_highres_consecutivedays, 
            aes(x = Longitude, y = Latitude, fill = consecutive_tmax35), colour = NA) + 
  geom_contour(data = cmip6_highres_consecutivedays, aes(x = Longitude, y = Latitude, 
                                                         z = consecutive_tmax35),
               breaks = seq(0, 100, 20), 
               colour = "white", alpha = 0.5, linewidth = 0.1) +
  geom_sf(data = southern_africa, fill = "white", color = "black",
          alpha = 0, lwd = 0.3) + 
  geom_sf_pattern(
    data = mkat_sf,
    aes(),
    pattern = "stripe",
    fill = 'white',
    alpha = 0,
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_size = 0.1,
    pattern_density = 0.01,  
    pattern_angle = 45) +
  geom_star(data = krr_coord, aes(x = Longitude, y = Latitude), fill = "black", size = 2) +
  annotate("rect", xmin= 12, xmax = 41, ymin = -14, ymax= -15.1,
           color="white", fill="white") + 
  facet_grid2(experiment ~ period, 
              strip = strip_themed(text_y = elem_list_text(colour = ssp_cols))) + 
  scale_fill_gradientn(colors = brewer_colors,  breaks = seq(0, 100, 20)) +
  #rescaler = ~ scales::rescale_mid(.x, mid = ))  + 
  scale_x_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), abs(x), "°")) +
  scale_y_continuous(labels = function(y) paste0(ifelse(y < 0, "-", ""), abs(y), "°")) + 
  labs(fill = "Consecutive\ndays ≥ 35°C") +  
  theme_minimal() +
  theme(panel.grid = element_line(linewidth = 0.3),
        axis.title = element_blank(), 
        axis.text = element_text(colour = "black", size = 7.5),
        strip.text.x.top =  element_text(size = 11.5),
        strip.text.y.right = element_text(size = 10.5),
        legend.title = element_text(size = 10.5, hjust = 0.5), 
        legend.text = element_text(size = 9.5)) + 
  coord_sf(expand = F)
p3   

#------------------

# n consecutive days over 37 degrees

p4 <- ggplot() +
  geom_tile(data = cmip6_highres_consecutivedays, 
            aes(x = Longitude, y = Latitude, fill = consecutive_tmax37), colour = NA) + 
  geom_contour(data = cmip6_highres_consecutivedays, aes(x = Longitude, y = Latitude, 
                                                         z = consecutive_tmax37),
               breaks = seq(0, 60, 10), 
               colour = "white", alpha = 0.5, linewidth = 0.1) +
  geom_sf(data = southern_africa, fill = "white", color = "black",
          alpha = 0, lwd = 0.3) + 
  geom_sf_pattern(
    data = mkat_sf,
    aes(),
    pattern = "stripe",
    fill = 'white',
    alpha = 0,
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_size = 0.1,
    pattern_density = 0.01,  
    pattern_angle = 45) +
  geom_star(data = krr_coord, aes(x = Longitude, y = Latitude), fill = "black", size = 2) +
  annotate("rect", xmin= 12, xmax = 41, ymin = -14, ymax= -15.1,
           color="white", fill="white") + 
  facet_grid2(experiment ~ period, 
              strip = strip_themed(text_y = elem_list_text(colour = ssp_cols))) + 
  scale_fill_gradientn(colors = brewer_colors,  breaks = seq(0, 60, 10)) +
  scale_x_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), abs(x), "°")) +
  scale_y_continuous(labels = function(y) paste0(ifelse(y < 0, "-", ""), abs(y), "°")) + 
  labs(fill = "Consecutive\ndays ≥ 37°C") +  
  theme_minimal() +
  theme(panel.grid = element_line(linewidth = 0.3),
        axis.title = element_blank(), 
        axis.text = element_text(colour = "black", size = 7.5),
        strip.text.x.top =  element_text(size = 11.5),
        strip.text.y.right = element_text(size = 10.5),
        legend.title = element_text(size = 10.5, hjust = 0.5), 
        legend.text = element_text(size = 9.5)) + 
  coord_sf(expand = F)
p4   


#-----------------

# Generate the final plot for the main text
 
# For the main plot I will use SSP2-4.5 as a medium-term scenario 
# Plot all the responses for the temperatures > 37 degrees
  cmip6_main <- filter(cmip6_highres_ndays, experiment %in% c("historical", "SSP2-4.5"))

# Put the coordinate for the reserve into a data frame for easier handling
krr_coord2 <- cmip6_main %>% 
  dplyr::select(experiment, period) %>% 
  distinct() %>% 
  mutate(Longitude = krr_lon, Latitude = krr_lat)

p5 <- ggplot() +
  geom_tile(data = cmip6_main, aes(x = Longitude, y = Latitude, fill = tmax37), colour = NA) + 
  geom_contour(data = cmip6_main, aes(x = Longitude, y = Latitude, z = tmax37),
               breaks = seq(0, 150, 25), 
               colour = "white", alpha = 0.2, linewidth = 0.1) +
  geom_sf(data = southern_africa, fill = "white", color = "black",
          alpha = 0, lwd = 0.3) + 
  geom_sf_pattern(
    data = mkat_sf,
    aes(),
    pattern = "stripe",
    fill = 'white',
    alpha = 0,
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_size = 0.1,
    pattern_density = 0.01,  
    pattern_angle = 45) +
  geom_star(data = krr_coord2, aes(x = Longitude, y = Latitude), fill = "black", size = 2) +
  annotate("rect", xmin= 12, xmax = 41, ymin = -14, ymax= -15.1,
           color="white", fill="white") + 
  facet_wrap(~period) + 
  scale_fill_gradientn(colors = brewer_colors,  breaks = seq(0, 100, 20))  + 
  scale_x_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), abs(x), "°")) +
  scale_y_continuous(labels = function(y) paste0(ifelse(y < 0, "-", ""), abs(y), "°")) + 
  labs(fill = "Days ≥ 37°C", tag = "A") +  
  theme_minimal() +
  theme(panel.grid = element_line(linewidth = 0.3),
        axis.title = element_blank(), 
        axis.text = element_text(colour = "black", size = 7.5),
        strip.text =  element_text(size = 11),
        legend.title = element_text(size = 10.5, hjust = 0.5), 
        legend.text = element_text(size = 9.5), 
        plot.tag = element_text(size = 16)) + 
  coord_sf(expand = F)
p5

# Repeat for the number of consecutive days > 37
  cmip6_main2 <- filter(cmip6_highres_consecutivedays, experiment %in% c("historical", "SSP2-4.5"))

p6 <- ggplot() +
  geom_tile(data = cmip6_main2, 
            aes(x = Longitude, y = Latitude, fill = consecutive_tmax37), 
            colour = NA) + 
  geom_contour(data = cmip6_main2, 
               aes(x = Longitude, y = Latitude, z = consecutive_tmax37),
               breaks = seq(0, 40, 10), 
               colour = "white", alpha = 0.2, linewidth = 0.1) +
  geom_sf(data = southern_africa, fill = "white", color = "black",
          alpha = 0, lwd = 0.3) + 
  geom_sf_pattern(
    data = mkat_sf,
    aes(),
    pattern = "stripe",
    fill = 'white',
    alpha = 0,
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_size = 0.1,
    pattern_density = 0.01,  
    pattern_angle = 45) +
  geom_star(data = krr_coord2, aes(x = Longitude, y = Latitude), fill = "black", size = 2) +
  annotate("rect", xmin= 12, xmax = 41, ymin = -14, ymax= -15.1,
           color="white", fill="white") + 
  facet_wrap(~period) + 
  scale_fill_gradientn(colors = brewer_colors,  breaks = seq(0, 40, 10))  + 
  scale_x_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), abs(x), "°")) +
  scale_y_continuous(labels = function(y) paste0(ifelse(y < 0, "-", ""), abs(y), "°")) + 
  labs(fill = "Consecutive\ndays ≥ 37°C", tag = "B") +  
  theme_minimal() +
  theme(panel.grid = element_line(linewidth = 0.3),
        axis.title = element_blank(), 
        axis.text = element_text(colour = "black", size = 7.5),
        strip.text =  element_text(size = 11),
        legend.title = element_text(size = 10.5, hjust = 0.5), 
        legend.text = element_text(size = 9.5), 
        plot.tag = element_text(size = 16)) + 
  coord_sf(expand = F)
p6

# Export the two plots as a pdf for the Manuscript. 
p5/p6

################### END OF SCRIPT   ##############################