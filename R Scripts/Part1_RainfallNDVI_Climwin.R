###########################################################

# LINKING CLIMATE VARIABILITY TO DEMOGRAPHY IN ARID ENVIRONMENTS: A MECHANISTIC STUDY IN KALAHARI MEERKATS 

# AUTHORS: Thorley, Duncan, ... and Clutton-Brock (2024) 

# PART 1: SLIDING WINDOW ANALYSIS OF THE EFFEFCT OF RAINFALL ON NDVI

###########################################################

# set up the working environment 
lapply(c("climwin", "ggplot2", "patchwork", "tidyverse"), FUN = library, character.only = TRUE)

# Set working directory
setwd("INSERT FILE PATH HERE")

# Load in the data sets: 
  # NDVI 
ndvi <- read.csv("Data\\MeerkatNDVI.csv", header = T) %>% 
  mutate(date = as.Date(date))

  # Rainfall (measured on-site)
rain <- read.csv("Data\\MeerkatRainReserve.csv", header = TRUE) %>% 
  mutate(date = as.Date(date)) %>% 
  rename(rain_mm = rain_reserve)

# Set up the linear model with all the climate windows
ndvi_slider <- slidingwin(xvar = list(rain = rain$rain_mm),
                          cdate = rain$date, 
                          bdate = ndvi$date, 
                          baseline = lm(mean_ndvi ~ 1, data = ndvi),
                          range = c(180, 0), 
                          type = "relative",
                          stat = c("sum"),  
                          func = c("lin"),
                          cinterval = "day", 
                          cmissing = "method1") 

# Get the top 10 best windows to first get some idea
ndvi_slider[[1]]$Dataset[1:10,]     # top window was 0-67 days
# Plot the information criteria comparisons for all windows
plotdelta(dataset = ndvi_slider[[1]]$Dataset)
# Look at the model output for the 'best' top window 
summary(ndvi_slider[[1]]$BestModel)
# Identify the median window of those within the 95% confidence set 
medwin(ndvi_slider[[1]]$Dataset) # median window 0.5 to 66 days 
# How many windows fell within the 95% 'confidence set'
dataset <- ndvi_slider[[1]]$Dataset
confidenceset <- dataset[which(cumsum(dataset$ModWeight) <= 0.95),]

# Run randomisations to confirm the signal did not occur by chance due to multiple testing.
job::job( { ndvi_rand <- randwin(repeats = 10, 
                                 xvar =  list(rain = rain$rain_mm),
                                 cdate = rain$date, 
                                 bdate = ndvi$date, 
                                 baseline = lm(mean_ndvi ~ 1, data = ndvi),
                                 range = c(180, 0), 
                                 type = "relative",
                                 stat = c("sum"), 
                                 func = c("lin"),
                                 cinterval = "day", 
                                 cmissing = "method1") })  
summary(ndvi_rand)

# save relevant outputs and plot useful graphics
output_rain <- ndvi_slider[[1]]$Dataset
rain_rand <- ndvi_rand[[1]]
rain_plots <- plotall(dataset = output_rain,
                      datasetrand = rain_rand,
                      bestmodel =   ndvi_slider[[1]]$BestModel, 
                      bestmodeldata = ndvi_slider[[1]]$BestModelData)

# From this, we can cherry pick some of the most informative to plot in paper SI
  # The heat map of AICc
  climwin_plot1 <- plotdelta(dataset = ndvi_slider[[1]]$Dataset) + 
    labs(title = " ") +
    theme_bw() + 
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.position.inside = c(0.8, 0.4), 
          legend.text = element_text(size = 10), 
          panel.grid = element_blank(), 
          legend.background = element_blank()) + 
    geom_point(aes(x = 0, y = 66), pch = 1, size = 3, stroke= 1)

  # The correlation for the best window (67-0) across the raw data 
  climwin_plot2 <- plotbest(dataset = output_rain,
                            bestmodel =  ndvi_slider[[1]]$BestModel, 
                            bestmodeldata =  ndvi_slider[[1]]$BestModelData) + 
    labs(title = " ", x = "Total rainfall (mm) \n67 to 0 days prior", 
         y = "NDVI") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 12),
          panel.grid = element_blank())

  # Lastly, visually show how lagged rainfall tracks NDVI in the raw time series
  # To do so, I will get the rolling cumulative rainfall 
  rain_lag <- rain %>% 
    dplyr::select(date, rain_mm) %>% 
    mutate(rain.lag = lag(rain_mm, n = 1)) %>%
    mutate(rain.67.previous = zoo::rollapply(data = rain.lag, 
                                             width = 67, 
                                             FUN = sum, 
                                             align = "right", 
                                             fill = NA, 
                                             na.rm = T))

  # Plot NDVI and cumulative rainfall time series
  climwin_plot3 <- ggplot(filter(rain_lag, date  > min(ndvi$date)), 
                          aes(x = date, y = rain.67.previous)) +
    geom_path(colour = "blue", alpha = 0.5, linewidth = 0.6) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 11, colour = "black", 
                                     angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 12),
        panel.grid = element_blank()) + 
    labs(y = "Total rainfall (mm) \n67 to 0 days prior", x = "Year") + 
    scale_x_date(date_breaks = "1 year", date_labels = "%Y")  + 
    scale_y_continuous(breaks = seq(0, 240, 40)) + 
    geom_line(data = ndvi,
              aes(y = mean_ndvi*1300), colour = "forestgreen", alpha = 0.5, linewidth = 1.05) + 
    scale_y_continuous(sec.axis = sec_axis(~./1300, name = "Mean NDVI"))

  # Combine all the plots using package 'patchwork'
  climwin_final <- (climwin_plot1 + climwin_plot2)/climwin_plot3

######################## END #############################