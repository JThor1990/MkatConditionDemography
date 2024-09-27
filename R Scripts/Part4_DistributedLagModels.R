###########################################################

# LINKING CLIMATE VARIABILITY TO DEMOGRAPHY IN ARID ENVIRONMENTS: A MECHANISTIC STUDY IN KALAHARI MEERKATS 

# AUTHORS: Thorley, Duncan, ... and Clutton-Brock (2024) 

# PART 4: USING DISTRIBUTED LAG MODELS TO MODEL THE LAGGED INTERACTION BETWEEN TEMPERATURE AND NDVI

###########################################################

# Load required packages
lapply(c("tidyverse", "ggplot2", "lubridate", "patchwork", "mgcv", "gratia", "ggh4x", "sf"),  FUN = library, character.only = TRUE)

# Set working directory 
setwd("INSERT FILE PATH HERE")

# Load the data sets: 
# temperature 
temp <- read.csv("Data\\MeerkatDailyTemp.csv", header = TRUE) %>% 
     mutate(date = as.Date(date), 
            tempmax = if_else(tempmax_noaa < 0, as.numeric(NA), tempmax_noaa ), 
            tempmin = if_else(tempmin_noaa  < -15, as.numeric(NA), tempmin_noaa),
            temp_mid = tempmax - (tempmax - tempmin)/2)

# NDVI 
ndvi_masked_agg <- read.csv("Data\\MeerkatNDVImasked.csv", header = T) %>% 
  mutate(date = as.Date(date)) %>% 
  rename(mean_ndvi = mean_ndvi_masked)
  
# Body condition (average deviation in adult body mass from their age-predicted mass)
  # The weekly average body condition of adults 
  weekly_mass_resid <- read.csv("Data\\MeerkatWeeklyConditionAllAdults.csv", header = TRUE) %>% 
    mutate(date = as.Date(strptime(date, format = "%d/%m/%Y")))
  # note that the data here refers to the end of the week
  
# Generate a plot theme
climate_theme <- theme_bw() + 
  theme(axis.title = element_text(size = 14, colour = "black"), 
        axis.text = element_text(size = 12, colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#--------------------------------------
# Create lagged variables for modelling 

lagard <- function(x,n.lag=16) {
  n <- length(x); X <- matrix(NA,n,n.lag)
  for (i in 1:n.lag) X[i:n,i] <- x[i:n-i+1]
  X
}

# For this we need the end of each week
#weekly_mass_resid$date <- as.Date(paste0(weekly_mass_resid$year, "-01-01")) + 
#  weekly_mass_resid$week*7 - 1

# Weekly temperatures 
weekly_tmax <- temp %>% 
  mutate(week = week(date), 
         week = ifelse(week == 53, 52, week)) %>% 
  group_by(week, year) %>% 
  summarise(meantemp = mean(tempmax_noaa, na.rm = TRUE), 
            meanmintemp = mean(tempmin_noaa , na.rm = TRUE), 
            meanavgtemp = mean(tempmid_noaa, na.rm = TRUE)) %>% 
  arrange(year, week) %>% 
  data.frame()

# and add the date again for merging 
weekly_tmax$date <- as.Date(paste0(weekly_tmax$year, "-01-01")) + 
  weekly_tmax$week*7 - 1  # get the end of each week

# Weekly NDVI 
weekly_ndvi <- data.frame(date = weekly_tmax$date)
max_date <- max(ndvi_masked_agg$date) 
weekly_ndvi$mean_ndvi <- sapply(weekly_ndvi$date, FUN = function(x) { 
  val <- ndvi_masked_agg$mean_ndvi[max(which(ndvi_masked_agg$date < x))]
  return(val)
})
weekly_ndvi <- filter(weekly_ndvi, !is.na(mean_ndvi), date < max_date)
tail(weekly_ndvi)

# Now set up the lagged temperature and ndvi matrices for the modelling
lagperiod <- c(min(weekly_mass_resid$date) - 16*7, max(weekly_mass_resid$date))

# lagged temperature data 
tempsub <- filter(weekly_tmax, between(date, lagperiod[1], lagperiod[2]))
templags <- data.frame(lagard(tempsub$meantemp, n.lag= 16))  
names(templags) <- paste0("tmax.", 1:16)
templags$date <- tempsub$date
templags <- templags[,c(17, 1:16)]

# lagged ndvi
ndvisub <- filter(weekly_ndvi, between(date, lagperiod[1], lagperiod[2]))
ndvilags <- data.frame(lagard(ndvisub$mean_ndvi, n.lag=16))  
names(ndvilags) <- paste0("ndvi.", 1:16)
ndvilags$date <- ndvisub$date
ndvilags <- ndvilags[,c(17, 1:16)]

# Now set up the final data set 
df_model <- left_join(weekly_mass_resid, templags) %>% 
  left_join(ndvilags) %>% 
  distinct()

df_model <- na.omit(df_model) 
df_model$yday <- yday(df_model$date)
range(df_model$date) 

#---------------------------------------------
# Run the full model across the complete time series 

tmax_weekly <- as.matrix(dplyr::select(df_model, tmax.1:tmax.16))
ndvi_weekly <- as.matrix(dplyr::select(df_model, ndvi.1:ndvi.16))
lags_weekly <-  matrix(rep(1:16, times = nrow(df_model)), ncol = 16, byrow = TRUE)
#df_model$year <- factor(df_model$year)

# Fit the model
dlm_full <- gam(weeklyresid ~ 
                  s(yday, k = 20, bs = "cc") +
                  s(GroupSize, k = 3, bs = "cr") +
                  te(ndvi_weekly, tmax_weekly, lags_weekly,
                     bs = "cr", k = c(5,5,5)),
                data = df_model,  
                method = "REML")
summary(dlm_full)  

layout(matrix(1:4, ncol = 2))
gam.check(dlm_full)
k.check(dlm_full)
layout(1)
mgcv::plot.gam(dlm_full, select = 1, main = NA)
mgcv::plot.gam(dlm_full, select = 2, main = NA)
mgcv::plot.gam(dlm_full, select = 3, too.far = 0.04, main = NA)

# check the concurvity: 
concurvity(dlm_full)
 # The values are high for the yday and the tensor smooth, which is expected given the highly seasonal nature of temperature and NDVI variation. 
 # For this reason, it is difficult to separate the effects of season from those of the specific environmental variables. However, I still think it is useful to include both terms, in particular because the day of the year smooth makes clear that the October increase occurs even in years with low early season rainfall.
 # Our approach will therefore be very transparent when describing the patterns. 
 # That said, the exclusion of the yday term to isolate the specific environmental effects produces similar inferences (run blocked out code to see this)

#dlm_full2 <- gam(weeklyresid ~ 
#                  s(GroupSize, k = 3, bs = "cr") +
#                  te(ndvi_weekly, tmax_weekly, lags_weekly,
#                    bs = "cr", k = c(5,5,5)),
#                data = df_model,  
                #method = "REML")
#summary(dlm_full2)  # quick to run
#mgcv::plot.gam(dlm_full2, select = 2, too.far = 0.04, main = NA)

# Other measures of fit for the full model: RMSE & MAE 
model_fit <- predict(dlm_full, type = "response") %>% 
  tibble(method = "GAM", pred = .) %>% 
  bind_cols(obs = df_model$weeklyresid)

model_fit %>%
  summarize(MAD = mean(abs(obs - pred), na.rm = TRUE),
            RMSE = sqrt(mean((obs - pred) ^ 2)))

# Predict the model on the raw data
# We can do this directly on the time series, which gives data weekly
df_model$fitted <- fitted(dlm_full)

# But to provide a time series without breaks across the panels where at either end of the year, it's better to set up the prediction for each day (with confidence intervals)
pdat <- data.frame(dplyr::select(df_model, date, year, yday, GroupSize))

pdat$lags_weekly <- lags_weekly  
pdat$tmax_weekly <- tmax_weekly  
pdat$ndvi_weekly <- ndvi_weekly  

### new bit ----->>>>>>>

# Get all the dates and back fill the gaps for the week 
#full_dates <- data.frame(date = seq(min(df_model$date) - 6, max(df_model$date), by = "day"))

#pdat <- full_dates %>%
#  left_join(pdat, by = "date") %>% 
#  fill(GroupSize, starts_with("lags"),  starts_with("tmax"),  starts_with("ndvi"), #.direction = "up") %>% 
#  mutate(year = year(date), yday = yday(date))

#tmax_weekly <- as.matrix(dplyr::select(df_model, tmax.1:tmax.16)#)
#ndvi_weekly <- as.matrix(dplyr::select(df_model, ndvi.1:ndvi.16)#)
#lags_weekly <-  matrix(rep(1:16, times = nrow(df_model)), ncol = 16, byrow = TRUE)#

pred <- predict(dlm_full, newdata = pdat, se.fit = TRUE)
crit <- qt(0.975, df = df.residual(dlm_full)) # ~95% interval critical
pdat <- transform(pdat, fitted = pred$fit, 
                  se = pred$se.fit)
pdat <- transform(pdat,
                  upper = fitted + (crit * se),
                  lower = fitted - (crit * se))

# get the days of the first month of each year
datelabs <- data.frame(date = seq.Date(as.Date("2010-01-01"), 
                                       as.Date("2010-12-31"), 1)) %>% 
  mutate(month = month(date), 
         yday = yday(date), 
         mday = mday(date)) %>% 
  filter(mday == 1, month %in% seq(2, 12, 2))

# Lastly, subtract 3 days from yday so that the raw data and predictions fall in the middle of each week. Means the continuity of the plot across panels is slightly improved

# make the plot
p <- ggplot(df_model, aes(x = yday - 3, y = weeklyresid)) +
  geom_point(data = mass_to_plot, aes(x = yday, y = mass_resid, 
                                      colour = mass_resid), 
            alpha = 0.06, size = 0.5) + 
  geom_hline(yintercept = 0, col = "darkblue", linetype = 2) +
  geom_line(linewidth = 1.01, col = "black") + 
  geom_ribbon(data = pdat,
              aes(ymin = lower, ymax = upper, y = fitted), fill = "red", 
              colour = NA, alpha = 0.3) +
  geom_line(data = pdat, 
            aes(y = fitted), col = "red", linewidth = 1.01, alpha = 0.7) + 
  facet_wrap(~year, ncol = 5) +
  theme_bw() + 
  theme(axis.title = element_text(size = 13), 
        axis.text = element_text(size = 10, colour = "black"), 
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(size = 10.5),
        legend.position = "none", 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(-0.01, "inch")) +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks = datelabs$yday, labels = datelabs$month) +
  scale_y_continuous(breaks = seq(-150, 150, 50), 
                     labels = seq(-150, 150, 50), 
                     limits = c(-175, 175)) +
  labs(x = NULL, y = "Residual body mass (g)") + 
  scale_colour_viridis_c()

dat_text <- data.frame(
  yday = c(100, 150),
  weeklyresid = c(75, 25),
  label = c("Data", "Predicted"), 
  year = c("2002", "2002"))

p <- p + geom_text(
  data    = dat_text,
  mapping = aes(label = label), 
  col = c("black", "red")) 

#---------------------------------------------
# Plot all of the smoothers from the full model 

# (1) Day of year effect: 
plot_g_a <- draw(dlm_full, select = 1) +
  geom_hline(yintercept = 0, linetype = 2, colour = "darkblue") +
  labs(title = NULL, caption = NULL,
       x = "Day of the year", y = "Body condition partial effect (g)", tag = "A") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        plot.tag = element_text(size = 15)) + 
  scale_y_continuous(breaks = seq(-80, 50, 10),  
                     labels = c("-80", "", "-60", "", "-40", 
                                "", "-20", "", "0", "", "20", "", "40", ""), 
                     limits = c(-85, 55))

# (2) Average group size effect: 
plot_g_b <- draw(dlm_full, select = 2) +
  geom_hline(yintercept = 0, linetype = 2, colour = "darkblue") +
  labs(title = NULL, caption = NULL,
       x = "Average group size", y = "Body condition partial effect (g)", 
       tag = "B") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        panel.grid = element_blank(), 
        plot.tag = element_text(size = 15)) + 
  scale_y_continuous(breaks = seq(-20, 35, 10), limits = c(-20, 35))

# plot them together
cowplot::plot_grid(plot_g_a, plot_g_b)

# (3) Lagged interaction between temperature and NDVI

# create the lag names
lagnames <- sapply(1:16, function(x) paste("Lag =", x))
names(lagnames) <- as.character(1:16)

# define jet colormap
jet.colors <- colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan", 
                                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))

# The problem with this is that the dist function doesn't work on anything above 2-D smooths. So we could need and plot this manually. 
sm <- smooth_estimates(dlm_full, n = 300, smooth = 3, n_3d = 16) 
# Only plot the prediction surface that falls within the 95% ellipse
d <- dplyr::select(df_model, tmax.1 , ndvi.1) %>% 
  rename(x = tmax.1, y = ndvi.1)

# density function
kd <- ks::kde(d, compute.cont=TRUE, h=0.2)

# extract the 95% polygon (5% probability mass) 
get_contour <- function(kd_out=kd, prob="5%") {
  contour_95 <- with(kd_out, 
                     contourLines(x=eval.points[[1]], 
                                  y=eval.points[[2]],
                                  z=estimate, 
                                  levels=cont[prob])[[1]])
  as_tibble(contour_95) %>% 
    mutate(prob = prob)
}
dat_out <- map_dfr(c("5%"), ~ get_contour(kd, .)) %>% 
  group_by(prob) %>% 
  mutate(n_val = 1:n()) %>% 
  ungroup() %>% 
  data.frame()

# check it's done a good job
ggplot(data = dat_out, aes(x, y)) +
  geom_point(data = d, alpha = 0.1) +
  geom_path(aes(group = prob), colour = "red", linewidth = 1) 

# Create a matrix of coordinates for the polygon
polygon_coords <- matrix(c(dat_out$x, dat_out$y), ncol = 2, byrow = FALSE)
# Close the polygon (add the first point at the end to close the loop)
polygon_coords <- rbind(polygon_coords, polygon_coords[1, ])
# Create an sf polygon object using st_polygon
polygon_sf <- st_polygon(list(polygon_coords))
# Create a matrix of coordinates for the points
points_coords <- matrix(c(sm$tmax_weekly, sm$ndvi_weekly), ncol = 2, byrow = FALSE)
# Create an sf geometry list of points using st_sfc
points_sf <- st_as_sf(data.frame(points_coords), coords = c(1,2))
# Which points in the prediction surface fall within the within the 95% ellipse
sm$inside <- st_intersects(points_sf, polygon_sf, sparse = T) %>% 
  as.numeric() %>% 
  replace_na(0)
intersection_pts <- sm$inside # store for later
# filter out only those points falling within the surface
sm <- filter(sm, inside == 1)

# plot at a smaller range of lags
plot_te <- ggplot(filter(sm, lags_weekly %in% c(1,2,4,8,12,16)), 
                   aes(x = tmax_weekly, y = ndvi_weekly)) + 
  geom_raster(aes(fill = .estimate)) + 
  geom_contour(aes(z = .estimate), breaks = seq(-80, 80, 5), 
               colour = "black", alpha = 0.2) +
  geom_contour(aes(z = .estimate), breaks = 0, 
               colour = "black", alpha = 0.5, linetype = 2) +
  labs(title = NULL, caption = NULL,
       y = "Weekly NDVI", x = "Mean weekly\n maximum temperature (°C)", 
       fill = "Body\ncondition\npartial\neffect\n(g)") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 9, colour = "black"), 
        strip.text = element_text(size = 10, colour = "black"), 
        axis.title = element_text(size = 10.5, colour = "black"), 
        legend.title = element_text(size = 10, hjust = 0.5),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "lightblue"), 
        panel.grid = element_blank()) + 
  facet_wrap(~lags_weekly, ncol = 3, labeller = as_labeller(lagnames)) +
  coord_cartesian(expand = FALSE) +
  scale_fill_gradientn(colors = jet.colors(10),
                       rescaler = ~ scales::rescale_mid(.x, mid = -12))  + 
  xlim(c(16, 43)) + 
  ylim(c(0.135, 0.32))

#---------------------------------------------
# Refit the lagged model on the different classes of individual 

  # Load in the data sets for sub females, sub males, dom females, and dom males
weekly_mass_resid_subfem <- read.csv("Data\\MeerkatWeeklyConditionSubFemales.csv", header = TRUE)  %>% 
  mutate(date = as.Date(strptime(date, format = "%d/%m/%Y")))

weekly_mass_resid_submale <- read.csv("Data\\MeerkatWeeklyConditionSubMales.csv", header = TRUE) %>% 
  mutate(date = as.Date(strptime(date, format = "%d/%m/%Y")))

weekly_mass_resid_domfem  <- read.csv("Data\\MeerkatWeeklyConditionDomFemales.csv", header = TRUE) %>% 
  mutate(date = as.Date(strptime(date, format = "%d/%m/%Y")))

weekly_mass_resid_dommales <- read.csv("Data\\MeerkatWeeklyConditionDomMales.csv", header = TRUE) %>% 
  mutate(date = as.Date(strptime(date, format = "%d/%m/%Y")))

#---------------------------------------------

# Fit the subordinate female model
df_model_subfem <- left_join(weekly_mass_resid_subfem, templags) %>% 
  left_join(ndvilags) %>% 
  distinct()

df_model_subfem <- na.omit(df_model_subfem) 
df_model_subfem$yday <- yday(df_model_subfem$date)

dlm_subfemale <- gam(weeklyresid ~ 
                       s(yday, k = 20, bs = "cc") +
                       s(GroupSize, k = 3, bs = "cr") +
                       te(ndvi_weekly, tmax_weekly, lags_weekly,
                          bs = "cr", k = c(5,5,5)),    
                     data=df_model_subfem,  
                     method = "REML")
summary(dlm_subfemale)  

layout(matrix(1:4, ncol = 2))
gam.check(dlm_subfemale) # misses some of the extreme residuals but its not bad
k.check(dlm_subfemale)
layout(1)
mgcv::plot.gam(dlm_subfemale, select = 1, main = NA)
mgcv::plot.gam(dlm_subfemale, select = 2, main = NA)
mgcv::plot.gam(dlm_subfemale, select = 3, too.far = 0.05, main = NA)

subfemale_model_fit <- predict(dlm_subfemale, type = "response") %>% 
  tibble(method = "GAM", pred = .) %>% 
  bind_cols(obs = df_model_subfem$weeklyresid)

subfemale_model_fit %>%
  summarize(MAD = mean(abs(obs - pred), na.rm = TRUE),
            RMSE = sqrt(mean((obs - pred) ^ 2)))

# Fit the subordinate male model
df_model_submale <- left_join(weekly_mass_resid_submale, templags) %>% 
  left_join(ndvilags) %>% 
  distinct()

df_model_submale <- na.omit(df_model_submale) 
df_model_submale$yday <- yday(df_model_submale$date)

dlm_submale <- gam(weeklyresid ~ 
                     s(yday, k = 20, bs = "cc") +
                     s(GroupSize, k = 3, bs = "cr") +
                     te(ndvi_weekly, tmax_weekly, lags_weekly,
                        bs = "cr", k = c(5,5,5)),    
                   data=df_model_submale,  
                   method = "REML")
summary(dlm_submale)  # quick to run

layout(matrix(1:4, ncol = 2))
gam.check(dlm_submale) # misses some of the extreme residuals but its not bad
k.check(dlm_submale)
layout(1)
mgcv::plot.gam(dlm_submale, select = 1, main = NA)
mgcv::plot.gam(dlm_submale, select = 2, main = NA)
mgcv::plot.gam(dlm_submale, select = 3, too.far = 0.05, main = NA)

submale_model_fit <- predict(dlm_submale, type = "response") %>% 
  tibble(method = "GAM", pred = .) %>% 
  bind_cols(obs = df_model_submale$weeklyresid)

submale_model_fit %>%
  summarize(MAD = mean(abs(obs - pred), na.rm = TRUE),
            RMSE = sqrt(mean((obs - pred) ^ 2)))

# Fit the dominant females model: 
df_model_domfem <- left_join(weekly_mass_resid_domfem, templags) %>% 
  left_join(ndvilags) %>% 
  distinct()

df_model_domfem <- na.omit(df_model_domfem) 
df_model_domfem$yday <- yday(df_model_domfem$date)

tmax_weekly2 <- as.matrix(dplyr::select(df_model_domfem, tmax.1:tmax.16))
ndvi_weekly2 <- as.matrix(dplyr::select(df_model_domfem, ndvi.1:ndvi.16))
lags_weekly2 <-  matrix(rep(1:16, times = nrow(df_model_domfem)), ncol = 16, byrow = TRUE)

dlm_domfemale <- gam(weeklyresid ~ 
                       s(yday, k = 20, bs = "cc") +
                       s(GroupSize, k = 3, bs = "cr") +
                       te(ndvi_weekly2, tmax_weekly2, lags_weekly2,
                          bs = "cr", k = c(5,5,5)),    
                     data=df_model_domfem,  
                     method = "REML")
summary(dlm_domfemale) 

layout(matrix(1:4, ncol = 2))
gam.check(dlm_domfemale) # misses some of the extreme residuals but its not bad
k.check(dlm_domfemale)
layout(1)
mgcv::plot.gam(dlm_domfemale, select = 1, main = NA)
mgcv::plot.gam(dlm_domfemale, select = 2, main = NA)
mgcv::plot.gam(dlm_domfemale, select = 3, too.far = 0.05, main = NA)

domfemale_model_fit <- predict(dlm_domfemale, type = "response") %>% 
  tibble(method = "GAM", pred = .) %>% 
  bind_cols(obs = df_model_domfem$weeklyresid)

domfemale_model_fit %>%
  summarize(MAD = mean(abs(obs - pred), na.rm = TRUE),
            RMSE = sqrt(mean((obs - pred) ^ 2)))

# Fit the dominant male model: 
df_model_dommale <- left_join(weekly_mass_resid_dommales, templags) %>% 
  left_join(ndvilags) %>% 
  distinct()

df_model_dommale <- na.omit(df_model_dommale) 
df_model_dommale$yday <- yday(df_model_dommale$date)

tmax_weekly3 <- as.matrix(dplyr::select(df_model_dommale, tmax.1:tmax.16))
ndvi_weekly3 <- as.matrix(dplyr::select(df_model_dommale, ndvi.1:ndvi.16))
lags_weekly3 <-  matrix(rep(1:16, times = nrow(df_model_dommale)), ncol = 16, byrow = TRUE)

dlm_dommale <- gam(weeklyresid ~ 
                     s(yday, k = 20, bs = "cc") +
                     s(GroupSize, k = 3, bs = "cr") +
                     te(ndvi_weekly3, tmax_weekly3, lags_weekly3,
                        bs = "cr", k = c(5,5,5)),    
                   data=df_model_dommale,  
                   method = "REML")
summary(dlm_dommale)  # quick to run

layout(matrix(1:4, ncol = 2))
gam.check(dlm_dommale) # misses some of the extreme residuals but its not bad
k.check(dlm_dommale)
layout(1)
mgcv::plot.gam(dlm_dommale, select = 1, main = NA)
mgcv::plot.gam(dlm_dommale, select = 2, main = NA)
mgcv::plot.gam(dlm_dommale, select = 3, too.far = 0.05, main = NA)

dommale_model_fit <- predict(dlm_dommale, type = "response") %>% 
  tibble(method = "GAM", pred = .) %>% 
  bind_cols(obs = df_model_dommale$weeklyresid)

dommale_model_fit %>%
  summarize(MAD = mean(abs(obs - pred), na.rm = TRUE),
            RMSE = sqrt(mean((obs - pred) ^ 2)))

# Predict all the model outputs 
df_model_subfem$fitted <- fitted(dlm_subfemale)
df_model_submale$fitted <- fitted(dlm_submale)
df_model_domfem$fitted <- fitted(dlm_domfemale)
df_model_dommale$fitted <- fitted(dlm_dommale)

# Predict the historic body condition time series for each class of individual

# subordinate females
pdat_subfem <- data.frame(dplyr::select(df_model_subfem, date, year, yday, GroupSize))
pdat_subfem$lags_weekly <- lags_weekly  
pdat_subfem$tmax_weekly <- tmax_weekly  
pdat_subfem$ndvi_weekly <- ndvi_weekly  

pred_subfem <- predict(dlm_subfemale, newdata = pdat_subfem, se.fit = TRUE)
crit_subfem <- qt(0.975, df = df.residual(dlm_subfemale)) # ~95% interval critical
pdat_subfem <- transform(pdat_subfem, fitted = pred_subfem$fit, 
                         se = pred_subfem$se.fit)
pdat_subfem <- transform(pdat_subfem,
                         upper = fitted + (crit * se),
                         lower = fitted - (crit * se)) %>% 
  dplyr::select(date:GroupSize, fitted:lower) %>% 
  mutate(DomStatus = "Subordinate females")

# subordinate males
pdat_submale <- data.frame(dplyr::select(df_model_submale, date, year, yday, GroupSize))
pdat_submale$lags_weekly <- lags_weekly  
pdat_submale$tmax_weekly <- tmax_weekly  
pdat_submale$ndvi_weekly <- ndvi_weekly  

pred_submale <- predict(dlm_submale, newdata = pdat_submale, se.fit = TRUE)
crit_submale <- qt(0.975, df = df.residual(dlm_submale)) # ~95% interval critical
pdat_submale <- transform(pdat_submale, fitted = pred_submale$fit, 
                          se = pred_submale$se.fit)
pdat_submale <- transform(pdat_submale,
                          upper = fitted + (crit * se),
                          lower = fitted - (crit * se)) %>% 
  dplyr::select(date:GroupSize, fitted:lower) %>% 
  mutate(DomStatus = "Subordinate males")

# dominant females
pdat_domfem <- data.frame(dplyr::select(df_model_domfem, date, year, yday, GroupSize))
pdat_domfem$lags_weekly2 <- lags_weekly2
pdat_domfem$tmax_weekly2 <- tmax_weekly2  
pdat_domfem$ndvi_weekly2 <- ndvi_weekly2  

pred_domfem <- predict(dlm_domfemale, newdata = pdat_domfem, se.fit = TRUE)
crit_domfem <- qt(0.975, df = df.residual(dlm_domfemale)) # ~95% interval critical
pdat_domfem <- transform(pdat_domfem, fitted = pred_domfem$fit, 
                         se = pred_domfem$se.fit)
pdat_domfem <- transform(pdat_domfem,
                         upper = fitted + (crit * se),
                         lower = fitted - (crit * se)) %>% 
  dplyr::select(date:GroupSize, fitted:lower) %>% 
  mutate(DomStatus = "Dominant females")

# dominant males
pdat_dommales <- data.frame(dplyr::select(df_model_dommale, date, year, yday, GroupSize))
pdat_dommales$lags_weekly3 <- lags_weekly3
pdat_dommales$tmax_weekly3 <- tmax_weekly3  
pdat_dommales$ndvi_weekly3 <- ndvi_weekly3  

pred_dommales <- predict(dlm_dommale, newdata = pdat_dommales, se.fit = TRUE)
crit_dommales <- qt(0.975, df = df.residual(dlm_dommale)) # ~95% interval critical
pdat_dommales <- transform(pdat_dommales, fitted = pred_dommales$fit, 
                           se = pred_dommales$se.fit)
pdat_dommales <- transform(pdat_dommales,
                           upper = fitted + (crit * se),
                           lower = fitted - (crit * se)) %>% 
  dplyr::select(date:GroupSize, fitted:lower) %>% 
  mutate(DomStatus = "Dominant males")

# Generate the plot as above for the predictions
pdat_allcat <- bind_rows(pdat_subfem, pdat_submale, pdat_domfem, pdat_dommales)

p_cats <- ggplot() + 
  # Population-level trend
  geom_hline(yintercept = 0, col = "darkblue", linetype = 1) +
  geom_line(data = df_model, aes(x = yday, y = weeklyresid), 
            linewidth = 1.01, col = "black") + # pop trend
  # status specific trend
  geom_line(data = pdat_allcat, 
            aes(x = yday, y = fitted, group = DomStatus, colour = DomStatus), 
            linewidth = 0.5, alpha = 1) + 
  facet_wrap(~year, ncol = 8) +
  theme_bw() + 
  theme(axis.title = element_text(size = 13), 
        axis.text = element_text(size = 10, colour = "black"), 
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(size = 11),
        legend.position = "right",
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 13), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.1, "lines")) +
  scale_x_continuous(breaks = datelabs$yday, labels = datelabs$month) +
  scale_y_continuous(breaks = seq(-125, 125, 50), 
                     labels = seq(-125, 125, 50), 
                     limits = c(-125, 100)) +
  labs(x = NULL, y = "Residual body mass (g)", colour = "Model predictions") + 
  scale_fill_manual(values = c("gold", "red", "turquoise", "purple")) + 
  scale_colour_manual(values = c("gold", "red", "turquoise", "purple"))

# Plot the 1-8 week lag for each of these models
sm_subfem <- smooth_estimates(dlm_subfemale, n = 300, smooth = 3, n_3d = 16) %>% 
  mutate(DomStatus = "Subordinate female")
sm_submale <- smooth_estimates(dlm_subfemale, n = 300, smooth = 3, n_3d = 16) %>% 
  mutate(DomStatus = "Subordinate male")
sm_domfem <- smooth_estimates(dlm_domfemale, n = 300, smooth = 3, n_3d = 16) %>% 
  mutate(DomStatus = "Dominant female") %>% 
  rename(ndvi_weekly = ndvi_weekly2, tmax_weekly = tmax_weekly2, lags_weekly = lags_weekly2) 
sm_dommale <- smooth_estimates(dlm_dommale, n = 300, smooth = 3, n_3d = 16) %>% 
  mutate(DomStatus = "Dominant male") %>% 
  rename(ndvi_weekly = ndvi_weekly3, tmax_weekly = tmax_weekly3, lags_weekly = lags_weekly3) 

sm_subfem$inside <- intersection_pts
sm_submale$inside <- intersection_pts
sm_domfem$inside <- intersection_pts
sm_dommale$inside <- intersection_pts
# filter out only those points falling within the surface
sm_subfem <- filter(sm_subfem, inside == 1)
sm_submale <- filter(sm_submale, inside == 1)
sm_domfem <- filter(sm_domfem, inside == 1)
sm_dommale <- filter(sm_dommale, inside == 1)

sm_allcat <- bind_rows(sm_subfem, sm_submale, sm_domfem, sm_dommale)
sm_allcat <- filter(sm_allcat, lags_weekly %in% c(1,2,4,8))

lagnames2 <- lagnames <- c(
  `1` = "Lag = 1",
  `2` = "Lag = 2",
  `4` = "Lag = 4", 
  `8` = "Lag = 8", 
  `Dominant female` = "Dominant females", 
  `Dominant male` = "Dominant males", 
  `Subordinate female` = "Subordinate females", 
  `Subordinate male` = "Subordinate males")

# create colours for the strips of the facets
status_cols <-  adjustcolor(c("gold", "red", "turquoise", "purple"), alpha.f = 0.5)

strip <- strip_themed(background_x = elem_list_rect(fill = rep("lightblue", 3)), 
                      background_y = elem_list_rect(fill =  status_cols))

plot_te2 <- ggplot(sm_allcat, aes(x = tmax_weekly, y = ndvi_weekly)) + 
  geom_raster(aes(fill = .estimate)) +
  geom_contour(aes(z = .estimate), breaks = seq(-80, 80, 5), 
               colour = "black", alpha = 0.2) +
  geom_contour(aes(z = .estimate), breaks = 0, 
               colour = "black", linetype = 2) +
  labs(title = NULL, caption = NULL,
       y = "Weekly NDVI", x = "Mean weekly\n maximum temperature (°C)", 
       fill = "Body\ncondition\npartial\neffect\n(g)") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 9, colour = "black"), 
        strip.text = element_text(size = 10, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.title.align = 0.5, 
        legend.text = element_text(size = 10),
        strip.background = element_rect(fill = "lightblue"), 
        panel.grid = element_blank()) + 
  facet_grid2(DomStatus~lags_weekly, labeller = as_labeller(lagnames2), strip = strip) +
  coord_cartesian(expand = FALSE) +
  scale_fill_gradientn(colors = jet.colors(10),  breaks = seq(-80, 40, 20), 
                       rescaler = ~ scales::rescale_mid(.x, mid = -15))  + 
  xlim(c(16, 43)) + 
  ylim(c(0.135, 0.32))

#-------------------------------

# Lastly, combine the estimates from the different smooths to give the predicted body condition in each week of the year, on average 
# I will base this on using the long-term median NDVI/Temp values in for each month 
pred_df <- df_model %>% 
  # remember ndvi/tempmax are lagged by 1 week so we want to first move them forward 1 week
  mutate(ndvi = lead(ndvi.1), tmax = lead(tmax.1)) %>% 
  dplyr::select(week, year, ndvi, tmax) %>% 
  group_by(week) %>% 
  summarise(ndvi = median(ndvi, na.rm = T), 
            tmax = median(tmax, na.rm = T)) %>% 
  mutate(yday = seq(7, 7*52, 7), 
         GroupSize = mean(df_model$GroupSize)) %>% 
  data.frame()

# plot the raw data for NDVI and temperature to see it 
p_rawndvi <- ggplot(pred_df, aes(x = yday, y = ndvi)) + 
  geom_vline(xintercept = yday(as.Date("2002-12-01")), lwd = 0.3, linetype = 3, colour  = "grey") + 
  geom_vline(xintercept = yday(as.Date("2002-05-01")), lwd = 0.3, linetype = 3, colour  = "grey") + 
  geom_vline(xintercept = yday(as.Date("2002-09-01")), lwd = 0.3, linetype = 3, colour  = "grey") +
  annotate("text", x = 60, y = 0.25, label = "MS", size = 3.5) +
  annotate("text", x = 180, y = 0.25, label = "W", size = 3.5) +
  annotate("text", x = 290, y = 0.25, label = "ES", size = 3.5) +
  geom_smooth(method = gam, formula =   y ~ s(x, bs = "cc", fx = TRUE, k = 6), 
              linewidth = 1, alpha = 0, colour = "forestgreen") + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 9, colour = "black"), 
        axis.title.y = element_text(size = 10.5, colour = "black"),
        plot.margin = margin(t = 0)) +
  scale_x_continuous(breaks = seq(0, 350, 50)) + 
  ylim(c(0.16, 0.253)) + 
  labs(y = "NDVI")

p_rawtemp <- ggplot(pred_df, aes(x = yday, y = tmax)) + 
  geom_vline(xintercept = yday(as.Date("2002-12-01")), lwd = 0.3, linetype = 3, colour  = "grey") + 
  geom_vline(xintercept = yday(as.Date("2002-05-01")), lwd = 0.3, linetype = 3, colour  = "grey") + 
  geom_vline(xintercept = yday(as.Date("2002-09-01")), lwd = 0.3, linetype = 3, colour  = "grey") +
  geom_smooth(method = gam, formula =   y ~ s(x, bs = "cc", fx = TRUE, k = 6), 
              linewidth = 1, alpha = 0, colour = "red") + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 9, colour = "black"), 
        axis.title.y = element_text(size = 10.5, colour = "black"),
        plot.margin = margin(t = 0)) +
  labs(y = "Temperature\n(°C)") + 
  scale_y_continuous(breaks = seq(23, 35, 3)) + 
  scale_x_continuous(breaks = seq(0, 350, 50)) 

p_rawndvi / p_rawtemp

# add another 16 rows to the bottom so that week 1 can also have a lag
pred_df <- bind_rows(pred_df, pred_df[1:16,])

pred_df <- pred_df %>% 
  mutate(ndvi.1 = lag(ndvi, 1), ndvi.2 = lag(ndvi, 2), ndvi.3 = lag(ndvi, 3),
         ndvi.4 = lag(ndvi, 4), ndvi.5 = lag(ndvi, 5), ndvi.6 = lag(ndvi, 6),
         ndvi.7 = lag(ndvi, 7), ndvi.8 = lag(ndvi, 8), ndvi.9 = lag(ndvi, 9),
         ndvi.10 = lag(ndvi, 10), ndvi.11 = lag(ndvi, 11), ndvi.12 = lag(ndvi, 12), 
         ndvi.13 = lag(ndvi, 13), ndvi.14 = lag(ndvi, 14), ndvi.15 = lag(ndvi, 15),
         ndvi.16 = lag(ndvi, 16),
         tmax.1 = lag(tmax, 1), tmax.2 = lag(tmax, 2), tmax.3 = lag(tmax, 3),
         tmax.4 = lag(tmax, 4), tmax.5 = lag(tmax, 5), tmax.6 = lag(tmax, 6),
         tmax.7 = lag(tmax, 7), tmax.8 = lag(tmax, 8), tmax.9 = lag(tmax, 9),
         tmax.10 = lag(tmax, 10), tmax.11 = lag(tmax, 11), tmax.12 = lag(tmax, 12), 
         tmax.13 = lag(tmax, 13), tmax.14 = lag(tmax, 14), tmax.15 = lag(tmax, 15),
         tmax.16 = lag(tmax, 16)) 

pred_df <- pred_df[-(1:16),] %>% 
  arrange(week)

tmax_weekly <- as.matrix(dplyr::select(pred_df, tmax.1:tmax.16))
ndvi_weekly <- as.matrix(dplyr::select(pred_df, ndvi.1:ndvi.16))
lags_weekly <-  matrix(rep(1:16, times = nrow(pred_df)), ncol = 16, byrow = TRUE)

pred_df$tmax_weekly <- tmax_weekly
pred_df$ndvi_weekly <- ndvi_weekly
pred_df$lags_weekly <- lags_weekly

pred_df$pred <- predict(dlm_full, newdata = pred_df)
pred_df$SE <- predict(dlm_full, newdata = pred_df, se.fit = T)$se.fit
pred_df$u95ci <- pred_df$pred  + 1.96*pred_df$SE
pred_df$l95ci <- pred_df$pred  - 1.96*pred_df$SE

# what's the average in the raw data 
raw_average <- df_model %>% 
  group_by(week) %>% 
  summarise(medianresid = median(weeklyresid, na.rm = TRUE), 
            madresid = stats::mad(weeklyresid, na.rm = TRUE), 
            meanresid = mean(weeklyresid, na.rm = TRUE),
            sdresid = sd(weeklyresid, na.rm = TRUE)) %>% 
  mutate(yday = seq(7, 7*52, 7))

#p_average <- ggplot(pred_df, aes(x = yday, y = pred)) + 
#  geom_hline(yintercept = 0, linetype = 2) +
#  geom_point(size = 1, pch  = 22, stroke = 1, fill = "black", alpha = 0.5) + 
#  geom_errorbar(aes(ymin = l95ci, ymax = u95ci), width = 0) + 
#geom_point(data = raw_average, aes(y = medianresid), size = 1.5, 
#           stroke = 1.2, pch = 1, colour = "red") + 
#  climate_theme + 
#  labs(y = "Average body condition (g)", 
#      x = "Day of year") +
#  scale_y_continuous() + 
#  theme(axis.text = element_text(size = 11), 
#        axis.title = element_text(size = 11.5)) + 
#  scale_x_continuous(breaks = seq(0, 300, 50)) + 
#  scale_y_continuous(breaks = seq(-40, 20, 10))
#p_average

# smooth out the blocky signal with a high dimension smoother
p_average <- ggplot(pred_df, aes(x = yday, y = pred)) + 
  geom_hline(yintercept = 0, colour = "lightgray") +
  geom_vline(xintercept = yday(as.Date("2002-12-01")), lwd = 0.3, linetype = 3, colour  = "grey") + 
  geom_vline(xintercept = yday(as.Date("2002-05-01")), lwd = 0.3, linetype = 3, colour  = "grey") + 
  geom_vline(xintercept = yday(as.Date("2002-09-01")), lwd = 0.3, linetype = 3, colour  = "grey") +
  geom_smooth(aes(y = u95ci), 
              method = gam, formula =   y ~ s(x, bs = "tp", fx = TRUE, k = 20), 
              linewidth = 0.5, alpha = 0, linetype = 2, colour = "black") +
  geom_smooth(aes(y = l95ci), 
              method = gam, formula =   y ~ s(x, bs = "tp", fx = TRUE, k = 20), 
              linewidth = 0.5, alpha = 0, linetype = 2, colour = "black") +
  geom_smooth(method = gam, formula =   y ~ s(x, bs = "tp", fx = TRUE, k = 20), 
              linewidth = 0.8, alpha = 0, colour = "black") +
  theme_classic() + 
  labs(y = "Predicted body\ncondition (g)", 
       x = "Day of year") +
  scale_y_continuous() + 
  theme(axis.text = element_text(size = 9, colour = "black"), 
        axis.title = element_text(size = 10.5, colour = "black")) + 
  scale_x_continuous(breaks = seq(0, 350, 50)) + 
  scale_y_continuous(breaks = seq(-40, 30, 10))
p_average

# Perhaps Add NDVI and temperature in a standard year on top of this year?
p_rawndvi <- p_rawndvi + 
  labs(tag = "B") + 
  theme(plot.tag = element_text(size = 15))
p_output <- p_rawndvi / p_rawtemp / p_average
p_output <- p_output + 
  plot_layout(heights =  c(0.25, 0.25, 0.5)) 

# add this below the other plot to give the final figure
# first need to add a bit of space to the one to make them line up more nicely
# and tweak the plot sizes
plot_te <- plot_te + 
  labs(tag = "A") + 
  theme(plot.tag = element_text(size = 15))

plot_te / p_rawndvi / p_rawtemp / p_average + 
  plot_layout(heights = c(0.48, 0.125, 0.125, 0.2))

#########################  END ###################################