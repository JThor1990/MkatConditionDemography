###########################################################

# LINKING CLIMATE VARIABILITY TO DEMOGRAPHY IN ARID ENVIRONMENTS: A MECHANISTIC STUDY IN KALAHARI MEERKATS 

# AUTHORS: Thorley, Duncan, ... and Clutton-Brock (2024) 

# PART 3: MODELLING DAY TO DAY CHANGES IN WEIGHT GAIN AND WEIGHT LOSS

###########################################################

# Load required packages
lapply(c("tidyverse", "ggplot2", "job", "ks", "lubridate", "mgcv", "gamm4",
         "gratia", "patchwork", "performance", "glmmTMB", "sf"), 
       FUN = library, character.only = TRUE)

# set up the working environment 
setwd("C:\\Users\\Jack\\Dropbox\\MeerkatEnvironmentWork\\Climate Paper 2\\Github scripts\\")

# create an inverse logit function for when backtransforming from binomial models
inv.logit <- function(x) { exp(x)/(1+exp(x)) }

# load up the data sets
weekly_dwg <- read.csv("Data\\MeerkatWeeklyPopulationLevelDWG.csv", header = TRUE)
oml <- read.csv("Data\\MeerkatOvernightMassLoss.csv", header = TRUE)
dwg <- read.csv("Data\\MeerkatIndividualDailyWeightGain.csv", header = TRUE)

#-------------------------------------------
# Correlation between average daily weight gain and the change in the average body condition of   all adults per week. 
 
# Average daily body mass gain per day (absolute) ~ condition 
  cor.test(weekly_dwg$mean_absmassgain , weekly_dwg$ConditionChange_t_t1) # 0.32

# Average daily body mass gain per day (%) ~ condition 
  cor.test(weekly_dwg$mean_percentmassgain, weekly_dwg$ConditionChange_t_t1) # 0.37

#--------------------------------------
# Linking daily weight gain to subsequent overnight mass loss 

length(unique(oml$IndividID)) # n = 1682 individuals
  
# Plot the monthly relationship between overnight mass loss and daily weight gain 
  p1 <- ggplot(oml, aes(x = PriorDailyMassGain_abs, 
                              y = OvernightMassLoss_g)) + 
    geom_vline(xintercept = 30, colour= "forestgreen") +
    geom_vline(xintercept = 0, colour= "grey") +
    geom_hline(yintercept = 0, colour= "grey") +
    geom_point(shape = ".", alpha = 0.3, colour = "dodgerblue") + 
    geom_abline(intercept = 0, slope = 1, colour= "black",  alpha = 0.8,
                linetype = 2) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3), colour = "red", 
                fill = "red", se = F) + 
    facet_wrap(~monthabb) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 11, colour = "black"), 
          axis.title = element_text(size = 12, colour = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_rect(fill = "lightblue"), 
          strip.text = element_text(size = 10.5), 
          legend.title = element_blank(), 
          legend.position = c(0.2, 0.85), 
          legend.text = element_text(size = 11), 
          legend.background = element_blank()) + 
    scale_x_continuous(breaks = seq(0, 100, 10),
                       limits = c(0, 100)) + 
    scale_y_continuous(limits = c(0, 100)) + 
    labs(y = "Overnight body mass loss to day 2 (g)",
         x = "Daily body mass gain on day 1 (g)")
  # the relationship is largely linear for areas of main data density

# Model the relationship between overnight mass loss and daily weight gain, allowing for the mass loss to be condition dependent (i.e. mass-dependent mass loss)

# For the purposes of modelling we will separate into three seasons to reflect the different environmental conditions that typically occur at three times of the year. 

# Repeat by season for the main text for a more direct comparison
  oml <- oml %>% 
    mutate(Season = as.factor(case_when(monthabb %in% month.abb[5:8] ~ "Winter", 
                                        monthabb %in% month.abb[9:11] ~ "EarlySummer", 
                                        monthabb %in% month.abb[c(1:4, 12)] ~ "LateSummer")), 
           breeding_season = if_else(month %in% 1:7, year(Date) - 1, year(Date))) %>% 
    mutate(breeding_season = paste0(breeding_season, "/", breeding_season + 1))
  
# set up the fixed and random effects  
  oml$IndividID <- as.factor(oml$IndividID)
  oml$GroupRef <- as.factor(oml$GroupRef)
  oml$breeding_season <- as.factor(oml$breeding_season)
  oml$PriorDailyMassGain_abs_z <- as.numeric(scale(oml$PriorDailyMassGain_abs))
  oml$mass_resid_z <- as.numeric(scale(oml$mass_resid))
 
  # fit the model 
  m1 <- glmmTMB(OvernightMassLoss_g ~ PriorDailyMassGain_abs_z*mass_resid_z*Season + 
             (1|IndividID) + (1|breeding_season), 
             data = oml)
  summary(m1)
  anova(m1, update(m1,~.-PriorDailyMassGain_abs_z:mass_resid_z:Season))
  r2(m1) # same name as adjusted R2

# Extract other measures of model predictive ability: RMSE & MAE 
  model_fit <- predict(m1, type = "response") %>% 
    tibble(method = "lm", pred = .) %>% 
    bind_cols(obs = oml$OvernightMassLoss_g)

  model_fit %>%
    summarize(MAD = mean(abs(obs - pred), na.rm = TRUE),
              RMSE = sqrt(mean((obs - pred) ^ 2)))
  # MAE = 7.59, RMSE = 10.1

# Predict at values of interest (need their condition)
  newdf <- expand.grid(PriorDailyMassGain_abs = seq(0, 60, 0.1), 
                        mass_resid = c(-50, 0, 50), 
                        Season = unique(oml$Season), 
                        IndividID = as.factor(37), 
                        breeding_season = as.factor("2010/2011")) %>% 
    mutate(PriorDailyMassGain_abs_z = (PriorDailyMassGain_abs - mean(oml$PriorDailyMassGain_abs))/
           sd(oml$PriorDailyMassGain_abs), 
           mass_resid_z = (mass_resid - mean(oml$mass_resid))/
           sd(oml$mass_resid))

  newdf$pred <- predict(m1, newdata = newdf, re.form = NA)
  newdf$SE <- predict(m1, newdata = newdf, re.form = NA, se.fit = TRUE)$se.fit
  newdf$l95ci <- newdf$pred - 1.96*newdf$SE 
  newdf$u95ci <- newdf$pred + 1.96*newdf$SE 
  newdf <- rename(newdf, `Body condition` = mass_resid) %>% 
    mutate(`Body condition` = case_when(`Body condition`  == -50 ~ "-50g", 
                                        `Body condition`  == 0 ~ "0g", 
                                        `Body condition`  == 50 ~ "+50g")) %>% 
    mutate(`Body condition` = factor(`Body condition`, 
                                     levels = c("-50g", "0g", "+50g")), 
           Season = case_when(Season  == "EarlySummer" ~ "Early summer", 
                              Season  == "LateSummer" ~ "Mid to late summer", 
                              Season  == "Winter" ~ "Winter"))

# Create the plot  
  p2 <- ggplot() + 
      geom_vline(xintercept = 0, colour= "grey") +
      geom_hline(yintercept = 0, colour= "grey") +
      geom_abline(intercept = 0, slope = 1, colour= "black", 
                  linetype = 2) +
      geom_line(data = newdf, aes(x = PriorDailyMassGain_abs, y = pred, 
                           group = `Body condition`, 
                           colour =  `Body condition`), linewidth = 1) +
      facet_wrap(~Season) + 
      theme_bw() + 
      theme(axis.text = element_text(size = 10.5, colour = "black"), 
            axis.title = element_text(size = 11.5, colour = "black"), 
            panel.grid = element_blank(), 
            strip.background = element_rect(fill = "lightblue"), 
            strip.text = element_text(size = 11), 
            legend.title = element_text(size = 10.5, hjust = 0.5), 
            legend.text = element_text(size = 10), 
            legend.background = element_blank(),
            plot.tag = element_text(size = 14)) + 
      labs(y = "Overnight body mass\nloss to day 2 (g)",
           x = "Daily body mass gain on day 1 (g)", 
           tag = "A") + 
      scale_x_continuous(breaks = seq(0, 60, 10), 
                         labels = c("", "10", "", "30", "", "50", "")) + 
    scale_y_continuous(breaks = seq(0, 50, 10)) + 
    scale_colour_manual(values = viridis::viridis(9)[c(2,5,8)]) + 
    scale_fill_manual(values = viridis::viridis(9)[c(2,5,8)])

# get the cut-points of the line where DMG = OML  
## (when mass gained in the day is offset by mass loss overnight)
  m1_cut <- newdf %>% 
    group_by(`Body condition`, Season) %>% 
    mutate(deviation = abs(pred - PriorDailyMassGain_abs)) %>% 
    filter(deviation == min(deviation)) %>% 
    mutate(label = paste0(PriorDailyMassGain_abs, "g")) %>% 
    data.frame()

  p2 <- p2 + 
    geom_segment(data = m1_cut, aes(x = PriorDailyMassGain_abs, 
                                    xend = PriorDailyMassGain_abs , 
                                    y = pred, 
                                    yend = 0, colour = Body.condition), 
                 alpha = 0.8, linetype = 2) + 
  geom_point(data = m1_cut, aes(x = PriorDailyMassGain_abs, y = pred, 
                                colour = Body.condition),   alpha = 0.8) +
  geom_text(data = filter(m1_cut, Body.condition == "0g"), 
            aes(x = PriorDailyMassGain_abs, y = -2, 
                label = paste0(signif(pred, 3), "g")), 
            colour = viridis::viridis(9)[5], size  = 2.7) 

# What does this work out as in percentages for the different body condition scores
  dwg %>%
      mutate(Season = case_when(month(Date) %in% 5:8 ~ "Winter",
                                month(Date)  %in% 9:11 ~ "EarlySummer", 
                                month(Date)  %in% c(12, 1:4) ~ "LateSummer")) %>% 
    filter(agecat == "Adult", between(mass_resid, -10, 10), 
           PregnantLactating != "Pregnant") %>% 
    group_by(Season) %>% 
    summarise(mean(MorningWeight)) %>% 
    data.frame()

  (657.5611 + 26.1)/657.5611 # early summer cut @ 26g  -- > 4.0%
  (647.1172 + 27.5)/647.1172 # late summer cut @ 27.5g -- > 4.3%
  (654.3306 + 28.9)/654.3306 # winter cut @ 28.9g  -- > 4.4%

# Calculate the mean condition, mass and daily mass gain in each period using the larger data set
dwg %>% 
   group_by(Season) %>% 
   summarise(mean(MorningWeight), 
             mean(mass_resid, na.rm = T), 
             mean(WeightGain_g)) %>% 
   data.frame() # This implies that most of the time meerkats are able to meet their daily mass needs across the year

#------------------------------------------------
# Estimate the probability of reaching the required threshold needed to maintain or gain body mass over consecutive days given the conditions in that month
# to do so we can ask whether DMG exceeds the pre-identified average threshold above (using the oml or the dwg data set)

#oml <- oml %>% 
#  mutate(neutralovertwodays = if_else(PriorDailyMassGain_abs >= OvernightMassLoss_g, 1, 0), 
#         exceedsthresholdmass = case_when(Season == "EarlySummer" & PriorDailyMassGain_abs >= 26 ~ 1,
#                                       Season == "EarlySummer" & PriorDailyMassGain_abs < 26 ~ 0, 
#                                       Season == "LateSummer" & PriorDailyMassGain_abs >= 27.5 ~ 1,
#                                       Season == "LateSummer" & PriorDailyMassGain_abs < 27.5  ~ 0,
#                                       Season == "Winter" & PriorDailyMassGain_abs >= 28.9 ~ 1,
#                                       Season == "Winter" & PriorDailyMassGain_abs < 28.9 ~ 0))

dwg <- dwg %>% 
  mutate(exceedsthresholdmass = case_when(Season == "EarlySummer" & WeightGain_g  >= 26.1 ~ 1,
                                          Season == "EarlySummer" & WeightGain_g  < 26.1 ~ 0, 
                                          Season == "LateSummer" & WeightGain_g  >= 27.5 ~ 1,
                                          Season == "LateSummer" & WeightGain_g  < 27.5  ~ 0,
                                          Season == "Winter" & WeightGain_g  >= 28.9 ~ 1,
                                          Season == "Winter" & WeightGain_g  < 28.9 ~ 0), 
         breeding_season = if_else(month %in% 1:7, year(Date) - 1, year(Date))) %>% 
  mutate(breeding_season = paste0(breeding_season, "/", breeding_season + 1))


# Set up the model with a season-specific smoother
oml$Season <- factor(oml$Season) ; dwg$Season <- factor(dwg$Season) # needs to be a factor to run 
oml$GroupRef <- factor(oml$GroupRef) ; dwg$GroupRef <- factor(dwg$GroupRef)

oml <- oml %>% 
  left_join(dplyr::select(dwg, GroupRef, Date, GroupSize_1moavg) %>% 
              distinct(), by = c("GroupRef", "Date")) %>% 
  filter(!is.na(ndvi_group)) %>% 
  mutate()

 # WAS DMG >= threshold mass #----------------------------

  # We are running as a gamm so as to control for individual and group-level variation in the form of conventional random effects
 # This can be carried out via mgcv::gamm or gamm4::gamm4. We use the former as it runs in a more reasonable time
  # We also fit via job to free up the console while it runs
  job( { 
    m2_full <- gamm(exceedsthresholdmass ~ Season + 
                 te(tempmax, ndvi_group, by = Season, k = c(5,5), bs = "cr") + 
                 s(GroupSize_1moavg, k = 3, bs = "cr"), 
               random = list(IndividID = ~ 1, breeding_season = ~ 1), 
               data = dwg, 
               method = "REML", 
               family = "binomial") 
    })

  saveRDS(m2_full, "ThresholdmassModel.rds")
  m2_full <- readRDS("ThresholdmassModel.rds")
  summary(m2_full$gam)
  summary(m2_full$lme)

  # diagnostics
  layout(matrix(1:4, ncol = 2))
  k.check(m2_full$gam)
  layout(1)

  # Quick look at the estimated smooths before predicting manually 
  # early summer (Sep-Nov)
  mgcv::plot.gam(m2_full$gam, select = 1, scheme = 2, 
                 shift = 0.5683,
                 trans = inv.logit,
                 main = NA, too.far = 0.03)

  # mid to late summer (Dec-Apr)
  mgcv::plot.gam(m2_full$gam, select = 2, scheme = 2, 
                 shift = 0.5683  + 0.1554,
                 trans = inv.logit,
                 main = NA, too.far = 0.03)

  # Winter (Jun-Aug)
  mgcv::plot.gam(m2_full$gam, select = 3, scheme = 2, 
                 shift = 0.5683 - 0.3257,
                 trans = inv.logit,
                 main = NA, too.far = 0.03)
  
  # group size
  mgcv::plot.gam(m2_full$gam, select = 4, trans = inv.logit)
  
  
  # Alternatively one can parameterise as gamm4 which uses lmer rather than nlme
  # Note that this will take a very long time to run
  #job( { 
  #  m2_full_gamm4 <- gamm4::gamm4(exceedsthresholdmass ~ Season + 
  #                 t2(tempmax, ndvi_group, by = Season, k = c(5,5), bs = "cr") + 
  #                 s(GroupSize_1moavg, k = 5, bs = "cr"),    # need slightly higher k to run
  #                 random = ~ (1|IndividID) + (1|breeding_season), 
  #                 data = dwg, 
  #                REML = TRUE, 
  #               family = "binomial") 
  #})
  

# -------------------

 # Predict for each month the probability of exceeded the daily threshold
  # We can predict either by taking it from the smooths and later adding the intercepts, 
  # or, we can predict directly from the model across a new data frame (see later code). 
  
  # Predicted from the smooths: First get the predictions for the three different periods
  # Note that these predictions are on the [logit] link scale
  sm_earlysummer <- smooth_estimates(m2_full, n = 300, select = 1)
  sm_latesummer <- smooth_estimates(m2_full, n = 300, select = 2)
  sm_winter <- smooth_estimates(m2_full, n = 300, select = 3)
 
  # Now we need to identify the region of the prediction surface to plot for each season, 
  # so that we can restrict the smoothers to where the data is. 
  # I will do this by isolating the 95% ellipse for the temp-ndvi combination in each season
  d_earlysummer <- dwg %>% 
    filter(Season == "EarlySummer") %>% 
    dplyr::select(Date, GroupRef, tempmax, ndvi_group) %>% 
    distinct() %>% 
    dplyr::select(tempmax, ndvi_group)
  
  d_latesummer <- dwg %>% 
    filter(Season == "LateSummer") %>% 
    dplyr::select(Date, GroupRef, tempmax, ndvi_group) %>% 
    distinct() %>% 
    dplyr::select(tempmax, ndvi_group)
  
  d_winter <- dwg %>% 
    filter(Season == "Winter" & tempmax > 10) %>% 
    dplyr::select(Date, GroupRef, tempmax, ndvi_group) %>% 
    distinct() %>% 
    dplyr::select(tempmax, ndvi_group)
  
  # density functions
  kd_earlysummer <- ks::kde(d_earlysummer, compute.cont=TRUE, h=0.2)
  kd_latesummer <- ks::kde(d_latesummer, compute.cont=TRUE, h=0.2)
  kd_winter <- ks::kde(d_winter, compute.cont=TRUE, h=0.2)
  
  # extract the 95% polygon (5% probability mass) 
  get_contour <- function(kd_out=kd, prob="5%") {
    contour_95 <- with(kd_out, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont[prob])[[1]])
    as_tibble(contour_95) %>% 
      mutate(prob = prob)
  }
  
  dat_out_earlysummer <- map_dfr(c("5%"), ~ get_contour(kd_earlysummer, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup() %>% 
    data.frame()
  
  dat_out_latesummer <- map_dfr(c("5%"), ~ get_contour(kd_latesummer, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup() %>% 
    data.frame()
  
  dat_out_winter <- map_dfr(c("5%"), ~ get_contour(kd_winter, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup() %>% 
    data.frame()
  
  # check it's done a good job
  ggplot(data = dat_out_earlysummer, aes(x, y)) +
    geom_point(data = d_earlysummer, aes(x = tempmax, y = ndvi_group), alpha = 0.1) + 
    geom_path(aes(group = prob), colour = "red", linewidth = 1) + 
    labs(title = "Early summer (Sep to Nov)")
  
  ggplot(data = dat_out_latesummer, aes(x, y)) +
    geom_point(data = d_latesummer, aes(x = tempmax, y = ndvi_group), alpha = 0.1) + 
    geom_path(aes(group = prob), colour = "red", linewidth = 1) + 
    labs(title = "Mid to late summer (Dec to Apr)")
  
  ggplot(data = dat_out_winter, aes(x, y)) +
    geom_point(data = d_winter, aes(x = tempmax, y = ndvi_group), alpha = 0.1) + 
    geom_path(aes(group = prob), colour = "red", linewidth = 1) + 
    labs(title = "Winter (May to Aug)")
  
  # Now I want to only consider the smoother predictions falling inside the 95% ellipse for each period
  # Create a matrix of coordinates for the polygon
  polygon_coords_1 <- matrix(c(dat_out_earlysummer$x, dat_out_earlysummer$y), ncol = 2, byrow = FALSE)
  polygon_coords_2 <- matrix(c(dat_out_latesummer$x, dat_out_latesummer$y), ncol = 2, byrow = FALSE)
  polygon_coords_3 <- matrix(c(dat_out_winter$x, dat_out_winter$y), ncol = 2, byrow = FALSE)
  # Close each polygon (add the first point at the end to close the loop)
  polygon_coords_1 <- rbind(polygon_coords_1, polygon_coords_1[1, ])
  polygon_coords_2 <- rbind(polygon_coords_2, polygon_coords_2[1, ])
  polygon_coords_3 <- rbind(polygon_coords_3, polygon_coords_3[1, ])
  # Create an sf polygon object using st_polygon
  polygon_sf_1 <- st_polygon(list(polygon_coords_1))
  polygon_sf_2 <- st_polygon(list(polygon_coords_2))
  polygon_sf_3 <- st_polygon(list(polygon_coords_3))
  # Create a matrix of coordinates for the points 
  points_coords_1 <- matrix(c(sm_earlysummer$tempmax, sm_earlysummer$ndvi_group), ncol = 2, byrow = FALSE)
  points_coords_2 <- matrix(c(sm_latesummer$tempmax, sm_latesummer$ndvi_group), ncol = 2, byrow = FALSE)
  points_coords_3 <- matrix(c(sm_winter$tempmax, sm_winter$ndvi_group), ncol = 2, byrow = FALSE)
  # Create an sf geometry list of points using st_sfc
  points_sf_1 <- st_as_sf(data.frame(points_coords_1), coords = c(1,2))
  points_sf_2 <- st_as_sf(data.frame(points_coords_2), coords = c(1,2))
  points_sf_3 <- st_as_sf(data.frame(points_coords_3), coords = c(1,2))

  # Which points in the prediction surface fall within the within the 95% ellipse
  sm_earlysummer$inside <- st_intersects(points_sf_1, polygon_sf_1, sparse = T) %>% 
    as.numeric() %>% 
    replace_na(0) 
  sm_earlysummer <- filter(sm_earlysummer, inside == 1)
  sm_latesummer$inside <- st_intersects(points_sf_2, polygon_sf_2, sparse = T) %>% 
    as.numeric() %>% 
    replace_na(0) 
  sm_latesummer <- filter(sm_latesummer, inside == 1)
  sm_winter$inside <- st_intersects(points_sf_3, polygon_sf_3, sparse = T) %>% 
    as.numeric() %>% 
    replace_na(0) 
  sm_winter <- filter(sm_winter, inside == 1)

  # bind them  together before plotting
  sm <- bind_rows(sm_earlysummer, sm_latesummer, sm_winter)
  
  # The smoother data as it sits are the partial effects. 
  # I need to add the model intercept to each period to get the 
  # absolute predicted probability
  summary(m2_full$gam)
  sm$pred <- NA
  sm$pred[sm$Season == "EarlySummer"] <- sm$.estimate[sm$Season == "EarlySummer"] + 0.5683
  sm$pred[sm$Season == "LateSummer"] <- sm$.estimate[sm$Season == "LateSummer"] + 0.5683 + 0.1554     
  sm$pred[sm$Season == "Winter"] <- sm$.estimate[sm$Season == "Winter"] + 0.5683 - 0.3257           
  
  # Need to now convert the smoothers back to the response scale
  hist(inv.logit(sm$pred))
  sm$pred <- inv.logit(sm$pred)
  
  # Define Jet colour ramp
  jet.colors <- colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan", 
                                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
  
  # Rename the factor levels for the plot
  sm$Season <- case_when(sm$Season == "EarlySummer" ~ "Early summer",
                         sm$Season == "LateSummer" ~ "Mid to late summer",
                         sm$Season == "Winter" ~ "Winter")
  
  # produce plot from the smoothers
  p3 <- ggplot(sm, aes(x = tempmax, y = ndvi_group)) + 
    geom_raster(aes(fill = pred)) +
    geom_contour(aes(z = pred),  breaks = seq(0.1, 1, 0.1), 
                 colour = "black", linewidth = 0.3, alpha = 0.2) +
    labs(title = NULL, caption = NULL,
         y = "NDVI", x = "Daily maximum temperature (°C)",   
         fill = "DMG > threshold", 
         tag = "B") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 10.5, colour = "black"), 
          strip.text = element_text(size = 11, colour = "black"), 
          axis.title = element_text(size = 11.5, colour = "black"), 
          legend.title.align = 0.5, 
          legend.title = element_text(size = 10.5),
          legend.text = element_text(size = 10),
          strip.background = element_rect(fill = "lightblue"), 
          panel.grid = element_blank(), 
          plot.tag = element_text(size = 14)) + 
    facet_wrap(~Season) +
    coord_cartesian(expand = FALSE) +
    scale_fill_gradientn(colors = jet.colors(20), 
                         breaks = seq(0.1, 1, 0.1), 
                         rescaler = ~ scales::rescale_mid(.x, mid = 0.55)) + 
    scale_x_continuous(breaks = seq(10, 45, 5), 
                       labels = c("10", "15", "20", "25", 
                                  "30", "35", "40", ""),
                       limits = c(10, 45)) + 
    scale_y_continuous(breaks = seq(0.15, 0.5, 0.05), 
                       labels = c("", "0.2", "", "0.3", 
                                  "", "0.4", "", "0.5"),
                       limits = c(0.13, 0.41))

  # put the two plots together
  p2/p3
  
  #N.B. If we instead wanted to predict the response directly this code is sufficient. 
  # The resulting plot is identical to going via smooth estimates, but is included to demonstrate this point
  # Create prediction data frame
    #range(dwg$tempmax) ; range(dwg$ndvi_group)
    #pred_df <- expand.grid(tempmax = seq(1.96, 44.36, length.out = 300), 
    #                       ndvi_group =  seq(0.1356, 0.543, length.out = 300), 
    #                       Season = levels(dwg$Season), 
    #                       GroupSize_1moavg = mean(dwg$GroupSize_1moavg))
  
  # predict the response
    #pred_df$pred <- predict(m2_full$gam, newdata = pred_df, type = "response")
    #pred_earlysummer <- filter(pred_df, Season == "EarlySummer")
    #pred_latesummer <- filter(pred_df, Season == "LateSummer")
    #pred_winter <- filter(pred_df, Season == "Winter")
  
  # create a matrix of coordinates for the points (using either the smooths or the raw data)
    #points_coords_1_pred <- matrix(c(pred_earlysummer$tempmax, pred_earlysummer$ndvi_group), 
     # ncol = 2, byrow = FALSE)
    #points_coords_2_pred <- matrix(c(pred_latesummer$tempmax, pred_latesummer$ndvi_group), 
     # ncol = 2, byrow = FALSE)
    #points_coords_3_pred <- matrix(c(pred_winter$tempmax, pred_winter$ndvi_group), 
     # ncol = 2, byrow = FALSE)
  
  # predicted response
    # points_sf_1_pred <- st_as_sf(data.frame(points_coords_1_pred), coords = c(1,2))
    # points_sf_2_pred <- st_as_sf(data.frame(points_coords_2_pred), coords = c(1,2))
    # points_sf_3_pred <- st_as_sf(data.frame(points_coords_3_pred), coords = c(1,2))
  
  # Which points in the prediction surface fall within the within the 95% ellipse
  #pred_earlysummer$inside <- st_intersects(points_sf_1_pred, polygon_sf_1, sparse = T) %>% 
  #  as.numeric() %>% 
  #  replace_na(0) 
  #pred_earlysummer <- filter(pred_earlysummer, inside == 1)
  #pred_latesummer$inside <- st_intersects(points_sf_2_pred, polygon_sf_2, sparse = T) %>% 
  #  as.numeric() %>% 
  #  replace_na(0) 
  #pred_latesummer <- filter(pred_latesummer, inside == 1)
  #pred_winter$inside <- st_intersects(points_sf_3_pred, polygon_sf_3, sparse = T) %>% 
  #  as.numeric() %>% 
  #  replace_na(0) 
  #pred_winter <- filter(pred_winter, inside == 1)
  
  # bind them  together before plotting
  #pred_df <- bind_rows(pred_earlysummer, pred_latesummer, pred_winter)
  
  # Rename the factor levels for the plot
  # pred_df$Season <- case_when(pred_df$Season == "EarlySummer" ~ "Early summer",
  #                            pred_df$Season == "LateSummer" ~ "Mid to late summer",
  #                            pred_df$Season == "Winter" ~ "Winter")
  
  # p3_alt <- ggplot(pred_df, aes(x = tempmax, y = ndvi_group)) + 
  #  geom_raster(aes(fill = pred)) +
  #  geom_contour(aes(z = pred),  breaks = seq(0.1, 1, 0.1), 
  #               colour = "black", linewidth = 0.3, alpha = 0.2) +
  #  labs(title = NULL, caption = NULL,
  #       y = "NDVI", x = "Daily maximum temperature (°C)",   
  #       fill = "DMG > threshold", 
  #       tag = "B") + 
  #  theme_bw() + 
  #  theme(axis.text = element_text(size = 10.5, colour = "black"), 
  #        strip.text = element_text(size = 11, colour = "black"), 
  #        axis.title = element_text(size = 11.5, colour = "black"), 
  #        legend.title.align = 0.5, 
  #        legend.title = element_text(size = 10.5),
  #        legend.text = element_text(size = 10),
  #        strip.background = element_rect(fill = "lightblue"), 
  #        panel.grid = element_blank(), 
  #        plot.tag = element_text(size = 14)) + 
  #  facet_wrap(~Season) +
  #  coord_cartesian(expand = FALSE) +
  #  scale_fill_gradientn(colors = jet.colors(20), 
  #                       breaks = seq(0.1, 1, 0.1), 
  #                       rescaler = ~ scales::rescale_mid(.x, mid = 0.55)) + 
  #  scale_x_continuous(breaks = seq(10, 45, 5), 
  #                     labels = c("10", "15", "20", "25", 
  #                               "30", "35", "40", ""),
  #                     limits = c(10, 45)) + 
  #  scale_y_continuous(breaks = seq(0.15, 0.5, 0.05), 
  #                     labels = c("", "0.2", "", "0.3", 
  #                                "", "0.4", "", "0.5"),
  #                     limits = c(0.13, 0.41))
    
  # p3_alt  # produces the same plot as p_3
  
#-------------------------------------------
  
  # We can also model the DWG directly (without controlling for condition)
  job( { 
    m3 <- gamm(WeightGain_g ~ Season + 
                 te(tempmax, ndvi_group, by = Season, k = c(5,5), bs = "cr") + 
                 s(GroupSize_1moavg, k = 3, bs = "cr"), 
               random = list(IndividID = ~ 1, breeding_season = ~ 1), 
               data = dwg, 
               method = "REML", 
               family = "gaussian") 
  })

  
  #saveRDS(m3, "DMGModel.rds")
  m3 <- readRDS("DMGModel.rds")
  summary(m3$gam)
  summary(m3$lme)

    # diagnostics
  layout(matrix(1:4, ncol = 2))
  k.check(m3$gam)
  layout(1)
  
  # Quick look at the estimated smooths before predicting manually 
  par(mfrow = c(1,3))
  # early summer (Sep-Nov)
  mgcv::plot.gam(m3$gam, select = 1, scheme = 2, 
                 shift = 35.273,
                 main = NA, too.far = 0.03)
  
  # mid to late summer (Dec-Apr)
  mgcv::plot.gam(m3$gam, select = 2, scheme = 2, 
                 shift = 35.273  + 1.065,
                 main = NA, too.far = 0.03)
  
  # Winter (Jun-Aug)
  mgcv::plot.gam(m3$gam, select = 3, scheme = 2, 
                 shift = 35.273 - -3.489,
                 main = NA, too.far = 0.03)
  
  # group size
  mgcv::plot.gam(m3$gam, select = 4,  shift = 35.273)
  
  # Predict for each month the probability of exceeded the daily threshold
  #First get the predictions for the three different periods
  #Note that these predictions are on the [logit] link scale
  sm2_earlysummer <- smooth_estimates(m3, n = 300, select = 1)
  sm2_latesummer <- smooth_estimates(m3, n = 300, select = 2)
  sm2_winter <- smooth_estimates(m3, n = 300, select = 3)
  
  # Not identify the region of the prediction surface to plot for each season 
  # Only plot the prediction surface that falls within the 95% ellipse
  d_earlysummer2 <- dwg %>% 
    filter(Season == "EarlySummer") %>% 
    dplyr::select(Date, GroupRef, tempmax, ndvi_group) %>% 
    distinct() %>% 
    dplyr::select(tempmax, ndvi_group)
  
  d_latesummer2 <- dwg %>% 
    filter(Season == "LateSummer") %>% 
    dplyr::select(Date, GroupRef, tempmax, ndvi_group) %>% 
    distinct() %>% 
    dplyr::select(tempmax, ndvi_group)
  
  d_winter2 <- dwg %>% 
    filter(Season == "Winter" & tempmax > 10) %>% 
    dplyr::select(Date, GroupRef, tempmax, ndvi_group) %>% 
    distinct() %>% 
    dplyr::select(tempmax, ndvi_group)
  
  # density functions
  kd_earlysummer2 <- ks::kde(d_earlysummer2, compute.cont=TRUE, h=0.2)
  kd_latesummer2 <- ks::kde(d_latesummer2, compute.cont=TRUE, h=0.2)
  kd_winter2 <- ks::kde(d_winter2, compute.cont=TRUE, h=0.2)
  
  # extract the 95% polygon (5% probability mass) 
  dat_out_earlysummer2 <- map_dfr(c("5%"), ~ get_contour(kd_earlysummer2, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup() %>% 
    data.frame()
  
  dat_out_latesummer2 <- map_dfr(c("5%"), ~ get_contour(kd_latesummer2, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup() %>% 
    data.frame()
  
  dat_out_winter2 <- map_dfr(c("5%"), ~ get_contour(kd_winter2, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup() %>% 
    data.frame()
  
  # Now I want to only consider the smoother predictions falling inside the 95% ellipse for each period
  # Create a matrix of coordinates for the polygon
  polygon_coords2_1 <- matrix(c(dat_out_earlysummer2$x, dat_out_earlysummer2$y), ncol = 2, byrow = FALSE)
  polygon_coords2_2 <- matrix(c(dat_out_latesummer2$x, dat_out_latesummer2$y), ncol = 2, byrow = FALSE)
  polygon_coords2_3 <- matrix(c(dat_out_winter2$x, dat_out_winter2$y), ncol = 2, byrow = FALSE)
  # Close each polygon (add the first point at the end to close the loop)
  polygon_coords2_1 <- rbind(polygon_coords2_1, polygon_coords2_1[1, ])
  polygon_coords2_2 <- rbind(polygon_coords2_2, polygon_coords2_2[1, ])
  polygon_coords2_3 <- rbind(polygon_coords2_3, polygon_coords2_3[1, ])
  # Create an sf polygon object using st_polygon
  polygon_sf2_1 <- st_polygon(list(polygon_coords2_1))
  polygon_sf2_2 <- st_polygon(list(polygon_coords2_2))
  polygon_sf2_3 <- st_polygon(list(polygon_coords2_3))
  # Create a matrix of coordinates for the points
  points_coords2_1 <- matrix(c(sm2_earlysummer$tempmax, sm2_earlysummer$ndvi_group), ncol = 2, byrow = FALSE)
  points_coords2_2 <- matrix(c(sm2_latesummer$tempmax, sm2_latesummer$ndvi_group), ncol = 2, byrow = FALSE)
  points_coords2_3 <- matrix(c(sm2_winter$tempmax, sm2_winter$ndvi_group), ncol = 2, byrow = FALSE)
  # Create an sf geometry list of points using st_sfc
  points_sf2_1 <- st_as_sf(data.frame(points_coords2_1), coords = c(1,2))
  points_sf2_2 <- st_as_sf(data.frame(points_coords2_2), coords = c(1,2))
  points_sf2_3 <- st_as_sf(data.frame(points_coords2_3), coords = c(1,2))
  # Which points in the prediction surface fall within the within the 95% ellipse
  sm2_earlysummer$inside <- st_intersects(points_sf2_1, polygon_sf2_1, sparse = T) %>% 
    as.numeric() %>% 
    replace_na(0) 
  sm2_earlysummer <- filter(sm2_earlysummer, inside == 1)
  sm2_latesummer$inside <- st_intersects(points_sf2_2, polygon_sf2_2, sparse = T) %>% 
    as.numeric() %>% 
    replace_na(0) 
  sm2_latesummer <- filter(sm2_latesummer, inside == 1)
  sm2_winter$inside <- st_intersects(points_sf2_3, polygon_sf2_3, sparse = T) %>% 
    as.numeric() %>% 
    replace_na(0) 
  sm2_winter <- filter(sm2_winter, inside == 1)
  
  # bind them  together before plotting
  sm2 <- bind_rows(sm2_earlysummer, sm2_latesummer, sm2_winter)
  
  # The data as it sits are the partial effects. 
  # I need to add the model intercept to each period to get the 
  # absolute predicted probability
  summary(m3$gam)
  sm2$pred <- NA
  sm2$pred[sm2$Season == "EarlySummer"] <- sm2$.estimate[sm2$Season == "EarlySummer"] + 35.273
  sm2$pred[sm2$Season == "LateSummer"] <- sm2$.estimate[sm2$Season == "LateSummer"] + 35.273 + 1.065     
  sm2$pred[sm2$Season == "Winter"] <- sm2$.estimate[sm2$Season == "Winter"] + 35.273 - 3.489           
  
  sm2$Season <- case_when(sm2$Season == "EarlySummer" ~ "Early summer",
                          sm2$Season == "LateSummer" ~ "Mid to late summer",
                          sm2$Season == "Winter" ~ "Winter")
  
  # produce plot
  p4 <- ggplot(sm2, aes(x = tempmax, y = ndvi_group)) + 
    geom_raster(aes(fill = pred)) +
    geom_contour(aes(z = pred),  breaks = seq(0, 100, 5), 
                 colour = "black", linewidth = 0.3, alpha = 0.2) +
    labs(title = NULL, caption = NULL,
         y = "NDVI", x = "Daily maximum temperature (°C)",   
         fill = "Daily body\nmass gain (g)") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 10.5, colour = "black"), 
          strip.text = element_text(size = 11, colour = "black"), 
          axis.title = element_text(size = 11.5, colour = "black"), 
          legend.title.align = 0.5, 
          legend.title = element_text(size = 10.5),
          legend.text = element_text(size = 10),
          strip.background = element_rect(fill = "lightblue"), 
          panel.grid = element_blank(), 
          plot.tag = element_text(size = 14)) + 
    facet_wrap(~Season) +
    coord_cartesian(expand = FALSE) +
    scale_fill_gradientn(colors = jet.colors(20), 
                         breaks = seq(10, 50, 5), 
                         rescaler = ~ scales::rescale_mid(.x, mid = 30)) + 
    scale_x_continuous(breaks = seq(10, 45, 5), 
                       labels = c("10", "15", "20", "25", 
                                  "30", "35", "40", ""),
                       limits = c(10, 45)) + 
    scale_y_continuous(breaks = seq(0.15, 0.5, 0.05), 
                       labels = c("", "0.2", "", "0.3", 
                                  "", "0.4", "", "0.5"),
                       limits = c(0.13, 0.41))
  
#######################   END   ###########################