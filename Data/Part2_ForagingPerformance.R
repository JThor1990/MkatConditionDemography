###########################################################

# LINKING CLIMATE VARIABILITY TO DEMOGRAPHY IN COOPERATIVELY BREEDING MEERKATS

# AUTHORS: Thorley, Duncan, Manser and Clutton-Brock (2025)
# PAPER: ECOLOGICAL MONOGRAPHS

# PART 2: MODELLING THE EFFECTS OF NDVI AND TEMPERATURE ON FORAGING PERFORMANCE

###########################################################

# set up the working environment and load required packages
lapply(c("ggplot2", "grDevices", "ks", "lubridate", "mgcv", "patchwork", "sf", "tidyverse"),  FUN = library, character.only = TRUE)

# Set working directory 
setwd("INSERT FILE PATH HERE")

# set up two colour palettes to use in the plots that are produced from here forward
jet.colours <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jet.colours.rev <- colorRampPalette(c("#7F0000", "red", "#FF7F00","yellow","#7FFF7F","cyan","#007FFF","blue","#00007F"))
  
# read master forage DF
forage_df <- read.csv("Data\\MeerkatForagingData.csv") %>% 
  filter(!is.na(GroupRef))

# Set up 3-way season 
forage_df <- forage_df %>%
  mutate(Season = factor(case_when(FocalMonth %in% c(12,1,2,3,4) ~ "LateSummer",
                          FocalMonth %in% c(5,6,7,8) ~ "Winter",
                          FocalMonth %in% c(9,10,11) ~ "EarlySummer"))) 

# And set up constraints on the temperature:ndvi surface for later plotting of the interaction
# This involves generating a 95% KDE polygon for each season
constraintPoly <- list(LateSummer=NA, Winter=NA, EarlySummer=NA)
for(i in c("LateSummer","Winter","EarlySummer")){
  d <- forage_df %>% filter(Season==i) %>%
    dplyr::select(numAvgHour,tempmax) %>%
    rename(x = numAvgHour, y = tempmax)
  # density function
  kd <- ks::kde(d, compute.cont=TRUE, h=0.2)
  # extract the 95% polygon (5% probability mass)
  get_contour <- function(kd_out=kd, prob="5%") {
    contour_95 <- with(kd_out, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont[prob])[[1]])
    as_tibble(contour_95) %>%
      mutate(prob = prob)
  }
  
  dat_out <- map_dfr(c("5%"), ~get_contour(kd, .)) %>%
    group_by(prob) %>%
    mutate(n_val = 1:n()) %>%
    ungroup() %>%
    data.frame()
  
  # Create a matrix of coordinates for the polygon
  polygon_coords <- matrix(c(dat_out$x, dat_out$y), ncol = 2, byrow = FALSE)
  # Close the polygon (add the first point at the end to close the loop)
  polygon_coords <- rbind(polygon_coords, polygon_coords[1, ])
  # Create an sf polygon object using st_polygon
  polygon_sf <- st_polygon(list(polygon_coords))
  
  constraintPoly[[i]] <- polygon_sf
}  # generate a list of three polygons, one for each season

# Before modelling, visualise the seasonal variation in temperature across the focal study period (1996-2001)
ggplot(forage_df, aes(x = as.Date(FocalDate), y = tempmax)) +
  geom_jitter(alpha = 0.3, width = 5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_hline(yintercept = 40, lty = 2, col = "red") +
  theme_bw()

tempMax <- forage_df %>% 
  group_by(Season) %>% 
  select(Season, FocalDate, tempmax) %>% 
  distinct()
summarise(tempMax, mean = mean(tempmax)) # seasonal mean temperature
tempMax %>% ungroup() %>% summarise(mean = mean(tempmax)) # annual mean temp 

# And also get some infromation on the distribution of focals across the period and across groups and individuals
ggplot(forage_df, aes(x = as.Date(FocalDate))) +
  geom_histogram(binwidth = 60)

nFocalMonth <- forage_df %>% 
  ungroup() %>%
  select(FocalYear, FocalMonth, GroupRef, FocalDate) %>%
  distinct() %>%
  group_by(FocalYear, FocalMonth, GroupRef) %>%
  summarise(n=n())
mean(nFocalMonth$n) ; sd(nFocalMonth$n) # mean +- sd focals / group / month

forage_df %>% ungroup() %>% summarise(mean = mean(numAvgHour)) # mean focal time-of-day

nFocalIndivid <- forage_df %>% 
  ungroup() %>% 
  group_by(IndividID) %>% 
  summarise(n=n())
mean(nFocalIndivid$n) ; sd(nFocalIndivid$n) # mean +- sd focals / individual
median(nFocalIndivid$n) ; range(nFocalIndivid$n)

length(unique(forage_df$Focal_ID)) # n = 6357 focals
length(unique(forage_df$IndividID)) # from 210 individuals

#---------------------------------------

# Model 1, FORAGING SUCCESS: Prey Acquisition Rates (GAM)

#---------------------------------------

forage_df <- mutate(forage_df, 
                    GroupRef = factor(GroupRef),
                    IndividID = factor(IndividID)) 
length(unique(forage_df$GroupRef))  # n = 13 groups
length(unique(forage_df$IndividID)) # n = 210 individuals

# Model of prey acquisition rate   
mod_forsuccess <- bam(SuccForageBoutCount ~ s(ndvi, bs = "cr", k=5) +
                        te(numAvgHour, tempmax, bs = "cr", by = Season, k = c(8,8)) + 
                        Season + 
                        offset(log(BehavObsDuration)) + 
                        s(IndividID,bs="re") + 
                        s(GroupRef,bs="re"),
               method = "REML",
               family = "nb",
               data = forage_df) 

# Model summary and basic smooth check
summary(mod_forsuccess)
anova(mod_forsuccess)
gam.check(mod_forsuccess)
plot(mod_forsuccess)

# If one instead wanted to fit via GAMM
# mod_forsuccess_GAMM <-  gamm(SuccForageBoutCount ~ s(ndvi, bs="cr", k=5) + 
#                               te(numAvgHour, tempmax, bs = "cr", k=c(8,8), by = Season) + 
#                               Season + 
#                               offset(log(BehavObsDuration)),
#                        random = list(IndividID = ~ 1,
#                                      GroupRef = ~ 1),
#                        data = forage_df,
#                        method = "REML",
#                        family=negbin(5.403))

#summary(mod_forsuccess_GAMM$gam)
#gam.check(mod_forsuccess_GAMM$gam)
    
#  Generate predictions from the model
pred_df <-  data.frame(Season = c("EarlySummer","LateSummer","Winter"),
                       ndvi = c(0.232, 0.274, 0.240)) %>%
  mutate(numAvgHour = list(seq(5.9, 13.2, 0.025))) %>% 
  unnest(col = numAvgHour) %>%
  mutate(tempmax = list(seq(5.5, 42.1, 0.01))) %>% 
  unnest(col = tempmax) %>%
  mutate(IndividID = 3, GroupRef = 7, BehavObsDuration = 1200, DiggingDur = 600) %>% 
  distinct()

points_coords <- matrix(c(pred_df$numAvgHour, pred_df$tempmax), ncol = 2, byrow = FALSE)
points_sf <- st_as_sf(data.frame(points_coords), coords = c(1,2))

# Restrict the predictions to the 95% KDE identified at the top of script
pred_df$insideES <- st_intersects(points_sf, constraintPoly[["EarlySummer"]], sparse = T) %>%
  as.numeric() %>%
  replace_na(0) 
pred_df$insideLS <- st_intersects(points_sf, constraintPoly[["LateSummer"]], sparse = T) %>%
  as.numeric() %>%
  replace_na(0) 
pred_df$insideW <- st_intersects(points_sf, constraintPoly[["Winter"]], sparse = T) %>%
  as.numeric() %>%
  replace_na(0) 

pred_df <- pred_df %>%
  mutate(InsideCombined = case_when(Season == "EarlySummer" ~ insideES,
                                    Season == "LateSummer" ~ insideLS,
                                    Season == "Winter" ~ insideW)) %>% 
  filter(InsideCombined == 1)

# Now everything is there for the prediction
forsucc_pred <- predict(mod_forsuccess, pred_df,
                        exclude = c("s(IndividID)", "s(GroupRef)"),
                        se.fit = TRUE, type = "response")

forsucc_pred <- cbind(pred_df, 
                      tibble(est = forsucc_pred$fit),
                      tibble(se = forsucc_pred$se.fit)) 

# Edit the names of the seasons to match other figures 
forsucc_pred <- forsucc_pred %>% 
  mutate(Season = factor(case_when(Season == "EarlySummer" ~ "Early summer", 
                                   Season == "LateSummer" ~ "Mid to late summer", 
                                   Season == "Winter" ~ "Winter"), 
                         levels = c("Mid to late summer", "Winter", "Early summer")),
         Strip = "Foraging success")

# Plot foraging success ~ temp at different times of the morning
si_plot1 <- ggplot(forsucc_pred, aes(x = numAvgHour, y = tempmax)) +
  geom_raster(aes(fill = est)) +
  geom_contour(aes(z = est), colour = "black", linewidth = 0.3, alpha = 0.2)  +
  scale_fill_gradientn(colours= jet.colours.rev(7)) +
  labs(x = "Hour of Day", y = "Daily maximum temperature (°C)", 
       title = "Foraging success: prey captures per focal", 
       fill = "Prey per foval", tag = "A") +
  facet_wrap(~Season) +
  theme_bw() +
  theme(axis.text = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 11.5, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"), 
        legend.title = element_text(size = 10.5, hjust = 0.5),
        legend.text = element_text(size = 10),
        strip.background = element_rect(fill = "lightblue"), 
        panel.grid = element_blank(), 
        plot.title = element_text(size = 12.5, hjust = 0.5), 
        plot.tag = element_text(size = 14))  +
  scale_x_continuous(breaks = 7:12) + 
  scale_y_continuous(breaks = seq(15, 40, 5))

# Plot a one panel version for the main text
p1 <- ggplot(filter(forsucc_pred, Season == "Mid to late summer"), 
             aes(x = numAvgHour, y = tempmax)) +
  geom_raster(aes(fill = est)) +
  geom_contour(aes(z = est), colour = "black", linewidth = 0.3, alpha = 0.2)  +
  scale_fill_gradientn(colours= jet.colours.rev(7), 
                       breaks = seq(6, 12, 1)) +
  labs(x = "Hour of Day", y = "Daily maximum temperature (°C)", 
       fill = "Prey captures\nper focal", tag = "B") +
  facet_wrap(~Strip) +
  theme_bw() +
  theme(axis.text = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 10.5, hjust = 0.5),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(), 
        strip.text = element_text(size = 12.5, colour = "black"), 
        strip.background = element_rect(fill = "lightblue"), 
        plot.title = element_text(size = 12), 
        plot.tag = element_text(size = 14))  +
  scale_x_continuous(breaks = 7:12) + 
  scale_y_continuous(breaks = seq(15, 40, 5))

# Plot the NDVI effect on foraging success
pred_df2 <- expand.grid(Season = c("LateSummer"), tempmax=32,
                        numAvgHour=9, ndvi=seq(0.20, 0.43, by = 0.01),
                        IndividID=3, GroupRef=7, BehavObsDuration=1200)

forsucc_pred2 <- predict(mod_forsuccess, pred_df2,
                         exclude = c("s(IndividID)","s(GroupRef)"),
                         se.fit = TRUE, type = "link")

forsucc_pred2 <- cbind(pred_df2,
                       tibble(est = forsucc_pred2$fit),
                       tibble(se = forsucc_pred2$se.fit)) 

forsucc_pred2 <- mutate(forsucc_pred2, 
                        l95ci = exp(est - 1.96*se),
                        u95ci = exp(est + 1.96*se), 
                        est = exp(est), 
                        Strip = "Foraging success")

p2 <- ggplot(forsucc_pred2, aes(x = ndvi,y = est)) +
  geom_ribbon(aes(x = ndvi, ymin = l95ci, ymax = u95ci), alpha=0.3) +
  geom_line(lwd=1) +
  labs(x = "NDVI", y = "Prey captures per focal", tag = "A") + 
  facet_wrap(~Strip) +
  theme_bw() +
  theme(axis.text = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 10.5, hjust = 0.5),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(), 
        strip.text = element_text(size = 12.5, colour = "black"), 
        strip.background = element_rect(fill = "lightblue"), 
        plot.title = element_text(size = 12), 
        plot.tag = element_text(size = 14)) 

# plot them next to each other
p2 | p1 

#---------------------------------------

# Model 2, FORAGING EFFORT: Proportion of time spent digging (GAM)

#---------------------------------------

# look at the distribution of the response
hist(sqrt(forage_df$DigProp))
  
mod_foreffort <- bam(sqrt(DigProp) ~ s(ndvi, bs="cr", k=5) +
                         te(numAvgHour, tempmax, bs="cr", by=Season, k=c(8,8)) + Season + 
                         s(IndividID, bs="re") + s(GroupRef, bs="re"),
                  method = "REML",
                  data = forage_df)

# Model summary and basic smooth check
summary(mod_foreffort)
anova(mod_foreffort)
gam.check(mod_foreffort)
plot(mod_foreffort)

# If one instead wanted to fit via GAMM
#  mod_foreffort_GAMM <-  gamm(sqrt(DigProp) ~ s(ndvi, bs="cr", k=5) + 
#                                te(numAvgHour, tempmax, bs="cr", by=Season, k=c(8,8)) +
#                                Season + 
#                                s(IndividID, bs="re") + 
#                                s(GroupRef,bs="re"),
#                              random = list(IndividID = ~ 1, GroupRef = ~ 1),
#                              data = forage_df,
#                              method = "REML",
#                              family=negbin(5.403))

# summary(mod_forsuccess_GAMM$gam)
# gam.check(mod_forsuccess_GAMM$gam)

# Plot foraging success ~ temp at different times of the morning
# (we can use the same prediction data frame as the models have similar structure)
foreffort_pred <- predict(mod_foreffort, pred_df,
                          exclude = c("s(IndividID)", "s(GroupRef)"),
                          se.fit = TRUE,
                          type = "response")

foreffort_pred <- cbind(pred_df,
                        tibble(est = foreffort_pred$fit),
                        tibble(se = foreffort_pred$se.fit))

# Edit the names of the seasons to match other figures 
foreffort_pred <- foreffort_pred %>% 
  mutate(Season = factor(case_when(Season == "EarlySummer" ~ "Early summer", 
                                   Season == "LateSummer" ~ "Mid to late summer", 
                                   Season == "Winter" ~ "Winter"), 
                         levels = c("Mid to late summer", "Winter", "Early summer")),
         Strip = "Foraging effort")

# Plot foraging success ~ temp at different times of the morning
invSq <- function(x) { x^2 }

si_plot2 <- ggplot(foreffort_pred, aes(x = numAvgHour, y = tempmax)) +
  geom_raster(aes(fill = invSq(est))) +
  geom_contour(aes(z = invSq(est)), colour = "black", linewidth = 0.3, alpha = 0.2)  +
  scale_fill_gradientn(colours= jet.colours(7), 
                       breaks = seq(0.24, 0.40, 0.04)) +
  labs(x = "Hour of Day", y = "Daily maximum temperature (°C)", 
       title = "Foraging effort: proportion of time spent digging per focal", 
       fill = "Proportion", tag = "B") +
  facet_wrap(~Season) +
  theme_bw() +
  theme(axis.text = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 11.5, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"), 
        legend.title = element_text(size = 10.5, hjust = 0.5),
        legend.text = element_text(size = 10),
        strip.background = element_rect(fill = "lightblue"), 
        panel.grid = element_blank(), 
        plot.title = element_text(size = 12.5, hjust = 0.5), 
        plot.tag = element_text(size = 14))  +
  scale_x_continuous(breaks = 7:12)

# We can now combine the two si_plots as they appear for the MS
si_plot1 / si_plot2

# Plot a one panel version of si_plot2 for the main text
p3 <- ggplot(filter(foreffort_pred, Season == "Mid to late summer"), 
             aes(x = numAvgHour, y = tempmax)) +
  geom_raster(aes(fill = invSq(est))) +
  geom_contour(aes(z = invSq(est)), colour = "black", linewidth = 0.3, alpha = 0.2)  +
  scale_fill_gradientn(colours= jet.colours(7), 
                       breaks = seq(0.24, 0.40, 0.02)) + 
  labs(x = "Hour of Day", y = "Daily maximum temperature (°C)", 
       fill = "Proportion of\ntime digging", tag = "D") +
  facet_wrap(~Strip) +
  theme_bw() +
  theme(axis.text = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 10.5, hjust = 0.5),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(), 
        strip.text = element_text(size = 12.5, colour = "black"), 
        strip.background = element_rect(fill = "lightblue"), 
        plot.title = element_text(size = 12), 
        plot.tag = element_text(size = 14))  +
  scale_x_continuous(breaks = 7:12)  + 
  scale_y_continuous(breaks = seq(15, 40, 5))


# Plot the NDVI effect on foraging effort
foreffort_pred2 <- predict(mod_foreffort, pred_df2,
                           exclude = c("s(IndividID)","s(GroupRef)"),
                           se.fit = TRUE, type = "response")

foreffort_pred2 <- cbind(pred_df2,
                         tibble(est = foreffort_pred2$fit),
                         tibble(se = foreffort_pred2$se.fit)) #%>%

foreffort_pred2 <- mutate(foreffort_pred2, 
                          l95ci = invSq(est - 1.96*se),
                          u95ci = invSq(est + 1.96*se), 
                          est = invSq(est), 
                          Strip = "Foraging effort")

p4 <- ggplot(foreffort_pred2, aes(x = ndvi,y = est)) +
  geom_ribbon(aes(x = ndvi, ymin = l95ci, ymax = u95ci), alpha=0.3) +
  geom_line(lwd=1) +
  labs(x = "NDVI", y = "Proportion of focal spent digging", tag = "C") + 
  facet_wrap(~Strip) +
  theme_bw() +
  theme(axis.text = element_text(size = 10.5, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 10.5, hjust = 0.5),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(), 
        strip.text = element_text(size = 12.5, colour = "black"), 
        strip.background = element_rect(fill = "lightblue"), 
        plot.title = element_text(size = 12), 
        plot.tag = element_text(size = 14)) + 
  scale_y_continuous(breaks = seq(0.2, 0.32, 0.02))

#-------------------------------------

# Combine all the plots together    

p2 + p1 + p4 + p3

####################################### END #######################################