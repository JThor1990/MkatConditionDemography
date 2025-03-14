###########################################################

# LINKING CLIMATE VARIABILITY TO DEMOGRAPHY IN COOPERATIVELY BREEDING MEERKATS

# AUTHORS: Thorley, Duncan, Manser and Clutton-Brock (2025)
# PAPER: ECOLOGICAL MONOGRAPHS

# PART 5: ANALYSING THE EFFECTS OF BODY CONDITION ON REPRODUCTION AND SURVIVAL

###########################################################

# Load required packages
lapply(c( "DHARMa", "gamm4","ggplot2", "job", "mgcv", "patchwork", "tidyverse"), FUN = library, character.only = TRUE)

# Set up the working directory
setwd("INSERT FILE PATH HERE")

# override MASS select conflict with the desired dplyr select
select <- dplyr::select

# create an inverse logit function for when backtransforming from binomial models
inv.logit <- function(x) { exp(x)/(1+exp(x)) }

# load in the data sets
LitterDF <- read.csv("Data/MeerkatLitterSuccess.csv")
ModDF <- read.csv("Data/MeerkatPregnancySurvival.csv")

# Final data filters for analyses --------

ModDF <- ModDF %>% 
  ungroup() %>% group_by(IndividID) %>% arrange(IndividID,MonthStart) %>%
  mutate(ImpCond=case_when(
    Died == 1 & is.na(avgCondition) == T ~ lag(avgCondition),
    Diss == 1 & is.na(avgCondition) == T ~ lag(avgCondition),
    TRUE ~ avgCondition)) 

ModDF <- ModDF %>% ungroup() %>% 
  mutate(Dominant = factor(Dominant),
         TB = factor(TB),
         IndividID = factor(IndividID),
         Died = as.numeric(Died),
         Diss = as.numeric(Diss))

LitterDF <- LitterDF %>%
  ungroup() %>% mutate(Dominant=factor(Dominant),
                       TB=factor(TB),
                       FemaleID=factor(FemaleID))

#---------------------------------------
# Modelling Pregnancy  
#----------------------------------------

ModDF %>% 
  ungroup() %>% 
  filter(is.na(avgCondition)==F) %>%
  summarise(Individ=n_distinct(IndividID),
            Preg=sum(Preg,na.rm=T),
            n=n())

# Base model with no interaction terms
gammPreg_Base <- gamm4(Preg ~ s(avgCondition, bs="cr", k=5) + TB + Dominant + s(Month,bs="cc",k=5),
                       random = ~(1|IndividID),
                       data = ModDF,
                       family = "binomial",REML=TRUE)
summary(gammPreg_Base$gam)
plot(gammPreg_Base$gam)

# Model with TB - condition interaction
gammPreg_TBint <- gamm4(Preg ~ s(avgCondition, bs="cr", k=5) + s(avgCondition, bs="cr", k=5, by=oTB) + 
                          TB + Dominant + s(Month, bs="cc", k=5),
                        random = ~(1|IndividID),
                        data = ModDF %>% 
                          ungroup() %>% 
                          mutate(oTB = ordered(TB, levels=c("N","Y"))),
                        family = "binomial",
                        REML = TRUE)
summary(gammPreg_TBint$gam)
plot(gammPreg_TBint$gam)

# Model with Dom - condition interaction
gammPreg_DOMint <- gamm4(Preg ~ s(avgCondition, bs="cr", k=5) + s(avgCondition, bs="cr", k=5, by=oDominant) + 
                           TB + Dominant + s(Month,bs="cc",k=5),
                         random = ~(1|IndividID),
                         data = ModDF %>% 
                           ungroup() %>% 
                           mutate(oDominant = ordered(Dominant,levels=c("N","Y"))),
                         family="binomial",
                         REML=TRUE)
summary(gammPreg_DOMint$gam)
plot(gammPreg_DOMint$gam)

# No interactions were significant, so we proceed with the base model
gamPreg_SummaryGAM <- summary(gammPreg_Base$gam)
gamPreg_SummaryMER <- summary(gammPreg_Base$mer)
print(gamPreg_SummaryGAM) ; print(gamPreg_SummaryMER)

# Basic visualisation of smooths
gam.check(gammPreg_Base$gam)
plot(gammPreg_Base$gam)
plot(ggpredict(gammPreg_Base,terms=c("ImpCond")))

# Model Diagnostics
sim_gamPreg <- simulateResiduals(gammPreg_Base$mer,n=1000)
testDispersion(sim_gamPreg) 
testOutliers(sim_gamPreg,type="bootstrap")
plot(sim_gamPreg)

# Generation predictions for paper plot
nd <- expand.grid(avgCondition = seq(-200, 200, 5), 
                  TB = c("Y","N"),
                  Dominant = c("Y","N"),
                  IndividID = 1, Month = 12)

Preg.preds <- cbind(nd,
                    est = predict(gammPreg_Base$gam, newdata = nd, exclude = c("IndividID")),
                    se = predict(gammPreg_Base$gam, newdata = nd, exclude = c("IndividID"), se.fit = TRUE)[[2]])

Preg.preds <- Preg.preds %>%
  mutate(l95ci = inv.logit(est-(1.96*se)), 
        u95ci = inv.logit(est + (1.96*se)), 
        est = inv.logit(est))

ggplot(Preg.preds, aes(x = avgCondition, y = est, group = Dominant, col = Dominant, fill = Dominant)) +
  geom_ribbon(aes(x = avgCondition, ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  theme_bw() +
  theme(panel.grid = element_blank())

p1 <- ggplot(Preg.preds %>% 
               filter(Dominant=="Y" & TB=="N"), aes(x = avgCondition, y = est)) +
  geom_ribbon(aes(ymin = l95ci, ymax= u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  labs(x = "Body condition (g)", y = "Pregnancy Probability") +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Plot the seasonal effect
nd2 <- expand.grid(avgCondition = 0, 
                   TB = "N", 
                   Dominant = c("Y","N"), 
                   IndividID = 1, 
                   Month = seq(1, 12, 0.1))

predMonth <- cbind(nd2,
                   est = predict(gammPreg_Base$gam, newdata = nd2, exclude = c("IndividID")),
                   se = predict(gammPreg_Base$gam, newdata = nd2, exclude = c("IndividID"), se.fit = TRUE)[[2]])

predMonth <- predMonth %>%
  mutate(l95ci = inv.logit(est-(1.96*se)), 
         u95ci = inv.logit(est + (1.96*se)), 
         est = inv.logit(est))

ggplot(predMonth, aes(x = Month, y = est, group = Dominant, col = Dominant, fill = Dominant)) +
  geom_ribbon(aes(ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  theme_bw() +
  theme(panel.grid = element_blank())

#---------------------------------------
# Modelling Litter Success  
#----------------------------------------

LitterDF <- LitterDF %>%  
  ungroup() %>%
  mutate(Dominant=factor(Dominant),
                       TB=factor(TB),
                       FemaleID=factor(FemaleID)) 

# Base model with no interaction terms
gammPregSucc_Base <- gamm4(LitterSuccess ~   s(Cond7Day, bs="cr",k=5) + Dominant + 
                             TB + s(Month, bs="cc", k=5),
                           random = ~(1|FemaleID),
                           data = LitterDF,
                           family = "binomial",
                           REML = TRUE)
summary(gammPregSucc_Base$gam)
plot(gammPregSucc_Base$gam)

# Model with Dom - condition interaction
gammPregSucc_DOMint <- gamm4(LitterSuccess ~ s(Cond7Day, bs="cr", k=5) +   s(Cond7Day, bs="cr", k=5, by=oDominant) + 
                               Dominant + TB + s(Month, bs="cc", k=5),
                             random = ~(1|FemaleID),
                             data = LitterDF %>% 
                               mutate(oDominant = ordered(Dominant, levels = c("N","Y"))),
                             family = "binomial",
                             REML = TRUE)
summary(gammPregSucc_DOMint$gam)
plot(gammPregSucc_DOMint$gam)

# Model with TB - condition interaction
gammPregSucc_TBint <- gamm4(LitterSuccess ~  s(Cond7Day, bs="cr", k=5) +  s(Cond7Day, bs="cr", k=5, by=oTB) + 
                              Dominant + TB + s(Month, bs="cc", k=5),
                            random = ~(1|FemaleID),
                            data = LitterDF   %>% 
                              mutate(oTB = ordered(TB, levels = c("N","Y"))),
                            family = "binomial",
                            REML=TRUE)
summary(gammPregSucc_TBint$gam)
plot(gammPregSucc_TBint$gam)

# No interactions were significant, so we proceed with the base model
gamLitter_SummaryGAM <- summary(gammPregSucc_Base$gam)
gamLitter_SummaryMER <- summary(gammPregSucc_Base$mer)
print(gamLitter_SummaryGAM) ; print(gamLitter_SummaryMER)

# Basic visualisation of smooths
gam.check(gammPregSucc_Base$gam)
plot(gammPregSucc_Base$gam)
plot(ggpredict(gammPregSucc_Base, terms=c("ImpCond")))

# Model Diagnostics
sim_gamLitter <- simulateResiduals(gammPregSucc_Base$mer,n=1000)
testDispersion(sim_gamLitter) 
testOutliers(sim_gamLitter,type="bootstrap")
plot(sim_gamLitter)

# Model predictions
nd <- expand.grid(TB = c("Y","N"), 
                  Dominant = c("Y","N"),
                  FemaleID = 1,
                  Month = 1,
                  Cond7Day = seq(-200, 200, 1))

predPrgSucc <- cbind(nd,
                     est = predict(gammPregSucc_Base$gam, newdata = nd, exclude = c("FemaleID")),
                     se = predict(gammPregSucc_Base$gam, newdata = nd, exclude = c("FemaleID"), se.fit=TRUE)[[2]])

predPrgSucc <- predPrgSucc %>%
  mutate(l95ci = inv.logit(est-(1.96*se)), 
         u95ci = inv.logit(est + (1.96*se)), 
         est = inv.logit(est))

p2 <- ggplot(predPrgSucc %>% 
               filter(Dominant=="Y" & TB=="N"), aes(x = Cond7Day, y = est)) +
  geom_ribbon(aes(ymin = l95ci, ymax= u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  labs(x = "Body condition (g)", y = "Litter Weaned Probability") +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Plot seasonal effect
nd2 <- expand.grid(TB = c("N"),
                   Dominant = c("Y","N"),
                   FemaleID = 1,
                   Month = seq(1,12,0.1),
                   Cond7Day = 0)

predMonth <- cbind(nd2,
                   est = predict(gammPregSucc_Base$gam, newdata = nd2, exclude = c("FemaleID")),
                   se = predict(gammPregSucc_Base$gam, newdata = nd2, exclude = c("FemaleID"), se.fit=TRUE)[[2]])

predMonth <- predMonth %>%
  mutate(l95ci = inv.logit(est-(1.96*se)), 
         u95ci = inv.logit(est + (1.96*se)), 
         est = inv.logit(est))

ggplot(predMonth, aes(x = Month, y = est, group = Dominant, col = Dominant, fill = Dominant)) +
  geom_ribbon(aes(ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  theme_bw() +
  theme(panel.grid = element_blank())

#---------------------------------------
# Modelling mortality 
#----------------------------------------

ModDF %>% 
  ungroup() %>% 
  filter(Diss==0 & is.na(ImpCond)==F) %>%
  summarise(Individ=n_distinct(IndividID),
                    Deaths=sum(Died),
                    n=n())

# Base model with no interactions
gammDeath_Base <- gamm4(Died ~ s(ImpCond, bs="cr", k=5) + TB + Dominant,
                         random=~(1|IndividID),
                         data = ModDF %>% 
                          ungroup() %>% 
                          filter(Diss==0),
                         family="binomial",
                        REML=TRUE)
summary(gammDeath_Base$gam)
plot(gammDeath_Base$gam)

# TB - Condition interaction
gammDeath_TBint <- gamm4(Died ~ s(ImpCond, bs="cr", k=5) + s(ImpCond, bs="cr", by=oTB, k=5) + TB + Dominant,
                        random = ~(1|IndividID),
                        data = ModDF %>% 
                          ungroup() %>% 
                          filter(Diss==0) %>% 
                          mutate(oTB = ordered(TB, levels = c("N","Y"))),
                        family="binomial",
                        REML=TRUE)
summary(gammDeath_TBint$gam)
plot(gammDeath_TBint$gam)

# Dominance - Condition interaction
gammDeath_DOMint <- gamm4(Died ~ s(ImpCond, bs="cr", k=5) + s(ImpCond, bs="cr", by=oDominant, k=5) + TB + Dominant,
                         random = ~(1|IndividID),
                         data = ModDF %>%
                           ungroup() %>% 
                           filter(Diss==0) %>% 
                           mutate(oDominant = ordered(Dominant, levels = c("N","Y"))),
                         family="binomial",
                         REML=TRUE)
summary(gammDeath_DOMint$gam)
plot(gammDeath_DOMint$gam)

# Only the TB interaction is significant, therefore we re-fit the model with the un-ordered factor to get level specific smooths
gammDeath_TBint_FINAL <- gamm4(Died ~ s(ImpCond, bs="cr", by=TB,k=5) + TB + Dominant,
                           random = ~(1|IndividID),
                           data = ModDF %>% 
                             ungroup() %>% 
                             filter(Diss == 0) %>% 
                             mutate(oTB = ordered(TB, levels = c("N","Y"))),
                           family = "binomial",
                           REML = TRUE)
  
# extract model objects to export table values
gamDeath_SummaryGAM <- summary(gammDeath_TBint_FINAL$gam)
gamDeath_SummaryMER <- summary(gammDeath_TBint_FINAL$mer)
print(gamDeath_SummaryGAM); print(gamDeath_SummaryMER)

# Basic visualisation of smooths
gam.check(gammDeath_TBint_FINAL$gam)
plot(gammDeath_TBint_FINAL$gam)
plot(ggpredict(gammDeath_TBint_FINAL, terms=c("ImpCond")))

# Model Diagnostics
sim_gamDeath <- simulateResiduals(gammDeath_TBint_FINAL$mer,n=1000)
testDispersion(sim_gamDeath) # this dispersion is due to the random effect structure, and disappears when Individ ID is dropped from the model
testOutliers(sim_gamDeath,type="bootstrap")
plot(sim_gamDeath)

# Generate new data for predictions
nd <- expand.grid(ImpCond = seq(-200, 200, 5),
                  TB = c("Y", "N"),
                  Dominant = c("Y", "N"),
                  IndividID = 0)

predDeath <- cbind(nd, 
                   est = predict(gammDeath_TBint_FINAL$gam, newdata = nd, exclude = c("IndividID")),
                   se = predict(gammDeath_TBint_FINAL$gam, newdata = nd, exclude = c("IndividID"), se.fit=TRUE)[[2]])

predDeath <- predDeath %>%
  mutate(l95ci = inv.logit(est-(1.96*se)), 
         u95ci = inv.logit(est + (1.96*se)), 
         est = inv.logit(est))

# Plot the predictions for dominants with TB and without
ggplot(predDeath %>% 
         filter(Dominant=="Y"), aes(x = ImpCond, y = est, group = TB, col = TB, fill = TB)) +
  geom_ribbon(aes(ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Final plot for paper (Dominant with no TB)
p3 <- ggplot(predDeath %>% 
               filter(Dominant=="Y" & TB=="N"), aes(x = ImpCond, y = est)) +
  geom_ribbon(aes(ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  labs(c = "Body condition (g)", y = "Mortality Probability") +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid = element_blank())

#---------------------------------------
# Modelling disappearance 
#----------------------------------------

ModDF %>% 
  ungroup() %>% 
  filter(Died == 0 & is.na(ImpCond) == F) %>%
  summarise(Individ = n_distinct(IndividID),
            Diss = sum(Diss),
            n = n())

# Base Model with no interactions
gammDiss_Base <- gamm4(Diss ~ s(ImpCond, bs="cr", k=5) + TB + Dominant + s(Month, bs="cc", k=5),
                       random = ~(1|IndividID),
                       data = ModDF %>% 
                         ungroup() %>% 
                         filter(Died == 0),
                       family = "binomial",
                       REML = TRUE)
summary(gammDiss_Base$gam)

# Model with TB - Condition interaction
gammDiss_TBint <- gamm4(Diss ~ s(ImpCond, bs="cr",k=5) + s(ImpCond, bs="cr", k=5, by=oTB) + 
                          TB + Dominant + s(Month,bs="cc",k=5),
                       random = ~(1|IndividID),
                       data = ModDF %>% 
                         ungroup() %>% 
                         filter(Died==0) %>% 
                         mutate(oTB = ordered(TB, levels = c("N","Y"))),
                       family="binomial",
                       REML=TRUE)
summary(gammDiss_TBint$gam)
plot(gammDiss_TBint$gam)

# Model with Dom - Condition interaction
gammDiss_DOMint <- gamm4(Diss ~ s(ImpCond, bs="cr", k=5) + s(ImpCond, bs="cr", k=5, by=oDominant) + TB + 
                           Dominant + s(Month, bs="cc", k=5),
                        random = ~(1|IndividID),
                        data = ModDF %>%
                          ungroup() %>% 
                          filter(Died==0) %>% 
                          mutate(oDominant = ordered(Dominant, levels = c("N","Y"))),
                        family="binomial",
                        REML=TRUE)
summary(gammDiss_DOMint$gam)
plot(gammDiss_DOMint$gam)

# Model with all interactions included
gammDiss_ALLint <- gamm4(Diss ~ s(ImpCond,bs="cr", by=DomTB, k=5) + Dominant + TB + s(Month, bs="cc", k=5),
                         random = ~(1|IndividID),
                         data = ModDF %>% 
                           ungroup() %>% 
                           filter(Died==0) %>%
                           mutate(DomTB=case_when(
                             Dominant=="Y" & TB=="Y" ~ "DomTB",
                             Dominant=="N" & TB=="Y" ~ "SubTB",
                             Dominant=="Y" & TB=="N" ~ "Dom",
                             Dominant=="N" & TB=="N" ~ "Sub")) %>% 
                           mutate(DomTB = factor(DomTB)) %>%
                           mutate(oDomTB = ordered(DomTB, levels = c("Sub","Dom","DomTB","SubTB"))),
                         family = "binomial",
                         REML = FALSE)

## extract model objects to export table values
gamDiss_SummaryGAM <- summary(gammDiss_ALLint$gam)
gamDiss_SummaryMER <- summary(gammDiss_ALLint$mer)
print(gamDiss_SummaryGAM) ; print(gamDiss_SummaryMER)

## Basic visualisation of smooths
gam.check(gammDiss_ALLint$gam)
plot(gammDiss_ALLint$gam)
plot(ggpredict(gammDiss_ALLint,terms=c("ImpCond")))

## Model Diagnostics
sim_gamDiss <- simulateResiduals(gammDiss_ALLint$mer,n=1000)
testDispersion(sim_gamDiss) ## this dispersion is due to the random effect structure, and disappears when Individ ID is dropped from the model
testOutliers(sim_gamDiss,type="bootstrap")
plot(sim_gamDiss)

# Generate predictions for paper plot
nd <- expand.grid(ImpCond = seq(-200,200,5),
                  TB = c("N","Y"),
                  Dominant = c("Y","N"),
                  IndividID = 1,
                  Month = 1) %>%
          mutate(DomTB = case_when(
            Dominant=="Y" & TB=="Y" ~ "DomTB",
            Dominant=="N" & TB=="Y" ~ "SubTB",
            Dominant=="Y" & TB=="N" ~ "Dom",
            Dominant=="N" & TB=="N" ~ "Sub"))

predDiss <- cbind(nd,
                  est = predict(gammDiss_ALLint$gam, newdata = nd, exclude = c("IndividID")),
                  se = predict(gammDiss_ALLint$gam, newdata = nd, exclude = c("IndividID"), se.fit=TRUE)[[2]])

predDiss <- predDiss %>%
  mutate(l95ci = inv.logit(est-(1.96*se)), 
         u95ci = inv.logit(est + (1.96*se)), 
         est = inv.logit(est))

ggplot(predDiss %>% 
         filter(TB == "N"), aes(x = ImpCond, y = est, group = Dominant, col = Dominant, fill = Dominant)) +
  geom_ribbon(aes(ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  theme_bw() +
  theme(panel.grid = element_blank())

p4 <- ggplot(predDiss %>% 
               filter(TB == "N" & Dominant == "Y"), aes(x=ImpCond, y = est)) +
  geom_ribbon(aes(ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd=1) +
  ylim(0,1) +
  labs(x = "Body condition (g)", y = "Disappearance Probability") +
  theme_bw() +
  theme(panel.grid = element_blank())

#-------------------------------------

# Finalise all the plots for the manuscript

# The plot for the main manuscript: combining p1, p2, p3 and p4 with different aesthetic
preds_final <- rbind(predDeath %>%
                       select(TB, Status = Dominant, Condition = ImpCond, est, l95ci, u95ci) %>%
                       mutate(Analysis = "Mortality", Status = ifelse(Status == "Y","DOM","SUB")),
                     predDiss %>% select(TB, Status = Dominant, Condition = ImpCond, est, l95ci, u95ci) %>%
                       mutate(Analysis = "Disappearance", Status = ifelse(Status == "Y","DOM","SUB")),
                     Preg.preds %>% select(TB, Status = Dominant, Condition = avgCondition, est, l95ci, u95ci) %>% 
                       mutate(Analysis = "Pregnant", Status = ifelse(Status == "Y","DOM","SUB")),
                     predPrgSucc %>% select(TB, Status = Dominant, Condition = Cond7Day, est, l95ci, u95ci) %>% 
                       mutate(Analysis = "Litter Weaned", Status = ifelse(Status == "Y","DOM","SUB"))) %>%
  mutate(Analysis = factor(Analysis, levels = c("Pregnant", "Litter Weaned", "Mortality", "Disappearance")))

# construct the plot
p_final <- ggplot(preds_final %>% filter(TB=="N" & Status=="DOM"),aes(x=Condition, y= est)) +
  geom_ribbon(aes(x=Condition, ymin= l95ci,ymax= u95ci), alpha=0.2) +
  geom_line(lwd=1) +
  xlim(-210,210) +
  labs(x = "Body condition (g)", y = "Probability") +
  scale_x_continuous(breaks=c(-160,-80,0,80,160)) +
  facet_wrap(~Analysis, nrow = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=11, colour = "black"),
        axis.title=element_text(size=13),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(size = 12))

# Make the supplementary plot. We plot it initially and then restructure it. 
preds_final_long <- rbind(preds_final %>% 
                            filter(TB=="N") %>%
                            pivot_longer(cols=c("Status")) %>% 
                            select(-TB),
                          preds_final %>% 
                            filter(Status=="DOM") %>%
                            pivot_longer(cols=c("TB")) %>% 
                            select(-Status)) %>% 
  mutate(value=factor(value,levels=c("DOM","SUB","Y","N")))

# Construct the plot 
si_final <- ggplot(preds_final_long, aes(x = Condition, y = est, group = value, col = value, fill = value)) +
  geom_ribbon(aes(x = Condition, ymin = l95ci, ymax = u95ci), alpha=0.2, col=alpha(0.2)) +
  geom_line(lwd=1) +
  xlim(-210,210) +
  labs(x = "Mean condition (g)", y = "Probability") +
  facet_grid(name ~ Analysis, scales = "fixed") +
  scale_x_continuous(breaks = c(-160,-80,0,80,160)) +
  scale_color_manual(values = c("steelblue3","orange2","mediumorchid2","seagreen3")) +
  scale_fill_manual(values = c("steelblue3","orange2","mediumorchid2","seagreen3")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11, colour= "black"),
        axis.title = element_text(size=13))

# We want to rework slightly to be more informative
si_pred <- preds_final %>% 
  mutate(FacetGroup = case_when(Analysis=="Mortality" & Status == "DOM" ~ "Mortality",
                                Analysis=="Pregnant" & TB == "N" ~ "Pregnant",
                                Analysis=="Litter Weaned" & TB == "N" ~ "Litter Weaned",
                                Analysis=="Disappearance" & Status == "DOM" ~ "Dominant Disappearance",
                                Analysis=="Disappearance" & Status == "SUB" ~ "Subordinate Disappearance",
                                TRUE ~ "Not used"),
         Legend = case_when(TB == "N" & Status == "DOM" ~ "Dominant",
                            TB == "Y" & Status == "DOM"~ "Dominant TB",
                            TB == "N" & Status == "SUB" ~ "Subordinate",
                            TB == "Y" & Status == "SUB" ~ "Subordinate TB")) %>%
  mutate(FacetGroup = ordered(FacetGroup, levels = c("Pregnant", "Litter Weaned", "Mortality",
                                                     "Dominant Disappearance","Subordinate Disappearance")))


# Remake the newer version
si_final <- ggplot(si_pred %>% filter(FacetGroup != "Not used"),
                   aes(x = Condition, y = est, group = Legend, col = Legend, fill = Legend)) +
  geom_ribbon(aes(x=Condition, ymin = l95ci, ymax = u95ci), alpha = 0.2, col = alpha(0.2)) +
  geom_line(lwd = 1) +
  xlim(-210, 210) +
  xlab("Mean Condition (grams)") +
  ylab("Probability") +
  labs(fill="Dominance/TB status",col="Dominance/TB status") +
  facet_wrap(~ FacetGroup,scales="free",ncol=3,nrow=2) +
  scale_x_continuous(breaks = c(-160,-80,0,80,160)) +
  scale_color_manual(values = c("steelblue3","orange2","mediumorchid2","seagreen3"))+
  scale_fill_manual(values = c("steelblue3","orange2","mediumorchid2","seagreen3"))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.title.x = element_text(size = 16, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        strip.background = element_rect(fill = "lightblue"),
        strip.text.x = element_text(size = 12),
        legend.position = c(0.82, 0.28),
        legend.background = element_rect(linetype = 1, colour = 1),
        legend.title=element_text(size = 14),
        legend.text=element_text(size = 12))

#######################   END   ###########################