# Linking climate variability to demography in meerkats

This repository contains R scripts for analyzing the relationship between climate variability and meerkat demography from our Ecological Monographs (2025) paper:
 *Linking climate variability to demography in cooperatively breeding meerkats*

The aim of our study was to develop a mechanistic understanding of how variation in rainfall and temperature affect meerkat demography. Specifically, we:
- Identified the critical period during which rainfall influences vegetation productivity at the field site.
- Investigated how vegetation productivity and temperature influence meerkat foraging performance.
- Examined how these climatic variables impact daily weight fluctuations and overall body condition.
- Analyzed conditions under which meerkats fail to accumulate enough mass during foraging to offset overnight losses.
- Explored the influence of body condition on reproduction and survival.
- Projected future temperature trends and their implications for meerkat populations up to the mid and end of the century.

The R code is placed within 6 scripts, each of which carries out a specific analysis or model contained within the paper:  

(1) `Part1_RainfallNDVI_Climwin.R` - Runs a climate window analysis to identify the lagged period over which rainfall affects NDVI

(2) `Part2_ForagingPerformance.R` - Analyses the effects of temperature and NDVI on foraging effort and foraging success, using a historical dataset of focal observations. 

(3) `Part3_DWGandOML.R` - Analyses the effects of temperature and NDVI on daily mass gain. By also exploring the relationship between daily mass gain and overnight mass loss, we also identidy environmental conditions under which the daily mass gains are often insufficient to offest overnight mass losses; at which point their body condition will decline over successive days. 

(4) `Part4_DistributedLagModels.R` - Analyses the lagged, interactive effects of temperature and NDVI on body condition (the relative mass of individuals compared to their esimated asymptotic mass). 

(5) `Part5_VitalRates.R` - Analyses the effects of changes in body condition on the reproductive performance and survival (or disappearance) of individuals. 

(6) `Part6_FutureWarmingKalahari.R` - Projects future temperatures across the meerkat distribution for the middle (2050-2060) and end of the century (2090-2100). Specifically, we estimate the mean number of days per year exceeding 35°C, and 37°C. Above 35°C it becomes increasingly difficult for adult meerkats to maintain their body condition, and at 37°C the probability that adult meerkats lose body condition over successive days is very high (daily mass gains < overnight mass losses). 

The scripts were produced by Dr Jack Thorley and Dr Christopher Duncan. All data sets are provided in the repository. 
