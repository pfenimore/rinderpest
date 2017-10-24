#EpiGrid global.R
# GPL v3

library(shiny)
library(xtable)
library(raster)
library(rgeos)  # Needed for Pakistan Roads
library(leaflet)
library(ggplot2)
library(geosphere)  # needed for distMeeus and distm (Pakistan roads)
library(shinyjs)
library(geojsonio)  # Used for as.json

releaseLoc <- "Public"        #******

diffeqs <- 1 # coupled diffeqs
animalCut <- 0.2 # If population in E is less than this don't progress.

DiseaseList <- c("Rinderpest")
SizeList <- c("Small: 0.8 x 0.8 degrees", "Medium: 2 x 2 degrees", "Large: 3 x 3 degrees", "Small Rect: 1 x 1.5 degrees", "Large Rect: 3 x 4 degrees")

# Parameters for starting the simulation  - used in server.r and ui.r
weeks_to_sim_BasicStart <- 5 
AggFact_BasicStart <- 1.0 
srmc_BasicStart <- 1.0  # mitigation that reduces transmission
lrmc_BasicStart <- 0.5  # Reduction in fraction of spread that is directional (long distance)
dl_BasicStart <- 1.3    # Decay of spread function

#Geographic Parameters
region_lengths_Basic <- 3
region_lengths_Medium <- 2
region_lengths_HighRes <- 0.8
lat_lengths_LargeRect <- 3
long_lengths_LargeRect <- 4
lat_lengths_SmallRect <- 1
long_lengths_SmallRect <- 1.5
longval <- 23.1   # 88.6 #                                                                                         #******# Kenya: 36.7  
longmin <- longval - region_lengths_HighRes/2
longmax <- longval + region_lengths_HighRes/2
latval <- 13.6    # 25.2#                                                                                       #****** # Kenya:-1.1   
latmin <- latval - region_lengths_HighRes/2
latmax <- latval + region_lengths_HighRes/2
corners <- list(longmin=longmin, longmax=longmax, latmin=latmin, latmax=latmax)
#...for Pakistan Rinderpest simulation
PakLoniC <- 74.56
PakLatC <- 35.87
PakHistCatRed <- 1.3  # Pakistan Historical Cattle Reduction                                                     #******
#FAO dairy report says Pakistans cattle population increased by ~44% between 1996 and 2006 (pg 1)           
#Using a slightly lower number for reduction, becaues Yakmos (cattle, Yak cross) were a significant part of the population

# Parameters the same for all diseases
beta_max <- 2.5  # All diseases.

################ Rinderpest disease Parameters -  rates are per day #############
kEI_Rind  <-  1/5.6  # Time to virus excretion

# Total time in I and H should be virus secretion, 6-8 days
kIH_Rind  <- 0.33  #1/3  There is a weird error with the sliders, have to reset slider to 0.33, when use 1/3 #               #******   
kIR_Rind  <- 0

kHD_Rind  <- 0.8/3 
kHR_Rind  <- 0.2/3

# Infectious period for untreated disease:
# 1/(kIH+kIR)*kIR/(kIH+kIR) + (1/(kIH+kIR)+1/(kHR+kHD))*kIH/(kIR+kIH)

# Mortality Fraction (i.e. fraction of seriously ill that die) for untreated disease (mortRate):
# kHD/(kHD+kHR)*kIH/(kIH+kIR)

Rr_IIt_Rind <- 0                                                                                                  #******
Rr_HHt_Rind <- 0

kItHt_Rind  <- 0  # not used
kItR_Rind   <- 0  # not used

kHtD_Rind <- 1 # not used, but sum must be greater than 0
kHtR_Rind <- 1 # not used

beta_Rind <- 1.1                                                                                                  #******


###################### Vaccine Parameters #####################
Vbase_BasicStart <- 0 # Baseline vaccination level
Vdelay_BasicStart <- 7  #days
vacradmin_BasicStart <- 5  #km
vacrad_BasicStart <- 14 #km
vac_TimeToImmunity_Rind <- 5 #7 #days

print("Loading cattle population data")
cattle <- raster("cattle/totcor/glbctd1t0503m")

print("Loading cattle population map")
load("Cmapadm.RData")