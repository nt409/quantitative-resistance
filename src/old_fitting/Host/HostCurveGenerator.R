# This script calls functions.R, which exports to Fitting/CSVs:
# - Host/{cultivarName}/FrameFull.csv;
# - Host/{cultivarName}/YearlyWorstSevs.csv
# - Host/allData.csv

library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)

setwd("c:/Users/user/Documents/Python/PhD/Polygenic/polygenic/Run/src/fitting/Host")
source('functions.R')
setwd("C:/Users/user/Documents/Python/PhD/Polygenic/polygenic/Run/data")

df_raw_data <- read.csv("01_raw/host_trials.csv", stringsAsFactors=F)


allData <- setup_data_frame(df_raw_data)

# OLD (but still right)
##

# minLocScore <- 5

# # 8 hosts from MoreThan8Years.csv
# hosts_to_test <- c(
#     'Mariboss',
#     'Hereford',
#     'Jensen',
#     'Tuareg',
#     'Ambition',
#     'JBAsano',
#     'Frument',
#     'KWSDacanto'
# )


# for(host_use in hosts_to_test){

#     allData <- setup_data_frame(df_raw_data)

#     # filter so that high pressure disease region
#     highLocScoreOnly = UseLocScoreInfo(allData, minLocScore)

#     ByYearCultLoc <- highLocScoreOnly %>%
#                 group_by(year, location, cultivar) %>% # LocScore
#                 summarise(
#                         MeanSeverity   = mean(stb),
#                         Count = n(),
#                         ) %>%
#                 arrange(year,MeanSeverity,Count)

#     CreateDataOutput(ByYearCultLoc, host_use)

# }
