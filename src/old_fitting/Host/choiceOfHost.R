# - exports Host/whichHost/MoreThan8Years
# - says which host in dataset frequently, and also a bit about protection levels

library(dplyr)


setwd("c:/Users/user/Documents/Python/Polygenic/polygenic/Run/src/fitting")


source('Host/functions.R')
# get LinearModel

setwd("C:/Users/user/Documents/Python/Polygenic/polygenic/Run/data")





highLSData = read.csv('02_processed/Host/DataOutput/highPressureLocations.csv')


byCult <- highLSData %>%
  group_by(cultivar) %>%
  summarise(
    dataPoints=n(),
    minYear  = min(year),
    maxYear  = max(year),
    numYears  = length(unique(year)),
    rSquared  = as.numeric(LinearModel(highLSData,cultivar)[1]),
    gdt       = as.numeric(LinearModel(highLSData,cultivar)[2]),
    pVal      = as.numeric(LinearModel(highLSData,cultivar)[3]),
    yMin= as.numeric(LinearModel(highLSData,cultivar)[4]),
    yMax= as.numeric(LinearModel(highLSData,cultivar)[5])
  ) %>%
  arrange(desc(numYears))


dataUse <- filter(byCult, numYears>=8, maxYear>=2010) %>% 
    arrange(yMin)



write.csv(dataUse,  file='04_justification/Host/whichHost/MoreThan8Years.csv')
