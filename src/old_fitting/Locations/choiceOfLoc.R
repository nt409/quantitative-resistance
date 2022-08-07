# exports Locations/whichLoc/MoreThan1000Points.csv

##############################################
library(dplyr)

setwd("C:/Users/user/Documents/Python/Polygenic/polygenic/Run/data")

allData = read.csv('02_processed/Host/allData.csv') # just tidied rawdata


mainLocs = allData %>% 
    group_by(location) %>%
    summarise(count=n(), 
        meanSev = mean(stb),
        nYears=length(unique(year)),
        perYearAv=count/nYears
    ) %>% 
    filter(count>1000) %>%
    arrange(desc(meanSev))


mainLocs

write.csv(mainLocs,
    file='04_justification/Locations/whichLoc/MoreThan1000Points.csv'
)

# pick highest, lowest, and one in middle:
# Soenderborg 41.6
# Trige 7.28
# Karise 23.0
