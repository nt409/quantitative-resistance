# uses functions.R to export Locations/{locUSE}/WorstCults.csv

library(dplyr)

##############################################
setwd("c:/Users/user/Documents/Python/Polygenic/polygenic/Run/src/data_processing")
source('Locations/functions.R')
# get split_best_worst_cults
setwd("C:/Users/user/Documents/Python/Polygenic/polygenic/Run/data")

##############################################
# script config - how to filter data, what to plot

allData = read.csv('02_processed/Host/allData.csv')

allData <- allData %>% 
      filter(!(cultivar %in% c('Blanding,vi-hved'," ","")))
allData$stb[allData$stb==0] = 0.001
allData$stb[allData$stb==100] = 99.999     



##############################################
# run functions

top_proportion <- 6

for(LocUSE in c('Trige', 'Karise', 'Soenderborg')){
  split_data = split_best_worst_cults(allData,top_proportion,LocUSE)
  WorstCults_SingleLoc = split_data$worst
  BestCults_SingleLoc = split_data$best
}








##############################################
# plots


plotSplitByCult = F
plotSevsLocUseByYear = F



if(plotSplitByCult){
  par(mfrow=c(3,1))
  plot(stb~jitter(year),BestCults_SingleLoc,col='green',pch=1)
  
  plot(stb~jitter(year),WorstCults_SingleLoc,col='black',pch=1)
  
  plot(stb~jitter(year),BestCults_SingleLoc,col='green',pch=1)
  points(stb~jitter(year),WorstCults_SingleLoc,col='black',pch=1)
}


if(plotSevsLocUseByYear){ # sevs in LocUSE by year
  par(mfrow=c(1,1))
  myLocData <- allData %>% filter(location==LocUSE)
  numPoints <- myLocData %>%
    group_by(year) %>%
    summarise(count=n())
  numPoints
  # 
  plot(stb~jitter(year),col=cultivar,data=myLocData,pch=20)
}

