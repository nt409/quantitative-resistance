# these functions are called from HostCurveGenerator, choiceOfHost

# exports to Fitting/CSVs:
# Host/{cultivarName}_FrameFull.csv;
# Host/{cultivarName}_YearlyWorstSevs.csv
# AllHostData/allData.csv

library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)

setwd("C:/Users/user/Documents/Python/PhD/Polygenic/polygenic/Run/data")


##--------------------------------------------------------------------------
# setup
setup_data_frame <- function(df_r){
  # Takes relevant columns for untreated, without cultivar mixtures

  df_untreated <- df_r[df_r$Behandling_faktor2=="Ingen bladsvampe-bek",] 

  # df_untreated$stb[df_untreated$stb==100] = 99.9 # values of 100 don't allow beta calculation
  # df_untreated$stb[df_untreated$stb==0] = 0.01  # values of 0 don't allow beta calculation

  relevant_columns <- data.frame(year=df_untreated$HOSTAAR,
                             location=df_untreated$postby,
                             post_code=df_untreated$postnr,
                             stb=df_untreated$stb,
                             cultivar=df_untreated$cultivar)

  output <- relevant_columns[complete.cases(relevant_columns),] %>%
    # exclude mixture (others to exclude?)
    filter(cultivar!='Blanding,vi-hved') %>%  
    # location!='Holeby')
    # exclude Holeby since it is an Island with unusually
    # low disease pressure and unusual climate?
    # then order by year
    arrange(year)

  write.csv(output, file='02_processed/Host/allData.csv')

  return(output)
}
##




UseLocScoreInfo = function(data_in, minLocScore){
  dfLoc <- read.csv("01_raw/Host/LocationScores.csv",stringsAsFactors = F)

  LocScoredf <- subset(dfLoc, select=c('postal.code', 'score_return'))
  colnames(LocScoredf) <- c("postal.code", "LocScore")
  out <- merge(LocScoredf, data_in, by.x="postal.code", by.y="post_code")

  highPressureLocations = out %>% filter(LocScore>minLocScore)

  write.csv(highPressureLocations, file='02_processed/Host/highPressureLocations.csv')

  return(highPressureLocations)

}







# used in CreateDataOutput
ExtractVal <- function(vect,index){
  if(length(index)>0){
    output <- vect[index]
  }else{
    output <- NA
  }
  return(output)
}



CreateDataOutput <- function(ByYearCultLoc, cultivarName){
  # exports:
  # - cultivarName,'/FrameFull.csv
  # - cultivarName,'/YearlyWorstSevs.csv

  OutputDataAll <- ByYearCultLoc %>%
      group_by(year, location) %>% # LocScore
      summarise(
        
        # my Cult info
        # "ExtractVal(vect,index)" adds an NA if cultivarName wasn't grown in this loc/year
        myCultSeverity   = ExtractVal(MeanSeverity, which(cultivar==cultivarName)),
        myCultCount      = ExtractVal(Count, which(cultivar==cultivarName)),

        # Worst Cult info
        WorstSeverity = max(MeanSeverity),

        # out of the worst, use one with maximum reps
        WorstCountAll = paste(Count[which(MeanSeverity == max(MeanSeverity))],collapse=', '),
        WorstCount = max(Count[which(MeanSeverity == max(MeanSeverity))]), 
        
        WorstCultivars = paste(cultivar[intersect( which(MeanSeverity == max(MeanSeverity)), which(Count==WorstCount))],collapse=', '),
        
        # how many cultivars have same severity as Worst Cult
        MultipleCultsWorst = length(which(MeanSeverity == max(MeanSeverity))),
        
        # min of data points in this year/loc for: either my cult or worst cult
        minNum = min(WorstCount,myCultCount)

      ) %>%
      filter(!is.na(myCultSeverity)) %>%
      mutate(Control = 100*(1-myCultSeverity/WorstSeverity))

  dataFiltered <- OutputDataAll %>% filter(WorstSeverity>=5)
  
  notUsed <- OutputDataAll %>% filter(WorstSeverity<5)
  
  print(paste(length(notUsed$Control), 
    "data points filtered out, leaving",
    length(dataFiltered$Control)))

  YearlyData <- dataFiltered %>%
                group_by(year) %>%
                summarise(MeanWorstSeverity = WorstSeverity %*% WorstCount/sum(WorstCount) )
  
  write.csv(dataFiltered, file=paste('03_model_inputs/Host/Varieties/',cultivarName,'/FrameFull.csv',sep=""))
  write.csv(YearlyData,file=paste('03_model_inputs/Host/Varieties/',cultivarName,'/YearlyWorstSevs.csv',sep=""))
  
  return(list("fullFrame"=OutputDataAll, "WorstSevs"=YearlyData))
}











LinearModel <- function(data_input, CultivarToTest){
  filtered_by_cultivar <- data_input %>%
    filter(cultivar %in% CultivarToTest)
  
  if(length(unique(filtered_by_cultivar$year))>1){
    LinModel  <- lm(stb~year,data = filtered_by_cultivar)
    lmSummary <- summary(LinModel)
    rSquared  <- lmSummary$r.squared
    gdt       <- lmSummary$coefficients[2,1]
    pVal      <- lmSummary$coefficients[2,4]
    
    TimeToPredict <- data.frame(
      year = c(min(filtered_by_cultivar$year), max(filtered_by_cultivar$year))
    )
    
    prediction <- predict(LinModel,TimeToPredict)
    
    # print(prediction)
    yMin <- prediction[1]
    yMax <- prediction[2]
    
  }else{
    rSquared <- 0
    gdt      <- 0
    pVal     <- 0
    yMin <- 0
    yMax <- 0
  }
  
  
  return(list(rSquared,gdt,pVal,yMin,yMax))
}
