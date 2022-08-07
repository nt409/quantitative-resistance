# used in other Locations files


library(kSamples)
library(dplyr)

setwd("C:/Users/user/Documents/Python/Polygenic/polygenic/Run/data")


##############################################



split_best_worst_cults <- function(DataIn, top_proportion, LocUSE){
    
    
  ByYearAndCultThisLoc <- DataIn %>%
    filter(location==LocUSE) %>%
    group_by(year,location,cultivar) %>%
    summarise(
      MeanSeverity   = mean(stb),
      count = n()
    ) %>%
    arrange(year,desc(MeanSeverity))
  
  
  dfTemp <- NULL
  dfTempBest <- NULL
  dfOutWorst  <- NULL
  dfOutBest  <- NULL
  for(y in unique(ByYearAndCultThisLoc$year)){
    
    thisYear <- filter(ByYearAndCultThisLoc,year==y)
  
    # how many bad cultivars?
    # top 1/Nth of dis scores
    # or just worst 2?
    num <- ceiling(length(thisYear$cultivar)/top_proportion)
    # num <- 4
    
    BadCult <- thisYear$cultivar[1:num]
    
    dfTemp  <- filter(DataIn,
                      year==y,
                      cultivar %in% BadCult,
                      location==LocUSE)
    
    dfOutWorst <- rbind(dfOutWorst,dfTemp)
    
    
    dfTempBest  <- filter(DataIn,
                      year==y,
                      !(cultivar %in% BadCult),
                      location==LocUSE)
    
    dfOutBest <- rbind(dfOutBest,dfTempBest)
  }
  
  WorstCultsThisLoc <- dfOutWorst %>%
    arrange(year)
  
  BestCultsThisLoc <- dfOutBest %>%
    arrange(year)

# export
write.csv(WorstCultsThisLoc,file=paste("03_model_inputs/Locations/",LocUSE,"/WorstCults_TopProp",top_proportion,".csv",sep=""))

out = list("best"=BestCultsThisLoc, "worst"=WorstCultsThisLoc)

return(out)
}

