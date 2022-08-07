  """
  This script gets data from FrankDataIn and transforms to 
  02_processed/Fungicide/ProthioconazoleControl.csv

  """
# writes Fungicide/FrankControl.csv

library(tidyr)
library(dplyr)
library(ggplot2)



setwd("C:/Users/user/Documents/Python/PhD/Polygenic/polygenic/Run/data")
df <- read.csv("01_raw/fungicide_frank.csv",stringsAsFactors=F)

RD2v <- NULL
k2v <- NULL
Srelv <- NULL

for(ii in 1:length(df$RD)){
  
  rd <- substr(df$RD[ii], 1,5)
  rd <- as.numeric(rd)
  RD2v <- c(RD2v,rd)
  
  len <- nchar(df$k[ii])
  kk <- substr(df$k[ii], 1, len-5)
  kk <- as.numeric(kk)
  k2v <- c(k2v,kk)
  
  xx <- 1 - rd + rd * exp(-kk)
  
  Srelv <- c(Srelv,xx)
  
}

df$RD2 <- RD2v
df$k2 <- k2v
df$Srel <- Srelv
df$Control <- 100 * (1-Srelv)

df$number_of_trials <- c(4,1,4,6,6,5,4,6,5,4,5,5,6,6,4,6,6,5)

df$minNum <- 3 * df$number_of_trials

plot(Control~year, df)

df_out <- subset(df, select=-c(RD,k))

write.csv(df_out, file="02_processed/proth_control_raw.csv")
