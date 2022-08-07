# exports YR/YieldRel.csv, YR/linearModelPars.csv

# This script finds our yield vs disease severity relationship

library(tidyr)
library(dplyr)
library(ggplot2)

setwd("C:/Users/user/Documents/Python/Polygenic/polygenic/Run/data")

df_raw_data <- read.csv("01_raw/YR/YR_in_sonderborg.csv",stringsAsFactors=F)

##
df <- df_raw_data
names(df)[9] <- 'STB' 
# check what time point this is for? ---
# 27.06 so around 75 (which was 25.06)



SoendDf_to_fit = data.frame("Yield" = df$yield/10, "STB"= df$STB/100)
SoendDf_to_fit = SoendDf_to_fit[complete.cases(SoendDf_to_fit),]

model <- lm(Yield~STB, SoendDf_to_fit)

model

gdt <- (model$coefficients[2])/model$coefficients[1] # Yield gradient
intercept <- model$coefficients[1]

linPars = data.frame("gdt" = gdt, "intercept"= intercept)

write.csv(SoendDf_to_fit,file = "03_model_inputs/YR/YieldRel.csv")
write.csv(linPars,file = "03_model_inputs/YR/linearModelPars.csv")
  
