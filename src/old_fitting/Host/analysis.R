# show that location, cultivar, year all significant effects

library(MASS)
library(dplyr)

setwd("C:/Users/user/Documents/Python/Polygenic/polygenic/Run/data")
allData = read.csv('02_processed/Host/DataOutput/allData.csv')
dataUse = read.csv('04_justification/Host/whichHost/MoreThan8Years.csv')





regression_df <- allData %>% filter(cultivar %in% dataUse$cultivar)


# linear model comparison

LM1 <- lm(formula = stb ~ year + location + cultivar,data = regression_df)
LM2 <- lm(formula = stb ~ year + cultivar, data = regression_df)
LM3 <- lm(formula = stb ~ year + location, data = regression_df)

summary(LM1)

anova(LM1,LM2)
anova(LM1,LM3)



# with interactions

LM_all_interactions  <- lm(formula = stb ~ year + location + cultivar
                   + year:location + cultivar:location + year:cultivar
                   + year:location:cultivar, data = regression_df)

par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(LM_all_interactions)

LM_2ndOrd_Interactions  <- lm(formula = stb ~ year + location + cultivar
                   + year:location + cultivar:location + year:cultivar
                   , data = regression_df)

plot(LM_2ndOrd_Interactions)





# AIC

step <- stepAIC(LM1, direction = 'both')
step$anova

step <- stepAIC(LM_all_interactions, direction = 'both')
step$anova

step <- stepAIC(LM_2ndOrd_Interactions, direction = 'both')
step$anova

