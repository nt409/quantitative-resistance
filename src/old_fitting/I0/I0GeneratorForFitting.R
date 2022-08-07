# writes: I0_Calculated.csv

# From Rose: NB:
# 001, 002 and 004 have treatments
# 
# A: untreated
# 
# B: Three treatments
# 
# C: Two treatments
# 
# D: One treatment
# 
# (i.e. column treat)
# 
# stb L2 vs75 28.06 = FINAL SEVERITY!
# 
# yield = AS SAYS
# 
# 4 replicates per treatment...
# 
# ==== IN 19350 ====
# 
# CULTIVAR 1 and 2 ARE MIXTURES; BENCHMARK IS 3 (I.E. ALMOST BROKEN); HEREFORD IS 5 (I.E. TOTALLY GONE); INFORMER IS 7 (I.E. BRAND NEW FULLY RESISTANT)
# 
# Note treatments are coded differently!!!!
# 
# Column 25.06.19 vs 75 stb l2 is the equivalent to the other files
# Column d39 is a sum of both leaf 1 and leaf 2 end pressure (i.e. 0 - 200)

library(tidyr)
library(dplyr)
library(ggplot2)

##--------------------------------------------------------------------------

setwd("C:/Users/user/Documents/Python/Polygenic/polygenic/Run/data")

# t_0 <- 1212
t_0 <- 1456
t_1 <- 2100 # 61-65 2066 to ?  2066 = GS61
t_2 <- 2350 # 71
t_3 <- 2500 # 75 < 2900 = 87   

RawData <- read.csv("01_raw/I0/I0_in.csv", stringsAsFactors=F)
# View(RawData)

fit <- lm(L2~treatment, RawData)

filteredRawdf <- RawData %>% filter(treatment=="a")

# 3: Benchmark almost broken, 5: Hereford
# A


I0_df <- data.frame(cultivar = filteredRawdf$cultivar,
                           I1 = 0.01*filteredRawdf$L2,
                           # I1: 06.06 vs 61-65 L3 (leaf 3)L2
                           I2 = 0.01*filteredRawdf$X17.06.Vs.71.stb.l2,
                           # I2: 17.06
                           I3 = 0.01*filteredRawdf$X25.06.19.vs.75.stb.l2
                           # I3: X25.06.19.vs.75.stb.l2
)


I0_df <- I0_df %>%
  filter(I1>0,
         I2>0,
         I3>0) 

I0_df = I0_df %>% mutate(I1logit = log(I1/(1-I1)),
                 I2logit = log(I2/(1-I2)),
                 I3logit = log(I3/(1-I3)))

I0_df$t0 <- t_0
I0_df$t1 <- t_1
I0_df$t2 <- t_2
I0_df$t3 <- t_3




find_logistic_params_diffBeta <- function(I1,I2,I3,t0,t1,t2,t3){
  
  I0logit = NULL
  beta = NULL
  I0s = NULL

  for(i in 1:length(I2)){

    # fitting to inverse logit transform
    new_df = data.frame(I = c(I1[i], I2[i], I3[i]), t = c(t1[i],t2[i],t3[i]))
    
    
    # only use 2 points:
    # new_df = data.frame(I = c(I1[i],I3[i]), t = c(t1[i],t3[i]))

    
    fit = lm(I~t,data=new_df)

    time_to_predict = data.frame(t=c(t0[1]))
    
    I0lgt = predict(fit,time_to_predict)

    I0 = exp(I0lgt)/(1+exp(I0lgt))
    
    I0s <- c(I0s,I0)
    I0logit <- c(I0logit,I0lgt)
    beta <- c(beta,fit$coefficients[2])

  }

  list_out <- list("beta" = beta, "I0logit" = I0logit, "I0s" = I0s)
  return(list_out)
}


fitDiffBetaEachRun <- function(I0_in){ # diff beta per run
  I0_out <- I0_in %>%
    mutate(
      beta   = find_logistic_params_diffBeta(I1logit,I2logit,I3logit,t0,t1,t2,t3)$beta,
      I0logit= find_logistic_params_diffBeta(I1logit,I2logit,I3logit,t0,t1,t2,t3)$I0logit,
      I0s    = find_logistic_params_diffBeta(I1logit,I2logit,I3logit,t0,t1,t2,t3)$I0s
    )
    #  %>% filter(I0<0.1) # rule out anomaly
    
  # use median (highly skewed data)
  I0_out$median_value <- median(I0_out$I0s)
  I0_out$mean_value <- mean(I0_out$I0s)


  # CSV
  write.csv(I0_out, file = "03_model_inputs/I0/I0_Calculated.csv")

  return(I0_out)
}

#--------------------------------------------------------------------
# RUN
I0_out <- fitDiffBetaEachRun(I0_df)








#--------------------------------------------------------------------
#plots

plot_i0_fit = function(){
  graphics.off()
  dev.new()
  plot(  I1logit~t1,data=I0_out,xlim=c(1000,2600),ylim=c(-17,5))
  points(I2logit~t2,data=I0_out)
  points(I3logit~t3,data=I0_out)
  points(I0logit~t0,data=I0_out)

  for(i in 1:length(I0_out$I0s)){
    t = c(1212,2500)
    y = I0_out$I0logit[i] + I0_out$beta[i]*(t-1212)
    lines(y~t)
  }
}










# I0 plot
if(F){
  ggplot(I0_df,aes(t1,I1,colour='I1')) +
    geom_point() +
    scale_x_continuous(name='t (degree-days)')+
    scale_y_continuous(name='Severity')+
    geom_point(data=I0_df,mapping=aes(t3,I3,colour='I3')) +
    geom_point(data=I0_df,mapping=aes(t2,I2,colour='I1a')) +
    geom_point(data=I0_df,mapping=aes(t0,I0,colour='I0')) +
    theme(legend.position='none')


  ggplot(data=I0_df, aes(x=1, y=I0)) + 
    geom_violin() + 
    stat_summary(fun.y=median, geom="point", size=2, color="red") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
}
