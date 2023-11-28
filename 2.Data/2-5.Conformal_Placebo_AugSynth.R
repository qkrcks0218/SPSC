##################################
# Functions for the SPSC approach
##################################
rm(list=ls())
source("0.Function_SPSC_Data.R")

ATT.Type <- "constant"
lengthb  <- 1  
Ypos     <- c(1,3)
gY.Bound <- c(6.5,7.1) 

##################################
# Library
##################################



for(gp in 1:6){
  
  ##################################
  # Hyperparameters
  ##################################
  
  mgrid <- c(24,24,24,24,10,48)
  m <- mgrid[gp]
  Num.Boot  <- 10
  Boot.valid.thr  <- 10000
  Boot.Scale <- c(0.8,0.9,1,1.1,1.2)
  
  gT.Bound <- c(-25,25)+c(0,167)
  stab.const <- 10^(-8)
  Boot.Factor <- 1
  bs.intercept <- F
  
  T0 <- 217
  mT <- round( (T0^(1/3)) ) 
  
  ##################################
  # Donor Choice: run 10-1-10-3 R files
  ##################################
  
  vc <- list()
  vc[[1]] <- c(4,10,11,12,14,17,18,24,26,34,43,45)
  vc[[2]] <- c(1,6,7,8,22,25,28,29,32,35,38,44)
  vc[[3]] <- c(2,3,16,19,20,21,27,31,36,37,40,41,46)
  vc[[4]] <- c(5,9,13,23,30,33,39,42,47,48,49)
  vc[[5]] <- c(1,2,5,9,18)
  vc[[6]] <- c(1,2,3,4,5,6,8,9,12,13,14,19,20,23,27,33,34,37,38,40,41,42,44,47)
  
  ##################################
  # Data Cleaning
  ##################################
  
  source("0.DataCleaning_Placebo.R")
  
  ##################################
  # SCPI
  ##################################
  
  library(augsynth)
  library(dplyr)
  
  Wmat.series <- Wmat.series[,vc[[gp]]]
  Wmat.Pre    <- Wmat.Pre   [,vc[[gp]]]
  Wmat.Post   <- Wmat.Post  [,vc[[gp]]]
  Tt <- T0 + T1
  Y1.series   <- Y1.series[1:Tt]
  
  N <- dim(Wmat.Pre)[2]
  
  
  ASC.Data.ATT <-
    data.frame(id=rep(c(1,1+1:N),each=(T0+T1)),
               time=rep(1:(T0+T1),N+1),
               outcome=c(Y1.series,as.vector(Wmat.series)),
               treatment=c(rep(c(0,1),c(T0,T1)),rep(0,(T0+T1)*N)))

  ASC.Fit <- augsynth(outcome~treatment,
                      unit=id,
                      time=time,
                      data=ASC.Data.ATT,
                      progfunc="Ridge",
                      scm=T)

  source("AugSynth_Fix.R")

  ASC.ATT <- summary.augsynth(ASC.Fit)
  
  save.image(file=sprintf("Conformal/Placebo_GP%d_Data_ASC.RData",gp))
  
  
  
}

