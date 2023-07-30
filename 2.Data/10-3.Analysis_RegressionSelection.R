rm(list=ls())

##################################
# Functions for the SPSC approach
##################################

source("0.Function_SPSC_Data.R")

ATT.Type <- "constant"
lengthb  <- 1  
Ypos     <- c(1,3)
gY.Bound <- c(6.5,7.1) 

##################################
# Library + Data Cleaning
##################################

library(splines)
library(MASS)
source("0.DataCleaning.R")

CHECK <- 0

while(CHECK==0){
  
  pvalue.work <- rep(0,dim(Wmat.series.work)[2])
  for(witer in 1:dim(Wmat.series.work)[2]){
    pvalue.work[witer] <- summary(lm(Wmat.series.work[1:T0,witer]~0+Y1.Pre+Wmat.series.work[1:T0,-witer]))$coefficients[1,4]
  }
  
  if(max(pvalue.work)>0.05 & length(pvalue.work)>2){
    Wmat.series.work <- Wmat.series.work[,-which.max(pvalue.work)]
    CHECK <- 0
  } else {
    CHECK <- 1
  }
}
sapply(1:24,function(vv){which( apply((Wmat.series.work[,vv] - Wmat.series)^2,2,mean)==0 )})
