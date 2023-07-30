##################################
# Functions for the SPSC approach
##################################

source("0.Function_SPSC_Data.R")

ATT.Type <- "constant"
lengthb  <- 1  
Ypos     <- c(1,3)
gY.Bound <- c(6.5,7.1) 

##################################
# Library
##################################

library(splines)
library(MASS)
library(Synth)

EFF <- matrix(0,6,2) # Effect Estimate

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
  
  source("0.DataCleaning.R")
  
  Wmat.series <- Wmat.series[,vc[[gp]]]
  Wmat.Pre    <- Wmat.Pre   [,vc[[gp]]]
  Wmat.Post   <- Wmat.Post  [,vc[[gp]]]
  
  N <- dim(Wmat.Pre)[2]
  m <- dim(gY.Pre)[2]
  if(ATT.Type=="constant"){
    lengthb <- 1
  } else {
    lengthb <- dim(Post.Time.Basis)[2]
  }
  
  ##################################
  # Make Synth Object
  ##################################
  
  Abadie.SC.Data <- data.frame(
    idvar = rep(c(0,1:N),each=(Tt)),
    idvarchr = as.character(rep(c(0,1:N),each=(Tt))),
    timevar = rep(1:Tt,N+1),
    outcomevar = c(Y1.series,as.vector(Wmat.series)),
    X1 = rep(1,(N+1)*Tt) + rnorm((N+1)*Tt)*0.00001
  )
  
  dataprep.out <- dataprep( foo=Abadie.SC.Data,
                            predictors="outcomevar",
                            predictors.op="mean",
                            dependent="outcomevar",
                            unit.variable="idvar",
                            time.variable="timevar",
                            treatment.identifier=0,
                            controls.identifier=1:N,
                            time.predictors.prior=1:T0,
                            unit.names.variable="idvarchr",
                            time.optimize.ssr=1:T0,
                            time.plot=1:Tt)
  
  synth.out <- synth(dataprep.out)
  path.plot(dataprep.res = dataprep.out,
            synth.res=synth.out)
  
  MSE <- function(WW){
    WW2 <- c(1-sum(WW),WW)
    sum( ((Y1.series - Wmat.series%*%WW2)[1:T0])^2 + sum(as.numeric(abs(WW2-0.5)>0.5)*100000) )
  }

  MSE2  <- rep(0,100)
  Wopt2 <- matrix(0,100,N)

  Wopt2[1,-1] <- optim(synth.out$solution.w[-1], MSE)$par
  Wopt2[1, 1] <- 1-sum(Wopt2[1,-1])
  MSE2[1]     <- MSE(Wopt2[1,-1])

  for(riter in 2:100){
    start <- runif(N)
    start <- start/sum(start)
    Wopt2[riter,-1] <- optim(start[-1], MSE)$par
    Wopt2[riter, 1] <- 1-sum(Wopt2[riter,-1])
    MSE2[riter]     <- MSE(Wopt2[riter,-1])
  }

  Wopt <- Wopt2[which.min(MSE2),]
  
  # Wopt <- synth.out$solution.w
  
  EFF[gp,1] <- mean( (Y1.series - Wmat.series%*%Wopt)[T0+1:T1] )
  
  Abadie.SC.Data.Placebo <- data.frame(
    idvar = rep(c(0,1:N),each=(T0)),
    idvarchr = as.character(rep(c(0,1:N),each=(T0))),
    timevar = rep(1:T0,N+1),
    outcomevar = c(Y1.series[1:T0],as.vector(Wmat.series[1:T0,])),
    X1 = rep(1,(N+1)*T0) + rnorm((N+1)*T0)*0.00001
  )
  
  T0.Placebo <- 181
  T1.Placebo <- 36
  
  dataprep.out.Placebo <- dataprep( foo=Abadie.SC.Data.Placebo,
                                    predictors="outcomevar",
                                    predictors.op="mean",
                                    dependent="outcomevar",
                                    unit.variable="idvar",
                                    time.variable="timevar",
                                    treatment.identifier=0,
                                    controls.identifier=1:N,
                                    time.predictors.prior=1:T0.Placebo,
                                    unit.names.variable="idvarchr",
                                    time.optimize.ssr=1:T0.Placebo,
                                    time.plot=1:T0)
  
  synth.out.Placebo <- synth(dataprep.out.Placebo)
  path.plot(dataprep.res = dataprep.out.Placebo,
            synth.res=synth.out.Placebo)
  
  
  MSE <- function(WW){
    WW2 <- c(1-sum(WW),WW)
    sum( ((Y1.series[T0.Placebo+1:T1.Placebo] -
             Wmat.series[T0.Placebo+1:T1.Placebo,]%*%WW2))^2 +
           sum(as.numeric(abs(WW2-0.5)>0.5)*100000) )
  }

  MSE2  <- rep(0,100)
  Wopt2 <- matrix(0,100,N)

  Wopt2[1,-1] <- optim(synth.out.Placebo$solution.w[-1], MSE)$par
  Wopt2[1, 1] <- 1-sum(Wopt2[1,-1])
  MSE2[1]     <- MSE(Wopt2[1,-1])

  for(riter in 2:100){
    start <- runif(N)
    start <- start/sum(start)
    Wopt2[riter,-1] <- optim(start[-1], MSE)$par
    Wopt2[riter, 1] <- 1-sum(Wopt2[riter,-1])
    MSE2[riter]     <- MSE(Wopt2[riter,-1])
  }
  
  # Wopt <- synth.out.Placebo$solution.w
  
  
  EFF[gp,2] <- mean( (Y1.series[T0.Placebo+1:T1.Placebo] - 
                        Wmat.series[T0.Placebo+1:T1.Placebo,]%*%Wopt) )
  
  
  
}

EFF 

# [1,] -0.1197775  0.771310895
# [2,] -1.0157451 -0.049089549
# [3,] -0.8568228  0.053757550
# [4,] -0.8627475 -0.008617287
# [5,]  0.7041317  1.576279078
# [6,] -0.8776551 -0.012510339
