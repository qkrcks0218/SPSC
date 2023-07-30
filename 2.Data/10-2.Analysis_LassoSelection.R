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


##################################
# Iteration over 48 donor candidates
##################################

Nm <- rep(0,48)

for(m.index in 1:48){
  
  source("0.DataCleaning.R")
  
  m <- m.index*2+2

  gY.Pre <- matrix(0,length(Y1.Pre),m)
  gY.Fit <- bs(Y1.Pre,df=m,
               intercept = F,
               Boundary.knots = gY.Bound)
  gY.Pre[,1:(m)]       <- gY.Fit
  
  overlap <- cbind(1:N,(sapply(1:N,
                               function(bb){
                                 RR1 <- range(Wmat.series[1:T0,bb])
                                 mean(c(as.numeric(RR1[1] <= Wmat.series[T0+1:T1,bb] & 
                                                     Wmat.series[T0+1:T1,bb] <= RR1[2] )))
                               })))
  
  
  N <- dim(Wmat.Pre)[2]
  m <- dim(gY.Pre)[2]
  lengthb <- 1
  
  GMM.Data <- cbind(rbind(gY.Pre,matrix(0,T1,m)),
                    rbind(Wmat.Pre,Wmat.Post),
                    c(rep(0,T0),rep(1,T1)),
                    c(Y1.Pre,Y1.Post))
  
  colnames(GMM.Data) <- (c(sprintf("G%0.4d",1:m),
                           sprintf("W%0.4d",1:N),
                           "A",
                           "Y"))
  
  T0 <- dim(Wmat.Pre)[1]
  T1 <- dim(Wmat.Post)[1]
  Tt <- T0 + T1
  m  <- dim(gY.Pre)[2]
  N  <- dim(Wmat.Pre)[2]
  A  <- rep(c(0,1),c(T0,T1))
  delta   <- 0
  
  ## GMM-1st
  Wmat <- rbind(Wmat.Pre,Wmat.Post)
  Y    <- c(Y1.Pre,Y1.Post)
  gY   <- rbind(gY.Pre, matrix(0,T1,m))
  
  GW <- t( sapply(1:dim(gY.Pre)[2],function(tt){
    apply( t(Wmat.Pre)*matrix(gY.Pre[,tt],N,T0,byrow=T), 1, mean)
  }) )
  
  GY <- ( sapply(1:dim(gY.Pre)[2],function(tt){
    mean(gY.Pre[,tt]*Y1.Pre)
  }) )
  
  CGLMNET <- glmnet::cv.glmnet(GW,GY,intercept=F, nfolds=T0)
  GLMNET  <- glmnet::glmnet(GW,GY,intercept=F,lambda=CGLMNET$lambda.min)
  
  Nm[m.index] <- sum(as.numeric(GLMNET$beta)!=0) # 1,2,3,9,44 with m = 24
  print(Nm[m.index])
}


MAR <- c(3.5, 3.5, 1, 0.5)

par(mar=MAR)

plot(seq(4,98,by=2), Nm,
     xlab="",ylab="",
     xlim=c(0,100),
     ylim=c(2,8),
     pch=19, cex=1)
abline(a=0,b=1/2,col=1,lty=2)

text(15.8, 7.8,labels="dim(g) = 2d",pos=4,col=1)
title(xlab="dim(g) = dimension of g(y)",
      ylab="d = Number of Selected Donors",
      line=2.5)

points(10,5,col=1,pch=19,cex=1.5)

graphics::arrows(12,3.5,
                 10+0.5,5-0.2,
                 length=0.1,col=1)

text(13,3.5,
     "Group 5",
     pos=1,
     col=1)

