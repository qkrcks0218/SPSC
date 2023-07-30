############################
# Recommend to implement parallel computing
# BATCH=1,...,504
# We only run for BATCH=1 for an illustration purpose
############################

BATCH <- 1
PARA  <- expand.grid(1:84,1:6)
gp    <- PARA[BATCH,2]
ITER  <- PARA[BATCH,1]

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

##################################
# Define GMM Data
##################################

gY.Pre <- matrix(0,length(Y1.Pre),m+mT)
gY.Fit <- bs(Y1.Pre,df=m,
             intercept = F,
             Boundary.knots = gY.Bound)
gY.Pre[,1:(m)]       <- gY.Fit
bT <- bs(1:T0,df=(mT),Boundary.knots = c(-2,T0+2))
gY.Pre[,m+1:mT] <- bT

Wmat.series <- Wmat.series[,vc[[gp]]]
Wmat.Pre    <- Wmat.Pre   [,vc[[gp]]]
Wmat.Post   <- Wmat.Post  [,vc[[gp]]]

N <- dim(Wmat.Pre)[2]

GMM.Data <- cbind(rbind(gY.Pre,matrix(0,T1,m+mT)),
                  rbind(Wmat.Pre,Wmat.Post),
                  c(rep(0,T0),rep(1,T1)),
                  c(Y1.Pre,Y1.Post))

colnames(GMM.Data) <- (c(sprintf("G%0.4d",1:(m+mT)),
                         sprintf("W%0.4d",1:N),
                         "A",
                         "Y"))

T0 <- dim(Wmat.Pre)[1]
T1 <- dim(Wmat.Post)[1]
Tt <- T0 + T1
N  <- dim(Wmat.Pre)[2]
A  <- rep(c(0,1),c(T0,T1))
delta   <- 0

##################################
# GMM without ridge regularization
# Point Estimate
##################################

Wmat <- rbind(Wmat.Pre,Wmat.Post)
Y    <- c(Y1.Pre,Y1.Post)
gY   <- rbind(gY.Pre, matrix(0,T1,m+mT))

GW <- t( sapply(1:dim(gY.Pre)[2],function(tt){
  apply( t(Wmat.Pre)*matrix(gY.Pre[,tt],N,T0,byrow=T), 1, mean)
}) )

GY <- ( sapply(1:dim(gY.Pre)[2],function(tt){
  mean(gY.Pre[,tt]*Y1.Pre)
}) )

GMM.gamma.naive <- my.inverse(A = t(GW)%*%(GW) )%*%(t(GW)%*%GY)
GMM.resid.naive <- as.numeric(Y1.Post - (Wmat.Post)%*%GMM.gamma.naive)


##################################
# GMM with ridge regularization
# Point Estimate
##################################

lambda.grid <- seq(-5,0,by=0.2)
lambda.min  <- lambda.grid[ which.min(sapply(lambda.grid,CV.Lambda)) ]

lambda.opt <- optimize(f=CV.Lambda,
                       lower= lambda.min - 2,
                       upper= lambda.min + 2)$minimum

GMM.gamma.naive.lambda <- my.inverse(A = t(GW)%*%(GW),
                                     stab.const = ((T.Pre+T.Post)/(T.Pre))^2*10^(lambda.opt),
                                     adjust = F)%*%(t(GW)%*%GY) 
GMM.resid.naive.lambda <- as.numeric(Y1.Post - (Wmat.Post)%*%GMM.gamma.naive.lambda)

##################################
# Conformal Inference
##################################

PPP <- matrix(1:(T1+1),84,2,byrow=T)
PPP[84,] <- c(166,167)

T.Window <- PPP[ITER,]

CP.CI <- matrix(0,2,7)
colnames(CP.CI) <- c("Time",
                     "SPSC",
                     "SPSC_Reg",
                     "SPSC_LB",
                     "SPSC_UB",
                     "SPSC_Reg_LB",
                     "SPSC_Reg_UB")
CP.CI[,1] <- T0+T.Window
CP.CI[,2] <- as.numeric(Y1.Post - (Wmat.Post)%*%GMM.gamma.naive)[T.Window]
CP.CI[,3] <- as.numeric(Y1.Post - (Wmat.Post)%*%GMM.gamma.naive.lambda)[T.Window]

CP.Beta.1 <- CP.PV.1 <- 
  CP.Beta.2 <- CP.PV.2 <- matrix(0,2,1001)

CPSPSC <- Conformal.Prediction.Data.Time(Wmat.Pre.Input  = Wmat.Pre,
                                         Wmat.Post.Input = Wmat.Post[T.Window,],
                                         Y1.Pre.Input    = Y1.Pre,
                                         Y1.Post.Input   = Y1.Post[T.Window],
                                         Xmat.Pre.Input  = NULL,
                                         Xmat.Post.Input = NULL,
                                         cov.ind.Input   = 0,
                                         center.Input    = CP.CI[,2],
                                         gY.Fit.Input    = gY.Fit,
                                         lambda.Input    = -Inf,
                                         constant.Input  = 0)

CPSPSC.Reg <- Conformal.Prediction.Data.Time(Wmat.Pre.Input  = Wmat.Pre,
                                             Wmat.Post.Input = Wmat.Post[T.Window,],
                                             Y1.Pre.Input    = Y1.Pre,
                                             Y1.Post.Input   = Y1.Post[T.Window],
                                             Xmat.Pre.Input  = NULL,
                                             Xmat.Post.Input = NULL,
                                             cov.ind.Input   = 0,
                                             center.Input    = CP.CI[,3],
                                             gY.Fit.Input    = gY.Fit,
                                             lambda.Input    = lambda.opt,
                                             constant.Input  = ((T.Pre+T.Post)/(T.Pre))^2)


CP.CI[,4:5]  <- CPSPSC$CI
CP.CI[,6:7]  <- CPSPSC.Reg$CI

CP.Beta.1[, 1:1001] <- t(CPSPSC$BG)
CP.Beta.1[, 1:1001] <- t(CPSPSC$BG)
CP.Beta.1[, 1:1001] <- t(CPSPSC$BG)
CP.Beta.2[, 1:1001] <- t(CPSPSC.Reg$BG)
CP.Beta.2[, 1:1001] <- t(CPSPSC.Reg$BG)
CP.Beta.2[, 1:1001] <- t(CPSPSC.Reg$BG)

CP.PV.1[, 1:1001] <- t(CPSPSC$PV)
CP.PV.1[, 1:1001] <- t(CPSPSC$PV)
CP.PV.1[, 1:1001] <- t(CPSPSC$PV)
CP.PV.2[, 1:1001] <- t(CPSPSC.Reg$PV)
CP.PV.2[, 1:1001] <- t(CPSPSC.Reg$PV)
CP.PV.2[, 1:1001] <- t(CPSPSC.Reg$PV)

write.csv(CP.CI,    sprintf("Conformal_Raw/Raw_Time_Group%s_CPCI_%0.3d.csv",           gp,ITER),row.names=F)
write.csv(CP.Beta.1,sprintf("Conformal_Raw/Raw_Time_Group%s_CPBeta_SPSC_%0.3d.csv",    gp,ITER),row.names=F)
write.csv(CP.Beta.2,sprintf("Conformal_Raw/Raw_Time_Group%s_CPBeta_SPSC_Reg_%0.3d.csv",gp,ITER),row.names=F)
write.csv(CP.PV.1,  sprintf("Conformal_Raw/Raw_Time_Group%s_CPPV_SPSC_%0.3d.csv",      gp,ITER),row.names=F)
write.csv(CP.PV.2,  sprintf("Conformal_Raw/Raw_Time_Group%s_CPPV_SPSC_Reg_%0.3d.csv",  gp,ITER),row.names=F)



