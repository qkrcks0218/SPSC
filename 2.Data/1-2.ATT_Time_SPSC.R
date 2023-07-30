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
  
  ##################################
  # GMM without ridge regularization
  # Point Estimate
  ##################################
  
  Wmat <- rbind(Wmat.Pre,Wmat.Post)
  Y    <- c(Y1.Pre,Y1.Post)
  gY   <- rbind(gY.Pre, matrix(0,T1,m))
  
  GW <- t( sapply(1:dim(gY.Pre)[2],function(tt){
    apply( t(Wmat.Pre)*matrix(gY.Pre[,tt],N,T0,byrow=T), 1, mean)
  }) )
  
  GY <- ( sapply(1:dim(gY.Pre)[2],function(tt){
    mean(gY.Pre[,tt]*Y1.Pre)
  }) )
  
  GMM.gamma.naive <- my.inverse(A = t(GW)%*%(GW) )%*%(t(GW)%*%GY)
  GMM.resid.naive <- as.numeric(Y1.Post - (Wmat.Post)%*%GMM.gamma.naive)
  
  GMM.beta.naive      <- mean(GMM.resid.naive)
  GMM.ATT.naive       <- Post.Time.Basis*GMM.beta.naive
  GMM.Simple.Coef <- c(GMM.beta.naive,GMM.gamma.naive)
  
  ##################################
  # GMM without ridge regularization
  # HAC variance
  ##################################
  
  GRAD <- GMM.Ft.Grad(GMM.Simple.Coef,GMM.Data)
  
  GMM.Simple.VAR.HAC <- Meat.HAC(Res.Mat = GMM.Ft(GMM.Simple.Coef, GMM.Data),
                                 GRAD = GRAD,
                                 beta.pos = 1:lengthb,
                                 bw.type="auto")
  
  SGRAD <- svd(GRAD)
  
  GMM.Simple.Var <- (SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)) %*% GMM.Simple.VAR.HAC$matrix %*% t((SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)))/(T0+T1)
  
  ##################################
  # GMM without ridge regularization
  # Block bootstrap
  ##################################
  
  GMM.Boot <- BOOT(round(GMM.Simple.VAR.HAC$bw.B*Boot.Scale),
                   Num.Boot,
                   T0,
                   T1,
                   Boot.Factor=Boot.Factor,
                   type="SPSC")
  
  boot.pos <- which.max(sapply(1:length(GMM.Boot),function(vv){
    
    NORM <- GMM.Boot[[vv]]^2
    valid.pos <- which( abs(log(NORM) - median(log(NORM))) <= Boot.valid.thr*IQR(log(NORM)) )
    return((var(GMM.Boot[[vv]][valid.pos])))
    
  } ))
  
  NORM <- GMM.Boot[[boot.pos]]^2
  valid.pos <- which( abs(log(NORM) - median(log(NORM))) <= Boot.valid.thr*IQR(log(NORM)) )
  GMM.Simple.Var.Boot <- matrix(var(GMM.Boot[[boot.pos]][valid.pos]),1,1)*(Boot.Factor)
  
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
  GMM.beta.naive.lambda      <- mean(GMM.resid.naive.lambda)
  GMM.ATT.naive.lambda       <- Post.Time.Basis*GMM.beta.naive.lambda
  
  ##################################
  # GMM with ridge regularization
  # Point Estimate
  # HAC variance
  # Block bootstrap
  ##################################
  
  GMM.Regular.Coef <- c(GMM.beta.naive.lambda,GMM.gamma.naive.lambda)
  
  GRAD1 <- GMM.Ft.Grad.lambda(GMM.Regular.Coef,GMM.Data,lambda.opt)
  GRAD2 <- GMM.Ft.Grad(GMM.Regular.Coef,GMM.Data)
  
  GMM.Regular.VAR.HAC <- Meat.HAC(Res.Mat = GMM.Ft(GMM.Regular.Coef, GMM.Data),
                                  GRAD = GRAD,
                                  beta.pos = 1:lengthb,
                                  bw.type="auto")
  
  SGRAD1 <- svd(GRAD1)
  SGRAD2 <- svd(GRAD2)
  
  GMM.Regular.Var <- (ginv(GRAD1) %*% t(ginv(GRAD1)) %*% t(GRAD2)) %*% 
    GMM.Regular.VAR.HAC$matrix %*% 
    t(ginv(GRAD1) %*% t(ginv(GRAD1)) %*% t(GRAD2))/(T0+T1)
  
  ##################################
  # GMM with ridge regularization
  # HAC variance
  ##################################
  
  GMM.Boot.lambda <- BOOT(round(GMM.Regular.VAR.HAC$bw.B*Boot.Scale),
                          Num.Boot,
                          T0,
                          T1,
                          Boot.Factor=Boot.Factor,
                          type="SPSC.Reg")
  
  boot.pos.lambda <- which.max(sapply(1:length(GMM.Boot.lambda),function(vv){
    
    NORM <- GMM.Boot.lambda[[vv]]^2
    valid.pos <- which( abs(log(NORM) - median(log(NORM))) <= Boot.valid.thr*IQR(log(NORM)) )
    return((var(GMM.Boot.lambda[[vv]][valid.pos])))
    
  } ))
  
  NORM <- GMM.Boot.lambda[[boot.pos.lambda]]^2
  valid.pos <- which( abs(log(NORM) - median(log(NORM))) <= Boot.valid.thr*IQR(log(NORM)) )
  GMM.Regular.Var.Boot <- matrix(var(GMM.Boot.lambda[[boot.pos.lambda]][valid.pos]),1,1)*(Boot.Factor)
  
  ##################################
  # OLS
  # Point Estimate
  ##################################
  
  OLS.Data <- cbind(rbind(Wmat.Pre,Wmat.Post),
                    c(rep(0,T.Pre),rep(1,T.Post)),
                    c(Y1.Pre,Y1.Post))
  
  colnames(OLS.Data) <- (c(sprintf("W%0.4d",1:N),
                           "A",
                           "Y"))
  
  OLS.gamma.naive <- my.inverse(A = t(Wmat.Pre)%*%(Wmat.Pre) )%*%(t(Wmat.Pre)%*%Y1.Pre) 
  OLS.resid.naive <- as.numeric(Y1.Post - (Wmat.Post)%*%OLS.gamma.naive)
  OLS.beta.naive     <- mean(OLS.resid.naive)
  OLS.ATT.naive      <- Post.Time.Basis*OLS.beta.naive
  
  ##################################
  # OLS
  # HAC variance
  ##################################
  
  OLS.Simple.Coef <- c(OLS.beta.naive,OLS.gamma.naive) 
  
  GRAD <- OLS.Ft.Grad(OLS.Simple.Coef,OLS.Data)
  OLS.Simple.VAR.HAC <- Meat.HAC(Res.Mat = OLS.Ft(OLS.Simple.Coef, OLS.Data),
                                 GRAD = GRAD,
                                 beta.pos = 1:lengthb,
                                 bw.type="auto")
  
  SGRAD <- svd(GRAD)
  
  OLS.Simple.Var <- (SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)) %*% OLS.Simple.VAR.HAC$matrix %*% t((SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)))/(T0+T1)
  
  ##################################
  # OLS
  # Block bootstrap
  ##################################
  
  OLS.Boot <- BOOT(round(OLS.Simple.VAR.HAC$bw.B*Boot.Scale),
                   Num.Boot,
                   T0,
                   T1,
                   Boot.Factor=Boot.Factor,
                   type="OLS")
  
  boot.pos <- which.max(sapply(1:length(GMM.Boot),function(vv){
    
    NORM <- OLS.Boot[[vv]]^2
    valid.pos <- which( abs(log(NORM) - median(log(NORM))) <= Boot.valid.thr*IQR(log(NORM)) )
    return((var(OLS.Boot[[vv]][valid.pos])))
  } )) 
  
  NORM <- OLS.Boot[[boot.pos]]^2
  valid.pos <- which( abs(log(NORM) - median(log(NORM))) <= Boot.valid.thr*IQR(log(NORM)) )
  OLS.Simple.Var.Boot <- matrix(var(OLS.Boot[[boot.pos]][valid.pos]),1,1)*(Boot.Factor)
  
  save.image(file=sprintf("Result/Time_GP%d_Data_SPSC_ATT.RData",gp))
  
}





