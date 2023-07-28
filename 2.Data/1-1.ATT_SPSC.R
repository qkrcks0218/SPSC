rm(list=ls())

# setwd("/Users/chanpark/Dropbox/Chan/Research/Postdoc2022/COCASC/SubmitCode/Data")
setwd("D:/Dropbox/Chan/Research/Postdoc2022/COCASC/SubmitCode/Data")
source("0.Function_SPSC_Data.R")

for(gp in 1:6){
  
  gY13 <- c(6.5,7.1)

  y.group <- 3
  ATT.Type <- "constant"
  lengthb  <- 1  
  
  mgrid <- c(24,24,24,24,10,48)
  m <- mgrid[gp]
  
  
  m.add <- 0
  Num.Boot  <- 10
  Boot.valid.thr  <- 10000
  Boot.Scale <- c(0.8,0.9,1,1.1,1.2)
  
  gY.Bound <- gY13
  Ypos <- c(1,3)
  
  gT.Bound <- c(-25,25)+c(0,167)
  stab.const <- 10^(-8)
  Boot.Factor <- 1
  bs.intercept <- F
  
  #########################################################
  
  library(splines)
  library(MASS)
  
  Data <- read.csv("Data.csv")
  Data$prc_log <- log( Data$mid_itp )
  Data$date <- as.numeric( as.Date(Data$date) )+23715
  
  # Firm ID 1 has missing time between date 665 (time 99) - 670 (time 101)
  # We impute this value
  
  (table(Data$ID, Data$date))[1,] 
  Data[Data$ID==1 & Data$date==665,]
  Data[Data$ID==1 & Data$date==670,]
  
  Data <- rbind(Data[1:99,],
                Data[Data$ID==1 & Data$date==665,],
                Data[100:dim(Data)[1],])
  Data$time <- rep(1:411,59)
  Data$date <- rep(sort(unique(Data$date)),59)
  
  # CUT
  Start <- as.numeric( as.Date("1905-12-31") )+23715
  END   <- as.numeric( as.Date("1909-01-02") )+23715
  Data  <- Data[Start <= Data$date & Data$date <= END, ]
  
  # 34, 37, 57 ; 37 = Lincoln
  
  Actual.Time <- unique(Data$date)
  Actual.Time.Date <- as.Date(Actual.Time-23715,origin="1970-01-01")
  
  Donor.Index <- setdiff(unique(Data$ID),unique( Data$ID[Data$treat_a+Data$treat_c==1] ))
  
  
  
  T.Post <- T1  <- sum(Data$treat_c)/3
  T.Pre  <- T0  <- dim(Data)[1]/length(unique(Data$ID)) - T1
  
  Post.Time.Basis     <- rep(1,T1)
  
  N <- length(Donor.Index)
  
  Wmat.series <- matrix(0,T0+T1,N)
  
  for(w.iter in 1:length(Donor.Index)){
    d.index <- Donor.Index[w.iter]
    Wmat.series[,w.iter] <- Data[Data$ID==d.index,]$prc_log
  }
  
  Wmat.series  <- Wmat.series
  
  Ymat.series <- rep(0,T0+T1)
  Ymat.series <- cbind(Data$prc_log[Data$ID==34],
                       Data$prc_log[Data$ID==37],
                       Data$prc_log[Data$ID==57])
  Y1.series <- apply(Ymat.series[,Ypos],1,mean)
  Ymat.Pre  <- Ymat.series[(1:T0),]
  Ymat.Post <- Ymat.series[T0+(1:T1),]
  
  ## Pre-treatment series 
  
  Wmat.Pre <- Wmat.series[(1:T0),] # + matrix(rnorm(T0*N),T0,N)*0.001
  Y1.Pre   <- Y1.series[(1:T0)]    # + rnorm(T0)*0.001
  
  Wmat.Post <- Wmat.series[T0+(1:T1),] # + matrix(rnorm(T1*N),T1,N)*0.001
  Y1.Post   <- Y1.series[T0+(1:T1)]    # + rnorm(T1)*0.001
  
  Wmat.series <- Wmat.series 
  Wmat.Pre    <- Wmat.Pre    
  Wmat.Post   <- Wmat.Post   
  
  
  gY.Pre <- matrix(0,length(Y1.Pre),m+m.add)
  gY.Fit <- bs(Y1.Pre,df=m+m.add,
               intercept = F,
               Boundary.knots = gY.Bound)
  gY.Pre[,1:(m+m.add)]       <- gY.Fit
  
  overlap <- cbind(1:N,(sapply(1:N,
                               function(bb){
                                 RR1 <- range(Wmat.series[1:T0,bb])
                                 mean(c(as.numeric(RR1[1] <= Wmat.series[T0+1:T1,bb] &
                                                     Wmat.series[T0+1:T1,bb] <= RR1[2] ) ))
                                 
                               })))
  
  vc <- list()
  
  SEQ <- c(12,24,36,48)+0.5
  
  vc[[1]] <- (1:N)[ (50-rank(overlap[,2]) < SEQ[1]) ]
  vc[[2]] <- (1:N)[ (SEQ[1] <= 50-rank(overlap[,2]) & 50-rank(overlap[,2]) < SEQ[2]) ]
  vc[[3]] <- (1:N)[ (SEQ[2] <= 50-rank(overlap[,2]) & 50-rank(overlap[,2]) < SEQ[3]) ]
  vc[[4]] <- (1:N)[ (SEQ[3] <= 50-rank(overlap[,2]) & 50-rank(overlap[,2]) < SEQ[4]) ]
  vc[[5]] <- c(1,2,5,9,18)
  vc[[6]] <- c(1,2,3,4,5,6,8,9,12,13,14,19,20,23,27,33,34,37,38,40,41,42,44,47)
  
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
  
  GMM.gamma.naive <- my.inverse(A = t(GW)%*%(GW) )%*%(t(GW)%*%GY)
  GMM.resid.naive <- as.numeric(Y1.Post - (Wmat.Post)%*%GMM.gamma.naive)
  
  GMM.beta.naive      <- mean(GMM.resid.naive)
  GMM.ATT.naive       <- Post.Time.Basis*GMM.beta.naive
  
  ##########################
  
  GMM.Simple.Coef <- c(GMM.beta.naive,GMM.gamma.naive)
  
  GRAD <- GMM.Ft.Grad(GMM.Simple.Coef,GMM.Data)
  
  GMM.Simple.VAR.HAC <- Meat.HAC(Res.Mat = GMM.Ft(GMM.Simple.Coef, GMM.Data),
                                 GRAD = GRAD,
                                 beta.pos = 1:lengthb,
                                 bw.type="auto")
  
  SGRAD <- svd(GRAD)
  
  GMM.Simple.Var <- (SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)) %*% GMM.Simple.VAR.HAC$matrix %*% t((SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)))/(T0+T1)
  
  #########################
  
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
  
  #####################
  
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
  
  
  
  ###############################
  
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
  
  ############################################
  
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
  
  OLS.Simple.Coef <- c(OLS.beta.naive,OLS.gamma.naive) 
  
  GRAD <- OLS.Ft.Grad(OLS.Simple.Coef,OLS.Data)
  OLS.Simple.VAR.HAC <- Meat.HAC(Res.Mat = OLS.Ft(OLS.Simple.Coef, OLS.Data),
                                 GRAD = GRAD,
                                 beta.pos = 1:lengthb,
                                 bw.type="auto")
  
  SGRAD <- svd(GRAD)
  
  OLS.Simple.Var <- (SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)) %*% OLS.Simple.VAR.HAC$matrix %*% t((SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)))/(T0+T1)
  
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
  
  save.image(file=sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  
}
