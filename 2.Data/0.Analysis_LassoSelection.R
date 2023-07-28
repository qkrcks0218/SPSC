rm(list=ls())

source("0.Function_SPSC_Data.R")

gp <- 5 ; y.group <- 3; att <- 1
ATT.Type <- "constant"; lengthb  <- 1  

Nm <- rep(0,48)

for(m.index in 1:48){
  
  m <- m.index*2+2
  
  m.add <- 0
  Num.Boot  <- 1000
  Boot.valid.thr  <- 10000
  Boot.Scale <- c(0.8,0.9,1,1.1,1.2)
  gY13 <- gY.Bound <- c(6.5,7.1)
  Ypos <- c(1,3)
  YF <- "Result_Unit13"
  
  
  gT.Bound <- c(-25,25)+c(0,167)
  stab.const <- 10^(-8)
  Boot.Factor <- 1
  bs.intercept <- F
  
  #########################################################
  
  # library(readstata13)
  library(splines)
  library(MASS)
  # library(gmm)
  
  # setwd("D:/Dropbox/Chan/Research/Postdoc2022/COCASC/Data/data")
  
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
  
  if(ATT.Type=="spline"){
    Post.Time.Basis.Fit <- bs(1:T1,
                              df=(lengthb),
                              intercept = bs.intercept,
                              Boundary.knots = gT.Bound)
    Post.Time.Basis     <- Post.Time.Basis.Fit
  } else if (ATT.Type=="constant") {
    Post.Time.Basis     <- rep(1,T1)
  } else if (ATT.Type=="exponential") {
    Post.Time.Basis     <- cbind(1,exp((1:T1-T1)/T1))
  }
  
  
  
  
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
                                                     Wmat.series[T0+1:T1,bb] <= RR1[2] )))
                               })))
  
  
  N <- dim(Wmat.Pre)[2]
  m <- dim(gY.Pre)[2]
  if(ATT.Type=="constant"){
    lengthb <- 1
  } else {
    lengthb <- dim(Post.Time.Basis)[2]
  }
  
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

