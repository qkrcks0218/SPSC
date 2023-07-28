############################
# Recommend to implement parallel computing
# BATCH=1,...,504
# We only run for BATCH=1 for an illustration purpose
############################

BATCH <- 1

source("0.Function_SPSC_Data.R")

gY13 <- c(6.5,7.1)

PARA <- expand.grid(1:84,1:6)

gp      <- PARA[BATCH,2]
BATCH   <- PARA[BATCH,1]
y.group <- 3

ATT.Type <- "constant"; lengthb  <- 1  

mgrid <- c(24,24,24,24,10,48)
m <- mgrid[gp]

m.add <- 0
T0 <- 217
mT <- round( (T0^(1/3)) )
Num.Boot  <- 10
Boot.valid.thr  <- 10000
Boot.Scale <- c(0.8,0.9,1,1.1,1.2)

gY.Bound <- gY13
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


gY.Pre <- matrix(0,length(Y1.Pre),m+mT)
gY.Fit <- bs(Y1.Pre,df=m,
             intercept = F,
             Boundary.knots = gY.Bound)
gY.Pre[,1:(m)]       <- gY.Fit
bT <- bs(1:T0,df=(mT),Boundary.knots = c(-2,T0+2))
gY.Pre[,m+1:mT] <- bT

overlap <- cbind(1:N,(sapply(1:N,
                             function(bb){
                               RR1 <- range(Wmat.series[1:T0,bb])
                               mean(c(as.numeric(RR1[1] <= Wmat.series[T0+1:T1,bb] &
                                                   Wmat.series[T0+1:T1,bb] <= RR1[2] ) ))
                               
                               # RR1 <- range(Wmat.series[1:T0,bb])
                               # RR2 <- range(Wmat.series[T0+1:T1,bb])
                               # mean(c(as.numeric(RR1[1] <= Wmat.series[T0+1:T1,bb] & 
                               #                     Wmat.series[T0+1:T1,bb] <= RR1[2]),
                               #        c(as.numeric(RR2[1] <= Wmat.series[1:T0,bb] &
                               #                       Wmat.series[1:T0,bb] <= RR2[2] ))))
                               
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

lengthb <- 1

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

## GMM-1st
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


##########################

lambda.grid <- seq(-5,0,by=0.2)
lambda.min  <- lambda.grid[ which.min(sapply(lambda.grid,CV.Lambda)) ]

lambda.opt <- optimize(f=CV.Lambda,
                       lower= lambda.min - 2,
                       upper= lambda.min + 2)$minimum

GMM.gamma.naive.lambda <- my.inverse(A = t(GW)%*%(GW),
                                     stab.const = ((T.Pre+T.Post)/(T.Pre))^2*10^(lambda.opt),
                                     adjust = F)%*%(t(GW)%*%GY)
GMM.resid.naive.lambda <- as.numeric(Y1.Post - (Wmat.Post)%*%GMM.gamma.naive.lambda)


############################################

PPP <- matrix(1:(T1+1),84,2,byrow=T)
PPP[84,] <- c(166,167)

T.Window <- PPP[BATCH,]

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

write.csv(CP.CI,    sprintf("Conformal_Raw/Raw_Time_Group%s_CPCI_%0.3d.csv",           gp,BATCH),row.names=F)
write.csv(CP.Beta.1,sprintf("Conformal_Raw/Raw_Time_Group%s_CPBeta_SPSC_%0.3d.csv",    gp,BATCH),row.names=F)
write.csv(CP.Beta.2,sprintf("Conformal_Raw/Raw_Time_Group%s_CPBeta_SPSC_Reg_%0.3d.csv",gp,BATCH),row.names=F)
write.csv(CP.PV.1,  sprintf("Conformal_Raw/Raw_Time_Group%s_CPPV_SPSC_%0.3d.csv",      gp,BATCH),row.names=F)
write.csv(CP.PV.2,  sprintf("Conformal_Raw/Raw_Time_Group%s_CPPV_SPSC_Reg_%0.3d.csv",  gp,BATCH),row.names=F)



