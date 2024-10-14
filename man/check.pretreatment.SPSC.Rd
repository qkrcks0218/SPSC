\name{check.pretreatment.SPSC}
\alias{check.pretreatment.SPSC}
\title{Checking the Pre-treatment Fit of a Single Proxy Synthetic Control Object}
\description{
  Checking the Pre-treatment Fit of a Single Proxy Synthetic Control Object}
\usage{
plot(spsc,                           # SPSC object
     ...)
}
\arguments{
  \item{spsc}{A SPSC object obtained from \code{SPSC}.}
}
\references{
  Chan Park & Eric Tchetgen Tchetgen (2023) \strong{Single Proxy Synthetic Control}, \emph{arXiv} \href{https://arxiv.org/abs/2307.16353}{[Link]}
}
\examples{

################################################################################
# Toy Example from Interactive Fixed Effect Models
################################################################################

## Parameters
T0         <- T1 <- 50                                       # length of pre- and post-treatment period
Tt         <- T1+T0                                          # length of total time series
N.Inv      <- N.Val   <- N.Donor <- 8                        # number of W=valid donors and V=invalid donors
rho.lambda <- 0.5                                            # AR coefficient of lantent factor
corr.Y0.W  <- 0.5                                            # corr(error of Y, error of W)
W.coef     <- rbind(8:1/4, rep(c(0.8,0.6,0.4,0.2),each=2))   # factor loadings of valid donors
V.coef     <- rbind(rep(1,8), rep(0.5,8))                    # factor loadings of invalid donors
Y0.coef    <- c(W.coef[,1:N.Donor]\%*\%rep(1/N.Donor,N.Donor)) # factor loadings of Y
SD         <- 0.1                                            # sd of errors
BT         <- 1+c(0,0,(1:Tt/T0))                             # baseline trend
True.ATT   <- 3                                              # true effect

## Latent factor of W
lambda.eps.series       <- matrix(0,2+Tt,2)
lambda.eps.series[1:2,] <- rnorm(4)*SD
for(time.index in 1:Tt){
  lambda.eps.series[time.index+2,] <-
    rho.lambda*lambda.eps.series[time.index+1,] +
    rho.lambda*lambda.eps.series[time.index,]/2 +
    rnorm(2)*SD
}
lambda.series <- lambda.eps.series + matrix(BT,2+Tt,2)

## Latent factor of V
zeta.eps.series       <- matrix(0,2+Tt,2)
zeta.eps.series[1:2,] <- rnorm(4)*SD
for(time.index in 1:Tt){
  zeta.eps.series[time.index+2,] <-
    rho.lambda*zeta.eps.series[time.index+1,] +
    rho.lambda*zeta.eps.series[time.index,]/2 +
    rnorm(2)*SD
}
zeta.series <- zeta.eps.series + matrix(BT,2+Tt,2)

## Generate Y,W,V
common.eps  <- rnorm(2+Tt)
Y0.idio.eps <- rnorm(2+Tt)
W.idio.eps  <- matrix(rnorm(N.Val*(2+Tt)),2+Tt,N.Val)
V.idio.eps  <- matrix(rnorm(N.Inv*(2+Tt)),2+Tt,N.Inv)
Y0.eps      <- SD*( corr.Y0.W*common.eps + sqrt(1-corr.Y0.W^2)*Y0.idio.eps )
W.eps       <- SD*( corr.Y0.W*matrix(common.eps,2+Tt,N.Val) + sqrt(1-corr.Y0.W^2)*W.idio.eps )
V.eps       <- SD*( V.idio.eps )

Y0.series <- c(lambda.series\%*\%Y0.coef) + Y0.eps
W.series  <- lambda.series\%*\%W.coef[,1:N.Val] + W.eps
V.series  <- zeta.series\%*\%V.coef[,1:N.Inv] + V.eps

## Generate error-prone treatment effect
beta.eps        <- c(rep(0,T0+2),rnorm(T1))
beta            <- c(rep(0,T0+2),rep(True.ATT,T1))
beta.with.noise <- beta + beta.eps*SD

## Post-treatment Y
Y1.series <- Y0.series + beta.with.noise

## Observed Y
Yobs.series                     <- rep(0,T0+T1+2)
Yobs.series[1:(2+T0)]           <- Y0.series[1:(2+T0)]
Yobs.series[(2+T0)+1:T1]        <- Y1.series[(2+T0)+1:T1]

## Pre-treatment series
Dmat.Pre     <- cbind(W.series,V.series)[2+(1:T0),]
Y1.Pre       <- Y0.series[2+(1:T0)]

## Post-treatment series
Dmat.Post    <- cbind(W.series,V.series)[2+T0+(1:T1),]
Y1.Post      <- Y1.series[2+T0+(1:T1)]
True.TT.Vec  <- beta.with.noise[2+T0+1:T1]

## SPSC
SPSC.Detrend <- SPSC(Y.Pre              = Y1.Pre,
                     Y.Post             = Y1.Post,
                     W.Pre              = Dmat.Pre,
                     W.Post             = Dmat.Post,
                     detrend            = TRUE,
                     detrend.ft         = function(t){Spline.Trend(t,T0,df=5)},
                     Y.basis            = function(y){matrix(c(y),1,1)},
                     att.ft             = function(t){matrix(c(1),1,1)},
                     lambda.type        = "cv",
                     lambda.value       = NULL,
                     lambda.grid        = seq(-6,2,by=0.5),
                     bootstrap.num      = 100,
                     conformal.period   = 1:T1,
                     conformal.cover    = TRUE,
                     true.effect        = True.TT.Vec,
                     conformal.interval = TRUE,
                     conformal.pvalue   = 0.05,
                     conformal.window   = 25)

# Average treatment effect and 95\% confidence interval
cbind(SPSC.Detrend$ATT - 1.96*SPSC.Detrend$ASE.ATT,
      SPSC.Detrend$ATT + 1.96*SPSC.Detrend$ASE.ATT)
True.ATT


par(mfrow=c(1,2))
## Graphical summary
plot.SPSC(SPSC.Detrend, caption.position="topleft")

## Pre-treatment Plot
check.pretreatment.SPSC(SPSC.Detrend)

################################################################################
# Toy Example of California Smoking Example
################################################################################

# install.packages(tidysynth)

library(tidysynth)
data("smoking")
State <- unique(smoking$state)
N <- length(unique(smoking$state))-1
Y <- smoking$cigsale[smoking$state=="California"]
Tt <- length(Y)
T0 <- 18
T1 <- Tt-T0
D <- matrix(0,Tt,N+1)
for(jj in 1:(N+1)){
  D[,jj] <- smoking$cigsale[smoking$state==State[jj]]
}
D <- D[,-which(State=="California")]

Y.Pre <- Y[1:T0]
Y.Post <- Y[T0+1:T1]
D.Pre <- D[1:T0,]
D.Post <- D[T0+1:T1,]

## SPSC
## smaller regularization parameters return numeric error due to small sample size
SPSC.Smoking <- SPSC(Y.Pre              = Y.Pre,
                     Y.Post             = Y.Post,
                     W.Pre              = D.Pre,
                     W.Post             = D.Post,
                     detrend            = TRUE,
                     detrend.ft         = function(t){Spline.Trend(t,T0,df=5)},
                     Y.basis            = function(y){matrix(c(y),1,1)},
                     att.ft             = function(t){matrix(c(1,t),1,2)},
                     lambda.type        = "cv",
                     lambda.value       = NULL,
                     lambda.grid        = seq(-1,2,by=0.5),
                     bootstrap.num      = 100,
                     conformal.period   = NULL,
                     conformal.cover    = FALSE,
                     true.effect        = NULL,
                     conformal.interval = FALSE,
                     conformal.pvalue   = 0.05,
                     conformal.window   = 25)

cbind(SPSC.Smoking$ATT - 1.96*SPSC.Smoking$ASE.ATT,
      SPSC.Smoking$ATT + 1.96*SPSC.Smoking$ASE.ATT)

cbind(SPSC.Smoking$ATT - 1.96*SPSC.Smoking$BSE.ATT,
      SPSC.Smoking$ATT + 1.96*SPSC.Smoking$BSE.ATT)

par(mfrow=c(1,2))
## Graphical summary
plot.SPSC(SPSC.Smoking, caption.position="bottomleft")

## Pre-treatment Plot
check.pretreatment.SPSC(SPSC.Smoking)
}
