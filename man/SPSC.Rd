\name{SPSC}
\alias{SPSC}
\title{Single Proxy Synthetic Control}
\description{
  Single Proxy Synthetic Control by Park and Tchetgen Tchetgen (2023).}
\usage{
SPSC(
    Y.Pre,                                               # vector of Y in the pre-trt period
    Y.Post,                                              # vector of Y in the post-trt period
    W.Pre,                                               # matrix of W (T0 x N) in the post-trt period
    W.Post,                                              # matrix of W (T1 x N) in the post-trt period
    detrend = TRUE,                                      # detrend Y0
    detrend.ft = function(t){Spline.Trend(t,T0,df=5)},   # detrend basis
    Y.basis = function(y){matrix(c(y),1,1)},             # basis function of Y0
    att.ft = function(t){matrix(c(1),1,1)},              # ATT basis
    lambda.type = "cv",                                  # how to select lambda
    lambda.value = NULL,                                 # if lambda.type="fix", specify lambda.value
    lambda.grid = seq(-6,2,by=0.5),                      # grid of lambda
    bootstrap.num = 0,                                   # bootstrap for SE
    conformal.period = NULL,                             # NULL or post-treatment time period
    conformal.cover = FALSE,                             # cover=0/1, interval=prediction interval
    true.effect = NULL,                                  # true.effect
    conformal.interval = TRUE,                           # cover=0/1, interval=prediction interval
    conformal.pvalue = 0.05,                             # level of conformal inference
    conformal.window = 25                                # range of grid search for prediction intervals
    )
}
\arguments{
  \item{Y.Pre}{A vector of the outcome in the pre-treatment period. The length is \eqn{T_0}, which is the length of the pre-treatment period.}
  \item{Y.Post}{A vector of the outcome in the post-treatment period. The length is \eqn{T_1}, which is the length of the post-treatment period.}
  \item{W.Pre}{A matrix of \eqn{N} donors in the pre-treatment period. The dimension is \eqn{T_0 \times N}.}
  \item{W.Post}{A matrix of \eqn{N} donors in the post-treatment period. The dimension is \eqn{T_1 \times N}.}
  \item{detrend}{\emph{(Optional)} A logical input \code{(TRUE/FALSE)} to detrend the time series. The default is \code{detrend = TRUE}.}
  \item{detrend.ft}{\emph{(Optional)} Basis functions used in detrending. The default is \code{detrend.ft = function(t){Spline.Trend(t,T0,df=5)}}, which is the cubic B-spline basis with 5 degrees of freedom. The function must return a row matrix, e.g., \code{detrend.ft = function(t){matrix(c(1,t/T0),1,2)}}.}
  \item{Y.basis}{\emph{(Optional)} Basis functions used in synthetic control. The default is \code{Y.basis = function(y){matrix(c(y),1,1)}}. The function must return a row matrix.}
  \item{att.ft}{\emph{(Optional)} Basis functions used in treatment effect estimation. The default is \code{att.ft = function(t){matrix(c(1),1,1)}}. The function must return a row matrix.}
  \item{lambda.type}{\emph{(Optional)} A option to how to determine the regularization parameter \eqn{\lambda}. If \code{lambda.type = "cv"}, \eqn{\lambda} is selected from bootstrap among the candidates specified in \code{lambda.grid}. If \code{lambda.type = "fix"}, \eqn{\lambda} is fixed as the value specified in \code{lambda.value}. The default is \code{lambda.type = "cv"}.}
  \item{lambda.value}{\emph{(Optional)}A user-specified \eqn{\lambda} value. If \code{lambda.type = "cv"}, this input will be ignored. The default is \code{lambda.value = NULL}.}
  \item{lambda.grid}{\emph{(Optional)}A vector of user-specified \eqn{\lambda} candidates used in bootstrap procedure. This grid represents exponents for \eqn{\lambda}, i.e. \eqn{\lambda = 10^\code{lambda.grid}}. If \code{lambda.type = "fix"}, this input will be ignored. The default is \code{lambda.grid = seq(-6,2,by=0.5)}.}
  \item{bootstrap.num}{\emph{(Optional)}The number of bootstrap for calculating the standard error of the average treatment effect. The default is \code{bootstrap.num = 0}.}
  \item{conformal.period}{\emph{(Optional)}A vector of post-treatment period to conduct inference based on prediction intervals. \code{When conformal.period = 1:length(Y.Post)}, prediction intervals over the entire post-treatment period are calculated. The default is \code{conformal.period = NULL}.}
  \item{conformal.cover}{\emph{(Optional)} A logical input \code{(TRUE/FALSE)} to obtain indicators of whether prediction intervals cover the true error-prone treatment effect. The default is \code{conformal.cover = FALSE}.}
  \item{true.effect}{\emph{(Optional)} A vector of the true error-prone treatment effect. The length must match to that of conformal.period. The default is \code{true.effect = NULL}.}
  \item{conformal.interval}{\emph{(Optional)} A logical input \code{(TRUE/FALSE)} to calculate prediction intervals. The default is \code{conformal.interval = TRUE}.}
  \item{conformal.pvalue}{\emph{(Optional)} A p-value used to calculate \eqn{100(1-conformal.pvalue)\%} prediction intervals. The default if \code{conformal.pvalue = 0.05}.}
  \item{conformal.window}{\emph{(Optional)} A numeric value to specify a range of grid search for prediction intervals. Prediction intervals will be searched over \eqn{Effect \pm \code{conformal.window}*SE} where \code{SE} is the asymptotic standard error of the average treatment effect. The default if \code{conformal.window = 25}.}
}
\value{
  \item{gamma}{A synthetic control weight of donors.}
  \item{SC}{A time series of synthetic control from period \eqn{1} to \eqn{T_0+T_1}.}
  \item{Y}{A time series of the outcome from period \eqn{1} to \eqn{T_0+T_1}.}
  \item{effect}{A time series of the differences between \code{Y} and \code{SC}.}
  \item{ATT}{An estimate of the average treatment effects from period \eqn{T_0+1} to \eqn{T_0+T_1}.}
  \item{ASE}{A standard error of \code{ATT} obtained from the HAC estimator.}
  \item{BSE}{A standard error of \code{ATT} obtained from block bootstrap.}
  \item{conformal.cover}{A time series of indicators of whether prediction intervals cover the true error-prone treatment effect from period \eqn{T_0+1} to \eqn{T_0+T_1}.}
  \item{conformal.interval}{A time series of prediction intervals from period \eqn{T_0+1} to \eqn{T_0+T_1}.}
  \item{lambda}{An exponent of the regularization parameter, i.e. \eqn{\code{lambda}=\log_{10}(\text{regularization parameter})}}
  \item{trend}{A time trend spanned by the \code{detrend.ft} function.}
  \item{detrend}{\code{detrend} input value}
  \item{detrend.ft}{\code{detrend.ft} input value}
  \item{Y.basis}{\code{Y.basis} input value}
  \item{att.ft}{\code{att.ft} input value}
  \item{bootstrap.num}{\code{bootstrap.num} input value}
  \item{conformal.period}{\code{conformal.period} input value}
  \item{conformal.pvalue}{\code{conformal.pvalue} input value}
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

## Graphical summary
plot.SPSC(SPSC.Detrend, caption.position="topleft")


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

## Graphical summary
plot.SPSC(SPSC.Smoking, caption.position="bottomleft")
}
