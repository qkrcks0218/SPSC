args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

library(MASS)
library(splines)
library(augsynth)
library(Synth)
library(scpi)

source("0.Function_SSC_Simulation.R")

beta.vec <- c(3,3) # true treatment effect
delta    <- expand.grid(0:1,1:2,c(2,5,9),1:10000)[BATCH,1] # delta=0: No X; delta=1: Use X
lengthb  <- expand.grid(0:1,1:2,c(2,5,9),1:10000)[BATCH,2] # b=1: constant ATT; b=2: linear ATT;
N.Donor  <- expand.grid(0:1,1:2,c(2,5,9),1:10000)[BATCH,3] # N=number of donors
Iter     <- expand.grid(0:1,1:2,c(2,5,9),1:10000)[BATCH,4] # random seed

Bootstrap          <- T     # bootstrap?
Num.Boot           <- 100   # number of bootstrap
Conformal          <- T     # conformal inference?
AugSynth.Ind       <- T

set.seed(BATCH)

for(ttt in 1:4){
  
  ## Time length
  
  T.Pre  <- c(50,100,250,500)[ttt]
  T.Post <- c(50,100,250,500)[ttt]
  Tt     <- T.Pre+T.Post
  Post.Conformal.Vec <- round(seq(0,T.Post,length=11)[-1])
  True.ATT <- beta.vec[1] + (1+T.Post)/T.Post/2*beta.vec[2]*(lengthb-1)
  
  ## Number of Donors
  N     <- N.Donor
  
  ## Dimension of g
  m  <- 2*N.Donor
  
  ## Correlations
  rho.Y0        <- 0.2
  rho.Y.eps     <- 0.2
  rho.W.eps     <- 0.2
  rho.beta      <- 0.2
  rho.cov       <- 0.2
  rho.cov.eps   <- 0.2
  SD            <- 1
  X.effect      <- 0.25
  sX            <- 0.5
  sY            <- 0.75
  sW            <- 0.75
  sB            <- 0.5
  
  ## Generate X.eps 
  X.eps.series <- matrix(0,T.Pre+T.Post+2,N+1)
  for(time.index in 1:(T.Pre+T.Post)){
    X.eps.series[time.index+2,] <- 
      rho.cov.eps*X.eps.series[time.index+1,] + rho.cov.eps*X.eps.series[time.index,]/2 + rnorm(N+1) * SD * sX
  }
  
  ## Generate Y.eps 
  Y0.eps <- rep(0,T.Pre+T.Post+2)
  for(time.index in 1:(T.Pre+T.Post)){
    Y0.eps[time.index+2] <-  
      rho.Y.eps*Y0.eps[time.index+1] + rho.Y.eps*Y0.eps[time.index]/2 + rnorm(1) * SD * sY
  }
  
  ## Generate beta.eps
  beta.eps <- rep(0,T.Pre+T.Post+2)
  for(time.index in (T.Pre+1):(T.Pre+T.Post)){
    beta.eps[time.index+2] <-
      rho.beta*beta.eps[time.index+1] + rho.beta*beta.eps[time.index]/2 + rnorm(1) * SD * sB
  }
  
  ## Effects
  if(lengthb==1){
    beta         <- c(0,0,rep(0,T.Pre),rep(beta.vec[1],T.Post))
  } else {
    beta         <- c(0,0,rep(0,T.Pre),
                      beta.vec[1] + beta.vec[2]*(1:T.Post)/T.Post)
  }
  
  ## Baseline Trend 
  BT <- c(0,0,1:(T.Pre+T.Post)/T.Pre)
  
  ## Generate X 
  X.series <- matrix(0,T.Pre+T.Post+2,N+1)
  for(time.index in 1:(T.Pre+T.Post)){
    X.series[time.index+2,] <- 
      rho.cov*X.series[time.index+1,] + rho.cov*X.series[time.index,]/2 + X.eps.series[time.index+2,]
  }
  
  X.series <- X.series*as.numeric(delta!=0)
  
  ## Generate Y
  Y0.series <- rep(0,T.Pre+T.Post+2)
  for(time.index in 1:(T.Pre+T.Post)){
    Y0.series[time.index+2] <- 
      rho.Y0*Y0.series[time.index+1] + rho.Y0*Y0.series[time.index]/2 + 
      Y0.eps[time.index+2] + BT[time.index+2] + delta*X.series[time.index+2,1]*X.effect
  }
  SP.Window <- c( min(Y0.series) - min(Y0.series[1:(2+T.Pre)]) - 0.5 , 
                  max(Y0.series) - max(Y0.series[1:(2+T.Pre)]) + 0.5 )
  Y1.series <- Y0.series + beta + beta.eps
  beta.with.noise <- beta + beta.eps
  Yobs.series <- rep(0,T.Pre+T.Post+2)
  Yobs.series[1:(2+T.Pre)] <- Y0.series[1:(2+T.Pre)]
  Yobs.series[(2+T.Pre)+1:T.Post] <- Y1.series[(2+T.Pre)+1:T.Post]
  
  ## Generate error
  eps.series <- matrix(0,T.Pre+T.Post+2,N)
  for(time.index in 1:(T.Pre+T.Post)){
    eps.series[time.index+2,] <- 
      rho.W.eps*eps.series[time.index+1,] + rho.W.eps*eps.series[time.index,]/2 + rnorm(N) * SD * sW
  }
  
  ## Generate W
  
  WEIGHT <- 1.5
  DEIGHT <- -1
  
  AAA <- matrix(0,N,N)
  
  diag(AAA) <- WEIGHT
  for(wmatindex in 1:(N-1)){
    AAA[wmatindex,wmatindex+1] <- DEIGHT
  }
  AAA[N,1] <- DEIGHT
  
  det(AAA)
  
  COEF <- rep(0,N)
  COEF[2] <- 1
  
  gamma <- ginv(AAA)%*%(COEF)
  gamma
  
  W.series      <- matrix(0,T.Pre+T.Post+2,N)
  Y0.original   <- Y0.series - delta*X.series[,1]*X.effect
  
  if(N==2){
    for(wit in 1:N){
      W.series[, wit] <- 
        AAA[1,wit] + 
        AAA[2,wit]*Y0.original
    }
  } else if (N==5) {
    for(wit in 1:N){
      W.series[, wit] <- 
        AAA[1,wit] + 
        AAA[2,wit]*Y0.original +
        AAA[3,wit]*(Y0.series^2)/2 + 
        AAA[4,wit]*as.numeric(Y0.series>3) +
        AAA[5,wit]*as.numeric(Y0.series<0)
    }
  } else if (N==9) {
    for(wit in 1:N){
      W.series[, wit] <- 
        AAA[1,wit] + 
        AAA[2,wit]*Y0.original +
        AAA[3,wit]*(Y0.series^2)/2 + 
        AAA[4,wit]*as.numeric(Y0.series>3) +
        AAA[5,wit]*as.numeric(Y0.series<0) + 
        AAA[6,wit]*as.numeric(0 <= Y0.series & Y0.series<1) + 
        AAA[7,wit]*as.numeric(1 <= Y0.series & Y0.series<2) + 
        AAA[8,wit]*exp((Y0.series-1.5)/2.5) +
        AAA[9,wit]*exp(-(Y0.series-1.5)/2.5)
      
    }
  }
  
  W.series      <- W.series + delta*X.series[,1+1:N]*X.effect + eps.series
  
  ## Pre-treatment series 
  Wmat.Pre <- W.series[2+(1:T.Pre),]
  Y1.Pre   <- Y0.series[2+(1:T.Pre)]
  Xmat.Pre <- X.series[2+(1:T.Pre),]
  
  ## Post-treatment series 
  Wmat.Post <- W.series[2+T.Pre+(1:T.Post),]
  Y1.Post   <- Y1.series[2+T.Pre+(1:T.Post)]
  Xmat.Post <- X.series[2+T.Pre+(1:T.Post),]
  
  ## Gy matrix
  
  if(delta==0){
    gY.Pre <- matrix(0,length(Y1.Pre),m) # each row = collection of g(Y)s at time t , each col = time series of each g(Y)
    
  } else {
    gY.Pre <- matrix(0,length(Y1.Pre),m+N+1) # each row = collection of g(Y)s at time t , each col = time series of each g(Y)
    
  }
  
  gY.Pre[,1:(m)] <- bs(Y1.Pre,df=(m),Boundary.knots = range(Y1.Pre))
  
  if(delta==1){
    gY.Pre[,m+1:(N+1)] <- Xmat.Pre # include covariates in g
  } 
  
  if(delta==0){
    GMM.Data <- cbind(rbind(gY.Pre,matrix(0,T.Post,m)),
                      rbind(Wmat.Pre,Wmat.Post),
                      c(rep(0,T.Pre),rep(1,T.Post)),
                      c(Y1.Pre,Y1.Post))
    
    colnames(GMM.Data) <- (c(sprintf("G%0.4d",1:m),
                             sprintf("W%0.4d",1:N),
                             "A",
                             "Y"))
  } else {
    GMM.Data <- cbind(rbind(gY.Pre,matrix(0,T.Post,m+N+1)),
                      rbind(Wmat.Pre,Wmat.Post),
                      c(rep(0,T.Pre),rep(1,T.Post)),
                      c(Y1.Pre,Y1.Post),
                      rbind(Xmat.Pre,Xmat.Post))
    
    colnames(GMM.Data) <- (c(sprintf("G%0.4d",1:(m+N+1)),
                             sprintf("W%0.4d",1:N),
                             "A",
                             "Y",
                             sprintf("X%0.4d",1:(N+1))))
  }
  
  
  
  ########################
  
  T0 <- dim(Wmat.Pre)[1]
  T1 <- dim(Wmat.Post)[1]
  Tt <- T0 + T1
  m  <- dim(gY.Pre)[2]
  N  <- dim(Wmat.Pre)[2]
  A  <- rep(c(0,1),c(T0,T1))
  
  ## GMM-1st
  Wmat <- rbind(Wmat.Pre,Wmat.Post)
  Y    <- c(Y1.Pre,Y1.Post)
  gY   <- rbind(gY.Pre, matrix(0,T1,m))
  
  if(delta==0){
    GW <- t( sapply(1:dim(gY.Pre)[2],function(tt){
      apply( t(Wmat.Pre)*matrix(gY.Pre[,tt],N,T.Pre,byrow=T), 1, mean)
    }) )
  } else {
    GW <- t( sapply(1:dim(gY.Pre)[2],function(tt){
      apply( t(cbind(Wmat.Pre,Xmat.Pre))*matrix(gY.Pre[,tt],N+N+1,T.Pre,byrow=T), 1, mean)
    }) )
  }
  
  
  GY <- ( sapply(1:dim(gY.Pre)[2],function(tt){
    mean(gY.Pre[,tt]*Y1.Pre)
  }) )
  
  gamma.naive <- c( ginv(t(GW)%*%GW)%*%(t(GW)%*%GY) )
  if(lengthb==1 & delta==0){
    beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%gamma.naive)~1)$coefficients )
  } else if(lengthb==2 & delta==0){
    beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%gamma.naive)~I((1:T1)/T1))$coefficients )
  } else if(lengthb==1 & delta==1){
    beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - cbind(Wmat.Post,Xmat.Post)%*%gamma.naive)~1)$coefficients )
  } else if(lengthb==2 & delta==1){
    beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - cbind(Wmat.Post,Xmat.Post)%*%gamma.naive)~I((1:T1)/T1))$coefficients )
  }
  
  
  
  ##########################
  
  GMM.Simple.Coef <- c(beta.naive,gamma.naive)
  
  GRAD <- GMM.Ft.Grad(GMM.Simple.Coef,GMM.Data)
  
  GMM.Simple.VAR.HAC <- Meat.HAC(Res.Mat = GMM.Ft(GMM.Simple.Coef, GMM.Data),
                                 GRAD = GRAD,
                                 beta.pos = 1:lengthb,
                                 bw.type="auto")
  
  SGRAD <- svd(GRAD)
  
  GMM.Simple.Var <- (SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)) %*% GMM.Simple.VAR.HAC$matrix %*% t((SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)))/(T0+T1)
  
  
  ###########################
  
  if(Bootstrap){
    
    GMM.Boot <- BOOT(round(GMM.Simple.VAR.HAC$bw.B*c(8,9,10,11,12)/10),
                     Num.Boot,
                     T0,
                     T1,
                     Boot.Factor=1,
                     type="SPSC")
    
    boot.pos <- which.max(sapply(1:length(GMM.Boot),function(vv){
      
      if(lengthb>1){
        return( norm(var(GMM.Boot[[vv]])) )
      } else {
        return((var(GMM.Boot[[vv]])))
      } 
    } ))  
    
    GMM.Simple.Var.Boot <- var(GMM.Boot[[boot.pos]])
  } else {
    GMM.Simple.Var.Boot <- matrix(100,lengthb,lengthb)
    
  }
  
  if(Conformal){
    
    if(delta==0){
      beta.series <- 
        Y1.Post[Post.Conformal.Vec] - Wmat.Post[Post.Conformal.Vec,]%*%gamma.naive
    } else if (delta==1) {
      beta.series <- 
        Y1.Post[Post.Conformal.Vec] - cbind(Wmat.Post,Xmat.Post)[Post.Conformal.Vec,]%*%gamma.naive
    }
    
    
    if(delta==0 & lengthb==1){
      CP.CI <- Conformal.Prediction(Wmat.Pre,
                                    Wmat.Post[Post.Conformal.Vec,],
                                    Y1.Pre,
                                    Y1.Post[Post.Conformal.Vec],
                                    Xmat.Pre=NULL,
                                    Xmat.Post=NULL,
                                    cov.ind=0,
                                    center=beta.series,
                                    bw=sd(beta.series))
    } else if(delta==0 & lengthb==2){
      CP.CI <- Conformal.Prediction(Wmat.Pre,
                                    Wmat.Post[Post.Conformal.Vec,],
                                    Y1.Pre,
                                    Y1.Post[Post.Conformal.Vec],
                                    Xmat.Pre=NULL,
                                    Xmat.Post=NULL,
                                    cov.ind=0,
                                    center=beta.series,
                                    bw=sd(beta.series))
    } else if(delta==1 & lengthb==1){
      CP.CI <- Conformal.Prediction(Wmat.Pre,
                                    Wmat.Post[Post.Conformal.Vec,],
                                    Y1.Pre,
                                    Y1.Post[Post.Conformal.Vec],
                                    Xmat.Pre,
                                    Xmat.Post[Post.Conformal.Vec,],
                                    cov.ind=1,
                                    center=beta.series,
                                    bw=sd(beta.series))
    } else if(delta==1 & lengthb==2){
      CP.CI <- Conformal.Prediction(Wmat.Pre,
                                    Wmat.Post[Post.Conformal.Vec,],
                                    Y1.Pre,
                                    Y1.Post[Post.Conformal.Vec],
                                    Xmat.Pre,
                                    Xmat.Post[Post.Conformal.Vec,],
                                    cov.ind=1,
                                    center=beta.series,
                                    bw=sd(beta.series))
    }
    
    
    if(delta==0 & lengthb==1){
      CP.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                         Wmat.Post[Post.Conformal.Vec,],
                                         Y1.Pre,
                                         Y1.Post[Post.Conformal.Vec],
                                         Xmat.Pre=NULL,
                                         Xmat.Post=NULL,
                                         cov.ind=0,
                                         effect=beta.with.noise[2+T0+Post.Conformal.Vec])
    } else if(delta==0 & lengthb==2){
      CP.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                         Wmat.Post[Post.Conformal.Vec,],
                                         Y1.Pre,
                                         Y1.Post[Post.Conformal.Vec],
                                         Xmat.Pre=NULL,
                                         Xmat.Post=NULL,
                                         cov.ind=0,
                                         effect=beta.with.noise[2+T0+Post.Conformal.Vec])
    } else if(delta==1 & lengthb==1){
      CP.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                         Wmat.Post[Post.Conformal.Vec,],
                                         Y1.Pre,
                                         Y1.Post[Post.Conformal.Vec],
                                         Xmat.Pre,
                                         Xmat.Post[Post.Conformal.Vec,],
                                         cov.ind=1,
                                         effect=beta.with.noise[2+T0+Post.Conformal.Vec])
    } else if(delta==1 & lengthb==2){
      CP.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                         Wmat.Post[Post.Conformal.Vec,],
                                         Y1.Pre,
                                         Y1.Post[Post.Conformal.Vec],
                                         Xmat.Pre,
                                         Xmat.Post[Post.Conformal.Vec,],
                                         cov.ind=1,
                                         effect=beta.with.noise[2+T0+Post.Conformal.Vec])
    }
    
    
  } else {
    CP.CI <- cbind(rep(-100,2),
                   rep( 100,2))
    CP.CV <- rep(1,T.Post)
  }
  
  
  
  ##########################
  
  lambda.grid <- seq(-5,0,by=0.2)
  lambda.opt  <- lambda.grid[ which.min(sapply(lambda.grid,CV.Lambda)) ]
  
  # lambda.opt <- optimize(f=CV.Lambda,
  #                        lower= lambda.min - 2,
  #                        upper= lambda.min + 2)$minimum
  
  gamma.naive.lambda <- ginv(t(GW)%*%(GW) + ((T.Pre+T.Post)/(T.Pre))^2*diag(rep(10^(lambda.opt),dim(GW)[2])) )%*%(t(GW)%*%(GY))
  
  if(lengthb==1 & delta==0){
    beta.naive.lambda  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%gamma.naive.lambda)~1)$coefficients )
  } else if(lengthb==2 & delta==0){
    beta.naive.lambda  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%gamma.naive.lambda)~I((1:T1)/T1))$coefficients )
  } else if(lengthb==1 & delta==1){
    beta.naive.lambda  <- as.numeric( lm(as.numeric(Y1.Post - cbind(Wmat.Post,Xmat.Post)%*%gamma.naive.lambda)~1)$coefficients )
  } else if(lengthb==2 & delta==1){
    beta.naive.lambda  <- as.numeric( lm(as.numeric(Y1.Post - cbind(Wmat.Post,Xmat.Post)%*%gamma.naive.lambda)~I((1:T1)/T1))$coefficients )
  }
  
  GMM.Regular.Coef <- c(beta.naive.lambda,gamma.naive.lambda)
  
  # GRAD1 <- GMM.Ft.Grad.lambda(GMM.Regular.Coef,GMM.Data,lambda.opt)
  GRAD2 <- GMM.Ft.Grad(GMM.Regular.Coef,GMM.Data)
  
  GMM.Regular.VAR.HAC <- Meat.HAC(Res.Mat = GMM.Ft(GMM.Regular.Coef, GMM.Data),
                                  GRAD = GRAD,
                                  beta.pos = 1:lengthb,
                                  bw.type="auto")
  
  # SGRAD1 <- svd(GRAD1)
  SGRAD2 <- svd(GRAD2)
  
  GMM.Regular.Var <- (ginv(GRAD2) %*% t(ginv(GRAD2)) %*% t(GRAD2)) %*% 
    GMM.Regular.VAR.HAC$matrix %*% 
    t(ginv(GRAD2) %*% t(ginv(GRAD2)) %*% t(GRAD2))/(T0+T1)
  
  
  ###########################
  
  
  if(Bootstrap){
    
    GMM.Boot <- BOOT(round(GMM.Regular.VAR.HAC$bw.B*c(8,9,10,11,12)/10),
                     Num.Boot,
                     T0,
                     T1,
                     Boot.Factor=1,
                     type="SPSC.Reg")
    
    boot.pos <- which.max(sapply(1:length(GMM.Boot),function(vv){
      
      if(lengthb>1){
        return( norm(var(GMM.Boot[[vv]])) )
      } else {
        return((var(GMM.Boot[[vv]])))
      } 
    } ))  
    
    GMM.Regular.Var.Boot <- var(GMM.Boot[[boot.pos]])
    
  } else {
    
    GMM.Regular.Var.Boot <- matrix(100,lengthb,lengthb)
    
  }
  
  if(Conformal){
    
    
    if(delta==0){
      beta.series <- 
        Y1.Post[Post.Conformal.Vec] - Wmat.Post[Post.Conformal.Vec,]%*%gamma.naive.lambda
    } else if (delta==1) {
      beta.series <- 
        Y1.Post[Post.Conformal.Vec] - cbind(Wmat.Post,Xmat.Post)[Post.Conformal.Vec,]%*%gamma.naive.lambda
    }
    
    if(delta==0 & lengthb==1){
      CP.Regular.CI <- Conformal.Prediction(Wmat.Pre,
                                            Wmat.Post[Post.Conformal.Vec,],
                                            Y1.Pre,
                                            Y1.Post[Post.Conformal.Vec],
                                            Xmat.Pre=NULL,
                                            Xmat.Post=NULL,
                                            cov.ind=0,
                                            center=beta.series,
                                            bw=sd(beta.series),
                                            lambda=lambda.opt,
                                            constant=4)
    } else if(delta==0 & lengthb==2){
      CP.Regular.CI <- Conformal.Prediction(Wmat.Pre,
                                            Wmat.Post[Post.Conformal.Vec,],
                                            Y1.Pre,
                                            Y1.Post[Post.Conformal.Vec],
                                            Xmat.Pre=NULL,
                                            Xmat.Post=NULL,
                                            cov.ind=0,
                                            center=beta.series,
                                            bw=sd(beta.series),
                                            lambda=lambda.opt,
                                            constant=4)
    } else if(delta==1 & lengthb==1){
      CP.Regular.CI <- Conformal.Prediction(Wmat.Pre,
                                            Wmat.Post[Post.Conformal.Vec,],
                                            Y1.Pre,
                                            Y1.Post[Post.Conformal.Vec],
                                            Xmat.Pre,
                                            Xmat.Post[Post.Conformal.Vec,],
                                            cov.ind=1,
                                            center=beta.series,
                                            bw=sd(beta.series),
                                            lambda=lambda.opt,
                                            constant=4)
    } else if(delta==1 & lengthb==2){
      CP.Regular.CI <- Conformal.Prediction(Wmat.Pre,
                                            Wmat.Post[Post.Conformal.Vec,],
                                            Y1.Pre,
                                            Y1.Post[Post.Conformal.Vec],
                                            Xmat.Pre,
                                            Xmat.Post[Post.Conformal.Vec,],
                                            cov.ind=1,
                                            center=beta.series,
                                            bw=sd(beta.series),
                                            lambda=lambda.opt,
                                            constant=4)
    }
    
    
    
    if(delta==0 & lengthb==1){
      CP.Regular.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                                 Wmat.Post[Post.Conformal.Vec,],
                                                 Y1.Pre,
                                                 Y1.Post[Post.Conformal.Vec],
                                                 Xmat.Pre=NULL,
                                                 Xmat.Post=NULL,
                                                 cov.ind=0,
                                                 effect=beta.with.noise[2+T0+Post.Conformal.Vec],
                                                 lambda=lambda.opt,
                                                 constant=4)
    } else if(delta==0 & lengthb==2){
      CP.Regular.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                                 Wmat.Post[Post.Conformal.Vec,],
                                                 Y1.Pre,
                                                 Y1.Post[Post.Conformal.Vec],
                                                 Xmat.Pre=NULL,
                                                 Xmat.Post=NULL,
                                                 cov.ind=0,
                                                 effect=beta.with.noise[2+T0+Post.Conformal.Vec],
                                                 lambda=lambda.opt,
                                                 constant=4)
    } else if(delta==1 & lengthb==1){
      CP.Regular.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                                 Wmat.Post[Post.Conformal.Vec,],
                                                 Y1.Pre,
                                                 Y1.Post[Post.Conformal.Vec],
                                                 Xmat.Pre,
                                                 Xmat.Post[Post.Conformal.Vec,],
                                                 cov.ind=1,
                                                 effect=beta.with.noise[2+T0+Post.Conformal.Vec],
                                                 lambda=lambda.opt,
                                                 constant=4)
    } else if(delta==1 & lengthb==2){
      CP.Regular.CV <- Conformal.Prediction.Fast(Wmat.Pre,
                                                 Wmat.Post[Post.Conformal.Vec,],
                                                 Y1.Pre,
                                                 Y1.Post[Post.Conformal.Vec],
                                                 Xmat.Pre,
                                                 Xmat.Post[Post.Conformal.Vec,],
                                                 cov.ind=1,
                                                 effect=beta.with.noise[2+T0+Post.Conformal.Vec],
                                                 lambda=lambda.opt,
                                                 constant=4)
    }
    
    
  } else {
    CP.Regular.CI <- cbind(rep(-100,2),
                           rep( 100,2))
    CP.Regular.CV <- rep(1,T.Post)
  }
  
  
  ##########################
  
  
  if(delta==0){
    OLS.Data <- cbind(rbind(Wmat.Pre,Wmat.Post),
                      c(rep(0,T.Pre),rep(1,T.Post)),
                      c(Y1.Pre,Y1.Post))
    
    colnames(OLS.Data) <- (c(sprintf("W%0.4d",1:N),
                             "A",
                             "Y"))
  } else {
    OLS.Data <- cbind(rbind(Wmat.Pre,Wmat.Post),
                      c(rep(0,T.Pre),rep(1,T.Post)),
                      c(Y1.Pre,Y1.Post),
                      rbind(Xmat.Pre,Xmat.Post))
    
    colnames(OLS.Data) <- (c(sprintf("W%0.4d",1:N),
                             "A",
                             "Y",
                             sprintf("X%0.4d",1:(N+1))))
  }
  
  
  
  
  if(delta==0){
    ols.gamma.naive <- ginv(t(Wmat.Pre)%*%(Wmat.Pre))%*%(t(Wmat.Pre)%*%Y1.Pre)
  } else if (delta==1){
    WXmat.Pre <- cbind(Wmat.Pre,Xmat.Pre)
    ols.gamma.naive <- ginv(t(WXmat.Pre)%*%(WXmat.Pre))%*%(t(WXmat.Pre)%*%Y1.Pre)
  }
  
  if(lengthb==1 & delta==0){
    ols.beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%ols.gamma.naive)~1)$coefficients )
  } else if(lengthb==2 & delta==0){
    ols.beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%ols.gamma.naive)~I((1:T1)/T1))$coefficients )
  } else if(lengthb==1 & delta==1){
    ols.beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - cbind(Wmat.Post,Xmat.Post)%*%ols.gamma.naive)~1)$coefficients )
  } else if(lengthb==2 & delta==1){
    ols.beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - cbind(Wmat.Post,Xmat.Post)%*%ols.gamma.naive)~I((1:T1)/T1))$coefficients )
  }
  
  OLS.Simple.Coef <- c(ols.beta.naive,ols.gamma.naive) 
  
  GRAD <- OLS.Ft.Grad(OLS.Simple.Coef,OLS.Data)
  OLS.Simple.VAR.HAC <- Meat.HAC(Res.Mat = OLS.Ft(OLS.Simple.Coef, OLS.Data),
                                 GRAD = GRAD,
                                 beta.pos = 1:lengthb,
                                 bw.type="auto")
  
  SGRAD <- svd(GRAD)
  
  OLS.Simple.Var <- (SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)) %*% OLS.Simple.VAR.HAC$matrix %*% 
    t((SGRAD$v%*%diag(1/SGRAD$d)%*%t(SGRAD$u)))/(T0+T1)
  
  
  if(Bootstrap){
    
    OLS.Boot <- BOOT(round(OLS.Simple.VAR.HAC$bw.B*c(8,9,10,11,12)/10),
                     Num.Boot,
                     T0,
                     T1,
                     Boot.Factor=1,
                     type="OLS")
    
    boot.pos <- which.max(sapply(1:length(OLS.Boot),function(vv){
      
      if(lengthb>1){
        return( norm(var(OLS.Boot[[vv]])) )
      } else {
        return((var(OLS.Boot[[vv]])))
      } 
    } ))  
    
    OLS.Simple.Var.Boot <- var(OLS.Boot[[boot.pos]])
    
  } else {
    
    if(lengthb==1){
      OLS.Simple.Var.Boot <- 100
    } else {
      OLS.Simple.Var.Boot <- c(100,100)
    }
    
  }
  
  
  ##########################
  
  
  if(delta==0){
    
    Abadie.SC.Data <- data.frame(
      idvar = rep(c(0,1:N),each=(Tt)),
      idvarchr = as.character(rep(c(0,1:N),each=(Tt))),
      timevar = rep(1:Tt,N+1),
      outcomevar = c(c(Y1.Pre,Y1.Post),
                     as.vector( rbind(Wmat.Pre,Wmat.Post) )),
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
    ATT.Abadie <- mean(Y1.Post - Wmat.Post%*%synth.out$solution.w)
  } else {
    Abadie.SC.Data <- data.frame(
      idvar = rep(c(0,1:N),each=(Tt)),
      idvarchr = as.character(rep(c(0,1:N),each=(Tt))),
      timevar = rep(1:Tt,N+1),
      outcomevar = c(c(Y1.Pre,Y1.Post),
                     as.vector( rbind(Wmat.Pre,Wmat.Post) )),
      X1 = c(as.vector( rbind(Xmat.Pre,Xmat.Post) ))
    )
    
    dataprep.out <- dataprep( foo=Abadie.SC.Data,
                              predictors=c("outcomevar","X1"),
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
    ATT.Abadie <- mean(Y1.Post - Wmat.Post%*%synth.out$solution.w)
  }
  
  ##########################
  
  if(Conformal){
    
    if(AugSynth.Ind){
      
      OLS.Data.DF <- data.frame(id=rep(c(1,1+1:N),each=(T0+T1)),
                                time=rep(1:(T0+T1),N+1),
                                outcome=c(c(Y1.Pre,Y1.Post),
                                          as.vector(rbind(Wmat.Pre,Wmat.Post))),
                                treatment=c(rep(c(0,1),c(T0,T1)),
                                            rep(0,(T0+T1)*N)),
                                covariate=as.vector(rbind(Xmat.Pre,Xmat.Post)))
      
      ASC.Time <- NULL
      ASC.Time <- c(ASC.Time,which(OLS.Data.DF$time<=T0))
      for(tt in 1:length(Post.Conformal.Vec)){
        ASC.Time <- c(ASC.Time,which(OLS.Data.DF$time==T0+Post.Conformal.Vec[tt]))
      }
      
      OLS.Data.DF.Short <- OLS.Data.DF[ASC.Time,]
      
      
      
      if(delta==0){
        ASC.L <- augsynth(outcome~treatment,
                          unit=id,
                          time=time,
                          data=OLS.Data.DF,
                          progfunc="Ridge",
                          scm=T)
        
        ASC.S <- augsynth(outcome~treatment,
                          unit=id,
                          time=time,
                          data=OLS.Data.DF.Short,
                          progfunc="Ridge",
                          scm=T)
      } else {
        ASC.L <- augsynth(outcome~treatment|covariate,
                          unit=id,
                          time=time,
                          data=OLS.Data.DF,
                          progfunc="Ridge",
                          scm=T)
        ASC.S <- augsynth(outcome~treatment|covariate,
                          unit=id,
                          time=time,
                          data=OLS.Data.DF.Short,
                          progfunc="Ridge",
                          scm=T)
      }
      
      ASC.L.Summary   <- summary(ASC.L,inf=F)
      ATT.ASC         <- as.numeric( ASC.L.Summary$average_att[1] )
      
      ASC.S.Summary   <- summary(ASC.S)
      ASC.Cover <- mean( ASC.S.Summary$att$lower_bound[T0+1:length(Post.Conformal.Vec)] <= 
                           beta.with.noise[2+T0+Post.Conformal.Vec] &
                           beta.with.noise[2+T0+Post.Conformal.Vec] <= 
                           ASC.S.Summary$att$upper_bound[T0+1:length(Post.Conformal.Vec)] )
      ASC.Length <- mean( ASC.S.Summary$att$upper_bound[T0+1:length(Post.Conformal.Vec)] 
                          - ASC.S.Summary$att$lower_bound[T0+1:length(Post.Conformal.Vec)] )
      
    } else {
      
      ATT.ASC <- -100
      ASC.Cover <- 0
      ASC.Length <- 100
      
    }
    
    
    
    if(delta==0){
      OLS.Data.DF <- data.frame(id=rep(c(0,1:N),each=(T0+T1)),
                                time=rep(1:(T0+T1),N+1),
                                outcome=c(c(Y1.Pre,Y1.Post),
                                          as.vector(rbind(Wmat.Pre,Wmat.Post))))
      
      SCD <- scdata(OLS.Data.DF,
                    id.var="id",
                    time.var="time",
                    outcome.var="outcome",
                    period.pre=1:T0,
                    period.post=T0+1:T1,
                    unit.tr=0,
                    unit.co=1:N,
                    constant=F)
      
      SCPI.Est <- scpi(SCD)
    } else {
      OLS.Data.DF <- data.frame(id=rep(c(0,1:N),each=(T0+T1)),
                                time=rep(1:(T0+T1),N+1),
                                outcome=c(c(Y1.Pre,Y1.Post),
                                          as.vector(rbind(Wmat.Pre,Wmat.Post))),
                                X=as.vector(X.series[2+1:(T0+T1),]))
      
      SCD <- scdata(OLS.Data.DF,
                    id.var="id",
                    time.var="time",
                    outcome.var="outcome",
                    period.pre=1:T0,
                    period.post=T0+1:T1,
                    unit.tr=0,
                    unit.co=1:N,
                    cov.adj=list(c("X")),
                    constant=F)
      
      SCPI.Est <- scpi(SCD)
    }
    
    
    SCPI.Effect <- mean(SCPI.Est$data$Y.post - SCPI.Est$est.results$Y.post.fit)
    SCPI.SE     <- mean(SCPI.Est$inference.results$CI.all.gaussian[,3])/2/qnorm(0.975)
    SCPI.Cover  <- mean( (SCPI.Est$data$Y.post - 
                            SCPI.Est$inference.results$CI.all.gaussian[,2] <= 
                            beta.with.noise[2+T0+1:T1] &
                            beta.with.noise[2+T0+1:T1] <= 
                            SCPI.Est$data$Y.post - 
                            SCPI.Est$inference.results$CI.all.gaussian[,1])[Post.Conformal.Vec] )
    SCPI.Cover.F  <- mean( (SCPI.Est$data$Y.post - 
                              SCPI.Est$inference.results$CI.all.gaussian[,2] <= 
                              beta.with.noise[2+T0+1:T1] &
                              beta.with.noise[2+T0+1:T1] <= 
                              SCPI.Est$data$Y.post - 
                              SCPI.Est$inference.results$CI.all.gaussian[,1]) )
    SCPI.Length <- mean( SCPI.Est$inference.results$CI.all.gaussian[Post.Conformal.Vec,3] )
    
  } else {
    
    ATT.ASC <- -100
    ASC.Cover <- 0
    ASC.Length <- 100

    SCPI.Effect <- -100
    SCPI.SE     <- 100
    SCPI.Cover  <- 0
    SCPI.Cover.F  <- 0
    SCPI.Length <- 100
  }
  
  ##########################
  
  ## Result summary
  
  if(lengthb==1){
    
    RRR <- c( OLS.Simple.Coef[1] - True.ATT, 
              sqrt(OLS.Simple.Var[1,1]),
              COVER(OLS.Simple.Coef[1] , sqrt(OLS.Simple.Var[1,1]) , True.ATT),
              sqrt(OLS.Simple.Var.Boot[1] ),
              COVER(OLS.Simple.Coef[1] , sqrt(OLS.Simple.Var.Boot[1]) , True.ATT),
              
              GMM.Simple.Coef[1] - True.ATT, 
              sqrt(GMM.Simple.Var[1,1]),
              COVER(GMM.Simple.Coef[1] , sqrt(GMM.Simple.Var[1,1]) , True.ATT),
              sqrt(GMM.Simple.Var.Boot[1] ),
              COVER(GMM.Simple.Coef[1] , sqrt(GMM.Simple.Var.Boot[1] ) , True.ATT),
              mean(CP.CI[,2]-CP.CI[,1]),
              mean(CP.CV),
              
              GMM.Regular.Coef[1] - True.ATT,
              sqrt(GMM.Regular.Var[1,1]),
              COVER(GMM.Regular.Coef[1] , sqrt(GMM.Regular.Var[1,1]) , True.ATT),
              sqrt(GMM.Regular.Var.Boot[1] ),
              COVER(GMM.Regular.Coef[1] , sqrt(GMM.Regular.Var.Boot[1] ) , True.ATT),
              mean(CP.Regular.CI[,2]-CP.Regular.CI[,1]),
              mean(CP.Regular.CV),
              
              lambda.opt,
              
              SCPI.Effect - True.ATT, SCPI.SE, SCPI.Cover, SCPI.Cover.F, SCPI.Length, 
              ATT.Abadie - True.ATT,  
              ATT.ASC - True.ATT, ASC.Cover, ASC.Length )
    
    
  } else {
    
    AVER <- function(EEE,VVV,est=T){
      v <- c(1,mean((1:T1)/(T1)))
      TE <- True.ATT
      EST <- as.numeric( t(v)%*%EEE[1:2] )
      SE <- sqrt( as.numeric(  t(v)%*%(VVV[1:2,1:2])%*%(v) ) )
      
      if(est==T){
        Result <- c( EST-TE, 
                     SE,
                     COVER(EST,SE,TE) )
      } else {
        Result <- c( SE,
                     COVER(EST,SE,TE) )
      }
      
      Result
    }
    
    RRR <- c( AVER(OLS.Simple.Coef,
                   OLS.Simple.Var),
              AVER(OLS.Simple.Coef,
                   OLS.Simple.Var.Boot,est=F),
              
              AVER(GMM.Simple.Coef,
                   GMM.Simple.Var),
              AVER(GMM.Simple.Coef,
                   GMM.Simple.Var.Boot,est=F),
              mean(CP.CI[,2]-CP.CI[,1]),
              mean(CP.CV),
              
              AVER(GMM.Regular.Coef,
                   GMM.Regular.Var),
              AVER(GMM.Regular.Coef,
                   GMM.Regular.Var.Boot,est=F),
              mean(CP.Regular.CI[,2]-CP.Regular.CI[,1]),
              mean(CP.Regular.CV),
              lambda.opt,
              
              SCPI.Effect - True.ATT, 
              SCPI.SE, SCPI.Cover, SCPI.Cover.F, SCPI.Length, 
              ATT.Abadie - True.ATT, 
              ATT.ASC - True.ATT, ASC.Cover, ASC.Length )
    
  }
  
  RRR <- matrix(RRR,1,length(RRR))
  colnames(RRR) <- c("OLSBias","OLSSE","OLSCover","OLSBSE","OLSBCover",
                     
                     "SSCBias","SSCSE","SSCCover","SSCBSE","SSCBCover",
                     "SSCLength","SSCCPCover",
                     
                     "SSCRegBias","SSCRegSE","SSCRegCover",
                     "SSCRegBSE","SSCRegBCover",
                     "SSCRegLength","SSCRegCPCover",
                     
                     "SSCRegLambda",
                     
                     "SCPIBias",
                     "SCPISE","SCPICover","SCPI.Cover.F","SCPILength",
                     "AbadieBias", 
                     "ASCBias", "ASCCover", "ASCLength")
  
  write.csv(RRR,
            sprintf("SPSC_T%0.5d_D%0.1d_B%0.1d_N%0.2d_Iter%0.3d.csv",
                    T.Pre,delta,lengthb,N,Iter),
            row.names=F)
   
}  
 
