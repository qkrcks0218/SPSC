library(MASS)

my.inverse <- function(A,
                       stab.const=NULL,
                       adjust=F){
  if( is.null(stab.const) ){
    return( ginv(A) )
  } else {
    return( solve( A + diag(rep(stab.const),dim(A)[2] ) ) - diag(rep(stab.const),dim(A)[2])*as.numeric(adjust) )
  }
}


Calculate.PValue <- function(Residual,
                             T0.CPV,
                             T1.CPV,
                             p=1){
  Residual.Extend <- c(Residual,Residual)
  if(p==1){
    S.Base <- mean(abs(Residual[T0.CPV+1:T1.CPV]))
    S.All  <- rep(0,T0.CPV+T1.CPV)
    for(bb in 1:(T0.CPV+T1.CPV)){
      S.All[bb] <- mean( abs(Residual.Extend[bb-1+1:T1.CPV]) )
    }
  } else {
    S.Base <- abs(mean(Residual[T0.CPV+1:T1.CPV]))
    S.All  <- rep(0,T0.CPV+T1.CPV)
    for(bb in 1:(T0.CPV+T1.CPV)){
      S.All[bb] <- abs(mean(Residual.Extend[bb-1+1:T1.CPV]) )
    }
  }
  
  return( 1-mean(S.All < S.Base) )
}

COVER <- function(a,b,center=3){
  as.numeric(abs(a-center)<=b*qnorm(0.975))
}

PValue.beta.SPSC <- function(ii,beta,
                             Wmat.Pre.Input.PV  ,
                             Wmat.Post.Input.PV ,
                             Y1.Pre.Input.PV    ,
                             Y1.Post.Input.PV   ,
                             Xmat.Pre.Input.PV  ,
                             Xmat.Post.Input.PV ,
                             cov.ind.Input.PV   ,
                             center.Input.PV    ,
                             gY.Fit.Input.PV    ,
                             lambda.Input.PV    ,
                             constant.Input.PV  ){
  
  T0.Input.PV <- dim(Wmat.Pre.Input.PV)[1]
  Tt.Input.PV <- T0.Input.PV + 1
  
  Wmat       <- rbind(Wmat.Pre.Input.PV,
                      Wmat.Post.Input.PV[ii,])
  if(cov.ind.Input.PV==1){
    Xmat       <- rbind(Xmat.Pre.Input.PV,
                        Xmat.Post.Input.PV[ii,])
  }
  Y          <- c(Y1.Pre.Input.PV,
                  Y1.Post.Input.PV[ii])
  
  beta.seq <- c(rep(0,T0.Input.PV),beta)
  
  Y.Pred.H0  <- Y - beta.seq
  
  if(cov.ind.Input.PV==0){
    gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m) 
  } else {
    gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m+N+1) 
  }
  
  gY.Pred.H0[1:T0.Input.PV,1:m] <- gY.Fit.Input.PV[1:T0.Input.PV,]
  gY.Pred.H0[Tt.Input.PV,1:m]   <- predict(gY.Fit.Input.PV, 
                                           Y.Pred.H0[Tt.Input.PV])
  
  Outside <- !(min(attributes(gY.Fit.Input.PV)$Boundary.knots) <= Y.Pred.H0[Tt.Input.PV] &
                 Y.Pred.H0[Tt.Input.PV] <= max(attributes(gY.Fit.Input.PV)$Boundary.knots))
  Notvalid <- sum(1-as.numeric(abs(gY.Pred.H0[Tt.Input.PV,1:m]-0.5)<=0.5)) > 0
  
  if(Outside | Notvalid){
    gY.Pred.H0[Tt.Input.PV,1:m] <- 0
  }
  
  if(cov.ind.Input.PV==1){
    gY.Pred.H0[,m+1:(N+1)] <- Xmat
  }
  
  if(cov.ind.Input.PV==0){
    GW.Input.PV <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
      apply( t(Wmat)*matrix(gY.Pred.H0[,tt],N,T0.Input.PV+1,byrow=T), 1, mean)
    }) )
  } else {
    GW.Input.PV <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
      apply( t(cbind(Wmat,Xmat))*matrix(gY.Pred.H0[,tt],N+N+1,T0.Input.PV+1,byrow=T), 1, mean)
    }) )
  }
  
  GY.Input.PV <- ( sapply(1:dim(gY.Pred.H0)[2],function(tt){
    mean(gY.Pred.H0[,tt]*Y.Pred.H0)
  }) )
  
  if(constant.Input.PV*10^(lambda.Input.PV)==0){
    gamma.naive <- my.inverse(A = t(GW.Input.PV)%*%GW.Input.PV )%*%(t(GW.Input.PV)%*%GY.Input.PV) 
  } else {
    gamma.naive <- my.inverse(A = t(GW.Input.PV)%*%GW.Input.PV,
                              stab.const = constant.Input.PV*10^(lambda.Input.PV),
                              adjust = F)%*%(t(GW.Input.PV)%*%GY.Input.PV) 
  }
  if(cov.ind.Input.PV==0){
    Residual <- c(Y.Pred.H0 - Wmat%*%gamma.naive)
  } else {
    Residual <- c(Y.Pred.H0 - cbind(Wmat,Xmat)%*%gamma.naive)
  }
  Calculate.PValue(Residual,T0.Input.PV,1)
  
}

Conformal.Prediction.Data <- function(Wmat.Pre.Input,
                                      Wmat.Post.Input,
                                      Y1.Pre.Input,
                                      Y1.Post.Input,
                                      Xmat.Pre.Input=NULL,
                                      Xmat.Post.Input=NULL,
                                      cov.ind.Input=0,
                                      center.Input,
                                      gY.Fit.Input,
                                      lambda.Input=-Inf,
                                      constant.Input=0){
  
  T0.Input <- dim(Wmat.Pre.Input)[1]
  T1.Input <- dim(Wmat.Post.Input)[1]
  Tt.Input <- T0.Input + T1.Input
  N.Input  <- dim(Wmat.Pre.Input)[2]
  A.Input  <- rep(c(0,1),c(T0.Input,T1.Input))
  
  ## GMM-1st
  
  CI <- matrix(0,T1.Input,2)
  BG.Mat <- matrix(0,1001,T1.Input)
  PV.Mat <- matrix(0,1001,T1.Input)
  
  for(ii in 1:T1.Input){
    
    beta.grid <- seq((center.Input[ii]-1),
                     (center.Input[ii]+1),
                     length=1001)
    
    diff.beta <- diff(beta.grid)[1]
    
    gsize <- diff(beta.grid)[1]
    
    PV <- sapply(beta.grid,
                 function(vv){PValue.beta.SPSC(ii,
                                               vv,
                                               Wmat.Pre.Input.PV  = Wmat.Pre.Input,
                                               Wmat.Post.Input.PV = Wmat.Post.Input,
                                               Y1.Pre.Input.PV    = Y1.Pre.Input,
                                               Y1.Post.Input.PV   = Y1.Post.Input,
                                               Xmat.Pre.Input.PV  = Xmat.Pre.Input,
                                               Xmat.Post.Input.PV = Xmat.Post.Input,
                                               cov.ind.Input.PV   = cov.ind.Input,
                                               center.Input.PV    = center.Input,
                                               gY.Fit.Input.PV    = gY.Fit.Input,
                                               lambda.Input.PV    = lambda.Input,
                                               constant.Input.PV  = constant.Input)})
    
    valid <- which(PV>=0.05)
    beta.valid <- beta.grid[valid[which(c(1,diff(valid))==1)]]
    PV.valid <- PV[valid[which(c(1,diff(valid))==1)]]
    
    CI[ii,] <- range(beta.valid)
    BG.Mat[,ii] <- beta.grid
    PV.Mat[,ii] <- PV
    
    print(ii)
  }
  
  Result <- list()
  
  Result$CI <- CI
  Result$BG <- BG.Mat
  Result$PV <- PV.Mat
  
  return(Result)
  
}




PValue.beta.SPSC.Time <- function(ii,beta,
                                  Wmat.Pre.Input.PV  ,
                                  Wmat.Post.Input.PV ,
                                  Y1.Pre.Input.PV    ,
                                  Y1.Post.Input.PV   ,
                                  Xmat.Pre.Input.PV  ,
                                  Xmat.Post.Input.PV ,
                                  cov.ind.Input.PV   ,
                                  center.Input.PV    ,
                                  gY.Fit.Input.PV    ,
                                  lambda.Input.PV    ,
                                  constant.Input.PV  ){
  
  T0.Input.PV <- dim(Wmat.Pre.Input.PV)[1]
  Tt.Input.PV <- T0.Input.PV + 1
  
  Wmat       <- rbind(Wmat.Pre.Input.PV,
                      Wmat.Post.Input.PV[ii,])
  if(cov.ind.Input.PV==1){
    Xmat       <- rbind(Xmat.Pre.Input.PV,
                        Xmat.Post.Input.PV[ii,])
  }
  Y          <- c(Y1.Pre.Input.PV,
                  Y1.Post.Input.PV[ii])
  
  beta.seq <- c(rep(0,T0.Input.PV),beta)
  
  Y.Pred.H0  <- Y - beta.seq
  
  if(cov.ind.Input.PV==0){
    gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m+mT) 
  } else {
    gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m+N+1) 
  }
  
  gY.Pred.H0[1:T0.Input.PV,1:m] <- gY.Fit.Input.PV[1:T0.Input.PV,]
  
  gY.Pred.H0[Tt.Input.PV,1:m]   <- predict(gY.Fit.Input.PV, 
                                           Y.Pred.H0[Tt.Input.PV])
  
  Outside <- !(min(attributes(gY.Fit.Input.PV)$Boundary.knots) <= Y.Pred.H0[Tt.Input.PV] &
                 Y.Pred.H0[Tt.Input.PV] <= max(attributes(gY.Fit.Input.PV)$Boundary.knots))
  Notvalid <- sum(1-as.numeric(abs(gY.Pred.H0[Tt.Input.PV,1:m]-0.5)<=0.5)) > 0
  
  if(Outside | Notvalid){
    gY.Pred.H0[Tt.Input.PV,1:m] <- 0
  }
  
  if(cov.ind.Input.PV==1){
    gY.Pred.H0[,m+1:(N+1)] <- Xmat
  }
  
  gY.Pred.H0[1:T0.Input.PV,m+1:mT] <- bT
  gY.Pred.H0[Tt.Input.PV,m+1:mT]   <- predict(bT, T0.Input.PV)
  
  if(cov.ind.Input.PV==0){
    GW.Input.PV <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
      apply( t(Wmat)*matrix(gY.Pred.H0[,tt],N,T0.Input.PV+1,byrow=T), 1, mean)
    }) )
  } else {
    GW.Input.PV <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
      apply( t(cbind(Wmat,Xmat))*matrix(gY.Pred.H0[,tt],N+N+1,T0.Input.PV+1,byrow=T), 1, mean)
    }) )
  }
  
  GY.Input.PV <- ( sapply(1:dim(gY.Pred.H0)[2],function(tt){
    mean(gY.Pred.H0[,tt]*Y.Pred.H0)
  }) )
  
  if(constant.Input.PV*10^(lambda.Input.PV)==0){
    gamma.naive <- my.inverse(A = t(GW.Input.PV)%*%GW.Input.PV )%*%(t(GW.Input.PV)%*%GY.Input.PV) 
  } else {
    gamma.naive <- my.inverse(A = t(GW.Input.PV)%*%GW.Input.PV,
                              stab.const = constant.Input.PV*10^(lambda.Input.PV),
                              adjust = F)%*%(t(GW.Input.PV)%*%GY.Input.PV) 
  }
  if(cov.ind.Input.PV==0){
    Residual <- c(Y.Pred.H0 - Wmat%*%gamma.naive)
  } else {
    Residual <- c(Y.Pred.H0 - cbind(Wmat,Xmat)%*%gamma.naive)
  }
  Calculate.PValue(Residual,T0.Input.PV,1)
  
}

Conformal.Prediction.Data.Time <- function(Wmat.Pre.Input,
                                           Wmat.Post.Input,
                                           Y1.Pre.Input,
                                           Y1.Post.Input,
                                           Xmat.Pre.Input=NULL,
                                           Xmat.Post.Input=NULL,
                                           cov.ind.Input=0,
                                           center.Input,
                                           gY.Fit.Input,
                                           lambda.Input=-Inf,
                                           constant.Input=0){
  
  T0.Input <- dim(Wmat.Pre.Input)[1]
  T1.Input <- dim(Wmat.Post.Input)[1]
  Tt.Input <- T0.Input + T1.Input
  N.Input  <- dim(Wmat.Pre.Input)[2]
  A.Input  <- rep(c(0,1),c(T0.Input,T1.Input))
  
  ## GMM-1st
  
  CI <- matrix(0,T1.Input,2)
  BG.Mat <- matrix(0,1001,T1.Input)
  PV.Mat <- matrix(0,1001,T1.Input)
  
  for(ii in 1:T1.Input){
    
    beta.grid <- seq((center.Input[ii]-1),
                     (center.Input[ii]+1),
                     length=1001)
    
    diff.beta <- diff(beta.grid)[1]
    
    gsize <- diff(beta.grid)[1]
    
    PV <- sapply(beta.grid,
                 function(vv){PValue.beta.SPSC.Time(ii,
                                                    vv,
                                                    Wmat.Pre.Input.PV  = Wmat.Pre.Input,
                                                    Wmat.Post.Input.PV = Wmat.Post.Input,
                                                    Y1.Pre.Input.PV    = Y1.Pre.Input,
                                                    Y1.Post.Input.PV   = Y1.Post.Input,
                                                    Xmat.Pre.Input.PV  = Xmat.Pre.Input,
                                                    Xmat.Post.Input.PV = Xmat.Post.Input,
                                                    cov.ind.Input.PV   = cov.ind.Input,
                                                    center.Input.PV    = center.Input,
                                                    gY.Fit.Input.PV    = gY.Fit.Input,
                                                    lambda.Input.PV    = lambda.Input,
                                                    constant.Input.PV  = constant.Input)})
    
    valid <- which(PV>=0.05)
    beta.valid <- beta.grid[valid[which(c(1,diff(valid))==1)]]
    PV.valid <- PV[valid[which(c(1,diff(valid))==1)]]
    
    CI[ii,] <- range(beta.valid)
    BG.Mat[,ii] <- beta.grid
    PV.Mat[,ii] <- PV
    
    print(ii)
  }
  
  Result <- list()
  
  Result$CI <- CI
  Result$BG <- BG.Mat
  Result$PV <- PV.Mat
  
  return(Result)
  
}


PValue.beta.OLS <-  function(ii,beta,
                             Wmat.Pre.Input.PV  ,
                             Wmat.Post.Input.PV ,
                             Y1.Pre.Input.PV    ,
                             Y1.Post.Input.PV   ,
                             Xmat.Pre.Input.PV  ,
                             Xmat.Post.Input.PV ,
                             cov.ind.Input.PV   ,
                             center.Input.PV    ){
  
  T0.Input.PV <- dim(Wmat.Pre.Input.PV)[1]
  Tt.Input.PV <- T0.Input.PV + 1
  
  Wmat       <- rbind(Wmat.Pre.Input.PV,
                      Wmat.Post.Input.PV[ii,])
  if(cov.ind.Input.PV==1){
    Xmat       <- rbind(Xmat.Pre.Input.PV,
                        Xmat.Post.Input.PV[ii,])
  }
  Y          <- c(Y1.Pre.Input.PV,
                  Y1.Post.Input.PV[ii])
  
  beta.seq <- c(rep(0,T0.Input.PV),beta)
  
  Y.Pred.H0  <- Y - beta.seq
  
  if(cov.ind.Input.PV==0){
    gamma.naive <- my.inverse(t(Wmat)%*%Wmat)%*%(t(Wmat)%*%Y.Pred.H0)
    Residual <- c(Y.Pred.H0 - Wmat%*%gamma.naive)
  } else {
    gamma.naive <- my.inverse(t(cbind(Wmat,Xmat))%*%(cbind(Wmat,Xmat)))%*%(t(cbind(Wmat,Xmat))%*%Y.Pred.H0)
    Residual <- c(Y.Pred.H0 - cbind(Wmat,Xmat)%*%gamma.naive)
  }
  
  Calculate.PValue(Residual,T0.Input.PV,1) 
  
}


Conformal.Prediction.Data.OLS <- function(Wmat.Pre.Input,
                                          Wmat.Post.Input,
                                          Y1.Pre.Input,
                                          Y1.Post.Input,
                                          Xmat.Pre.Input=NULL,
                                          Xmat.Post.Input=NULL,
                                          cov.ind.Input=0,
                                          center.Input){
  
  T0.Input <- dim(Wmat.Pre.Input)[1]
  T1.Input <- dim(Wmat.Post.Input)[1]
  Tt.Input <- T0.Input + T1.Input
  N.Input  <- dim(Wmat.Pre.Input)[2]
  A.Input  <- rep(c(0,1),c(T0.Input,T1.Input))
  
  ## GMM-1st
  
  CI <- matrix(0,T1.Input,2)
  
  BG.Mat <- matrix(0,1001,T1.Input)
  PV.Mat <- matrix(0,1001,T1.Input)
  
  for(ii in 1:T1.Input){
    
    beta.grid <- seq((center.Input[ii]-1),
                     (center.Input[ii]+1),
                     length=1001)
    
    diff.beta <- diff(beta.grid)[1]
    
    gsize <- diff(beta.grid)[1]
    
    PV <- sapply(beta.grid,
                 function(vv){PValue.beta.OLS(ii,
                                              vv,
                                              Wmat.Pre.Input.PV  = Wmat.Pre.Input,
                                              Wmat.Post.Input.PV = Wmat.Post.Input,
                                              Y1.Pre.Input.PV    = Y1.Pre.Input,
                                              Y1.Post.Input.PV   = Y1.Post.Input,
                                              Xmat.Pre.Input.PV  = Xmat.Pre.Input,
                                              Xmat.Post.Input.PV = Xmat.Post.Input,
                                              cov.ind.Input.PV   = cov.ind.Input,
                                              center.Input.PV    = center.Input)})
    
    valid <- which(PV>=0.05)
    beta.valid <- beta.grid[valid[which(c(1,diff(valid))==1)]]
    PV.valid <- PV[valid[which(c(1,diff(valid))==1)]]
    
    CI[ii,] <- range(beta.valid)
    BG.Mat[,ii] <- beta.grid
    PV.Mat[,ii] <- PV
    
    print(ii)
  }
  
  Result <- list()
  
  Result$CI <- CI
  Result$BG <- BG.Mat
  Result$PV <- PV.Mat
  
  return(Result)
  
}


SANDWICH <- function(SVD,MEAT,Tt){
  SVD$v%*%diag(1/SVD$d)%*%t(SVD$u)%*%MEAT%*%t(SVD$v%*%diag(1/SVD$d)%*%t(SVD$u))/Tt
}

GMM.Ft <- function(theta,x){
  
  # theta
  long.gamma <- theta[-c(1:lengthb)]
  beta  <- theta[1:lengthb]
  
  # gY
  gY <- x[,substr(colnames(x),1,1)=="G"]
  
  # donor
  W  <- x[,substr(colnames(x),1,1)=="W"]
  
  # Time
  A  <- x[,substr(colnames(x),1,1)=="A"]
  
  # Outcome
  Y  <- x[,substr(colnames(x),1,1)=="Y"]
  
  # Covariate
  X  <- x[,substr(colnames(x),1,1)=="X"]
  
  if(delta==1){
    # Covariate
    X  <- x[,substr(colnames(x),1,1)=="X"]
    Nx <- dim(X)[2]
    WX <- cbind(W,X)
    gX <- gY
  } else {
    gX <- gY
    WX <- W
  }
  
  T1 <- sum(A)
  T0 <- length(A) - T1
  Tt <- T1+T0
  
  Result <- matrix(0,Tt,dim(gX)[2]+lengthb)
  
  for(tt in 1:T0){
    Result[tt,1:dim(gX)[2]] <- gX[tt,]*(Y[tt] - sum(WX[tt,]*long.gamma))
  }
  
  # for(tt in T0+1:T1){
  #   RESIDUAL <- Y[tt] - sum(WX[tt,]*long.gamma) - ATT(beta,tt-T0)
  #   Result[tt,dim(gX)[2]+1:lengthb] <- RESIDUAL*GRAD.ATT(beta,tt-T0)
  #
  # }
  
  if(ATT.Type=="spline"){
    for(tt in T0+1:T1){
      RESIDUAL <- Y[tt] - sum(WX[tt,]*long.gamma) - sum(Post.Time.Basis[tt-T0,]*beta)
      Result[tt,dim(gX)[2]+1:lengthb] <- -RESIDUAL*Post.Time.Basis[tt-T0,]
      
    }
  }
  
  if(ATT.Type=="constant"){
    for(tt in T0+1:T1){
      RESIDUAL <- Y[tt] - sum(WX[tt,]*long.gamma) - beta
      Result[tt,dim(gX)[2]+1:lengthb] <- -RESIDUAL
    }
  }
  
  
  
  return(Result)
  
}

GMM.Ft.lambda <- function(theta,x,lambda){
  
  # theta
  long.gamma <- theta[-c(1:lengthb)]
  beta  <- theta[1:lengthb]
  
  # gY
  gY <- x[,substr(colnames(x),1,1)=="G"]
  
  # donor
  W  <- x[,substr(colnames(x),1,1)=="W"]
  
  # Time
  A  <- x[,substr(colnames(x),1,1)=="A"]
  
  # Outcome
  Y  <- x[,substr(colnames(x),1,1)=="Y"]
  
  if(delta==1){
    # Covariate
    X  <- x[,substr(colnames(x),1,1)=="X"]
    Nx <- dim(X)[2]
    WX <- cbind(W,X)
    gX <- gY
  } else {
    gX <- gY
    WX <- W
  }
  
  T1 <- sum(A)
  T0 <- length(A) - T1
  Tt <- T1+T0
  
  Result <- matrix(0,Tt,dim(gX)[2]+lengthb+length(long.gamma))
  # Result <- matrix(0,Tt,dim(gX)[2]+lengthb)
  
  for(tt in 1:T0){
    Result[tt,1:dim(gX)[2]] <- gX[tt,]*(Y[tt] - sum(WX[tt,]*long.gamma))
  }
  
  if(ATT.Type=="spline"){
    for(tt in T0+1:T1){
      RESIDUAL <- Y[tt] - sum(WX[tt,]*long.gamma) - sum(Post.Time.Basis[tt-T0,]*beta)
      Result[tt,dim(gX)[2]+1:lengthb] <- -RESIDUAL*Post.Time.Basis[tt-T0,]
      
    }
  }
  
  if(ATT.Type=="constant"){
    for(tt in T0+1:T1){
      RESIDUAL <- Y[tt] - sum(WX[tt,]*long.gamma) - beta
      Result[tt,dim(gX)[2]+1:lengthb] <- -RESIDUAL
    }
  }
  
  
  
  for(tt in 1:(T0+T1)){
    Result[tt,dim(gX)[2]+lengthb+1:length(long.gamma)] <- 10^(lambda/2)*(long.gamma)
  }
  
  return(Result)
  
}


OLS.Ft <- function(theta,x){
  
  # theta
  long.gamma <- theta[-c(1:lengthb)]
  beta  <- theta[1:lengthb]
  
  # donor
  W  <- x[,substr(colnames(x),1,1)=="W"]
  
  # Time
  A  <- x[,substr(colnames(x),1,1)=="A"]
  
  # Outcome
  Y  <- x[,substr(colnames(x),1,1)=="Y"]
  
  # Covariate
  X  <- x[,substr(colnames(x),1,1)=="X"]
  
  if(delta==1){
    # Covariate
    X  <- x[,substr(colnames(x),1,1)=="X"]
    Nx <- dim(X)[2]
    WX <- cbind(W,X)
  } else {
    WX <- W
  }
  
  T1 <- sum(A)
  T0 <- length(A) - T1
  Tt <- T1+T0
  
  Result <- matrix(0,Tt,dim(WX)[2]+lengthb)
  
  for(tt in 1:T0){
    Result[tt,1:dim(WX)[2]] <- WX[tt,]*(Y[tt] - sum(WX[tt,]*long.gamma))
  }
  
  if(ATT.Type=="spline"){
    for(tt in T0+1:T1){
      RESIDUAL <- Y[tt] - sum(WX[tt,]*long.gamma) - sum(Post.Time.Basis[tt-T0,]*beta)
      Result[tt,dim(WX)[2]+1:lengthb] <- -RESIDUAL*Post.Time.Basis[tt-T0,]
      
    }
  }
  
  if(ATT.Type=="constant"){
    for(tt in T0+1:T1){
      RESIDUAL <- Y[tt] - sum(WX[tt,]*long.gamma) - beta
      Result[tt,dim(WX)[2]+1:lengthb] <- -RESIDUAL
    }
  }
  
  return(Result)
  
}



GMM.Ft.Moment <- function(theta,x){
  sum(apply( GMM.Ft(theta,x),2,mean)^2)
}

GMM.Ft.Moment.lambda <- function(theta,x,lambda=0){
  sum(apply( GMM.Ft.lambda(theta,x,lambda),2,mean)^2)
  # sum(apply( GMM.Ft(theta,x),2,mean)^2) + 10^(lambda)*sum(theta[-1]^2)
}

GMM.Ft.Grad <- function(theta,x){
  
  dim.GMM <- length( apply(GMM.Ft(theta,x),2,mean) )
  EPS <- 10^(-6)
  GRAD <- matrix(0,dim.GMM,length(theta))
  
  for(grad.iter in 1:length(theta)){
    perturb <- rep(0,length(theta))
    perturb[grad.iter] <- EPS
    theta.p <- theta + perturb
    theta.n <- theta - perturb
    GRAD[,grad.iter] <- ( apply(GMM.Ft(theta+perturb,x),2,mean) - apply(GMM.Ft(theta-perturb,x),2,mean) )/(2*EPS)
  }
  return(GRAD)
}

GMM.Ft.Grad.lambda <- function(theta,x,lambda=0){
  
  dim.GMM <- length( apply(GMM.Ft.lambda(theta,x,lambda),2,mean) )
  EPS <- 10^(-6)
  GRAD <- matrix(0,dim.GMM,length(theta))
  
  for(grad.iter in 1:length(theta)){
    perturb <- rep(0,length(theta))
    perturb[grad.iter] <- EPS
    theta.p <- theta + perturb
    theta.n <- theta - perturb
    GRAD[,grad.iter] <- ( apply(GMM.Ft.lambda(theta+perturb,x,lambda),2,mean) - apply(GMM.Ft.lambda(theta-perturb,x,lambda),2,mean) )/(2*EPS)
  }
  return(GRAD)
}


CV.Lambda <- function(lambda){
  
  if(delta==0){
    posW <- which(substr(colnames(GMM.Data),1,1)=="W")
    posG <- which(substr(colnames(GMM.Data),1,1)=="G")
  } else {
    posW <- c( which(substr(colnames(GMM.Data),1,1)=="W"),
               which(substr(colnames(GMM.Data),1,1)=="X") )
    posG <- c( which(substr(colnames(GMM.Data),1,1)=="G"),
               which(substr(colnames(GMM.Data),1,1)=="X") )
  }
  
  CV.Residual <- rep(0,T.Pre)
  
  
  for(jj in 1:T.Pre){
    Temp.T <- (1:T.Pre)[-jj]
    X.Temp  <- GMM.Data[Temp.T,posW]
    Y.Temp  <- Y[Temp.T]
    gY.Temp <- GMM.Data[Temp.T,posG]
    GY <- c(apply(gY.Temp*Y.Temp,2,mean))
    GX <- matrix(0,length(posG),length(posW))
    for(loo in 1:length(Temp.T)){
      GX <- GX + gY.Temp[loo,]%*%t(X.Temp[loo,])
    }
    GX <- GX/(length(Temp.T))
    CV.Residual[jj] <-  (Y[jj] - as.numeric(GMM.Data[jj,posW]%*%
                                              my.inverse(A = t(GX)%*%(GX),
                                                         stab.const = ((T.Pre+T.Post)/(T.Pre))^2*10^(lambda),
                                                         adjust = F)%*%
                                              (t(GX)%*%(GY))))
    
  }
  
  # Prod.Res <- GMM.Data[,posG]*CV.Residual
  
  # return( sum(apply(Prod.Res,2,mean)^2) )
  
  return(log(mean(CV.Residual)^2))
  
}

CV.Lambda.Block <- function(lambda){
  
  if(delta==0){
    posW <- which(substr(colnames(GMM.Data),1,1)=="W")
    posG <- which(substr(colnames(GMM.Data),1,1)=="G")
  } else {
    posW <- c( which(substr(colnames(GMM.Data),1,1)=="W"),
               which(substr(colnames(GMM.Data),1,1)=="X") )
    posG <- c( which(substr(colnames(GMM.Data),1,1)=="G"),
               which(substr(colnames(GMM.Data),1,1)=="X") )
  }
  
  CV.Residual <- matrix(0,T.Pre,length(posG))
  
  MS <- 1:(T.Pre/2)
  AS <- (T.Pre/2)+1:(T.Pre/2)
  
  Temp.T <- MS
  
  X.Temp  <- GMM.Data[Temp.T,posW]
  Y.Temp  <- Y[Temp.T]
  gY.Temp <- GMM.Data[Temp.T,posG]
  GY <- c(apply(gY.Temp*Y.Temp,2,mean))
  GX <- matrix(0,length(posG),length(posW))
  for(loo in 1:length(Temp.T)){
    GX <- GX + gY.Temp[loo,]%*%t(X.Temp[loo,])
  }
  GX <- GX/(length(Temp.T))
  
  for(jj in AS){
    CV.Residual[jj,] <- GMM.Data[jj,posG]*
      (Y[jj] - as.numeric(GMM.Data[jj,posW]%*%
                            my.inverse(A = t(GX)%*%(GX),
                                       stab.const = ((T.Pre+T.Post)/(T.Pre))^2*10^(lambda),
                                       adjust = F)%*%
                            (t(GX)%*%(GY))))
  }
  
  Temp.T <- AS
  
  X.Temp  <- GMM.Data[Temp.T,posW]
  Y.Temp  <- Y[Temp.T]
  gY.Temp <- GMM.Data[Temp.T,posG]
  GY <- c(apply(gY.Temp*Y.Temp,2,mean))
  GX <- matrix(0,length(posG),length(posW))
  for(loo in 1:length(Temp.T)){
    GX <- GX + gY.Temp[loo,]%*%t(X.Temp[loo,])
  }
  GX <- GX/(length(Temp.T))
  
  for(jj in MS){
    CV.Residual[jj,] <- GMM.Data[jj,posG]*
      (Y[jj] - as.numeric(GMM.Data[jj,posW]%*%
                            my.inverse(A = t(GX)%*%(GX),
                                       stab.const = ((T.Pre+T.Post)/(T.Pre))^2*10^(lambda),
                                       adjust = F)%*%
                            (t(GX)%*%(GY))))
  }
  
  return( log( sum( apply(CV.Residual,2,mean)^2 ) ) )
  
}


OLS.Ft.Moment <- function(theta,x){
  sum(apply( OLS.Ft(theta,x),2,mean)^2)
}

OLS.Ft.Grad <- function(theta,x){
  
  dim.OLS <- length( apply(OLS.Ft(theta,x),2,mean) )
  EPS <- 10^(-6)
  GRAD <- matrix(0,dim.OLS,length(theta))
  
  for(grad.iter in 1:length(theta)){
    perturb <- rep(0,length(theta))
    perturb[grad.iter] <- EPS
    theta.p <- theta + perturb
    theta.n <- theta - perturb
    GRAD[,grad.iter] <- ( apply(OLS.Ft(theta+perturb,x),2,mean) - apply(OLS.Ft(theta-perturb,x),2,mean) )/(2*EPS)
  }
  return(GRAD)
}










Meat.HAC <- function(Res.Mat,
                     GRAD,
                     BREAD=NULL,
                     beta.pos, 
                     bw.type="auto"){
  
  CROSS <- list()
  for(ll in 0:(T0+T1-1)){
    CROSS[[ll+1]] <- matrix(0,dim(Res.Mat)[2],dim(Res.Mat)[2])
    for(tt in 1:(T0+T1-ll)){
      CROSS[[ll+1]] <- CROSS[[ll+1]]+
        (t(t(Res.Mat[tt,]))%*%(t(Res.Mat[tt+ll,]))) +
        +as.numeric(ll>0)*(t(t(Res.Mat[tt+ll,]))%*%(t(Res.Mat[tt,])))
    }
  }
  
  if(bw.type=="conservative"){
    
    ADR.VAR <- function(bw){
      
      BW <- 10^bw
      
      Meat.Weight.F <- rep(0,T0+T1)
      Meat.Weight.F[1] <- 1
      for(tt in 2:(T0+T1)){
        z <- 6*pi/5*(tt-1)/BW
        Meat.Weight.F[tt] <- 3/(z^2)*(sin(z)/z - cos(z))
      }
      
      Meat.Meat <- matrix(0,dim(Res.Mat)[2],dim(Res.Mat)[2])
      for(ll in 0:(T0+T1-1)){
        Meat.Meat <- Meat.Meat+CROSS[[ll+1]]*Meat.Weight.F[ll+1]
      }
      
      Meat.Meat <- Meat.Meat/(T0+T1)
      
      Meat.Simple <- (t(GRAD))%*%Meat.Meat%*%((GRAD))
      
      GMM.Simple.Var <- (BREAD)%*%(Meat.Simple)%*%t((BREAD))/(T0+T1)
      
      if(lengthb>1){
        return( - norm(GMM.Simple.Var[beta.pos,beta.pos]) )
      } else {
        return( - GMM.Simple.Var[beta.pos,beta.pos] )
      }
    }
    
    OPTIM <- optim(par=0,
                   fn=ADR.VAR)
    opt.BW <- 10^OPTIM$par
    
  } else if (bw.type=="auto"){
    
    rho.Vec    <- rep(0,lengthb) 
    sigma2.Vec <- rep(0,lengthb)
    
    for(beta.iter in 1:lengthb){
      ARIMA <- try( arima(x=Res.Mat[,dim(Res.Mat)[2]-lengthb+beta.iter],
                          order=c(1,0,0)),silent=T )
      if(class(ARIMA)=="try-error"){
        ARIMA <- arima(x=Res.Mat[,dim(Res.Mat)[2]-lengthb+beta.iter],
                       order=c(1,1,0))
      } 
      
      rho.Vec[beta.iter]    <- ARIMA$coef[1]
      sigma2.Vec[beta.iter] <- ARIMA$sigma2
      
      if(abs(ARIMA$coef[1])>0.95){
        rho.Vec[beta.iter]    <- 0.95*sign(ARIMA$coef[1])
        tempy <- Res.Mat[,dim(Res.Mat)[2]-lengthb+beta.iter]
        sigma2.Vec[beta.iter] <- var(tempy[-1] - rho.Vec[beta.iter]*tempy[-length(tempy)])
      }
      
      
    }
    
    alpha1 <- sum( 4*rho.Vec^2*sigma2.Vec^2/(1-rho.Vec)^6/(1+rho.Vec)^2 )/
      sum( sigma2.Vec^2/(1-rho.Vec)^4 )
    
    alpha2 <- sum( 4*rho.Vec^2*sigma2.Vec^2/(1-rho.Vec)^8 )/
      sum( sigma2.Vec^2/(1-rho.Vec)^4 )
    
    opt.BW <- 1.3221*(alpha2*(T1+T0))^(0.2)
    
    opt.Bartlett.BW <- 1.1447*(alpha1*(T1+T0))^(1/3)
    
  }
  
  
  Meat.Weight.F <- rep(0,T0+T1)
  Meat.Weight.F[1] <- 1
  for(tt in 2:(T0+T1)){
    z <- 6*pi/5*(tt-1)/opt.BW
    Meat.Weight.F[tt] <- 3/(z^2)*(sin(z)/z - cos(z))
  }
  
  Meat.Meat <- matrix(0,dim(Res.Mat)[2],dim(Res.Mat)[2])
  for(ll in 0:(T0+T1-1)){
    Meat.Meat <- Meat.Meat+CROSS[[ll+1]]*Meat.Weight.F[ll+1]
  }
  
  Meat.Meat <- Meat.Meat/(T0+T1)
  Meat.Simple <- Meat.Meat
  
  RESULT <- list()
  RESULT$matrix <- Meat.Simple
  RESULT$bw <- opt.BW
  RESULT$bw.B <- opt.Bartlett.BW
  return( RESULT  )
  
} 





BOOT <- function(Time.Window,
                 Num.Boot,
                 T0,
                 T1,
                 Boot.Factor=1,
                 type="SPSC"){
  
  Time.Window <- apply(cbind(1,Time.Window),1,max)
  Time.Window <- unique(Time.Window)
  
  if(type=="SPSC"){
    Data <- GMM.Data 
  } else if (type=="SPSC.Reg"){
    Data <- GMM.Data  
  } else if (type=="OLS"){
    Data <- OLS.Data 
  }
  
  Boot.Grid <- length(Time.Window)
  
  posG <- c( which(substr(colnames( Data),1,1)=="G") )
  posW <- c( which(substr(colnames( Data),1,1)=="W") )
  posY <- c( which(substr(colnames( Data),1,1)=="Y") )
  posX <- c( which(substr(colnames( Data),1,1)=="X") )
  
  if(length(posX)>0){
    posG <- c(posG,posX)
    posW <- c(posW,posX)
  }
  
  Boot <- list()
  for(tw.iter in 1:Boot.Grid){
    Boot[[tw.iter]] <- matrix(0,Num.Boot,lengthb)
  }
  
  for(tw.iter in 1:Boot.Grid){
    
    for(boot.iter in 1:Num.Boot){
      
      T0.Extend  <- c(1:T0,1:T0)
      tw         <- Time.Window[tw.iter]
      Boot.Start <- sample(1:T0,ceiling(Boot.Factor*T0/tw),replace=T)
      Boot.Pre   <- as.vector(sapply(Boot.Start,function(v){ T0.Extend[v:(v+tw-1)] }))[1:(Boot.Factor*T0)]
      
      T1.Extend  <- c(1:T1,1:T1)
      tw         <- Time.Window[tw.iter]
      Boot.Start <- sample(T1,ceiling(Boot.Factor*T1/tw),replace=T)
      Boot.Post  <- as.vector(sapply(Boot.Start,function(v){ T1.Extend[v:(v+tw-1)] }))[1:(Boot.Factor*T1)] + T0
      
      Boot.Data <- cbind(Data,1:(T0+T1))[c(Boot.Pre,Boot.Post),]
      
      if(type=="SPSC"){
        
        Wmat.Boot.Pre <- Boot.Data[1:(Boot.Factor*T0),posW]
        gY.Boot.Pre   <- Boot.Data[1:(Boot.Factor*T0),posG]
        Y1.Boot.Pre   <- Boot.Data[1:(Boot.Factor*T0),posY]
        
        Wmat.Boot.Post <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posW]
        gY.Boot.Post   <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posG]
        Y1.Boot.Post   <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posY]
        
        Time.Post      <- Boot.Data[Boot.Factor*T0+1:(Boot.Factor*T1),dim(Boot.Data)[2]]
        
        GW.Boot <- t( sapply(1:dim(gY.Boot.Pre)[2],function(tt){
          apply( t(Wmat.Boot.Pre)*matrix(gY.Boot.Pre[,tt],dim(Wmat.Boot.Pre)[2],Boot.Factor*T0,byrow=T), 1, mean)  }) )
        
        GY.Boot <- ( sapply(1:dim(gY.Boot.Pre)[2],function(tt){
          mean(gY.Boot.Pre[,tt]*Y1.Boot.Pre)
        }) )
        
        gamma.naive.Boot <- my.inverse(A = t(GW.Boot)%*%(GW.Boot) )%*%(t(GW.Boot)%*%(GY.Boot))
        
        resid.naive.Boot <- as.numeric(Y1.Boot.Post - (Wmat.Boot.Post)%*%gamma.naive.Boot)
        
        if(ATT.Type=="spline"){
          X.Temp <- predict(Post.Time.Basis.Fit, Time.Post-T0)
          beta.naive.Boot <- my.inverse(A = t(X.Temp)%*%(X.Temp) )%*%(t(X.Temp)%*%resid.naive.Boot)
          Boot[[tw.iter]][boot.iter,] <- beta.naive.Boot
        } else {
          Boot[[tw.iter]][boot.iter,] <- mean(resid.naive.Boot)
        }
        
      } else if (type=="SPSC.Reg"){
        
        Wmat.Boot.Pre <- Boot.Data[1:(Boot.Factor*T0),posW]
        gY.Boot.Pre   <- Boot.Data[1:(Boot.Factor*T0),posG]
        Y1.Boot.Pre   <- Boot.Data[1:(Boot.Factor*T0),posY]
        
        Wmat.Boot.Post <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posW]
        gY.Boot.Post   <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posG]
        Y1.Boot.Post   <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posY]
        
        Time.Post      <- Boot.Data[Boot.Factor*T0+1:(Boot.Factor*T1),dim(Boot.Data)[2]]
        
        GW.Boot <- t( sapply(1:dim(gY.Boot.Pre)[2],function(tt){
          apply( t(Wmat.Boot.Pre)*matrix(gY.Boot.Pre[,tt],dim(Wmat.Boot.Pre)[2],Boot.Factor*T0,byrow=T), 1, mean)  }) )
        
        GY.Boot <- ( sapply(1:dim(gY.Boot.Pre)[2],function(tt){
          mean(gY.Boot.Pre[,tt]*Y1.Boot.Pre)
        }) )
        
        gamma.naive.Boot <- my.inverse(A = t(GW.Boot)%*%(GW.Boot),
                                       stab.const = ((T.Pre+T.Post)/(T.Pre))^2*10^(lambda.opt),
                                       adjust = F)%*%(t(GW.Boot)%*%(GY.Boot))
        
        resid.naive.Boot <- as.numeric(Y1.Boot.Post - (Wmat.Boot.Post)%*%gamma.naive.Boot)
        
        if(ATT.Type=="spline"){
          X.Temp <- predict(Post.Time.Basis.Fit, Time.Post-T0)
          beta.naive.Boot <- my.inverse(A = t(X.Temp)%*%(X.Temp) )%*%(t(X.Temp)%*%resid.naive.Boot)
          Boot[[tw.iter]][boot.iter,] <- beta.naive.Boot
        } else {
          Boot[[tw.iter]][boot.iter,] <- mean(resid.naive.Boot)
        }
        
      } else if (type=="OLS"){
        
        Wmat.Boot.Pre <- Boot.Data[1:(Boot.Factor*T0),posW]
        Y1.Boot.Pre   <- Boot.Data[1:(Boot.Factor*T0),posY]
        
        Wmat.Boot.Post <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posW]
        Y1.Boot.Post   <- Boot.Data[(Boot.Factor*T0)+1:(Boot.Factor*T1),posY]
        
        Time.Post      <- Boot.Data[Boot.Factor*T0+1:(Boot.Factor*T1),dim(Boot.Data)[2]]
        
        gamma.naive.Boot <- my.inverse(A = t(Wmat.Boot.Pre)%*%(Wmat.Boot.Pre) )%*%(t(Wmat.Boot.Pre)%*%(Y1.Boot.Pre))
        
        resid.naive.Boot <- as.numeric(Y1.Boot.Post - (Wmat.Boot.Post)%*%gamma.naive.Boot)
        
        if(ATT.Type=="spline"){
          X.Temp <- predict(Post.Time.Basis.Fit, Time.Post-T0)
          beta.naive.Boot <- my.inverse(A = t(X.Temp)%*%(X.Temp) )%*%(t(X.Temp)%*%resid.naive.Boot)
          Boot[[tw.iter]][boot.iter,] <- beta.naive.Boot
        } else {
          Boot[[tw.iter]][boot.iter,] <- mean(resid.naive.Boot)
        }
        
      }
      
      
    }
    
  }
  
  return(Boot)
  
}





