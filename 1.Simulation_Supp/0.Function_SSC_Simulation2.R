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
                             T0,
                             T1,
                             p=1){
  Residual.Extend <- c(Residual,Residual)
  if(p==1){
    S.Base <- mean(abs(Residual[T0+1:T1]))
    S.All  <- rep(0,T0+T1)
    for(bb in 1:(T0+T1)){
      S.All[bb] <- mean( abs(Residual.Extend[bb-1+1:T1]) )
    }
  } else {
    S.Base <- abs(mean(Residual[T0+1:T1]))
    S.All  <- rep(0,T0+T1)
    for(bb in 1:(T0+T1)){
      S.All[bb] <- abs(mean(Residual.Extend[bb-1+1:T1]) )
    }
  }
  
  return( 1-mean(S.All < S.Base) )
}

COVER <- function(a,b,center=3){
  as.numeric(abs(a-center)<=b*qnorm(0.975))
}

Conformal.Prediction <- function(Wmat.Pre,
                                 Wmat.Post,
                                 Y1.Pre,
                                 Y1.Post,
                                 Xmat.Pre=NULL,
                                 Xmat.Post=NULL,
                                 cov.ind=0,
                                 center,
                                 bw,
                                 lambda=-Inf,
                                 constant=0){
  
  T0 <- dim(Wmat.Pre)[1]
  T1 <- dim(Wmat.Post)[1]
  Tt <- T0 + T1
  N  <- dim(Wmat.Pre)[2]
  A  <- rep(c(0,1),c(T0,T1))
  
  ## GMM-1st
  
  PValue.beta <- function(ii,beta){
    
    Wmat       <- rbind(Wmat.Pre,Wmat.Post[ii,])
    if(cov.ind==1){
      Xmat       <- rbind(Xmat.Pre,Xmat.Post[ii,])
    }
    Y          <- c(Y1.Pre,Y1.Post[ii])
    
    beta.seq <- c(rep(0,T0),beta)
    
    Y.Pred.H0  <- Y - beta.seq
    
    if(cov.ind==0){
      gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m) 
    } else {
      gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m+N+1) 
    }
    
    BS.gY <- bs(Y.Pred.H0[1:T0],df=m,
                Boundary.knots = range(Y.Pred.H0[1:T0])+c(-2,2))
    
    gY.Pred.H0[1:T0,1:m] <- BS.gY
    gY.Pred.H0[T0+1,1:m] <- predict(BS.gY,Y.Pred.H0[1+T0])
    
    Outside <- !(min(attributes(BS.gY)$Boundary.knots) <= Y.Pred.H0[T0+1] &
                   Y.Pred.H0[T0+1] <= max(attributes(BS.gY)$Boundary.knots))
    Notvalid <- sum(1-as.numeric(abs(gY.Pred.H0[T0+1,1:m]-0.5)<=0.5)) > 0
    
    
    if(Outside | Notvalid){
      gY.Pred.H0[T0+1,1:m] <- 0
    } 
    
    if(cov.ind==1){
      gY.Pred.H0[,m+1:(N+1)] <- Xmat
    }
    
    if(cov.ind==0){
      GW <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
        apply( t(Wmat)*matrix(gY.Pred.H0[,tt],N,T0+1,byrow=T), 1, mean)
      }) )
    } else {
      GW <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
        apply( t(cbind(Wmat,Xmat))*matrix(gY.Pred.H0[,tt],N+N+1,T0+1,byrow=T), 1, mean)
      }) )
    }
    
    GY <- ( sapply(1:dim(gY.Pred.H0)[2],function(tt){
      mean(gY.Pred.H0[,tt]*Y.Pred.H0)
    }) )
    
    gamma.naive <- c( ginv(t(GW)%*%GW + constant*10^(lambda)*diag(rep(1,dim(GW)[2])))%*%(t(GW)%*%GY) )
    if(cov.ind==0){
      Residual <- c(Y.Pred.H0 - Wmat%*%gamma.naive)
    } else {
      Residual <- c(Y.Pred.H0 - cbind(Wmat,Xmat)%*%gamma.naive)
    }
    Calculate.PValue(Residual,T0,1)
  }
  
  CI <- matrix(0,T1,2)
  for(ii in 1:T1){
    beta.grid <- seq(-10,20,length=3001)
    gsize <- diff(beta.grid)[1]
    PV <- sapply(beta.grid,
                 function(vv){PValue.beta(ii,vv)})
    valid <- which(PV>=0.05)
    beta.valid <- beta.grid[valid[which(c(1,diff(valid))==1)]]
    PV.valid <- PV[valid[which(c(1,diff(valid))==1)]]
    
    BLOCKS <- list()
    JUMP <- c(0,which(diff(beta.valid)>1.1*gsize),length(beta.valid))
    for(jindex in 1:(length(JUMP)-1)){
      BLOCKS[[jindex]] <- (JUMP[jindex]+1):JUMP[jindex+1]
    }
    
    max.block <- which.max(sapply(1:length(BLOCKS),
                                  function(vv){ max(PV.valid[BLOCKS[[vv]]]) }))
    
    CI[ii,] <- range(beta.valid[BLOCKS[[max.block]]])
  }
  
  
  return(CI)
  
}




Conformal.Prediction.Time <- function(Wmat.Pre,
                                      Wmat.Post,
                                      Y1.Pre,
                                      Y1.Post,
                                      Xmat.Pre=NULL,
                                      Xmat.Post=NULL,
                                      cov.ind=0,
                                      center,
                                      bw,
                                      lambda=-Inf,
                                      constant=0){
  
  T0 <- dim(Wmat.Pre)[1]
  T1 <- dim(Wmat.Post)[1]
  Tt <- T0 + T1
  N  <- dim(Wmat.Pre)[2]
  A  <- rep(c(0,1),c(T0,T1))
  
  ## GMM-1st
  
  PValue.beta <- function(ii,beta){
    
    Wmat       <- rbind(Wmat.Pre,Wmat.Post[ii,])
    if(cov.ind==1){
      Xmat       <- rbind(Xmat.Pre,Xmat.Post[ii,])
    }
    Y          <- c(Y1.Pre,Y1.Post[ii])
    
    beta.seq <- c(rep(0,T0),beta)
    
    Y.Pred.H0  <- Y - beta.seq
    
    gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m+mT) 
    
    BS.gY <- bs(Y.Pred.H0[1:T0],df=m,
                Boundary.knots = range(Y.Pred.H0[1:T0])+c(-2,2))
    
    gY.Pred.H0[1:T0,1:(m)] <- BS.gY
    gY.Pred.H0[T0+1,1:(m)] <- predict(BS.gY,Y.Pred.H0[1+T0])
    
    if(mT>0){ 
      gY.Pred.H0[1:T0,(m)+1:(mT)] <- bT
      if(bT.Type=="bs"){
        gY.Pred.H0[T0+1,(m)+1:(mT)] <- predict(bT,T0+1)
      } else if (bT.Type=="Constant"){
        gY.Pred.H0[T0+1,(m)+1:(mT)] <- 1
      } else if (bT.Type=="Linear"){
        gY.Pred.H0[T0+1,(m)+1:(mT)] <- c(1,T0+1)
      }
    }
    
    
    
    Outside <- !(min(attributes(BS.gY)$Boundary.knots) <= Y.Pred.H0[T0+1] &
                   Y.Pred.H0[T0+1] <= max(attributes(BS.gY)$Boundary.knots))
    Notvalid <- sum(1-as.numeric(abs(gY.Pred.H0[T0+1,1:m]-0.5)<=0.5)) > 0
    
    
    if(Outside | Notvalid){
      gY.Pred.H0[T0+1,1:m] <- 0
    } 
    
    if(cov.ind==1){
      gY.Pred.H0[,m+1:(N+1)] <- Xmat
    }
    
    if(cov.ind==0){
      GW <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
        apply( t(Wmat)*matrix(gY.Pred.H0[,tt],N,T0+1,byrow=T), 1, mean)
      }) )
    } else {
      GW <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
        apply( t(cbind(Wmat,Xmat))*matrix(gY.Pred.H0[,tt],N+N+1,T0+1,byrow=T), 1, mean)
      }) )
    }
    
    GY <- ( sapply(1:dim(gY.Pred.H0)[2],function(tt){
      mean(gY.Pred.H0[,tt]*Y.Pred.H0)
    }) )
    
    gamma.naive <- c( ginv(t(GW)%*%GW + constant*10^(lambda)*diag(rep(1,dim(GW)[2])))%*%(t(GW)%*%GY) )
    if(cov.ind==0){
      Residual <- c(Y.Pred.H0 - Wmat%*%gamma.naive)
    } else {
      Residual <- c(Y.Pred.H0 - cbind(Wmat,Xmat)%*%gamma.naive)
    }
    Calculate.PValue(Residual,T0,1)
  }
  
  CI <- matrix(0,T1,2)
  for(ii in 1:T1){
    beta.grid <- seq(-10,20,length=3001)
    gsize <- diff(beta.grid)[1]
    PV <- sapply(beta.grid,
                 function(vv){PValue.beta(ii,vv)})
    valid <- which(PV>=0.05)
    beta.valid <- beta.grid[valid[which(c(1,diff(valid))==1)]]
    PV.valid <- PV[valid[which(c(1,diff(valid))==1)]]
    
    BLOCKS <- list()
    JUMP <- c(0,which(diff(beta.valid)>1.1*gsize),length(beta.valid))
    for(jindex in 1:(length(JUMP)-1)){
      BLOCKS[[jindex]] <- (JUMP[jindex]+1):JUMP[jindex+1]
    }
    
    max.block <- which.max(sapply(1:length(BLOCKS),
                                  function(vv){ max(PV.valid[BLOCKS[[vv]]]) }))
    
    CI[ii,] <- range(beta.valid[BLOCKS[[max.block]]])
  }
  
  
  return(CI)
  
}



Conformal.Prediction.Fast <- function(Wmat.Pre,
                                      Wmat.Post,
                                      Y1.Pre,
                                      Y1.Post,
                                      Xmat.Pre=NULL,
                                      Xmat.Post=NULL,
                                      cov.ind=0,
                                      effect,
                                      lambda=-Inf,
                                      constant=0){
  
  T0 <- dim(Wmat.Pre)[1]
  T1 <- dim(Wmat.Post)[1]
  Tt <- T0 + T1
  N  <- dim(Wmat.Pre)[2]
  A  <- rep(c(0,1),c(T0,T1))
  
  ## GMM-1st
  
  PValue.beta <- function(ii,beta){
    
    Wmat       <- rbind(Wmat.Pre,Wmat.Post[ii,])
    if(cov.ind==1){
      Xmat       <- rbind(Xmat.Pre,Xmat.Post[ii,])
    }
    Y          <- c(Y1.Pre,Y1.Post[ii])
    
    beta.seq <- c(rep(0,T0),beta)
    
    Y.Pred.H0  <- Y - beta.seq
    
    if(cov.ind==0){
      gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m) 
    } else {
      gY.Pred.H0 <- matrix(0,length(Y.Pred.H0),m+N+1) 
    }
    
    BS.gY <- bs(Y.Pred.H0[1:T0],df=m,
                Boundary.knots = range(Y.Pred.H0[1:T0])+c(-2,2))
    
    gY.Pred.H0[1:T0,1:m] <- BS.gY
    gY.Pred.H0[T0+1,1:m] <- predict(BS.gY,Y.Pred.H0[1+T0])
    
    Outside <- !(min(attributes(BS.gY)$Boundary.knots) <= Y.Pred.H0[T0+1] &
                   Y.Pred.H0[T0+1] <= max(attributes(BS.gY)$Boundary.knots))
    Notvalid <- sum(1-as.numeric(abs(gY.Pred.H0[T0+1,1:m]-0.5)<=0.5)) > 0
    
    if(Outside | Notvalid){
      gY.Pred.H0[T0+1,1:m] <- 0
    } 
    
    if(cov.ind==1){
      gY.Pred.H0[,m+1:(N+1)] <- Xmat
    }
    
    if(cov.ind==0){
      GW <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
        apply( t(Wmat)*matrix(gY.Pred.H0[,tt],N,T0+1,byrow=T), 1, mean)
      }) )
    } else {
      GW <- t( sapply(1:dim(gY.Pred.H0)[2],function(tt){
        apply( t(cbind(Wmat,Xmat))*matrix(gY.Pred.H0[,tt],N+N+1,T0+1,byrow=T), 1, mean)
      }) )
    }
    
    GY <- ( sapply(1:dim(gY.Pred.H0)[2],function(tt){
      mean(gY.Pred.H0[,tt]*Y.Pred.H0)
    }) )
    
    gamma.naive <- c( ginv(t(GW)%*%GW + constant*10^(lambda)*diag(rep(1,dim(GW)[2])))%*%(t(GW)%*%GY) )
    if(cov.ind==0){
      Residual <- c(Y.Pred.H0 - Wmat%*%gamma.naive)
    } else {
      Residual <- c(Y.Pred.H0 - cbind(Wmat,Xmat)%*%gamma.naive)
    }
    Calculate.PValue(Residual,T0,1)
  }
  
  COVER <- rep(0,T1)
  for(ii in 1:T1){
    COVER[ii] <- as.numeric(PValue.beta(ii,effect[ii])>=0.05)
  }
  
  return(COVER)
  
}



# 
# 
# Conformal.Prediction.OLS <- function(Wmat.Pre,
#                                      Wmat.Post,
#                                      Y1.Pre,
#                                      Y1.Post,
#                                      Xmat.Pre=NULL,
#                                      Xmat.Post=NULL,
#                                      cov.ind=0,
#                                      center){
#   
#   T0 <- dim(Wmat.Pre)[1]
#   T1 <- dim(Wmat.Post)[1]
#   Tt <- T0 + T1
#   N  <- dim(Wmat.Pre)[2]
#   A  <- rep(c(0,1),c(T0,T1))
#   
#   ## GMM-1st
#   
#   PValue.beta <- function(ii,beta){
#     
#     Wmat       <- rbind(Wmat.Pre,Wmat.Post[ii,])
#     if(cov.ind==1){
#       Xmat       <- rbind(Xmat.Pre,Xmat.Post[ii,])
#     }
#     Y          <- c(Y1.Pre,Y1.Post[ii])
#     
#     beta.seq <- c(rep(0,T0),beta)
#     
#     Y.Pred.H0  <- Y - beta.seq
#     
#     
#     if(cov.ind==0){
#       gamma.naive <- lm(Y.Pred.H0~0+Wmat)$coefficients
#       Residual <- c(Y.Pred.H0 - Wmat%*%gamma.naive)
#     } else {
#       gamma.naive <- lm(Y.Pred.H0~0+Wmat+Xmat)$coefficients
#       Residual <- c(Y.Pred.H0 - cbind(Wmat,Xmat)%*%gamma.naive)
#     }
#     Calculate.PValue(Residual,T0,1)
#   }
#   
#   BOUND <- function(ii,
#                     center){
#     
#     GRID <- seq( center - 10,
#                  center + 10,
#                  length=11)
#     
#     pv <- sapply(GRID,function(vv){PValue.beta(ii,vv)})
#     
#     CEN <- GRID[which.max(pv)]
#     
#     LB <- optimize(f=function(beta){
#       (PValue.beta(ii, beta)-0.05)^2
#     },
#     interval=c(center-10,CEN))$minimum
#     
#     UB <- optimize(f=function(beta){
#       (PValue.beta(ii, beta)-0.05)^2
#     },
#     interval=c(CEN,center+10))$minimum
#     
#     return(c(LB,UB))
#   }
#   
#   
#   CI <- matrix(0,T1,2)
#   for(ii in 1:T1){
#     CI[ii,] <- BOUND(ii,center[ii])
#   }
#   
#   return(CI)
#   
# }



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
  
  if(delta==1){
    # Covariate
    X  <- x[,substr(colnames(x),1,1)=="X"]
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
  
  for(tt in T0+1:T1){
    if(lengthb==1){
      Result[tt,dim(gX)[2]+1] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1])
    } else if(lengthb==2){
      Result[tt,dim(gX)[2]+1] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1]-beta[2]*(tt-T0)/T1)
      Result[tt,dim(gX)[2]+2] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1]-beta[2]*(tt-T0)/T1)*((tt-T0)/T1)
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
  
  for(tt in T0+1:T1){
    if(lengthb==1){
      Result[tt,dim(gX)[2]+1] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1])
    } else if(lengthb==2){
      Result[tt,dim(gX)[2]+1] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1]-beta[2]*(tt-T0)/T1)
      Result[tt,dim(gX)[2]+2] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1]-beta[2]*(tt-T0)/T1)*((tt-T0)/T1)
    }
  }
  
  for(tt in 1:(T0+T1)){
    if(lengthb==1){
      # constant term returs error, so add a very tiny noise
      Result[tt,dim(gX)[2]+1+1:length(long.gamma)] <- 10^(lambda/2)*(long.gamma)
    } else if(lengthb==2){
      # constant term returs error, so add a very tiny noise
      Result[tt,dim(gX)[2]+2+1:length(long.gamma)] <- 10^(lambda/2)*(long.gamma)
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
  EPS <- 10^(-8)
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
  EPS <- 10^(-8)
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
                                              (solve( t(GX)%*%(GX) + ((T.Pre+T.Post)/(T.Pre))^2*10^(lambda)*diag(rep(1,dim(GX)[2])) )%*%(t(GX)%*%(GY)))))
    
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
                            (solve( t(GX)%*%(GX) + ((T.Pre+T.Post)/(T.Pre/2))^2*10^(lambda)*diag(rep(1,dim(GX)[2])) )%*%(t(GX)%*%(GY)))))
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
                            (solve( t(GX)%*%(GX) + ((T.Pre+T.Post)/(T.Pre/2))^2*10^(lambda)*diag(rep(1,dim(GX)[2])) )%*%(t(GX)%*%(GY)))))
  }
  
  return( log( sum( apply(CV.Residual,2,mean)^2 ) ) )
  
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
  
  if(delta==1){
    # Covariate
    X  <- x[,substr(colnames(x),1,1)=="X"]
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
  
  for(tt in T0+1:T1){
    if(lengthb==1){
      Result[tt,dim(WX)[2]+1] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1])
    } else if(lengthb==2){
      Result[tt,dim(WX)[2]+1] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1]-beta[2]*(tt-T0)/T1)
      Result[tt,dim(WX)[2]+2] <- -(Y[tt] - sum(WX[tt,]*long.gamma)-beta[1]-beta[2]*(tt-T0)/T1)*((tt-T0)/T1)
    }
  }
  
  return(Result)
  
}


OLS.Ft.Moment <- function(theta,x){
  sum(apply( OLS.Ft(theta,x),2,mean)^2)
}

OLS.Ft.Grad <- function(theta,x){
  
  dim.OLS <- length( apply(OLS.Ft(theta,x),2,mean) )
  EPS <- 10^(-8)
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
      
      GMM.Simple.Var <- ginv(BREAD)%*%(Meat.Simple)%*%t(ginv(BREAD))/(T0+T1) 
      
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
      ARIMA <- arima(x=Res.Mat[,dim(Res.Mat)[2]-lengthb+beta.iter],
                     order=c(1,0,0))
      rho.Vec[beta.iter]    <- ARIMA$coef[1]
      sigma2.Vec[beta.iter] <- ARIMA$sigma2
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
        
        if(lengthb>1){
          X.Temp <- cbind(1,(1:T1)/T1)
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
        
        if(lengthb>1){
          X.Temp <- cbind(1,(1:T1)/T1)
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
        
        if(lengthb>1){
          X.Temp <- cbind(1,(1:T1)/T1)
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

