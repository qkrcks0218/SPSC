SPSC <- function(
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
    conformal.period = 1,                                # NULL or post-treatment time period
    conformal.cover = TRUE,                              # cover=0/1, interval=prediction interval
    true.effect = c(3),                                  # true.effect
    conformal.interval = TRUE,                           # cover=0/1, interval=prediction interval
    conformal.pvalue = 0.05,                             # level of conformal inference
    conformal.window = 25                                # range of grid search for prediction intervals
){
  
  T0 <- length(Y.Pre)
  T1 <- length(Y.Post)
  Tt <- T0+T1
  
  ## check detrend
  
  if(detrend){
    if(is.null(detrend.ft)){
      warning("No function is specified for detrending, linear function function(t){c(1,t)} is used")
      detrend.ft = function(t){c(1,t)}
    }
    detrend.basis <- make.matrix( t(sapply(1:Tt,detrend.ft)), Tt )
    detrend.basis.Pre <- make.matrix( t(sapply(1:T0,detrend.ft)), T0 )
  }
  
  ## check lambda
  
  if(lambda.type=="fix"&is.null(lambda.value)){
    warning("No lambda value is specified, 10^lambda is selected as 10^(-6)")
    lambda.value <- 10^(-6)
  }
  
  att.basis.Post <- make.matrix( t(sapply(1:T1,att.ft)), T1 )
  
  ## check conformal.cover
  if(conformal.cover){
    if(is.null(true.effect)){
      warning("No true effect is specified despite conformal.cover=TRUE\n Coverage cannot be calculated")
      conformal.period <- NULL
    }
    if(length(conformal.period)!=length(true.effect)){
      warning("length(true effect)!=length(conformal.period)")
      nl <- min(length(true.effect),length(conformal.period))
      true.effect <- true.effect[1:nl]
      conformal.period <- conformal.period[1:nl]
    }
  }
  
  ## check conformal.interval
  if(conformal.interval){
    if(length(conformal.period)==0){
      warning("No period is specified despite conformal.interval=TRUE")
      conformal.period <- 1
    }
  }
  
  ###### GMM Data
  
  if(detrend){
    GMM.Data <- cbind( detrend.basis,
                       rbind( matrix(0,T0,ncol(att.basis.Post)), att.basis.Post ),
                       rbind(W.Pre, W.Post),
                       rep(c(0,1),c(T0,T1)),
                       c(Y.Pre, Y.Post) )
  } else {
    GMM.Data <- cbind( rbind( matrix(0,T0,ncol(att.basis.Post)), att.basis.Post ),
                       rbind(W.Pre, W.Post),
                       rep(c(0,1),c(T0,T1)),
                       c(Y.Pre, Y.Post) )
  }
  
  GMM.Data <- as.matrix(GMM.Data)
  
  if(detrend){
    d <- ncol(detrend.basis)
  }
  b <- ncol(att.basis.Post)
  N <- ncol(W.Pre)
  
  if(detrend){
    colnames(GMM.Data) <-   c(sprintf("D%0.4d",1:d),
                              sprintf("B%0.4d",1:b),
                              sprintf("W%0.4d",1:N),
                              "A",
                              "Y")
  } else {
    colnames(GMM.Data) <-   c(sprintf("B%0.4d",1:b),
                              sprintf("W%0.4d",1:N),
                              "A",
                              "Y")
  }
  
  GMM.Data <- data.frame(GMM.Data)
  
  for(jj in 1:ncol(GMM.Data)){
    GMM.Data[,jj] <- as.numeric( GMM.Data[,jj] )
  }
  
  pos.D <- which(substr(colnames(GMM.Data),1,1)=="D")
  pos.B <- which(substr(colnames(GMM.Data),1,1)=="B")
  pos.W <- which(substr(colnames(GMM.Data),1,1)=="W")
  pos.Y <- which(substr(colnames(GMM.Data),1,1)=="Y")
  
  if(detrend){
    Factor <- mean(apply(GMM.Data[1:T0,pos.W],2,sd)) / mean(apply(GMM.Data[1:T0,pos.D],2,sd))
    GMM.Data[,pos.D] <- Factor * GMM.Data[,pos.D]
  }
  
  ###### effect calculate
  
  theta.estimate <- SPSC.Effect(GMM.Data,
                                Y.basis,
                                lambda.type,
                                lambda.value,
                                lambda.grid,
                                calculate.beta=TRUE)
  
  if(detrend){
    dim.detrend <- length(theta.estimate$detrend)
  }
  dim.gamma <- length(theta.estimate$gamma)
  dim.beta <- length(theta.estimate$beta)
  lambda.opt <- theta.estimate$lambda
  
  ###### avar
  
  Grad.Psi <- Grad.Psi.Aver.Ft(theta.estimate, GMM.Data, Y.basis)
  Psi <- Psi.Ft(theta.estimate, GMM.Data, Y.basis)
  
  if(detrend){
    Grad.lambda <- matrix(0,dim.detrend+dim.gamma+dim.beta,dim.detrend+dim.gamma+dim.beta)
    Grad.lambda[dim.detrend+1:dim.gamma,
                dim.detrend+1:dim.gamma] <- diag(rep(1,dim.gamma))
  } else {
    Grad.lambda <- matrix(0,
                          dim.gamma+dim.beta,dim.gamma+dim.beta)
    Grad.lambda[1:dim.gamma,
                1:dim.gamma] <- diag(rep(1,dim.gamma))
  }
  
  S1 <- t(Grad.Psi)%*%(Grad.Psi) + (Tt/T0)^2*10^(theta.estimate$lambda) * Grad.lambda
  if(detrend){
    Sigma <- HAC.Meat(Psi,dim.detrend+dim.gamma+1:dim.beta,T0,T1)
  } else {
    Sigma <- HAC.Meat(Psi,dim.gamma+1:dim.beta,T0,T1)
  }
  
  S2 <- t(Grad.Psi)%*%(Sigma$matrix)%*%(Grad.Psi)
  
  AVAR <- MASS::ginv(S1)%*%S2%*%t(MASS::ginv(S1))/Tt
  if(detrend){
    AVAR.beta <- AVAR[dim.detrend+dim.gamma+1:dim.beta,dim.detrend+dim.gamma+1:dim.beta]
  } else {
    AVAR.beta <- AVAR[dim.gamma+1:dim.beta,dim.gamma+1:dim.beta]
  }
  
  ATT <- as.numeric( matrix(as.matrix(GMM.Data[T0+1:T1,pos.B]),T1,dim.beta)%*%theta.estimate$beta )
  AVAR.ATT <- diag( (matrix(as.matrix(GMM.Data[T0+1:T1,pos.B]),T1,dim.beta))%*%AVAR.beta%*%t((matrix(as.matrix(GMM.Data[T0+1:T1,pos.B]),T1,dim.beta))) )
  ASE.ATT <- sqrt(AVAR.ATT)
  
  RESULT <- list()
  
  RESULT$gamma <- theta.estimate$gamma
  RESULT$SC <- as.numeric(as.matrix(GMM.Data[1:Tt,pos.W])%*%theta.estimate$gamma)
  RESULT$Y <- as.numeric(GMM.Data$Y[1:Tt])
  
  RESULT$effect <- as.numeric(RESULT$Y-RESULT$SC)
  RESULT$ATT <- ATT
  RESULT$ASE.ATT <- ASE.ATT
  
  ###### bootstrap
  
  if(bootstrap.num>1){
    
    T0.Extend  <- c(1:T0,1:T0)
    T1.Extend  <- c(1:T1,1:T1)
    tw <- round(max(Sigma$bw.B,1))
    
    boot.mat <- matrix(0,bootstrap.num,dim.beta)
    
    for(boot.iter in 1:bootstrap.num){
      
      Boot.Start <- sample(T0,ceiling(T0/tw),replace=T)
      Boot.Pre   <- as.vector(sapply(Boot.Start,function(v){ T0.Extend[v:(v+tw-1)] }))[1:(T0)]
      
      Boot.Start <- sample(T1,ceiling(T1/tw),replace=T)
      Boot.Post  <- as.vector(sapply(Boot.Start,function(v){ T1.Extend[v:(v+tw-1)] }))[1:(T1)] + T0
      
      GMM.Data.Boot <- GMM.Data[c(Boot.Pre,Boot.Post),]
      
      boot.mat[boot.iter,] <- SPSC.Effect(GMM.Data.Boot,
                                          Y.basis,
                                          lambda.type="fix",
                                          lambda.value=theta.estimate$lambda,
                                          lambda.grid=NULL,
                                          calculate.beta=TRUE)$beta
    }
    
    RESULT$BSE.ATT <- apply(matrix(as.matrix(GMM.Data[T0+1:T1,pos.B]),T1,dim.beta)%*%t(boot.mat),1,sd)
  }
  
  if(!is.null(conformal.period)){
    
    valid.pvalue <- max(((1:T0+1)/(T0+1))[(1:T0+1)/(T0+1)<=conformal.pvalue])
    
    if(conformal.cover){
      
      GMM.Data.Conformal <- GMM.Data
      GMM.Data.Conformal[T0+conformal.period,]$Y <- GMM.Data.Conformal[T0+conformal.period,]$Y - true.effect
      GMM.Data.Conformal$A <- 0
      
      conformal.cover.vec <- rep(0,length(conformal.period))
      
      for(conformal.iter in 1:length(conformal.period)){
        GMM.Data.Conformal.Iter <- GMM.Data.Conformal[c(1:T0,T0+conformal.period[conformal.iter]),]
        conformal.gamma <- SPSC.Effect(GMM.Data.Conformal.Iter,
                                       Y.basis,
                                       lambda.type="fix",
                                       lambda.value=theta.estimate$lambda,
                                       calculate.beta=FALSE)$gamma
        
        conformal.residual <-
          as.numeric(GMM.Data.Conformal.Iter$Y - as.matrix(GMM.Data.Conformal.Iter[,pos.W])%*%conformal.gamma)
        
        conformal.cover.vec[conformal.iter] <- Calculate.PValue(conformal.residual,T0,1,p=1)
        
      }
      
      RESULT$conformal.cover <- as.numeric(valid.pvalue<=conformal.cover.vec)
      
    }
    
    if(conformal.interval){
      
      GMM.Data.Conformal <- GMM.Data
      GMM.Data.Conformal$A <- 0
      
      conformal.interval.vec <- matrix(0,length(conformal.period),2)
      
      for(conformal.iter in 1:length(conformal.period)){
        
        beta.point <- GMM.Data.Conformal[T0+conformal.period[conformal.iter],pos.Y] -
          sum(as.numeric(GMM.Data.Conformal[T0+conformal.period[conformal.iter],pos.W])*theta.estimate$gamma)
        
        window.unit <- ASE.ATT[conformal.iter]
        
        if(is.na(window.unit)){
          window.unit <- sd(RESULT$Y - RESULT$SC)
        }
        
        beta.grid <- beta.point + seq(-conformal.window,conformal.window,length=101)*window.unit
        
        pvalue.vec <-
          sapply(1:length(beta.grid),
                 function(beta.index){
                   GMM.Data.Conformal.Iter <- GMM.Data.Conformal[c(1:T0,T0+conformal.period[conformal.iter]),]
                   GMM.Data.Conformal.Iter$Y[T0+1] <-
                     GMM.Data.Conformal.Iter$Y[T0+1] - beta.grid[beta.index]
                   
                   conformal.gamma <- SPSC.Effect(GMM.Data.Conformal.Iter,
                                                  Y.basis,
                                                  lambda.type="fix",
                                                  lambda.value=theta.estimate$lambda,
                                                  calculate.beta=FALSE)$gamma
                   
                   conformal.residual <-
                     as.numeric(GMM.Data.Conformal.Iter$Y - as.matrix(GMM.Data.Conformal.Iter[,pos.W])%*%conformal.gamma)
                   Calculate.PValue(conformal.residual,T0,1,p=1)
                 })
        
        position <- which(pvalue.vec>=valid.pvalue)
        d.position <- diff(position)
        breaks <- c(0, which(d.position != 1), length(position))
        consecutive_chunks <- lapply(seq_along(breaks[-1]),
                                     function(i) position[(breaks[i] + 1):breaks[i + 1]])
        valid.gp <- which(sapply(1:length(consecutive_chunks),
                                 function(tt){
                                   51%in%consecutive_chunks[[tt]]
                                 }))
        
        real.beta.grid <- consecutive_chunks[[valid.gp]]
        
        blu <- beta.grid[min(real.beta.grid)]
        bul <- beta.grid[max(real.beta.grid)]
        
        if( min(real.beta.grid)>1 ){
          bll <- beta.grid[which(beta.grid==blu)-1]
        } else {
          bll <- blu - 10*window.unit
        }
        
        if( max(real.beta.grid)<101 ){
          buu <- beta.grid[which(beta.grid==bul)+1]
        } else {
          buu <- bul + 10*window.unit
        }
        
        beta.grid.narrow <- c(seq(bll,blu,length=51),
                              seq(bul,buu,length=51))
        
        pvalue.vec.narrow <-
          sapply(1:length(beta.grid.narrow),
                 function(beta.index){
                   GMM.Data.Conformal.Iter <- GMM.Data.Conformal[c(1:T0,T0+conformal.period[conformal.iter]),]
                   GMM.Data.Conformal.Iter$Y[T0+1] <-
                     GMM.Data.Conformal.Iter$Y[T0+1] - beta.grid.narrow[beta.index]
                   
                   conformal.gamma <- SPSC.Effect(GMM.Data.Conformal.Iter,
                                                  Y.basis,
                                                  lambda.type="fix",
                                                  lambda.value=theta.estimate$lambda,
                                                  calculate.beta=FALSE)$gamma
                   
                   conformal.residual <-
                     as.numeric(GMM.Data.Conformal.Iter$Y - as.matrix(GMM.Data.Conformal.Iter[,pos.W])%*%conformal.gamma)
                   Calculate.PValue(conformal.residual,T0,1,p=1)
                 })
        
        conformal.interval.vec[conformal.iter,] <- range(beta.grid.narrow[pvalue.vec.narrow>=valid.pvalue])
      }
      
      RESULT$conformal.interval <- conformal.interval.vec
      colnames(RESULT$conformal.interval) <- c("LB","UB")
    }
    
  }
  
  RESULT$lambda <- theta.estimate$lambda
  
  RESULT$detrend <- detrend
  RESULT$detrend.ft <- detrend.ft
  if(detrend){
    RESULT$trend <- as.numeric(as.matrix(GMM.Data[,pos.D])%*%theta.estimate$detrend)
  }
  RESULT$Y.basis <- Y.basis
  RESULT$att.ft <- att.ft
  RESULT$bootstrap.num <- bootstrap.num
  RESULT$conformal.period <- conformal.period
  RESULT$conformal.pvalue <- conformal.pvalue
  
  class(RESULT) <- c("list","SPSC")
  
  return(RESULT)
}

SPSC.Effect <- function(GMM.Data,
                        Y.basis,
                        lambda.type,
                        lambda.value,
                        lambda.grid,
                        calculate.beta=TRUE){
  
  T0 <- sum(1-GMM.Data$A)
  T1 <- sum(GMM.Data$A)
  Tt <- T0+T1
  
  pos.D <- which(substr(colnames(GMM.Data),1,1)=="D")
  pos.B <- which(substr(colnames(GMM.Data),1,1)=="B")
  pos.W <- which(substr(colnames(GMM.Data),1,1)=="W")
  pos.Y <- which(substr(colnames(GMM.Data),1,1)=="Y")
  
  Y.Pre <- GMM.Data$Y[GMM.Data$A==0]
  Y.Post <- GMM.Data$Y[GMM.Data$A==1]
  W.Pre <- GMM.Data[GMM.Data$A==0,pos.W]
  W.Post <- GMM.Data[GMM.Data$A==1,pos.W]
  
  detrend <- length(pos.D)>0
  
  ###### Detrend
  
  if(detrend){
    
    detrend.formula <- as.formula(paste(c("Y~0+",
                                          paste(colnames(GMM.Data)[pos.D],collapse="+")),
                                        collapse=""))
    LM <- lm(detrend.formula,data=GMM.Data[GMM.Data$A==0,])
    Y.detrend.Pre <- as.numeric(LM$residuals)
    Y.detrend.estimate <- as.numeric(LM$coefficients)
    Y.detrend.basis.Pre <- make.matrix( t(sapply(1:T0,function(t){Y.basis(Y.detrend.Pre[t])})), T0 )
    gY.Pre <- cbind(GMM.Data[GMM.Data$A==0,pos.D],
                    Y.detrend.basis.Pre)
    gY.Post <- matrix(0,T1,ncol(gY.Pre))
    
  } else {
    gY.Pre <- matrix(Y.Pre,T0,1)
    gY.Post <- matrix(0,T1,ncol(gY.Pre))
  }
  
  ###### GMM detrend Data
  
  GMM.detrend.Data <- cbind(rbind(as.matrix(gY.Pre), as.matrix(gY.Post)),
                            rbind(W.Pre,  W.Post),
                            rep(c(0,1),   c(T0,T1)),
                            c(Y.Pre,      Y.Post))
  GMM.detrend.Data <- as.matrix(GMM.detrend.Data)
  
  m <- ncol(gY.Pre)
  N <- ncol(W.Pre)
  
  colnames(GMM.detrend.Data) <-   c(sprintf("G%0.4d",1:m),
                                    sprintf("W%0.4d",1:N),
                                    "A", "Y")
  
  pos.detrend.W <- which(substr(colnames(GMM.detrend.Data),1,1)=="W")
  pos.detrend.G <- which(substr(colnames(GMM.detrend.Data),1,1)=="G")
  
  GMM.detrend.Data <- data.frame(GMM.detrend.Data)
  for(jj in 1:ncol(GMM.detrend.Data)){
    GMM.detrend.Data[,jj] <- as.numeric( GMM.detrend.Data[,jj] )
  }
  
  ###### select lambda.opt
  
  if(lambda.type=="cv"){
    LOO.error <- sapply( lambda.grid, function(lambda){CV.lambda(GMM.detrend.Data, lambda)} )
    lambda.opt <- lambda.grid[which.min(LOO.error)]
  } else {
    lambda.opt <- lambda.value
  }
  
  ###### gamma
  
  if(length(pos.detrend.G)>1){
    GY <- apply(GMM.detrend.Data$Y[1:T0] * GMM.detrend.Data[1:T0,pos.detrend.G],2,mean)
    GW <- as.matrix(sapply(1:N,function(nn){apply(GMM.detrend.Data[1:T0,pos.detrend.W[nn]]*GMM.detrend.Data[1:T0,pos.detrend.G],2,mean)}))
    gamma.estimate <- as.numeric( MASS::ginv(t(GW)%*%(GW) + (Tt/T0)^2*diag(rep(10^(lambda.opt),N)))%*%(t(GW)%*%GY) )
  } else {
    GY <- mean(GMM.detrend.Data$Y[1:T0] * GMM.detrend.Data[1:T0,pos.detrend.G])
    GW <- as.numeric( apply(GMM.detrend.Data[1:T0,pos.detrend.W] * GMM.detrend.Data[1:T0,pos.detrend.G],2,mean) )
    GW <- matrix(GW,1,length(GW))
    gamma.estimate <- as.numeric( MASS::ginv(t(GW)%*%(GW) + (Tt/T0)^2*diag(rep(10^(lambda.opt),N)))%*%(t(GW)%*%GY) )
  }
  
  ###### beta
  
  if(calculate.beta){
    att.basis.Post <- as.matrix(GMM.Data[GMM.Data$A==1,pos.B])
    Y.effect.Post <- as.numeric(Y.Post - as.matrix(W.Post)%*%gamma.estimate)
    beta.estimate <- as.numeric( lm( Y.effect.Post~0+att.basis.Post )$coefficients )
  }
  
  ###### effect estimate
  
  theta.estimate <- list()
  if(detrend){
    theta.estimate$detrend <- Y.detrend.estimate
  } else {
    theta.estimate$detrend <- NULL
  }
  theta.estimate$gamma <- gamma.estimate
  if(calculate.beta){
    theta.estimate$beta <- beta.estimate
  }
  theta.estimate$lambda <- lambda.opt
  
  return(theta.estimate)
}

HAC.Meat <- function(Psi,
                     beta.pos,
                     T0,T1){
  
  lengthb <- length(beta.pos)
  
  CROSS <- list()
  for(ll in 0:(T0+T1-1)){
    CROSS[[ll+1]] <- matrix(0,dim(Psi)[2],dim(Psi)[2])
    for(tt in 1:(T0+T1-ll)){
      CROSS[[ll+1]] <- CROSS[[ll+1]]+
        (t(t(Psi[tt,]))%*%(t(Psi[tt+ll,]))) +
        +as.numeric(ll>0)*(t(t(Psi[tt+ll,]))%*%(t(Psi[tt,])))
    }
  }
  
  rho.Vec    <- rep(0,lengthb)
  sigma2.Vec <- rep(0,lengthb)
  
  for(beta.iter in 1:lengthb){
    ARIMA <- try( arima(x=Psi[,dim(Psi)[2]-lengthb+beta.iter],
                        order=c(1,0,0)),
                  silent=T)
    if(class(ARIMA)!="try-error"){
      rho.Vec[beta.iter]    <- ARIMA$coef[1]
      sigma2.Vec[beta.iter] <- ARIMA$sigma2
    } else {
      ARIMA <- try( arima(x=Psi[,dim(Psi)[2]-lengthb+beta.iter],
                          order=c(1,0,0),method="ML"),
                    silent=T)
      if(class(ARIMA)!="try-error"){
        rho.Vec[beta.iter]    <- ARIMA$coef[1]
        sigma2.Vec[beta.iter] <- ARIMA$sigma2
      }
    }
    
  }
  
  alpha1 <- sum( 4*rho.Vec^2*sigma2.Vec^2/(1-rho.Vec)^6/(1+rho.Vec)^2 )/
    sum( sigma2.Vec^2/(1-rho.Vec)^4 )
  
  alpha2 <- sum( 4*rho.Vec^2*sigma2.Vec^2/(1-rho.Vec)^8 )/
    sum( sigma2.Vec^2/(1-rho.Vec)^4 )
  
  opt.BW <- 1.3221*(alpha2*(T1+T0))^(0.2)
  
  opt.Bartlett.BW <- 1.1447*(alpha1*(T1+T0))^(1/3)
  
  Meat.Weight.F <- rep(0,T0+T1)
  Meat.Weight.F[1] <- 1
  for(tt in 2:(T0+T1)){
    z <- 6*pi/5*(tt-1)/opt.BW
    Meat.Weight.F[tt] <- 3/(z^2)*(sin(z)/z - cos(z))
  }
  
  Meat.Meat <- matrix(0,dim(Psi)[2],dim(Psi)[2])
  for(ll in 0:(T0+T1-1)){
    Meat.Meat <- Meat.Meat+CROSS[[ll+1]]*Meat.Weight.F[ll+1]
  }
  
  Meat.Meat <- Meat.Meat/(T0+T1)
  Meat.Simple <- Meat.Meat
  
  RESULT <- list()
  RESULT$matrix <- Meat.Simple
  RESULT$bw <- opt.BW
  RESULT$bw.B <- opt.Bartlett.BW
  return( RESULT )
  
}

Psi.Ft <- function(theta,
                   GMM.Data,
                   Y.basis){
  
  pos.D <- which(substr(colnames(GMM.Data),1,1)=="D")
  pos.B <- which(substr(colnames(GMM.Data),1,1)=="B")
  pos.W <- which(substr(colnames(GMM.Data),1,1)=="W")
  pos.Y <- which(substr(colnames(GMM.Data),1,1)=="Y")
  
  T0 <- sum(1-GMM.Data$A)
  T1 <- sum(GMM.Data$A)
  Tt <- T0+T1
  
  detrend <- !is.null(theta$detrend)
  
  if(detrend){
    Psi.detrend <- as.matrix( (1-GMM.Data$A)*GMM.Data[,pos.D]*(GMM.Data$Y - as.matrix(GMM.Data[,pos.D])%*%theta$detrend) )
    Psi.gamma.detrend <- as.matrix( (1-GMM.Data$A)*GMM.Data[,pos.D]*(GMM.Data$Y - as.matrix(GMM.Data[,pos.W])%*%theta$gamma) )
    Y.detrend <- (GMM.Data$Y - as.matrix(GMM.Data[,pos.D])%*%theta$detrend)
    Y.detrend.basis <- make.matrix( t(sapply( 1:Tt, function(jj){Y.basis(GMM.Data$Y[jj])} )), Tt )
    Psi.gamma.SPSC <- as.matrix( (1-GMM.Data$A)*Y.detrend.basis*c(GMM.Data$Y - as.matrix(GMM.Data[,pos.W])%*%theta$gamma) )
  } else {
    Psi.detrend <- NULL
    Psi.gamma.detrend <- NULL
    Y.detrend.basis <- make.matrix( t(sapply( 1:Tt, function(jj){Y.basis(GMM.Data$Y[jj])} )), Tt )
    Psi.gamma.SPSC <- as.matrix( (1-GMM.Data$A)*Y.detrend.basis*c(GMM.Data$Y - as.matrix(GMM.Data[,pos.W])%*%theta$gamma) )
  }
  
  Psi.beta <- as.matrix( GMM.Data$A*GMM.Data[,pos.B]*(GMM.Data$Y - as.matrix(GMM.Data[,pos.W])%*%theta$gamma - as.matrix(GMM.Data[,pos.B])%*%theta$beta) )
  
  Psi.merge <- as.matrix( cbind(Psi.detrend,
                                Psi.gamma.detrend,
                                Psi.gamma.SPSC,
                                Psi.beta) )
  colnames(Psi.merge) <- NULL
  
  return( Psi.merge )
}

Grad.Psi.Aver.Ft <- function(theta,
                             GMM.Data,
                             Y.basis){
  
  eps <- 10^(-6)
  
  Psi.ncol <- ncol( Psi.Ft(theta,GMM.Data,Y.basis) )
  T0 <- sum(1-GMM.Data$A)
  T1 <- sum(GMM.Data$A)
  Tt <- T0+T1
  
  if(!is.null(theta$detrend)){
    
    num.detrend <- length(theta$detrend)
    num.gamma <- length(theta$gamma)
    num.beta <- length(theta$beta)
    
    Grad.Matrix <- matrix(0,Psi.ncol,num.detrend+num.gamma+num.beta)
    
    theta.vec <- c(theta$detrend,
                   theta$gamma,
                   theta$beta)
    
    for(jj in 1:length(theta.vec)){
      
      theta.vec.plus <- theta.vec.minus <- theta.vec
      theta.vec.plus[jj] <- theta.vec.plus[jj] + eps
      theta.vec.minus[jj] <- theta.vec.minus[jj] - eps
      
      theta.plus <- theta.minus <- list()
      
      theta.plus$detrend <- theta.vec.plus[1:num.detrend]
      theta.plus$gamma <- theta.vec.plus[num.detrend+1:num.gamma]
      theta.plus$beta <- theta.vec.plus[num.detrend+num.gamma+1:num.beta]
      
      theta.minus$detrend <- theta.vec.minus[1:num.detrend]
      theta.minus$gamma <- theta.vec.minus[num.detrend+1:num.gamma]
      theta.minus$beta <- theta.vec.minus[num.detrend+num.gamma+1:num.beta]
      
      Grad.Matrix[,jj] <-
        ( apply( Psi.Ft(theta.plus,GMM.Data,Y.basis), 2, mean ) - apply( Psi.Ft(theta.minus,GMM.Data,Y.basis), 2, mean ) ) / (2*eps)
      
    }
    
  } else {
    
    num.gamma <- length(theta$gamma)
    num.beta <- length(theta$beta)
    
    Grad.Matrix <- matrix(0,Psi.ncol,num.gamma+num.beta)
    
    theta.vec <- c(theta$gamma,
                   theta$beta)
    
    for(jj in 1:length(theta.vec)){
      
      theta.vec.plus <- theta.vec.minus <- theta.vec
      theta.vec.plus[jj] <- theta.vec.plus[jj] + eps
      theta.vec.minus[jj] <- theta.vec.minus[jj] - eps
      
      theta.plus <- theta.minus <- list()
      
      theta.plus$gamma <- theta.vec.plus[1:num.gamma]
      theta.plus$beta <- theta.vec.plus[num.gamma+1:num.beta]
      
      theta.minus$gamma <- theta.vec.minus[1:num.gamma]
      theta.minus$beta <- theta.vec.minus[num.gamma+1:num.beta]
      
      Grad.Matrix[,jj] <-
        ( apply( Psi.Ft(theta.plus,GMM.Data,Y.basis), 2, mean ) - apply( Psi.Ft(theta.minus,GMM.Data,Y.basis), 2, mean ) ) / (2*eps)
      
    }
    
  }
  
  return(Grad.Matrix)
  
}


CV.lambda <- function(GMM.detrend.Data,
                      lambda){
  
  posW <- which(substr(colnames(GMM.detrend.Data),1,1)=="W")
  posG <- which(substr(colnames(GMM.detrend.Data),1,1)=="G")
  Y <- GMM.detrend.Data$Y
  T1 <- sum(GMM.detrend.Data$A)
  T0 <- sum(1-GMM.detrend.Data$A)
  Tt <- T1+T0
  CV.Residual <- rep(0,T0)
  for(jj in 1:T0){
    Temp.T <- (1:T0)[-jj]
    X.Temp  <- GMM.detrend.Data[Temp.T,posW]
    Y.Temp  <- Y[Temp.T]
    gY.Temp <- GMM.detrend.Data[Temp.T,posG]
    if(length(posG)==1){
      gY.Temp <- matrix(gY.Temp,length(gY.Temp),1)
    }
    N <- ncol(X.Temp)
    m <- ncol(gY.Temp)
    GY <- c(apply(gY.Temp*Y.Temp,2,mean))
    GX <- as.matrix(sapply(1:N,function(nn){apply(gY.Temp*X.Temp[,nn],2,mean)}))
    if(length(posG)==1){
      GX <- matrix(GX,1,length(GX))
    }
    gamma.loo <- (MASS::ginv( t(GX)%*%(GX) + (Tt/T0)^2*10^(lambda)*diag(rep(1,N)) )%*%(t(GX)%*%(GY)))
    CV.Residual[jj] <-  as.numeric((Y[jj] - (t(as.numeric(GMM.detrend.Data[jj,posW]))%*%gamma.loo)))
  }
  
  return(log(mean(CV.Residual^2)))
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

make.matrix <- function(vv,TT){
  matrix( vv, TT, ncol(vv)*nrow(vv)/TT )
}

Spline.Trend <- function(t,T0,df=5){
  MAT <- splines::bs(1:(T0+1),df=df,intercept=T)
  return(matrix(MAT[min(t,T0+1),],1,ncol(MAT)))
}


plot.SPSC <- function(spsc,                           # SPSC object
                      COL=c(1,2),                     # Y/SC color
                      LTY=c(1,2),                     # Y/SC line type
                      PI=T,                           # draw prediction interval
                      COL.PI=rgb(1,0,0,0.2),          # prediction interval color
                      caption=T,                      # Denote variables
                      ...
){
  Tt <- length(spsc$Y)
  T1 <- length(spsc$ATT)
  T0 <- Tt-T1
  Yvec <- c(spsc$Y, spsc$SC)
  YL <- range(Yvec)+c(-0.1,0.1)*sd(Yvec)
  
  plot(1:Tt,
       spsc$Y,
       type='l',
       col=COL[1],
       lty=LTY[1],
       xlab="time",
       ylab="outcome",
       xlim=c(0,Tt+1),
       ylim=YL) # observed Y
  points(1:Tt,
         spsc$SC,
         type='l',
         col=COL[2],
         lty=LTY[2]) # SC = synthetic control
  
  if(PI & !is.null(spsc$conformal.interval)){
    T1.PI <- nrow(spsc$conformal.interval)
    
    PI.UB <- spsc$Y[T0+spsc$conformal.period] - spsc$conformal.interval[,1]
    PI.LB <- spsc$Y[T0+spsc$conformal.period] - spsc$conformal.interval[,2]
    POLYGON <- cbind(c(T0,T0+spsc$conformal.period,
                       T0+spsc$conformal.period[T1.PI:1],T0),
                     c(spsc$SC[T0],
                       PI.UB[1:T1.PI],
                       PI.LB[T1.PI:1],
                       spsc$SC[T0]))
    polygon(POLYGON,
            col=COL.PI,
            border=NA)
  }
  
  if(caption){
    if(spsc$Y[Tt]>spsc$SC[Tt]){
      text(Tt+1,spsc$Y[Tt],"Y",col=COL[1],pos=3)
      text(Tt+1,spsc$SC[Tt],"SC",col=COL[2],pos=1)
    } else {
      text(Tt+1,spsc$Y[Tt],"Y",col=COL[1],pos=1)
      text(Tt+1,spsc$SC[Tt],"SC",col=COL[2],pos=3)
    }
  }
  
  
  
}
