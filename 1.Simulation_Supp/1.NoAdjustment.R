BATCH <- 1

source("0.Function_SSC_Simulation2.R")
library(scpi)
library(MASS)
library(splines)

RRR <- list()

for(rho.index in 1:3){
  
  rho <- c(0,0.5,1)[rho.index]
  RRR[[rho.index]] <- list()
  
  for(constant.c.index in 1:5){
    
    constant.c <- c(-1,-0.5,0,0.5,1)[constant.c.index]
    
    set.seed(BATCH)
    
    M  <- 1
    J  <- N <- 10
    K  <- 1
    T0 <- T.Pre <- 100
    T1 <- T.Post <- 1
    Tt <- 101
    
    bjt <- matrix(0,Tt,N)
    
    for(jj in 1:N){
      for(tt in 1:Tt){
        
        bjt[tt,jj] <- ifelse(tt>1, 
                             rho*bjt[tt-1,jj]+rnorm(1),
                             rnorm(1))
      }
      
    }
    
    bjt[Tt,1] <- bjt[Tt,1] + constant.c*sd(bjt[1:T0,1])
    weight <- c(3,4,3,rep(0,N-3))/10
    at <- as.numeric( bjt%*%weight ) + rnorm(Tt)*0.5
    
    ##########################################
    # Our Terminology
    ##########################################
    
    lengthb <- 1
    delta   <- 0
    m <- 2*N
     
    W.series  <- bjt
    Wmat.Pre  <- W.series[1:T0,]
    Wmat.Post <- matrix(W.series[T0+1:T1,],1,N)
    
    Y1.series <- at
    Y1.Pre    <- Y1.series[1:T0]
    Y1.Post   <- Y1.series[T0+1:T1]
    
    ## Gy matrix
    
    gY.Pre <- matrix(0,length(Y1.Pre),m) # each row = collection of g(Y)s at time t , each col = time series of each g(Y)
    gY.Pre[,1:(m)] <- bs(Y1.Pre,df=(m),Boundary.knots = range(Y1.Pre)+c(-2,2))
    
    
    GMM.Data <- cbind(rbind(gY.Pre,matrix(0,T.Post,m)),
                      rbind(Wmat.Pre,Wmat.Post),
                      c(rep(0,T.Pre),rep(1,T.Post)),
                      c(Y1.Pre,Y1.Post))
    
    colnames(GMM.Data) <- (c(sprintf("G%0.4d",1:m),
                             sprintf("W%0.4d",1:N),
                             "A",
                             "Y"))
    
    
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
    
    GW <- t( sapply(1:dim(gY.Pre)[2],function(tt){
      apply( t(Wmat.Pre)*matrix(gY.Pre[,tt],N,T.Pre,byrow=T), 1, mean)
    }) )
    
    GY <- ( sapply(1:dim(gY.Pre)[2],function(tt){
      mean(gY.Pre[,tt]*Y1.Pre)
    }) )
    
    gamma.naive <- c( ginv(t(GW)%*%GW)%*%(t(GW)%*%GY) )
    beta.naive  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%gamma.naive)~1)$coefficients )
    
    ##########################
    
    GMM.Simple.Coef <- c(beta.naive,gamma.naive)
    
    CP.CI <- Conformal.Prediction(Wmat.Pre,
                                  matrix(Wmat.Post[1,],1,N),
                                  Y1.Pre,
                                  Y1.Post[1],
                                  Xmat.Pre=NULL,
                                  Xmat.Post=NULL,
                                  cov.ind=0,
                                  center=beta.naive,
                                  bw=3*sd(Y1.Pre))
    
    ##########################
    
    lambda.grid <- seq(-10,0,by=0.1)
    lambda.min  <- lambda.grid[ which.min(sapply(lambda.grid,CV.Lambda)) ]
    
    lambda.opt <- optimize(f=CV.Lambda,
                           lower= lambda.min - 2,
                           upper= lambda.min + 2)$minimum
    
    
    gamma.naive.lambda <- ginv(t(GW)%*%(GW) + ((T.Pre+T.Post)/(T.Pre))^2*diag(rep(10^(lambda.opt),dim(GW)[2])) )%*%(t(GW)%*%(GY))
    beta.naive.lambda  <- as.numeric( lm(as.numeric(Y1.Post - (Wmat.Post)%*%gamma.naive.lambda)~1)$coefficients )
    GMM.Regular.Coef <- c(beta.naive.lambda,gamma.naive.lambda)
    
    CP.Regular.CI <- Conformal.Prediction(Wmat.Pre,
                                          matrix(Wmat.Post[1,],1,N),
                                          Y1.Pre,
                                          Y1.Post[1],
                                          Xmat.Pre=NULL,
                                          Xmat.Post=NULL,
                                          cov.ind=0,
                                          center=beta.naive.lambda,
                                          bw=3*sd(Y1.Pre),
                                          lambda=lambda.opt,
                                          constant=((T.Pre+T.Post)/(T.Pre))^2)
    
    
    ############################
    
    OLS.Data <- cbind(rbind(Wmat.Pre,Wmat.Post),
                      c(rep(0,T.Pre),rep(1,T.Post)),
                      c(Y1.Pre,Y1.Post))
    
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
                  constant=T)
    
    SCPI.Est <- scpi(SCD)
    
    CENTER <- Y1.Post
    
    RRR[[rho.index]][[constant.c.index]] <- 
      c(as.numeric(median(c(CP.CI,0))==0),
        as.numeric(median(c(CP.Regular.CI,0))==0),
        as.numeric(median(c(SCPI.Est$inference.results$CI.all.gaussian[1:2]-CENTER,0))==0),
        diff(c(CP.CI)),
        diff(c(CP.Regular.CI)),
        SCPI.Est$inference.results$CI.all.gaussian[3],
        beta.naive,
        beta.naive.lambda,
        as.numeric( Y1.Post - SCPI.Est$est.results$Y.post.fit )
      )
    
    
  }
  
}


RESULT <- c( RRR[[1]][[1]],
             RRR[[1]][[2]],
             RRR[[1]][[3]],
             RRR[[1]][[4]],
             RRR[[1]][[5]],
             
             RRR[[2]][[1]],
             RRR[[2]][[2]],
             RRR[[2]][[3]],
             RRR[[2]][[4]],
             RRR[[2]][[5]],
             
             RRR[[3]][[1]],
             RRR[[3]][[2]],
             RRR[[3]][[3]],
             RRR[[3]][[4]],
             RRR[[3]][[5]] )

RESULT <- matrix(RESULT,
                 1,length(RESULT))

TYPE <- c("SPSC","SPSCReg","SCPI")

colnames(RESULT) <- 
  c( sprintf("Rho00_Constant1_Cover%s",TYPE), sprintf("Rho00_Constant1_Length%s",TYPE),  sprintf("Rho00_Constant1_Bias%s",TYPE),
     sprintf("Rho00_Constant2_Cover%s",TYPE), sprintf("Rho00_Constant2_Length%s",TYPE),  sprintf("Rho00_Constant2_Bias%s",TYPE),
     sprintf("Rho00_Constant3_Cover%s",TYPE), sprintf("Rho00_Constant3_Length%s",TYPE),  sprintf("Rho00_Constant3_Bias%s",TYPE),
     sprintf("Rho00_Constant4_Cover%s",TYPE), sprintf("Rho00_Constant4_Length%s",TYPE),  sprintf("Rho00_Constant4_Bias%s",TYPE),
     sprintf("Rho00_Constant5_Cover%s",TYPE), sprintf("Rho00_Constant5_Length%s",TYPE),  sprintf("Rho00_Constant5_Bias%s",TYPE),
     
     sprintf("Rho05_Constant1_Cover%s",TYPE), sprintf("Rho05_Constant1_Length%s",TYPE),  sprintf("Rho05_Constant1_Bias%s",TYPE),
     sprintf("Rho05_Constant2_Cover%s",TYPE), sprintf("Rho05_Constant2_Length%s",TYPE),  sprintf("Rho05_Constant2_Bias%s",TYPE),
     sprintf("Rho05_Constant3_Cover%s",TYPE), sprintf("Rho05_Constant3_Length%s",TYPE),  sprintf("Rho05_Constant3_Bias%s",TYPE),
     sprintf("Rho05_Constant4_Cover%s",TYPE), sprintf("Rho05_Constant4_Length%s",TYPE),  sprintf("Rho05_Constant4_Bias%s",TYPE),
     sprintf("Rho05_Constant5_Cover%s",TYPE), sprintf("Rho05_Constant5_Length%s",TYPE),  sprintf("Rho05_Constant5_Bias%s",TYPE),
     
     sprintf("Rho10_Constant1_Cover%s",TYPE), sprintf("Rho10_Constant1_Length%s",TYPE),  sprintf("Rho10_Constant1_Bias%s",TYPE),
     sprintf("Rho10_Constant2_Cover%s",TYPE), sprintf("Rho10_Constant2_Length%s",TYPE),  sprintf("Rho10_Constant2_Bias%s",TYPE),
     sprintf("Rho10_Constant3_Cover%s",TYPE), sprintf("Rho10_Constant3_Length%s",TYPE),  sprintf("Rho10_Constant3_Bias%s",TYPE),
     sprintf("Rho10_Constant4_Cover%s",TYPE), sprintf("Rho10_Constant4_Length%s",TYPE),  sprintf("Rho10_Constant4_Bias%s",TYPE),
     sprintf("Rho10_Constant5_Cover%s",TYPE), sprintf("Rho10_Constant5_Length%s",TYPE),  sprintf("Rho10_Constant5_Bias%s",TYPE) )


write.csv(RESULT,
          sprintf("Result_Raw/Result_BATCH%0.4d.csv",BATCH),
          row.names=F)


