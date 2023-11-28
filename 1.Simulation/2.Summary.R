TT  <- c(50,100,250,500)
DD  <- c(0,1)
BB  <- c(1,2)
NN  <- c(2,5,9)

CSV <- list()
for(tt in 1:4){
  CSV[[tt]] <- list()
  for(dd in 1:2){
    CSV[[tt]][[dd]] <- list()
    for(bb in 1:2){
      CSV[[tt]][[dd]][[bb]] <- list()
      for(nn in 1:3){
        FILE <- sprintf("Result/Merge_SPSC_T%0.5d_D%0.1d_B%0.1d_N%0.2d.csv",
                        TT[tt],DD[dd],BB[bb],NN[nn])
        FFFF <- read.csv(FILE)
        CSV[[tt]][[dd]][[bb]][[nn]] <- cbind(FFFF,NA)
        colnames( CSV[[tt]][[dd]][[bb]][[nn]] ) <- c(colnames(FFFF),"NA")
      }
    }
  }
}

CN <- colnames( CSV[[1]][[1]][[1]][[1]] )
col.EST <- c(which(CN=="OLSBias"),
             which(CN=="AbadieBias"),
             which(CN=="ASCBias"),
             which(CN=="SCPIBias"),
             which(CN=="SSCBias"),
             which(CN=="SSCRegBias"))

Estimator.Name <- c("OLS-NoReg","OLS-Standard","ASC","SCPI","SPSC-NoReg","SPSC-Ridge")

################################################################################
# Plot Function
################################################################################

PLOT <- function(dd,bb){
  
  BIAS <- function(dd,bb,nn){
    Bias.Raw <- list()
    Bias.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST]  
    
    MM <- max(dim(Bias.Raw[[1]])[1],
              dim(Bias.Raw[[2]])[1],
              dim(Bias.Raw[[3]])[1],
              dim(Bias.Raw[[4]])[1])
    
    Bias.Merge <- matrix(NA,MM,4*length(col.EST))
    for(bb in 1:length(col.EST)){
      Bias.Merge[ 1:dim(Bias.Raw[[1]])[1],4*bb-3] <- Bias.Raw[[1]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[2]])[1],4*bb-2] <- Bias.Raw[[2]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[3]])[1],4*bb-1] <- Bias.Raw[[3]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[4]])[1],4*bb-0] <- Bias.Raw[[4]][,bb]
    }
    
    Bias.Merge
  }
  
  BB1 <- BIAS(dd,bb,1)
  BB2 <- BIAS(dd,bb,2)
  BB3 <- BIAS(dd,bb,3)
  
  MD <- max(dim(BB1)[1], dim(BB2)[1], dim(BB3)[1])
  if(dim(BB1)[1]<MD){
    BB1 <- rbind(BB1,
                 matrix(NA,MD-dim(BB1)[1],dim(BB1)[2]))
  }
  if(dim(BB2)[1]<MD){
    BB2 <- rbind(BB2,
                 matrix(NA,MD-dim(BB2)[1],dim(BB2)[2]))
  }
  if(dim(BB3)[1]<MD){
    BB3 <- rbind(BB3,
                 matrix(NA,MD-dim(BB3)[1],dim(BB3)[2]))
  }
  
  
  
  Bias.Merge <- cbind( BB1, BB2, BB3 )
  
  T.t   <- length(TT)
  Xpos  <- 1:(dim(Bias.Merge)[2]) + 
    rep(c(1:length(col.EST)-1, 
          1:length(col.EST)+length(col.EST), 
          1:length(col.EST)+2*length(col.EST)+1),
        each=T.t)
  
  COL <- COL.Line <- rep(rep( c("gray80","gray80", "gray60","gray60", "black","black"), each=T.t),length(col.EST))
  
  PCH <- rep(c(15,17,18,19)[1:T.t],length(col.EST)*length(NN))
  CEX <- rep(c(1,1,1.25,1)[1:T.t],length(col.EST)*length(NN))
  LTY <- rep(rep(c(2,1,2,1,2,1),each=T.t),length(NN))
  
  Bias <- apply(Bias.Merge,2,function(v){mean(as.numeric(na.omit(v)))})
  UB   <- apply(Bias.Merge,2,function(v){quantile(as.numeric(na.omit(v)),0.975)})
  LB   <- apply(Bias.Merge,2,function(v){quantile(as.numeric(na.omit(v)),0.025)})
  
  # YL <- c(min(floor(LB)),max(ceiling(UB)))
  YL <- c(-1.5,3.5)
  
  MAR <- c(0.5,3.5,0.5,0.5)
  
  layout(matrix(c(1,2),1,2,byrow=T),
         widths=c(5,1),
         heights=c(0.5,5))
  
  par(mar=MAR)
  plot.new()
  plot.window(xlim=range(Xpos)+c(-0.5,0.5),
              ylim=YL)
  
  LINE <- function(jj){
    points(Xpos[jj],Bias[jj],col=COL[jj],pch=PCH[jj],cex=CEX[jj])
    segments(Xpos[jj],LB[jj],Xpos[jj],UB[jj],lty=LTY[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,LB[jj],Xpos[jj]+0.2,LB[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,UB[jj],Xpos[jj]+0.2,UB[jj],col=COL[jj])
    
  }
  abline(h=0,lty=3,col="gray75")
  sapply(1:length(Xpos),LINE)
  abline(v=mean(Xpos[c(length(col.EST)*T.t,
                       length(col.EST)*T.t+1)]),
         col=1,lty=3)
  abline(v=mean(Xpos[c(2*length(col.EST)*T.t,
                       2*length(col.EST)*T.t+1)]),
         col=1,lty=3)
  
  axis(2,at=seq(min(YL),max(YL),by=0.5))
  axis(2,at=0,
       labels="Bias",
       tick=F,
       line=1.25,
       cex.axis=1.25)
  
  text(mean( Xpos[  1:(length(col.EST)*T.t)] ),
       (YL)[2], "d=2",cex=1.25)
  text(mean( Xpos[length(col.EST)*T.t+1:(length(col.EST)*T.t)] ),  
       (YL)[2], "d=5",cex=1.25)
  text(mean( Xpos[2*length(col.EST)*T.t+1:(length(col.EST)*T.t)] ), 
       (YL)[2], "d=9",cex=1.25)
  
  
  ##
  
  par(mar=MAR*c(1,0,1,0))
  plot.new()
  
  NR <- 1+1+1+length(col.EST)+4
  
  YT <- (seq(95,5,length=NR)/100)[1:NR]
  XT <- c(0,0.05,0.25,0.3,
          0,0.15,0.3)
  
  text(XT[1],YT[1],"Estimator",pos=4,cex=1.25)
  for(bb in 1:length(col.EST)){
    segments(XT[2],YT[bb+1],XT[3],YT[bb+1],
             col=COL.Line[4*bb-3],lty=LTY[4*bb-3],lwd=3)
    text(XT[4],YT[bb+1],
         Estimator.Name[bb],
         pos=4)
  }
  
  bb <- length(col.EST)+3
  text(XT[5],YT[bb],expression(T["0"]),pos=4,cex=1.25)
  points(XT[6],YT[bb+1],pch=PCH[1],cex=CEX[1])
  text(XT[7],YT[bb+1],expression(T["0"]*"=50"),pos=4)
  points(XT[6],YT[bb+2],pch=PCH[2],cex=CEX[2])
  text(XT[7],YT[bb+2],expression(T["0"]*"=100"),pos=4)
  points(XT[6],YT[bb+3],pch=PCH[3],cex=CEX[3])
  text(XT[7],YT[bb+3],expression(T["0"]*"=250"),pos=4)
  points(XT[6],YT[bb+4],pch=PCH[4],cex=CEX[4])
  text(XT[7],YT[bb+4],expression(T["0"]*"=500"),pos=4)
  
}

HT <- 3.25 ; WD <- 14

png("Plot_Cov_No_ATT_0.png",height=HT,width=WD,unit="in",res=500)
PLOT(1,1)  # No Cov + Constant Trend
dev.off()
png("Plot_Cov_Yes_ATT_0.png",height=HT,width=WD,unit="in",res=500)
PLOT(2,1)  # Cov + Constant Trend
dev.off()
png("Plot_Cov_No_ATT_1.png",height=HT,width=WD,unit="in",res=500)
PLOT(1,2)  # No Cov + Linear Trend
dev.off()
png("Plot_Cov_Yes_ATT_1.png",height=HT,width=WD,unit="in",res=500)
PLOT(2,2)  # Cov + Linear Trend
dev.off()


################################################################################
# Plot For Paper
################################################################################

PLOT.Paper <- function(dd){
  
  BIAS <- function(dd,bb,nn){
    Bias.Raw <- list()
    Bias.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST]  
    
    MM <- max(dim(Bias.Raw[[1]])[1],
              dim(Bias.Raw[[2]])[1],
              dim(Bias.Raw[[3]])[1],
              dim(Bias.Raw[[4]])[1])
    
    Bias.Merge <- matrix(NA,MM,4*length(col.EST))
    for(bb in 1:length(col.EST)){
      Bias.Merge[ 1:dim(Bias.Raw[[1]])[1],4*bb-3] <- Bias.Raw[[1]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[2]])[1],4*bb-2] <- Bias.Raw[[2]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[3]])[1],4*bb-1] <- Bias.Raw[[3]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[4]])[1],4*bb-0] <- Bias.Raw[[4]][,bb]
    }
    
    Bias.Merge
  }
  
  bb <- 1
  
  BB1 <- BIAS(dd,bb,1)
  BB2 <- BIAS(dd,bb,2)
  BB3 <- BIAS(dd,bb,3)
  
  MD <- max(dim(BB1)[1], dim(BB2)[1], dim(BB3)[1])
  if(dim(BB1)[1]<MD){
    BB1 <- rbind(BB1,
                 matrix(NA,MD-dim(BB1)[1],dim(BB1)[2]))
  }
  if(dim(BB2)[1]<MD){
    BB2 <- rbind(BB2,
                 matrix(NA,MD-dim(BB2)[1],dim(BB2)[2]))
  }
  if(dim(BB3)[1]<MD){
    BB3 <- rbind(BB3,
                 matrix(NA,MD-dim(BB3)[1],dim(BB3)[2]))
  }
  
  
  
  Bias.Merge <- cbind( BB1, BB2, BB3 )
  
  T.t   <- length(TT)
  Xpos  <- 1:(dim(Bias.Merge)[2]) + 
    rep(c(1:length(col.EST)-1, 
          1:length(col.EST)+length(col.EST), 
          1:length(col.EST)+2*length(col.EST)+1),
        each=T.t)
  
  Bias <- apply(Bias.Merge,2,function(v){mean(as.numeric(na.omit(v)))})
  UB   <- apply(Bias.Merge,2,function(v){quantile(as.numeric(na.omit(v)),0.975)})
  LB   <- apply(Bias.Merge,2,function(v){quantile(as.numeric(na.omit(v)),0.025)})
  
  
  PCH <- rep(c(15,17,18,19)[1:T.t],length(col.EST)*length(NN))
  CEX <- rep((c(1,1,1.25,1)[1:T.t])*1.5,length(col.EST)*length(NN))
  LTY <- rep(rep(c(2,1,2,1,2,1),each=T.t),length(NN))
  
  COL <- COL.Line <- rep(rep( c("gray75","gray75", "gray55","gray55", "black","black"), each=T.t),length(col.EST))
  # YL <- c(min(floor(LB)),max(ceiling(UB)))
  YL <- c(-1.5,3.5)
  
  MAR <- c(0.5,0.5,0.5,0.5)
  
  layout(cbind(6,matrix(c(1,7,2,3,4,3,5,3),4,2,byrow=T)),
         widths=c(0.4,5,1.2),
         heights=c(0.5,5,0.5,5))
  
  par(mar=MAR*c(0,1,0,1))
  plot.new()
  text(0.5,0.5,"Constant Effect",cex=1.75)
  
  plot.new()
  plot.window(xlim=range(Xpos)+c(-0.5,0.5),
              ylim=YL)
  
  LINE <- function(jj){
    points(Xpos[jj],Bias[jj],col=COL[jj],pch=PCH[jj],cex=CEX[jj])
    segments(Xpos[jj],LB[jj],Xpos[jj],UB[jj],lty=LTY[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,LB[jj],Xpos[jj]+0.2,LB[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,UB[jj],Xpos[jj]+0.2,UB[jj],col=COL[jj])
    
  }
  abline(h=0,lty=3,col="gray75")
  sapply(1:length(Xpos),LINE)
  abline(v=mean(Xpos[c(length(col.EST)*T.t,
                       length(col.EST)*T.t+1)]),
         col=1,lty=3)
  abline(v=mean(Xpos[c(2*length(col.EST)*T.t,
                       2*length(col.EST)*T.t+1)]),
         col=1,lty=3)
  
  axis(2,at=seq(min(YL),max(YL),by=0.5),
       labels=c("-1.5","","-0.5","","0.5","","1.5","","2.5","","3.5"))
  axis(2,at=0,
       labels="",
       tick=F,
       line=1.25,
       cex.axis=1.25)
  
  text(mean( Xpos[  1:(length(col.EST)*T.t)] ),
       (YL)[2], "d=2",cex=1.5)
  text(mean( Xpos[length(col.EST)*T.t+1:(length(col.EST)*T.t)] ),  
       (YL)[2], "d=5",cex=1.5)
  text(mean( Xpos[2*length(col.EST)*T.t+1:(length(col.EST)*T.t)] ), 
       (YL)[2], "d=9",cex=1.5)
  
  
  ##
  
  par(mar=MAR*c(1,0,1,0))
  plot.new()
  
  NR <- 1+1+1+length(col.EST)+4
  
  YT <- (seq(95,5,length=NR)/100)[1:NR]
  XT <- c(0,0.05,0.45,0.5,
          0,0.25,0.5)
  
  text(XT[1],YT[1],"Estimator",pos=4,cex=1.4)
  for(bb in 1:length(col.EST)){
    segments(XT[2],YT[bb+1],XT[3],YT[bb+1],
             col=COL.Line[4*bb-3],lty=LTY[4*bb-3],lwd=3)
    text(XT[4],YT[bb+1],
         Estimator.Name[bb],
         cex=1.25,
         pos=4)
  }
  
  bb <- length(col.EST)+3
  text(XT[5],YT[bb],expression(T["0"]),pos=4,cex=1.4)
  points(XT[6],YT[bb+1],pch=PCH[1],cex=CEX[1])
  text(XT[7],YT[bb+1],expression(T["0"]*"=50"),pos=4,cex=1.25)
  points(XT[6],YT[bb+2],pch=PCH[2],cex=CEX[2])
  text(XT[7],YT[bb+2],expression(T["0"]*"=100"),pos=4,cex=1.25)
  points(XT[6],YT[bb+3],pch=PCH[3],cex=CEX[3])
  text(XT[7],YT[bb+3],expression(T["0"]*"=250"),pos=4,cex=1.25)
  points(XT[6],YT[bb+4],pch=PCH[4],cex=CEX[4])
  text(XT[7],YT[bb+4],expression(T["0"]*"=500"),pos=4,cex=1.25)
  
  ##
  
  
  bb <- 2
  
  BB1 <- BIAS(dd,bb,1)
  BB2 <- BIAS(dd,bb,2)
  BB3 <- BIAS(dd,bb,3)
  
  MD <- max(dim(BB1)[1], dim(BB2)[1], dim(BB3)[1])
  if(dim(BB1)[1]<MD){
    BB1 <- rbind(BB1,
                 matrix(NA,MD-dim(BB1)[1],dim(BB1)[2]))
  }
  if(dim(BB2)[1]<MD){
    BB2 <- rbind(BB2,
                 matrix(NA,MD-dim(BB2)[1],dim(BB2)[2]))
  }
  if(dim(BB3)[1]<MD){
    BB3 <- rbind(BB3,
                 matrix(NA,MD-dim(BB3)[1],dim(BB3)[2]))
  }
  
  
  
  Bias.Merge <- cbind( BB1, BB2, BB3 )
  
  T.t   <- length(TT)
  Xpos  <- 1:(dim(Bias.Merge)[2]) + 
    rep(c(1:length(col.EST)-1, 
          1:length(col.EST)+length(col.EST), 
          1:length(col.EST)+2*length(col.EST)+1),
        each=T.t)
  
  Bias <- apply(Bias.Merge,2,function(v){mean(as.numeric(na.omit(v)))})
  UB   <- apply(Bias.Merge,2,function(v){quantile(as.numeric(na.omit(v)),0.975)})
  LB   <- apply(Bias.Merge,2,function(v){quantile(as.numeric(na.omit(v)),0.025)})
  
  
  par(mar=MAR*c(0,1,0,1))
  plot.new()
  text(0.5,0.5,"Linear Effect",cex=1.75)
  
  plot.new()
  plot.window(xlim=range(Xpos)+c(-0.5,0.5),
              ylim=YL)
  
  LINE <- function(jj){
    points(Xpos[jj],Bias[jj],col=COL[jj],pch=PCH[jj],cex=CEX[jj])
    segments(Xpos[jj],LB[jj],Xpos[jj],UB[jj],lty=LTY[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,LB[jj],Xpos[jj]+0.2,LB[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,UB[jj],Xpos[jj]+0.2,UB[jj],col=COL[jj])
    
  }
  abline(h=0,lty=3,col="gray75")
  sapply(1:length(Xpos),LINE)
  abline(v=mean(Xpos[c(length(col.EST)*T.t,
                       length(col.EST)*T.t+1)]),
         col=1,lty=3)
  abline(v=mean(Xpos[c(2*length(col.EST)*T.t,
                       2*length(col.EST)*T.t+1)]),
         col=1,lty=3)
  
  axis(2,at=seq(min(YL),max(YL),by=0.5),
       labels=c("-1.5","","-0.5","","0.5","","1.5","","2.5","","3.5"))
  axis(2,at=0,
       labels="",
       tick=F,
       line=1.25,
       cex.axis=1.25)
  
  text(mean( Xpos[  1:(length(col.EST)*T.t)] ),
       (YL)[2], "d=2",cex=1.5)
  text(mean( Xpos[length(col.EST)*T.t+1:(length(col.EST)*T.t)] ),  
       (YL)[2], "d=5",cex=1.5)
  text(mean( Xpos[2*length(col.EST)*T.t+1:(length(col.EST)*T.t)] ), 
       (YL)[2], "d=9",cex=1.5)
  
  
  ##
  
  par(mar=MAR*c(1,0,1,0))
  plot.new()
  
  text(0.5,0.5,"Bias",srt=90,cex=1.75)
  
}

png("Plot_Cov_No.png",height=6,width=12,unit="in",res=500)
PLOT.Paper(1)
dev.off()

################################################################################
# Conformal Inference
################################################################################

col.EST.Cover <- c(which(CN=="ASCCover"),
                   which(CN=="ASCLength"),
                   which(CN=="SCPICover"),
                   which(CN=="SCPILength"),
                   which(CN=="SSCCPCover"),
                   which(CN=="SSCLength"),
                   which(CN=="SSCRegCPCover"),
                   which(CN=="SSCRegLength") )

Cover <- function(dd,bb,nn){
  Cover.Raw <- list()
  Cover.Raw[[1]] <- apply( CSV[[1]][[dd]][[bb]][[nn]][,col.EST.Cover], 2, mean )
  # print( dim(CSV[[1]][[dd]][[bb]][[nn]])[1] )
  Cover.Raw[[2]] <- apply( CSV[[2]][[dd]][[bb]][[nn]][,col.EST.Cover], 2, mean )
  # print( dim(CSV[[2]][[dd]][[bb]][[nn]])[1] )
  Cover.Raw[[3]] <- apply( CSV[[3]][[dd]][[bb]][[nn]][,col.EST.Cover], 2, mean )
  # print( dim(CSV[[3]][[dd]][[bb]][[nn]])[1] )
  Cover.Raw[[4]] <- apply( CSV[[4]][[dd]][[bb]][[nn]][,col.EST.Cover], 2, mean )
  # print( dim(CSV[[4]][[dd]][[bb]][[nn]])[1] )
  
  Cover.Merge <- c(Cover.Raw[[2]],
                   Cover.Raw[[4]]) 
  Cover.Merge
}

TABLE.Conformal <- function(dd,bb,add=0){
  
  
  
  CM <- list()
  for(nn in 1:3){
    CM[[nn]] <- Cover(dd,bb,nn)
  }
  
  RRR <- rbind(CM[[1]],CM[[2]],CM[[3]])[,c(1, 
                                           9, 
                                           3, 
                                           11,
                                           5, 
                                           13,
                                           7, 
                                           15)+add]
  
  RRR
}

TABLE.Conformal(1,1)
TABLE.Conformal(1,2)

TABLE.Conformal.Paper <- function(dd){
  
  T1 <- TABLE.Conformal(dd,1)
  T2 <- TABLE.Conformal(dd,2)
  TT <- rbind(T1,T2)
  
  C1 <- c("\\multicolumn{1}{|c|}{\\multirow{3}{*}{\\begin{tabular}[c]{@{}c@{}}Constant\\\\ Effect\\end{tabular}}}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{\\multirow{3}{*}{\\begin{tabular}[c]{@{}c@{}}Linear\\\\ Effect\\end{tabular}}}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}")
  
  AND <- c(" & ")
  
  Cd <- rep(c("$d=2$","$d=5$","$d=9$"),2)
  
  C2 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,1])
  C3 <- sprintf("%0.3f",TT[,2])
  C4 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,3])
  C5 <- sprintf("%0.3f",TT[,4])
  C6 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,5])
  C7 <- sprintf("%0.3f",TT[,6])
  C8 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,7])
  C9 <- sprintf("%0.3f",TT[,8])
  
  END <- c("\\\\ \\cline{2-10}",
           "\\\\ \\cline{2-10}",
           "\\\\ \\hline",
           "\\\\ \\cline{2-10}",
           "\\\\ \\cline{2-10}",
           "\\\\ \\hline")
  
  DDD <- apply(cbind(C1,AND,Cd,AND,C2,AND,C3,AND,C4,AND,C5,AND,C6,AND,C7,AND,C8,AND,C9,END),
               1,
               function(v){paste(v,collapse = "")})
  
  print(data.frame(DDD),row.names=F)
}

TABLE.Conformal.Paper(1)

TABLE.Conformal.Paper.Appendix <- function(add){
  
  T1 <- TABLE.Conformal(1,1,add)
  T2 <- TABLE.Conformal(2,1,add)
  T3 <- TABLE.Conformal(1,2,add)
  T4 <- TABLE.Conformal(2,2,add)
  TT <- rbind(T1,T2,T3,T4)
  
  C1 <- c("\\multicolumn{1}{|c|}{\\multirow{6}{*}{\\begin{tabular}[c]{@{}c@{}}Constant\\\\ Effect\\end{tabular}}}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{\\multirow{6}{*}{\\begin{tabular}[c]{@{}c@{}}Linear\\\\ Effect\\end{tabular}}}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}")
  
  AND <- c(" & ")
  
  Cdelta <- c("\\multicolumn{1}{c|}{\\multirow{3}{*}{$\\delta=0$}}",
              "\\multicolumn{1}{c|}{}",
              "\\multicolumn{1}{c|}{}",
              "\\multicolumn{1}{c|}{\\multirow{3}{*}{$\\delta=1$}}",
              "\\multicolumn{1}{c|}{}",
              "\\multicolumn{1}{c|}{}",
              "\\multicolumn{1}{c|}{\\multirow{3}{*}{$\\delta=0$}}",
              "\\multicolumn{1}{c|}{}",
              "\\multicolumn{1}{c|}{}",
              "\\multicolumn{1}{c|}{\\multirow{3}{*}{$\\delta=1$}}",
              "\\multicolumn{1}{c|}{}",
              "\\multicolumn{1}{c|}{}")
  
  Cd <- rep(c("$d=2$","$d=5$","$d=9$"),4)
  
  C2 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,1])
  C3 <- sprintf("%0.3f",TT[,2])
  C4 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,3])
  C5 <- sprintf("%0.3f",TT[,4])
  C6 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,5])
  C7 <- sprintf("%0.3f",TT[,6])
  C8 <- sprintf("\\multicolumn{1}{c|}{%0.3f}",TT[,7])
  C9 <- sprintf("%0.3f",TT[,8])
  
  END <- c("\\\\ \\cline{3-11}",
           "\\\\ \\cline{3-11}",
           "\\\\ \\cline{2-11}",
           "\\\\ \\cline{3-11}",
           "\\\\ \\cline{3-11}",
           "\\\\ \\hline",
           "\\\\ \\cline{3-11}",
           "\\\\ \\cline{3-11}",
           "\\\\ \\cline{2-11}",
           "\\\\ \\cline{3-11}",
           "\\\\ \\cline{3-11}",
           "\\\\ \\hline")
  
  DDD <- apply(cbind(C1,AND,Cdelta,AND,Cd,AND,C2,AND,C3,AND,C4,AND,C5,AND,C6,AND,C7,AND,C8,AND,C9,END),
               1,
               function(v){paste(v,collapse = "")})
  
  print(data.frame(DDD),row.names=F)
}

TABLE.Conformal.Paper.Appendix(0)
 



################################################################################
# Table Function
################################################################################

col.EST.ASE <- c(which(CN=="OLSSE"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="SSCSE"),
                 which(CN=="SSCRegSE"))

col.EST.Cover <- c(which(CN=="OLSCover"),
                   which(CN=="NA"),
                   which(CN=="NA"),
                   which(CN=="NA"),
                   which(CN=="SSCCover"),
                   which(CN=="SSCRegCover"))

col.EST.BSE <- c(which(CN=="OLSBSE"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="SSCBSE"),
                 which(CN=="SSCRegBSE"))

col.EST.BCover <- c(which(CN=="OLSBCover"),
                    which(CN=="NA"),
                    which(CN=="NA"),
                    which(CN=="NA"),
                    which(CN=="SSCBCover"),
                    which(CN=="SSCRegBCover"))

ASE <- function(dd,bb,nn){
  ASE.Raw <- list()
  ASE.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  ASE.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  ASE.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  ASE.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  
  ASE.Merge <- cbind(ASE.Raw[[1]][,1],
                     ASE.Raw[[2]][,1],
                     ASE.Raw[[3]][,1],
                     ASE.Raw[[4]][,1])
  for(bb in 2:length(col.EST.ASE)){
    ASE.Merge <- cbind(ASE.Merge,
                       ASE.Raw[[1]][,bb],
                       ASE.Raw[[2]][,bb],
                       ASE.Raw[[3]][,bb],
                       ASE.Raw[[4]][,bb])
  }
  ASE.Merge
}

Cover <- function(dd,bb,nn){
  Cover.Raw <- list()
  Cover.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  Cover.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  Cover.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  Cover.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  
  Cover.Merge <- cbind(Cover.Raw[[1]][,1],
                       Cover.Raw[[2]][,1],
                       Cover.Raw[[3]][,1],
                       Cover.Raw[[4]][,1])
  for(bb in 2:length(col.EST.Cover)){
    Cover.Merge <- cbind(Cover.Merge,
                         Cover.Raw[[1]][,bb],
                         Cover.Raw[[2]][,bb],
                         Cover.Raw[[3]][,bb],
                         Cover.Raw[[4]][,bb])
  }
  Cover.Merge
}

BSE <- function(dd,bb,nn){
  BSE.Raw <- list()
  BSE.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  BSE.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  BSE.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  BSE.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  
  BSE.Merge <- cbind(BSE.Raw[[1]][,1],
                     BSE.Raw[[2]][,1],
                     BSE.Raw[[3]][,1],
                     BSE.Raw[[4]][,1])
  for(bb in 2:length(col.EST.BSE)){
    BSE.Merge <- cbind(BSE.Merge,
                       BSE.Raw[[1]][,bb],
                       BSE.Raw[[2]][,bb],
                       BSE.Raw[[3]][,bb],
                       BSE.Raw[[4]][,bb])
  }
  BSE.Merge
}

BCover <- function(dd,bb,nn){
  Cover.Raw <- list()
  Cover.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  Cover.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  Cover.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  Cover.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  
  Cover.Merge <- cbind(Cover.Raw[[1]][,1],
                       Cover.Raw[[2]][,1],
                       Cover.Raw[[3]][,1],
                       Cover.Raw[[4]][,1])
  for(bb in 2:length(col.EST.Cover)){
    Cover.Merge <- cbind(Cover.Merge,
                         Cover.Raw[[1]][,bb],
                         Cover.Raw[[2]][,bb],
                         Cover.Raw[[3]][,bb],
                         Cover.Raw[[4]][,bb])
  }
  Cover.Merge
}

TABLE.Paper <- function(dd){
  
  bb <- 1
  
  Report.Column <- 48 + seq(1:12)*2
  
  BIAS <- function(dd,bb,nn){
    Bias.Raw <- list()
    Bias.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST]  
    
    MM <- max(dim(Bias.Raw[[1]])[1],
              dim(Bias.Raw[[2]])[1],
              dim(Bias.Raw[[3]])[1],
              dim(Bias.Raw[[4]])[1])
    
    Bias.Merge <- matrix(NA,MM,4*length(col.EST))
    for(bb in 1:length(col.EST)){
      Bias.Merge[ 1:dim(Bias.Raw[[1]])[1],4*bb-3] <- Bias.Raw[[1]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[2]])[1],4*bb-2] <- Bias.Raw[[2]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[3]])[1],4*bb-1] <- Bias.Raw[[3]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[4]])[1],4*bb-0] <- Bias.Raw[[4]][,bb]
    }
    
    Bias.Merge
  }
  
  MMM <- function(BB1,BB2,BB3){
    MD <- max(dim(BB1)[1], dim(BB2)[1], dim(BB3)[1])
    if(dim(BB1)[1]<MD){
      BB1 <- rbind(BB1,
                   matrix(NA,MD-dim(BB1)[1],dim(BB1)[2]))
    }
    if(dim(BB2)[1]<MD){
      BB2 <- rbind(BB2,
                   matrix(NA,MD-dim(BB2)[1],dim(BB2)[2]))
    }
    if(dim(BB3)[1]<MD){
      BB3 <- rbind(BB3,
                   matrix(NA,MD-dim(BB3)[1],dim(BB3)[2]))
    }
    cbind( BB1, BB2, BB3 )
  }
  
  Bias.Merge <- MMM(BIAS(dd,bb,1),
                    BIAS(dd,bb,2),
                    BIAS(dd,bb,3))
  ASE.Merge <- MMM(ASE(dd,bb,1),
                   ASE(dd,bb,2),
                   ASE(dd,bb,3))
  BSE.Merge <- MMM(BSE(dd,bb,1),
                   BSE(dd,bb,2),
                   BSE(dd,bb,3))
  
  Cover.Merge <- MMM(Cover(dd,bb,1),
                     Cover(dd,bb,2),
                     Cover(dd,bb,3))
  BCover.Merge <- MMM(BCover(dd,bb,1),
                      BCover(dd,bb,2),
                      BCover(dd,bb,3))
  
  Bias.Table.1 <- apply(Bias.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  MSE.Table.1  <- apply(Bias.Merge[,Report.Column]^2,2,function(vv){mean(na.omit(vv))})
  ESE.Table.1  <- apply(Bias.Merge[,Report.Column],2,function(vv){sd(na.omit(vv))})
  ASE.Table.1 <- apply(ASE.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  BSE.Table.1 <- apply(BSE.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  Cover.Table.1  <- apply(Cover.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  BCover.Table.1 <- apply(BCover.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  
  Final.Table.1 <- rbind(Bias.Table.1*10,
                         ASE.Table.1*10,
                         ESE.Table.1*10,
                         MSE.Table.1*100,
                         Cover.Table.1)
  for(tt in 1:dim(Final.Table.1)[1]){
    Final.Table.1[tt,is.nan(Final.Table.1[tt,])] <- NA
  }
  
  FT1 <- Final.Table.1
  
  
  
  bb <- 2
  
  Report.Column <- 48 + seq(1:12)*2
  
  BIAS <- function(dd,bb,nn){
    Bias.Raw <- list()
    Bias.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST]  
    
    MM <- max(dim(Bias.Raw[[1]])[1],
              dim(Bias.Raw[[2]])[1],
              dim(Bias.Raw[[3]])[1],
              dim(Bias.Raw[[4]])[1])
    
    Bias.Merge <- matrix(NA,MM,4*length(col.EST))
    for(bb in 1:length(col.EST)){
      Bias.Merge[ 1:dim(Bias.Raw[[1]])[1],4*bb-3] <- Bias.Raw[[1]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[2]])[1],4*bb-2] <- Bias.Raw[[2]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[3]])[1],4*bb-1] <- Bias.Raw[[3]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[4]])[1],4*bb-0] <- Bias.Raw[[4]][,bb]
    }
    
    Bias.Merge
  }
  
  MMM <- function(BB1,BB2,BB3){
    MD <- max(dim(BB1)[1], dim(BB2)[1], dim(BB3)[1])
    if(dim(BB1)[1]<MD){
      BB1 <- rbind(BB1,
                   matrix(NA,MD-dim(BB1)[1],dim(BB1)[2]))
    }
    if(dim(BB2)[1]<MD){
      BB2 <- rbind(BB2,
                   matrix(NA,MD-dim(BB2)[1],dim(BB2)[2]))
    }
    if(dim(BB3)[1]<MD){
      BB3 <- rbind(BB3,
                   matrix(NA,MD-dim(BB3)[1],dim(BB3)[2]))
    }
    cbind( BB1, BB2, BB3 )
  }
  
  Bias.Merge <- MMM(BIAS(dd,bb,1),
                    BIAS(dd,bb,2),
                    BIAS(dd,bb,3))
  ASE.Merge <- MMM(ASE(dd,bb,1),
                   ASE(dd,bb,2),
                   ASE(dd,bb,3))
  BSE.Merge <- MMM(BSE(dd,bb,1),
                   BSE(dd,bb,2),
                   BSE(dd,bb,3))
  
  Cover.Merge <- MMM(Cover(dd,bb,1),
                     Cover(dd,bb,2),
                     Cover(dd,bb,3))
  BCover.Merge <- MMM(BCover(dd,bb,1),
                      BCover(dd,bb,2),
                      BCover(dd,bb,3))
  
  Bias.Table.1 <- apply(Bias.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  MSE.Table.1  <- apply(Bias.Merge[,Report.Column]^2,2,function(vv){mean(na.omit(vv))})
  ESE.Table.1  <- apply(Bias.Merge[,Report.Column],2,function(vv){sd(na.omit(vv))})
  ASE.Table.1 <- apply(ASE.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  BSE.Table.1 <- apply(BSE.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  Cover.Table.1  <- apply(Cover.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  BCover.Table.1 <- apply(BCover.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  
  Final.Table.1 <- rbind(Bias.Table.1*10,
                         ASE.Table.1*10,
                         ESE.Table.1*10,
                         MSE.Table.1*100,
                         Cover.Table.1)
  for(tt in 1:dim(Final.Table.1)[1]){
    Final.Table.1[tt,is.nan(Final.Table.1[tt,])] <- NA
  }
  
  FT2 <- Final.Table.1
  
  FT <- rbind(FT1,FT2)
  
  
  AND <- " & "
  
  C1 <- c("\\multicolumn{1}{|c|}{\\multirow{5}{*}{\\begin{tabular}[c]{@{}c@{}}Constant\\\\ Effect\\end{tabular}}}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}", 
          "\\multicolumn{1}{|c|}{}", 
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{\\multirow{5}{*}{\\begin{tabular}[c]{@{}c@{}}Linear\\\\ Effect\\end{tabular}}}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}",
          "\\multicolumn{1}{|c|}{}")
  
  C2 <- rep(c("Bias {\\tiny $(\\times 10)$}",
              "ASE {\\tiny $(\\times 10)$}",
              "ESE {\\tiny $(\\times 10)$}",
              "MSE {\\tiny $(\\times 100)$}",
              "Coverage "),2)
  
  C3 <- sprintf("\\multicolumn{1}{c|}{%0.2f}",FT[,1])
  C4 <- sprintf("%0.2f",FT[,2])
  C5 <- sprintf("\\multicolumn{1}{c|}{%0.2f}",FT[,3])
  C6 <- sprintf("%0.2f",FT[,4])
  C7 <- sprintf("\\multicolumn{1}{c|}{%0.2f}",FT[,5])
  C8 <- sprintf("%0.2f",FT[,6])
  C9 <- sprintf("\\multicolumn{1}{c|}{%0.2f}",FT[,7])
  C10 <- sprintf("%0.2f",FT[,8])
  C11 <- sprintf("\\multicolumn{1}{c|}{%0.2f}",FT[,9])
  C12 <- sprintf("%0.2f",FT[,10])
  C13 <- sprintf("\\multicolumn{1}{c|}{%0.2f}",FT[,11])
  C14 <- sprintf("%0.2f",FT[,12])
  
  END <- c(" \\\\ \\cline{2-14}",
           " \\\\ \\cline{2-14}",
           " \\\\ \\cline{2-14}",
           " \\\\ \\cline{2-14}",
           " \\\\ \\hline",
           " \\\\ \\cline{2-14}",
           " \\\\ \\cline{2-14}",
           " \\\\ \\cline{2-14}",
           " \\\\ \\cline{2-14}",
           " \\\\ \\hline")
  
  DDD <- apply(cbind(C1,AND,C2,AND,C3,AND,C4,AND,C5,AND,C6,AND,C7,AND,C8,AND,C9,AND,C10,AND,C11,AND,C12,AND,C13,AND,C14,END),
               1,
               function(v){paste(v,collapse = "")})
  
  print(data.frame(DDD),row.names=F)
  
}


TABLE.Paper(1)












################################################################################
# Table Function
################################################################################

col.EST.ASE <- c(which(CN=="OLSSE"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="SSCSE"),
                 which(CN=="SSCRegSE"))

col.EST.Cover <- c(which(CN=="OLSCover"),
                   which(CN=="NA"),
                   which(CN=="NA"),
                   which(CN=="NA"),
                   which(CN=="SSCCover"),
                   which(CN=="SSCRegCover"))

col.EST.BSE <- c(which(CN=="OLSBSE"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="NA"),
                 which(CN=="SSCBSE"),
                 which(CN=="SSCRegBSE"))

col.EST.BCover <- c(which(CN=="OLSBCover"),
                    which(CN=="NA"),
                    which(CN=="NA"),
                    which(CN=="NA"),
                    which(CN=="SSCBCover"),
                    which(CN=="SSCRegBCover"))

ASE <- function(dd,bb,nn){
  ASE.Raw <- list()
  ASE.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  ASE.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  ASE.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  ASE.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.ASE]  
  
  ASE.Merge <- cbind(ASE.Raw[[1]][,1],
                     ASE.Raw[[2]][,1],
                     ASE.Raw[[3]][,1],
                     ASE.Raw[[4]][,1])
  for(bb in 2:length(col.EST.ASE)){
    ASE.Merge <- cbind(ASE.Merge,
                       ASE.Raw[[1]][,bb],
                       ASE.Raw[[2]][,bb],
                       ASE.Raw[[3]][,bb],
                       ASE.Raw[[4]][,bb])
  }
  ASE.Merge
}

Cover <- function(dd,bb,nn){
  Cover.Raw <- list()
  Cover.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  Cover.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  Cover.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  Cover.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.Cover]  
  
  Cover.Merge <- cbind(Cover.Raw[[1]][,1],
                       Cover.Raw[[2]][,1],
                       Cover.Raw[[3]][,1],
                       Cover.Raw[[4]][,1])
  for(bb in 2:length(col.EST.Cover)){
    Cover.Merge <- cbind(Cover.Merge,
                         Cover.Raw[[1]][,bb],
                         Cover.Raw[[2]][,bb],
                         Cover.Raw[[3]][,bb],
                         Cover.Raw[[4]][,bb])
  }
  Cover.Merge
}

BSE <- function(dd,bb,nn){
  BSE.Raw <- list()
  BSE.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  BSE.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  BSE.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  BSE.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.BSE]  
  
  BSE.Merge <- cbind(BSE.Raw[[1]][,1],
                     BSE.Raw[[2]][,1],
                     BSE.Raw[[3]][,1],
                     BSE.Raw[[4]][,1])
  for(bb in 2:length(col.EST.BSE)){
    BSE.Merge <- cbind(BSE.Merge,
                       BSE.Raw[[1]][,bb],
                       BSE.Raw[[2]][,bb],
                       BSE.Raw[[3]][,bb],
                       BSE.Raw[[4]][,bb])
  }
  BSE.Merge
}

BCover <- function(dd,bb,nn){
  Cover.Raw <- list()
  Cover.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  Cover.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  Cover.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  Cover.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST.BCover]  
  
  Cover.Merge <- cbind(Cover.Raw[[1]][,1],
                       Cover.Raw[[2]][,1],
                       Cover.Raw[[3]][,1],
                       Cover.Raw[[4]][,1])
  for(bb in 2:length(col.EST.Cover)){
    Cover.Merge <- cbind(Cover.Merge,
                         Cover.Raw[[1]][,bb],
                         Cover.Raw[[2]][,bb],
                         Cover.Raw[[3]][,bb],
                         Cover.Raw[[4]][,bb])
  }
  Cover.Merge
}

TABLE.No.Label <- function(dd,bb,Wind){
  
  Report.Column <- 24*(Wind-1) + seq(1:12)*2
  
  BIAS <- function(dd,bb,nn){
    Bias.Raw <- list()
    Bias.Raw[[1]] <- CSV[[1]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[2]] <- CSV[[2]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[3]] <- CSV[[3]][[dd]][[bb]][[nn]][,col.EST]  
    Bias.Raw[[4]] <- CSV[[4]][[dd]][[bb]][[nn]][,col.EST]  
    
    MM <- max(dim(Bias.Raw[[1]])[1],
              dim(Bias.Raw[[2]])[1],
              dim(Bias.Raw[[3]])[1],
              dim(Bias.Raw[[4]])[1])
    
    Bias.Merge <- matrix(NA,MM,4*length(col.EST))
    for(bb in 1:length(col.EST)){
      Bias.Merge[ 1:dim(Bias.Raw[[1]])[1],4*bb-3] <- Bias.Raw[[1]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[2]])[1],4*bb-2] <- Bias.Raw[[2]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[3]])[1],4*bb-1] <- Bias.Raw[[3]][,bb]
      Bias.Merge[ 1:dim(Bias.Raw[[4]])[1],4*bb-0] <- Bias.Raw[[4]][,bb]
    }
    
    Bias.Merge
  }
  
  MMM <- function(BB1,BB2,BB3){
    MD <- max(dim(BB1)[1], dim(BB2)[1], dim(BB3)[1])
    if(dim(BB1)[1]<MD){
      BB1 <- rbind(BB1,
                   matrix(NA,MD-dim(BB1)[1],dim(BB1)[2]))
    }
    if(dim(BB2)[1]<MD){
      BB2 <- rbind(BB2,
                   matrix(NA,MD-dim(BB2)[1],dim(BB2)[2]))
    }
    if(dim(BB3)[1]<MD){
      BB3 <- rbind(BB3,
                   matrix(NA,MD-dim(BB3)[1],dim(BB3)[2]))
    }
    cbind( BB1, BB2, BB3 )
  }
  
  Bias.Merge <- MMM(BIAS(dd,bb,1),
                    BIAS(dd,bb,2),
                    BIAS(dd,bb,3))
  ASE.Merge <- MMM(ASE(dd,bb,1),
                   ASE(dd,bb,2),
                   ASE(dd,bb,3))
  BSE.Merge <- MMM(BSE(dd,bb,1),
                   BSE(dd,bb,2),
                   BSE(dd,bb,3))
  
  Cover.Merge <- MMM(Cover(dd,bb,1),
                     Cover(dd,bb,2),
                     Cover(dd,bb,3))
  BCover.Merge <- MMM(BCover(dd,bb,1),
                      BCover(dd,bb,2),
                      BCover(dd,bb,3))
  
  Bias.Table.1   <- apply(Bias.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  MSE.Table.1    <- apply(Bias.Merge[,Report.Column]^2,2,function(vv){mean(na.omit(vv))})
  ESE.Table.1    <- apply(Bias.Merge[,Report.Column],2,function(vv){sd(na.omit(vv))})
  ASE.Table.1    <- apply(ASE.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  BSE.Table.1    <- apply(BSE.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  Cover.Table.1  <- apply(Cover.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  BCover.Table.1 <- apply(BCover.Merge[,Report.Column],2,function(vv){mean(na.omit(vv))})
  
  Final.Table.1 <- rbind(Bias.Table.1*10,
                         ASE.Table.1*10,
                         BSE.Table.1*10,
                         ESE.Table.1*10,
                         MSE.Table.1*100,
                         Cover.Table.1,
                         BCover.Table.1)
  for(tt in 1:dim(Final.Table.1)[1]){
    Final.Table.1[tt,is.nan(Final.Table.1[tt,])] <- NA
  }
  
  Final.Table.1
  
}

PRINT.WHOLE <- function(dd,bb){
  AND <- rep(" & ",21)
  T1 <- rbind( TABLE.No.Label(dd,bb,1) , 
               TABLE.No.Label(dd,bb,2) , 
               TABLE.No.Label(dd,bb,3) )
  PT1 <- function(T1){
    cT1 <- dim(T1)[2]
    RR  <- cbind(AND, sprintf("%0.3f",T1[,1]))
    for(tt in 2:cT1){
      RR <- cbind(RR, AND, sprintf("%0.3f",T1[,tt]))
    }
    cT1 <- dim(RR)[2]
    for(tt in 1:(cT1/4)){
      RR[,4*tt-2] <- 
        sprintf("\\multicolumn{1}{c|}{%s}",RR[,4*tt-2])
    }
    RR
  }
  
  END <- rep(c(rep( "\\\\ \\cline{3-15}",6),
               "\\\\ \\cline{2-15}"),3)
  END[21] <- "\\\\ \\hline"
  
  Col3 <- rep(c( "Bias $(\\times 10)$",
                 "ASE $(\\times 10)$",
                 "BSE $(\\times 10)$",
                 "ESE $(\\times 10)$",
                 "MSE $(\\times 100)$",
                 "Coverage (ASE)",
                 "Coverage (BSE)"),3)
  
  Col2 <- c("\\multicolumn{1}{c|}{\\multirow{7}{*}{2}}",
            rep("\\multicolumn{1}{c|}{}",6),
            "\\multicolumn{1}{c|}{\\multirow{7}{*}{5}}",
            rep("\\multicolumn{1}{c|}{}",6),
            "\\multicolumn{1}{c|}{\\multirow{7}{*}{9}}",
            rep("\\multicolumn{1}{c|}{}",6))
  
  Col1 <- c(sprintf("\\multicolumn{1}{|c|}{\\multirow{21}{*}{%0.1d}}",dd-1),
            rep("\\multicolumn{1}{|c|}{}",20))
  
  TT1 <- data.frame(cbind(Col1,AND,Col2,AND,Col3,PT1(T1),END))
  apply(TT1,1,function(v){paste(v,collapse="")})
  # TT1
}

print(data.frame( PRINT.WHOLE(1,1) ),row.names=F)
print(data.frame( PRINT.WHOLE(2,1) ),row.names=F)

print(data.frame( PRINT.WHOLE(1,2) ),row.names=F)
print(data.frame( PRINT.WHOLE(2,2) ),row.names=F)
