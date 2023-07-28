# Merge Result

# DD  <- c(0,1)
# BB  <- c(1,2)
# NN  <- c(2,5,9)
# TT  <- c(50,100,250,1000)
# 
# for(dd in DD){
#   for(bb in BB){
#     for(nn in NN){
#         Folder <- sprintf("Result_D%0.1d_B%0.1d_N%0.2d",dd,bb,nn)
#         LF     <- list.files(sprintf("%s/%s",getwd(),Folder))
#         for(tt in TT){
#           LFT <- LF[grepl(sprintf("T%0.5d",tt), LF)]
#           RRR <- read.csv(sprintf("%s/%s/%s",getwd(),Folder,LFT[1]))
#           if(length(LFT)>1){
#             for(lft in 2:length(LFT)){
#               RRR <- rbind(RRR,read.csv(sprintf("%s/%s/%s",getwd(),Folder,LFT[lft])))
#             }
#           }
#           write.csv(RRR,
#                     sprintf("Result/Result_D%0.1d_B%0.1d_N%0.2d_T%0.5d.csv",dd,bb,nn,tt),
#                     row.names=F)
#         }
#     }
#   }
# }

CSV <- list()
for(dd in 1:2){
  CSV[[dd]] <- list()
  for(bb in 1:2){
    CSV[[dd]][[bb]] <- list()
    for(nn in 1:3){
      CSV[[dd]][[bb]][[nn]] <- list()
      for(tt in 1:4){
        FILE <- sprintf("Result/Result_D%0.1d_B%0.1d_N%0.2d_T%0.5d.csv",DD[dd],BB[bb],NN[nn],TT[tt])
        CSV[[dd]][[bb]][[nn]][[tt]] <- read.csv(FILE)
      }
    }
  }
}

col.est   <- 
  col.ase   <- 
  col.bse   <- 
  col.asecv <- 
  col.bsecv <- 
  col.cfl   <- 
  col.cfcv  <- 
  col.reg   <- 
  col.scpi  <- list()

col.est   [[1]] <- c(1,6,13)
col.ase   [[1]] <- c(2,7,14)
col.bse   [[1]] <- c(4,9,16)
col.asecv [[1]] <- c(3,8,15)
col.bsecv [[1]] <- c(5,10,17)
col.cfl   [[1]] <- c(11,18,24)
col.cfcv  [[1]] <- c(12,19,23)
col.reg   [[1]] <- c(20) 
col.scpi  [[1]] <- c(21,22,23,24)

col.est   [[2]] <- c(1,6,13)+20
col.ase   [[2]] <- c(2,7,14)+20
col.bse   [[2]] <- c(4,9,16)+20
col.asecv [[2]] <- c(3,8,15)+20
col.bsecv [[2]] <- c(5,10,17)+20
col.cfl   [[2]] <- c(11,18,24)+20
col.cfcv  [[2]] <- c(12,19,23)+20
col.reg   [[2]] <- c(20)+20
col.scpi  [[2]] <- c(21,22,23,24)+20

MERGE <- function(dd,bb,VEC,bloc=1){
  Result.N2 <- cbind(
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]][1]], 
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]][2]], 
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]][3]],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]][3]],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]][3]],
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]][3]]) 
  
  Result.N5 <- cbind(
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]][1]], 
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]][2]], 
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]][3]],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]][3]],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]][3]],
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]][3]]) 
  
  Result.N9 <- cbind(
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]][1]],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]][1]], 
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]][2]],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]][2]], 
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]][3]],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]][3]],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]][3]],
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]][3]]) 
  
  
  Result <- cbind(Result.N2,
                  Result.N5,
                  Result.N9)
  return(Result)
}

Result <- MERGE(1,1,col.est)-3
Result.ASE <- MERGE(1,1,col.ase)
Result.BSE <- MERGE(1,1,col.bse)
Result.ASECV <- MERGE(1,1,col.asecv)
Result.BSECV <- MERGE(1,1,col.bsecv)

Result.X <- MERGE(2,1,col.est)-3
Result.X.ASE <- MERGE(2,1,col.ase)
Result.X.BSE <- MERGE(2,1,col.bse)
Result.X.ASECV <- MERGE(2,1,col.asecv)
Result.X.BSECV <- MERGE(2,1,col.bsecv)

################################
# Figure 1 of the main paper
# Constant ATT, no covariate
################################

T.t <- dim(Result)[2]/3/3
Xpos  <- 1:dim(Result)[2] + rep(c(0:2, 4:6, 8:10),each=T.t)


COL <- rep(rep( c(rgb(0,0,0,0.2),
                  rgb(0,0,0,0.5),
                  rgb(0,0,0,1)),
                each=T.t),3)

PCH <- rep(c(15,17,18,19)[1:T.t],9)
CEX <- rep(c(1,1,1.5,1)[1:T.t],9)
LTY <- rep(1,dim(Result)[2])

Bias <- apply(Result,2,mean)
UB   <- apply(Result,2,function(v){quantile(v,0.975)})
LB   <- apply(Result,2,function(v){quantile(v,0.025)})

YL <- c(floor(min(LB)*4)/4,ceiling(max(UB)*4)/4)

MAR <- c(0.5,3.5,0.5,0.5)

layout(matrix(c(1,2),1,2,byrow=T),
       widths=c(5,1),
       heights=c(0.5,5))

par(mar=MAR)
plot.new()
plot.window(xlim=range(Xpos)+c(-0.5,0.5),ylim=YL)

LINE <- function(jj){
  points(Xpos[jj],Bias[jj],col=COL[jj],pch=PCH[jj],cex=CEX[jj])
  segments(Xpos[jj],LB[jj],Xpos[jj],UB[jj],lty=LTY[jj],col=COL[jj])
  segments(Xpos[jj]-0.2,LB[jj],Xpos[jj]+0.2,LB[jj],col=COL[jj])
  segments(Xpos[jj]-0.2,UB[jj],Xpos[jj]+0.2,UB[jj],col=COL[jj])

}
abline(h=0,col=1,lty=2)
sapply(1:dim(Result)[2],LINE)
abline(v=mean(Xpos[c(3*T.t,3*T.t+1)]),col=1,lty=3)
abline(v=mean(Xpos[3*T.t+c(3*T.t,3*T.t+1)]),col=1,lty=3)
axis(2,at=seq(-1,1,by=0.5))
axis(2,at=0,labels="Bias",tick=F,line=1.25,cex.axis=1.25)

text(mean( Xpos[  1:(3*T.t)] ),     1.25, "d=2",cex=1.5)
text(mean( Xpos[3*T.t+1:(3*T.t)] ), 1.25, "d=5",cex=1.5)
text(mean( Xpos[6*T.t+1:(3*T.t)] ), 1.25, "d=9",cex=1.5)


##

par(mar=MAR*c(1,0,1,0))
plot.new()

YT <- (seq(95,5,by=-10)/100)[1:5]
XT <- c(0,0.05,0.25,0.3,
        0.6,0.65,0.7)

text(XT[1],YT[1],"Estimator",pos=4,cex=1.25)
segments(XT[2],YT[2],XT[3],YT[2],col=COL[1],lwd=1.5)
text(XT[4],YT[2],"OLS",pos=4)
segments(XT[2],YT[3],XT[3],YT[3],col=COL[5],lwd=1.5)
text(XT[4],YT[3],"SPSC",pos=4)
segments(XT[2],YT[4],XT[3],YT[4],col=COL[9],lwd=1.5)
text(XT[4],YT[4],"SPSC-Ridge",pos=4)

##

YT <- (seq(95,5,by=-10)/100)[6:10]
XT <- c(0,0.05,0.25,0.3,
        0,0.15,0.3)

text(XT[5],YT[1],expression(T["0"]),pos=4,cex=1.25)
points(XT[6],YT[2],pch=PCH[1],cex=CEX[1])
text(XT[7],YT[2],expression(T["0"]*"=50"),pos=4)
points(XT[6],YT[3],pch=PCH[2],cex=CEX[2])
text(XT[7],YT[3],expression(T["0"]*"=100"),pos=4)
points(XT[6],YT[4],pch=PCH[3],cex=CEX[3])
text(XT[7],YT[4],expression(T["0"]*"=250"),pos=4)
points(XT[6],YT[5],pch=PCH[4],cex=CEX[4])
text(XT[7],YT[5],expression(T["0"]*"=1000"),pos=4)

################################
# Figure S1 of the supplement
# Constant ATT
################################

## both delta=0, delta=1

T.t <- dim(Result)[2]/3/3
Xpos  <- 1:dim(Result)[2] + rep(c(0:2, 4:6, 8:10),each=T.t)


COL <- rep(rep( c(rgb(0,0,0,0.2),
                  rgb(0,0,0,0.5),
                  rgb(0,0,0,1)),
                each=T.t),3)

PCH <- rep(c(15,17,18,19)[1:T.t],9)
CEX <- rep(c(1,1,1.5,1)[1:T.t],9)
LTY <- rep(1,dim(Result)[2])

Bias <- apply(Result,2,mean)
UB   <- apply(Result,2,function(v){quantile(v,0.975)})
LB   <- apply(Result,2,function(v){quantile(v,0.025)})

YL <- c(floor(min(LB)*4)/4,ceiling(max(UB)*4)/4)

MAR <- c(0.5,3,0.5,0.5)

layout(matrix(c(8,5,9,
                6,1,3,
                7,2,4),3,3,byrow=T),
       widths=c(0.5,5,1),
       heights=c(0.5,5,5))

par(mar=MAR)
plot.new()
plot.window(xlim=range(Xpos)+c(-0.5,0.5),ylim=YL)

LINE <- function(jj){
  points(Xpos[jj],Bias[jj],col=COL[jj],pch=PCH[jj],cex=CEX[jj])
  segments(Xpos[jj],LB[jj],Xpos[jj],UB[jj],lty=LTY[jj],col=COL[jj])
  segments(Xpos[jj]-0.2,LB[jj],Xpos[jj]+0.2,LB[jj],col=COL[jj])
  segments(Xpos[jj]-0.2,UB[jj],Xpos[jj]+0.2,UB[jj],col=COL[jj])

}
abline(h=0,col=2,lty=2)
sapply(1:dim(Result)[2],LINE)
abline(v=mean(Xpos[c(3*T.t,3*T.t+1)]),col=1,lty=3)
abline(v=mean(Xpos[3*T.t+c(3*T.t,3*T.t+1)]),col=1,lty=3)
axis(2,at=seq(-1,1,by=0.5))
axis(2,at=0,labels="Bias",tick=F,line=1,cex.axis=1.25)

##


Bias <- apply(Result.X,2,mean)
UB   <- apply(Result.X,2,function(v){quantile(v,0.975)})
LB   <- apply(Result.X,2,function(v){quantile(v,0.025)})
YL <- c(floor(min(LB)*4)/4,ceiling(max(UB)*4)/4)

par(mar=MAR)
plot.new()
plot.window(xlim=range(Xpos)+c(-0.5,0.5),ylim=YL)

abline(h=0,col=2,lty=2)
sapply(1:dim(Result.X)[2],LINE)
abline(v=mean(Xpos[c(3*T.t,3*T.t+1)]),col=1,lty=3)
abline(v=mean(Xpos[3*T.t+c(3*T.t,3*T.t+1)]),col=1,lty=3)

axis(2,at=seq(-1,1,by=0.5))
axis(2,at=0,labels="Bias",tick=F,line=1,cex.axis=1.25)

##

par(mar=MAR*c(1,0,1,0))
plot.new()

YT <- c(9,7,5,3,1)/10
XT <- c(0,0.05,0.25,0.3,
        0.6,0.65,0.7)

text(XT[1],YT[1],"Estimator",pos=4,cex=1.5)
segments(XT[2],YT[2],XT[3],YT[2],col=COL[1],lwd=1.5)
text(XT[4],YT[2],"OLS",pos=4)
segments(XT[2],YT[3],XT[3],YT[3],col=COL[5],lwd=1.5)
text(XT[4],YT[3],"SPSC",pos=4)
segments(XT[2],YT[4],XT[3],YT[4],col=COL[9],lwd=1.5)
text(XT[4],YT[4],"SPSC-Reg",pos=4)

##

par(mar=MAR*c(1,0,1,0))
plot.new()

YT <- c(9,7,5,3,1)/10
XT <- c(0,0.05,0.25,0.3,
        0,0.15,0.3)

text(XT[5],YT[1],expression(T["0"]),pos=4,cex=1.5)
points(XT[6],YT[2],pch=PCH[1],cex=CEX[1])
text(XT[7],YT[2],expression(T["0"]*"=50"),pos=4)
points(XT[6],YT[3],pch=PCH[2],cex=CEX[2])
text(XT[7],YT[3],expression(T["0"]*"=100"),pos=4)
points(XT[6],YT[4],pch=PCH[3],cex=CEX[3])
text(XT[7],YT[4],expression(T["0"]*"=250"),pos=4)
points(XT[6],YT[5],pch=PCH[4],cex=CEX[4])
text(XT[7],YT[5],expression(T["0"]*"=1000"),pos=4)


##

par(mar=MAR*c(0,1,0,1))
plot.new()
plot.window(xlim=range(Xpos)+c(-0.5,0.5),ylim=c(0,1))
text(mean( Xpos[  1:(3*T.t)] ), 0.5, "d=2",cex=1.5)
text(mean( Xpos[3*T.t+1:(3*T.t)] ), 0.5, "d=5",cex=1.5)
text(mean( Xpos[6*T.t+1:(3*T.t)] ), 0.5, "d=9",cex=1.5)

##

par(mar=MAR*c(1,0,1,0))
plot.new()
text(0.33, 0.5, expression(delta*"=0"), srt=90, cex=1.5)
text(0.67, 0.5, "(without X)", srt=90, cex=1.2)

##

par(mar=MAR*c(1,0,1,0))
plot.new()
text(0.33, 0.5, expression(delta*"=1"), srt=90, cex=1.5)
text(0.67, 0.5, "(with X)", srt=90, cex=1.2)

################################
# Table 1 of the main paper
# Constant ATT, no covariate, d=9
################################

Col1 <- c( "\\multicolumn{1}{|c|}{\\multirow{21}{*}{0}}",
           rep("\\multicolumn{1}{|c|}{}",20),
           "\\multicolumn{1}{|c|}{\\multirow{21}{*}{1}}",
           rep("\\multicolumn{1}{|c|}{}",20) )

AND <- rep(" & ",rep=42)

Col2 <- rep( c( "\\multicolumn{1}{c|}{\\multirow{7}{*}{2}}",
                rep("\\multicolumn{1}{c|}{}",6),
                "\\multicolumn{1}{c|}{\\multirow{7}{*}{5}}",
                rep("\\multicolumn{1}{c|}{}",6),
                "\\multicolumn{1}{c|}{\\multirow{7}{*}{9}}",
                rep("\\multicolumn{1}{c|}{}",6) ), 2)

Col3 <- rep( c("Bias ($\\times$10)",
               "ASE ($\\times$10)",
               "BSE ($\\times$10)",
               "ESE ($\\times$10)",
               "MSE ($\\times$100)",
               "Cover (ASE)",
               "Cover (BSE)"),6 )

END <- c( rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\hline \\hline",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\hline")


Summary <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.ASE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.BSE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", sd(Result[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result[,jj]^2)*100) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.ASECV[,jj])) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.BSECV[,jj])) ) )
}

Summary.X <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.ASE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.BSE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", sd(Result.X[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X[,jj]^2)*100) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.ASECV[,jj])) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.BSECV[,jj])) ) )
}

Result.Constant <- 
  cbind( Col1, AND, 
         Col2, AND,
         Col3, 
         rbind( cbind(sapply(1:12,Summary)),
                cbind(sapply(12+1:12,Summary)),
                cbind(sapply(24+1:12,Summary)),
                cbind(sapply(1:12,Summary.X)),
                cbind(sapply(12+1:12,Summary.X)),
                cbind(sapply(24+1:12,Summary.X)) ),
         END )


print(data.frame( Result.Constant ),row.names=F)


print(data.frame(cbind(Col3[1:7], 
                       sapply(24+1:12,Summary),
                       rep("\\\\ \\hline", 7) )),
      row.names=F)

################################
# Figures S2 and S2 of the supplement
# Linear ATT
################################

for(bloc in 1:2){
  
  Result <- MERGE(1,2,col.est,bloc)-3
  Result.ASE <- MERGE(1,2,col.ase,bloc)
  Result.BSE <- MERGE(1,2,col.bse,bloc)
  Result.ASECV <- MERGE(1,2,col.asecv,bloc)
  Result.BSECV <- MERGE(1,2,col.bsecv,bloc)
  
  Result.X <- MERGE(2,2,col.est,bloc)-3
  Result.X.ASE <- MERGE(2,2,col.ase,bloc)
  Result.X.BSE <- MERGE(2,2,col.bse,bloc)
  Result.X.ASECV <- MERGE(2,2,col.asecv,bloc)
  Result.X.BSECV <- MERGE(2,2,col.bsecv,bloc)
  
  T.t <- dim(Result)[2]/3/3
  Xpos  <- 1:dim(Result)[2] + rep(c(0:2, 4:6, 8:10),each=T.t)
  
  
  COL <- rep(rep( c(rgb(0,0,0,0.2),
                    rgb(0,0,0,0.5),
                    rgb(0,0,0,1)),
                  each=T.t),3)
  
  PCH <- rep(c(15,17,18,19)[1:T.t],9)
  CEX <- rep(c(1,1,1.5,1)[1:T.t],9)
  LTY <- rep(1,dim(Result)[2])
  
  Bias <- apply(Result,2,mean)
  UB   <- apply(Result,2,function(v){quantile(v,0.975)})
  LB   <- apply(Result,2,function(v){quantile(v,0.025)})
  
  YL <- c(floor(min(LB)*4)/4,ceiling(max(UB)*4)/4)
  
  MAR <- c(0.5,3,0.5,0.5)
  
  layout(matrix(c(8,5,9,
                  6,1,3,
                  7,2,4),3,3,byrow=T),
         widths=c(0.5,5,1),
         heights=c(0.5,5,5))
  
  par(mar=MAR)
  plot.new()
  plot.window(xlim=range(Xpos)+c(-0.5,0.5),ylim=YL)
  
  LINE <- function(jj){
    points(Xpos[jj],Bias[jj],col=COL[jj],pch=PCH[jj],cex=CEX[jj])
    segments(Xpos[jj],LB[jj],Xpos[jj],UB[jj],lty=LTY[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,LB[jj],Xpos[jj]+0.2,LB[jj],col=COL[jj])
    segments(Xpos[jj]-0.2,UB[jj],Xpos[jj]+0.2,UB[jj],col=COL[jj])
    
  }
  abline(h=0,col=2,lty=2)
  sapply(1:dim(Result)[2],LINE) 
  abline(v=mean(Xpos[c(3*T.t,3*T.t+1)]),col=1,lty=3)
  abline(v=mean(Xpos[3*T.t+c(3*T.t,3*T.t+1)]),col=1,lty=3)
  axis(2,at=seq(-1,1,by=0.5))
  axis(2,at=0,labels="Bias",tick=F,line=1,cex.axis=1.25)
  
  ##
  
  
  Bias <- apply(Result.X,2,mean)
  UB   <- apply(Result.X,2,function(v){quantile(v,0.975)})
  LB   <- apply(Result.X,2,function(v){quantile(v,0.025)})
  YL <- c(floor(min(LB)*4)/4,ceiling(max(UB)*4)/4)
  
  par(mar=MAR)
  plot.new()
  plot.window(xlim=range(Xpos)+c(-0.5,0.5),ylim=YL)
  
  abline(h=0,col=2,lty=2)
  sapply(1:dim(Result.X)[2],LINE) 
  abline(v=mean(Xpos[c(3*T.t,3*T.t+1)]),col=1,lty=3)
  abline(v=mean(Xpos[3*T.t+c(3*T.t,3*T.t+1)]),col=1,lty=3)
  # axis(3, 
  #      at=c( mean( Xpos[  1:(3*T.t)] ),
  #            mean( Xpos[3*T.t+1:(3*T.t)] ),
  #            mean( Xpos[6*T.t+1:(3*T.t)] ) ),
  #      labels=c("N=2","N=5","N=9"),
  #      tick=F,
  #      line=-2)
  
  axis(2,at=seq(-1,1,by=0.5))
  axis(2,at=0,labels="Bias",tick=F,line=1,cex.axis=1.25)
  
  ##
  
  par(mar=MAR*c(1,0,1,0))
  plot.new()
  
  YT <- c(9,7,5,3,1)/10
  XT <- c(0,0.05,0.25,0.3,
          0.6,0.65,0.7)
  
  text(XT[1],YT[1],"Estimator",pos=4,cex=1.5)
  segments(XT[2],YT[2],XT[3],YT[2],col=COL[1],lwd=1.5)
  text(XT[4],YT[2],"OLS",pos=4)
  segments(XT[2],YT[3],XT[3],YT[3],col=COL[5],lwd=1.5)
  text(XT[4],YT[3],"SPSC",pos=4)
  segments(XT[2],YT[4],XT[3],YT[4],col=COL[9],lwd=1.5)
  text(XT[4],YT[4],"SPSC-Reg",pos=4)
  
  ##
  
  par(mar=MAR*c(1,0,1,0))
  plot.new()
  
  YT <- c(9,7,5,3,1)/10
  XT <- c(0,0.05,0.25,0.3,
          0,0.15,0.3)
  
  text(XT[5],YT[1],expression(T["0"]),pos=4,cex=1.5)
  points(XT[6],YT[2],pch=PCH[1],cex=CEX[1])
  text(XT[7],YT[2],expression(T["0"]*"=50"),pos=4)
  points(XT[6],YT[3],pch=PCH[2],cex=CEX[2])
  text(XT[7],YT[3],expression(T["0"]*"=100"),pos=4)
  points(XT[6],YT[4],pch=PCH[3],cex=CEX[3])
  text(XT[7],YT[4],expression(T["0"]*"=250"),pos=4)
  points(XT[6],YT[5],pch=PCH[4],cex=CEX[4])
  text(XT[7],YT[5],expression(T["0"]*"=1000"),pos=4)
  
  
  ##
  
  par(mar=MAR*c(0,1,0,1))
  plot.new()
  plot.window(xlim=range(Xpos)+c(-0.5,0.5),ylim=c(0,1))
  text(mean( Xpos[  1:(3*T.t)] ), 0.5, "d=2",cex=1.5)
  text(mean( Xpos[3*T.t+1:(3*T.t)] ), 0.5, "d=5",cex=1.5)
  text(mean( Xpos[6*T.t+1:(3*T.t)] ), 0.5, "d=9",cex=1.5)
  
  ##
  
  par(mar=MAR*c(1,0,1,0))
  plot.new()
  text(0.33, 0.5, expression(delta*"=0"), srt=90, cex=1.5)
  text(0.67, 0.5, "(without X)", srt=90, cex=1.2)
  
  ##
  
  par(mar=MAR*c(1,0,1,0))
  plot.new()
  text(0.33, 0.5, expression(delta*"=1"), srt=90, cex=1.5)
  text(0.67, 0.5, "(with X)", srt=90, cex=1.2)
  
}

################################
# Tables S1-S3 of the supplement
# Constant ATT
################################

# (lengthb,bloc)=(1,1),(2,1),(2,2) correspond to Tables S1-S3, respectively

lengthb <- 2 
bloc    <- 2 

Result <- MERGE(1,lengthb,col.est,bloc)-3
Result.ASE <- MERGE(1,lengthb,col.ase,bloc)
Result.BSE <- MERGE(1,lengthb,col.bse,bloc)
Result.ASECV <- MERGE(1,lengthb,col.asecv,bloc)
Result.BSECV <- MERGE(1,lengthb,col.bsecv,bloc)

Result.X <- MERGE(2,lengthb,col.est,bloc)-3
Result.X.ASE <- MERGE(2,lengthb,col.ase,bloc)
Result.X.BSE <- MERGE(2,lengthb,col.bse,bloc)
Result.X.ASECV <- MERGE(2,lengthb,col.asecv,bloc)
Result.X.BSECV <- MERGE(2,lengthb,col.bsecv,bloc)



Col1 <- c( "\\multicolumn{1}{|c|}{\\multirow{21}{*}{0}}",
           rep("\\multicolumn{1}{|c|}{}",20),
           "\\multicolumn{1}{|c|}{\\multirow{21}{*}{1}}",
           rep("\\multicolumn{1}{|c|}{}",20) )

AND <- rep(" & ",rep=42)

Col2 <- rep( c( "\\multicolumn{1}{c|}{\\multirow{7}{*}{2}}",
                rep("\\multicolumn{1}{c|}{}",6),
                "\\multicolumn{1}{c|}{\\multirow{7}{*}{5}}",
                rep("\\multicolumn{1}{c|}{}",6),
                "\\multicolumn{1}{c|}{\\multirow{7}{*}{9}}",
                rep("\\multicolumn{1}{c|}{}",6) ), 2)

Col3 <- rep( c("Bias ($\\times$10)",
               "ASE ($\\times$10)",
               "BSE ($\\times$10)",
               "ESE ($\\times$10)",
               "MSE ($\\times$100)",
               "Cover (ASE)",
               "Cover (BSE)"),6 )

END <- c( rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\hline \\hline",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",6),
          "\\\\ \\hline")


Summary <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.ASE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.BSE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", sd(Result[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result[,jj]^2)*100) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.ASECV[,jj])) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.BSECV[,jj])) ) )
}

Summary.X <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.ASE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.BSE[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", sd(Result.X[,jj])*10) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X[,jj]^2)*100) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.ASECV[,jj])) ),
     ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result.X.BSECV[,jj])) ) )
}

Result.Constant <- 
  cbind( Col1, AND, 
         Col2, AND,
         Col3, 
         rbind( cbind(sapply(1:12,Summary)),
                cbind(sapply(12+1:12,Summary)),
                cbind(sapply(24+1:12,Summary)),
                cbind(sapply(1:12,Summary.X)),
                cbind(sapply(12+1:12,Summary.X)),
                cbind(sapply(24+1:12,Summary.X)) ),
         END )


print(data.frame( Result.Constant ),row.names=F)

################################
# Table S4 of the supplement
# Conformal inference 
################################

MERGE.CF <- function(dd,bb,VEC,bloc=1){
  Result.N2 <- cbind(
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]]][1], 
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]]][2], 
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]]][3]) 
  
  Result.N5 <- cbind(
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]]][1], 
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]]][2], 
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]]][3]) 
  
  Result.N9 <- cbind(
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]]][1], 
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]]][3]) 
  
  
  Result <- cbind(Result.N2,
                  Result.N5,
                  Result.N9)
  return(Result)
}





bloc <- 1

Col1 <- c( "\\multicolumn{1}{|c|}{\\multirow{6}{*}{Constant}}",
           rep("\\multicolumn{1}{|c|}{}",5),
           "\\multicolumn{1}{|c|}{\\multirow{6}{*}{Linear}}",
           rep("\\multicolumn{1}{|c|}{}",5) )

AND <- rep(" & ",rep=12)

Col2 <- rep( c( "\\multicolumn{1}{c|}{\\multirow{3}{*}{0}}",
                rep("\\multicolumn{1}{c|}{}",2),
                "\\multicolumn{1}{c|}{\\multirow{3}{*}{1}}",
                rep("\\multicolumn{1}{c|}{}",2) ), 2)

Col3 <- rep( c( "\\multicolumn{1}{c|}{2}",
                "\\multicolumn{1}{c|}{5}",
                "\\multicolumn{1}{c|}{9}"),
             4)

END <- c( rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\hline \\hline",
          rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\hline")


Summary11 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result11[,jj])) ) )
}
Summary12 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result12[,jj])) ) )
}
Summary21 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result21[,jj])) ) )
}
Summary22 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result22[,jj])) ) )
}

Result11 <- MERGE.CF(1,1, col.cfcv, bloc)
Result12 <- MERGE.CF(1,2, col.cfcv,bloc)
Result21 <- MERGE.CF(2,1, col.cfcv,bloc)
Result22 <- MERGE.CF(2,2, col.cfcv,bloc)

RRR1 <- cbind(Col1,AND,Col2,AND,Col3,
              rbind( matrix(sapply(1:dim(Result11)[2],Summary11),3,12,byrow=T),
                     matrix(sapply(1:dim(Result12)[2],Summary21),3,12,byrow=T),
                     matrix(sapply(1:dim(Result21)[2],Summary12),3,12,byrow=T),
                     matrix(sapply(1:dim(Result22)[2],Summary22),3,12,byrow=T)),
              END)

print( data.frame( RRR1 ), row.names=F)


################################
# Table S5 of the supplement
# Conformal inference 
################################

MERGE.CF <- function(dd,bb,VEC,bloc=1){
  Result.N2 <- cbind(
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]]][1], 
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]]][2], 
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]]][3]) 
  
  Result.N5 <- cbind(
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]]][1], 
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]]][2], 
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]]][3]) 
  
  Result.N9 <- cbind(
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]]][1], 
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]]][2],   
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]]][3],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]]][3]) 
  
  
  Result <- cbind(Result.N2,
                  Result.N5,
                  Result.N9)
  return(Result)
}





bloc <- 1

Col1 <- c( "\\multicolumn{1}{|c|}{\\multirow{6}{*}{Constant}}",
           rep("\\multicolumn{1}{|c|}{}",5),
           "\\multicolumn{1}{|c|}{\\multirow{6}{*}{Linear}}",
           rep("\\multicolumn{1}{|c|}{}",5) )

AND <- rep(" & ",rep=12)

Col2 <- rep( c( "\\multicolumn{1}{c|}{\\multirow{3}{*}{0}}",
                rep("\\multicolumn{1}{c|}{}",2),
                "\\multicolumn{1}{c|}{\\multirow{3}{*}{1}}",
                rep("\\multicolumn{1}{c|}{}",2) ), 2)

Col3 <- rep( c( "\\multicolumn{1}{c|}{2}",
                "\\multicolumn{1}{c|}{5}",
                "\\multicolumn{1}{c|}{9}"),
             4)

END <- c( rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\hline \\hline",
          rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",2),
          "\\\\ \\hline")


Summary11 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result11[,jj])) ) )
}
Summary12 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result12[,jj])) ) )
}
Summary21 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result21[,jj])) ) )
}
Summary22 <- function(jj){
  c( ( sprintf(" & \\multicolumn{1}{c|}{$%0.3f$}", mean(Result22[,jj])) ) )
}

Result11 <- MERGE.CF(1,1, col.cfl, bloc)
Result12 <- MERGE.CF(1,2, col.cfl,bloc)
Result21 <- MERGE.CF(2,1, col.cfl,bloc)
Result22 <- MERGE.CF(2,2, col.cfl,bloc)

RRR1 <- cbind(Col1,AND,Col2,AND,Col3,
              rbind( matrix(sapply(1:dim(Result11)[2],Summary11),3,12,byrow=T),
                     matrix(sapply(1:dim(Result12)[2],Summary21),3,12,byrow=T),
                     matrix(sapply(1:dim(Result21)[2],Summary12),3,12,byrow=T),
                     matrix(sapply(1:dim(Result22)[2],Summary22),3,12,byrow=T)),
              END)

print( data.frame( RRR1 ), row.names=F)

################################
# Table S6 of the supplement
# Conformal inference 
################################

MERGE.SCPI <- function(dd,bb,VEC,bloc=1){
  Result.N2 <- cbind(
    CSV[[dd]][[bb]][[1]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[1]][[4]][, VEC[[bloc]]][1]) 
  
  Result.N5 <- cbind(
    CSV[[dd]][[bb]][[2]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[2]][[4]][, VEC[[bloc]]][1]) 
  
  Result.N9 <- cbind(
    CSV[[dd]][[bb]][[3]][[1]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[2]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[3]][, VEC[[bloc]]][1],   
    CSV[[dd]][[bb]][[3]][[4]][, VEC[[bloc]]][1]) 
  
  
  Result <- cbind(Result.N2,
                  Result.N5,
                  Result.N9)
  return(Result)
}

Result11 <- MERGE.SCPI(1,1, col.scpi, 1)

apply(Result11 ,2, mean)-3
apply(Result11 ,2, sd)



Result11 <- MERGE.SCPI(1,1, col.scpi, 1)
Result21 <- MERGE.SCPI(2,1, col.scpi, 1)

apply(Result21 ,2, mean)-3
apply(Result21 ,2, sd)



TEMP <- rbind( c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result,2,mean)[5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result11-3,2,mean)[1:4]*10) ),
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result,2,sd)[5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result11,2,sd)[1:4]*10) ),
               
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result,2,mean)[12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result11-3,2,mean)[4+1:4]*10) ),
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result,2,sd)[12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result11,2,sd)[4+1:4]*10) ),
               
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result,2,mean)[12+12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result11-3,2,mean)[4+4+1:4]*10) ),
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result,2,sd)[12+12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result11,2,sd)[4+4+1:4]*10) ),
               
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result.X,2,mean)[5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result21-3,2,mean)[1:4]*10) ),
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result.X,2,sd)[5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result21,2,sd)[1:4]*10) ),
               
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result.X,2,mean)[12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result21-3,2,mean)[4+1:4]*10) ),
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result.X,2,sd)[12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result21,2,sd)[4+1:4]*10) ),
               
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result.X,2,mean)[12+12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result21-3,2,mean)[4+4+1:4]*10) ),
               c( sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result.X,2,sd)[12+12+5:12]*10), 
                  sprintf("& \\multicolumn{1}{c|}{%0.3f}",apply(Result21,2,sd)[4+4+1:4]*10) ) )




bloc <- 1

AND <- rep(" & ",rep=12)

Col1 <- c( "\\multicolumn{1}{|c|}{\\multirow{6}{*}{0}}",
           rep("\\multicolumn{1}{|c|}{}",5),
           "\\multicolumn{1}{|c|}{\\multirow{6}{*}{1}}",
           rep("\\multicolumn{1}{|c|}{}",5) )

Col2 <- rep( c( "\\multicolumn{1}{c|}{\\multirow{2}{*}{2}}",
                "\\multicolumn{1}{c|}{}",
                "\\multicolumn{1}{c|}{\\multirow{2}{*}{5}}",
                "\\multicolumn{1}{c|}{}",
                "\\multicolumn{1}{c|}{\\multirow{2}{*}{9}}",
                "\\multicolumn{1}{c|}{}"),
             2)

Col3 <- rep( c( "Bias ($\\times$10)",
                "ESE ($\\times$10)"),
             6)

END <- c( rep("\\\\ \\cline{3-15}",1),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",1),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",1),
          "\\\\ \\hline \\hline",
          rep("\\\\ \\cline{3-15}",1),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",1),
          "\\\\ \\cline{2-15}",
          rep("\\\\ \\cline{3-15}",1),
          "\\\\ \\hline")

print(data.frame(cbind(Col1, AND, Col2, AND, Col3, TEMP, END)),row.names=F)
