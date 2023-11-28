library(png)
caption1 <- readPNG("Caption1.png")
caption2 <- readPNG("Caption2.png")

addImg <- function(
    obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
    x = NULL, # mid x coordinate for image
    y = NULL, # mid y coordinate for image
    width = NULL, # width of image (in x coordinate units)
    interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}


ATT.Type <- "constant"
Abadie <- c(0.771310895,
            -0.049089549,
            0.053757550,
            -0.008617287,
            1.576279078,
            -0.012510339)

TABLE <- function(gp){
  
  ATT.Type="constant"
  T1 <- 36
  T0 <- 217 - T1
  
  one <- matrix(1/T1,T1,1)
  
  one <- matrix(1/T1,T1,1)
  
  load(sprintf("Result/Placebo_GP%d_Data_SPSC_ATT.RData",gp))
  
  Result.constant <- 
    rbind( c( as.numeric(t(one)%*%Post.Time.Basis%*%GMM.Simple.Coef[1:lengthb]), 
              as.numeric(t(one)%*%Post.Time.Basis%*%GMM.Regular.Coef[1:lengthb]),
              as.numeric(t(one)%*%Post.Time.Basis%*%OLS.Simple.Coef[1:lengthb]) ) ,
           
           c( as.numeric((t(one)%*%Post.Time.Basis)%*%
                           GMM.Simple.Var[1:lengthb,1:lengthb]%*%t(t(one)%*%Post.Time.Basis))^0.5,
              as.numeric((t(one)%*%Post.Time.Basis)%*%
                           GMM.Regular.Var[1:lengthb,1:lengthb]%*%t(t(one)%*%Post.Time.Basis))^0.5,
              as.numeric((t(one)%*%Post.Time.Basis)%*%
                           OLS.Simple.Var[1:lengthb,1:lengthb]%*%t(t(one)%*%Post.Time.Basis))^0.5 ),
           
           c( as.numeric((t(one)%*%Post.Time.Basis)%*%
                           GMM.Simple.Var.Boot[1:lengthb,1:lengthb]%*%t(t(one)%*%Post.Time.Basis))^0.5,
              as.numeric((t(one)%*%Post.Time.Basis)%*%
                           GMM.Regular.Var.Boot[1:lengthb,1:lengthb]%*%t(t(one)%*%Post.Time.Basis))^0.5,
              as.numeric((t(one)%*%Post.Time.Basis)%*%
                           OLS.Simple.Var.Boot[1:lengthb,1:lengthb]%*%t(t(one)%*%Post.Time.Basis))^0.5 ) )
  
  C1 <- c( sprintf("& %0.3f",Result.constant[1,3]),
           sprintf("& %0.3f",Result.constant[1,1]),
           sprintf("& %0.3f",Result.constant[1,2]) )
  
  C2 <- c( sprintf("& %0.3f",Result.constant[2,3]),
           sprintf("& %0.3f",Result.constant[2,1]),
           sprintf("& %0.3f",Result.constant[2,2]) )
  
  C3 <- c( sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,3]-qnorm(0.975)*Result.constant[2,3],
                   Result.constant[1,3]+qnorm(0.975)*Result.constant[2,3]),
           sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,1]-qnorm(0.975)*Result.constant[2,1],
                   Result.constant[1,1]+qnorm(0.975)*Result.constant[2,1]),
           sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,2]-qnorm(0.975)*Result.constant[2,2],
                   Result.constant[1,2]+qnorm(0.975)*Result.constant[2,2]) )
  
  load(sprintf("Conformal/Placebo_GP%d_Data_SCPI.RData",gp))
  SCPI.ATT.GP  <- mean( SCPI$data$Y.post - SCPI$est.results$Y.post.fit )
  
  load(sprintf("Conformal/Placebo_GP%d_Data_ASC.RData",gp))
  ASC.ATT.GP <- as.numeric(ASC.ATT$average_att[1])
  
  C1.Extend <- c(sprintf("& %0.3f",c(Abadie[gp],
                                   SCPI.ATT.GP,
                                   ASC.ATT.GP)),
                 C1)
  C2.Extend <- c("& -","& -","& -",C2)
  C3.Extend <- c("& -","& -","& -",C3)
  
  C11 <- c("Estimate","ASE","95\\% CI")
  cbind(C11,
        rbind(c(C1.Extend),c(C2.Extend),c(C3.Extend)),
        c("\\\\ \\cline{2-8}",
          "\\\\ \\cline{2-8}",
          "\\\\ \\hline"))
  
}

Group <- c("\\multicolumn{1}{|c|}{\\multirow{3}{*}{Group 1 (12)}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{Group 2 (12)}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{Group 3 (13)}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{Group 4 (12)}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{Group 5 (5)}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{Group 6 (24)}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}")

YFYF <- "Result_Unit13_Placebo"

## Table 12 of the supplement
print(data.frame(cbind(Group, " & ",
                       rbind(TABLE(1),
                             TABLE(2),
                             TABLE(3),
                             TABLE(4),
                             TABLE(5),
                             TABLE(6)))),row.names=F)





PPP <- function(gp,YL.True.P,head=F,TP=c(1,1,1)){
  
  
  load(sprintf("Result/Placebo_GP%d_Data_SPSC_ATT.RData",gp))
  
  ATT.Type <- "constant"
  LastC <- 0.4
  ATT.Type="constant"
  T1 <- 36
  T0 <- 181
  
  
  SHADE <- 0.5
  
  
  
  LWD <- c(1,2,1,1)
  COL <- c(1,
           rgb(0,0,0,0.5),
           rgb(0,0,0,1),
           1)
  COL.Overlap <- 1
  AREA.COL <- rgb(0,0,0,0.25)
  
  LTY <- c(1,3,1,1) 
  
  
  
  # png(sprintf("../plot_paper/%s_PI_GP%s.png",YFYF,gp),width=10,height=4.5,unit="in",res=500)
  
  MAR <- c(2,4,0.5,0)
  layout(rbind(rbind(5:8,matrix(c(1:4),1,4))), 
         widths=c(1,1,1,LastC), heights=c(1,8))
  par(mar=MAR,oma=c(0.5,0.5,1.5,0.5))
  
  PLOT2 <- function(Actual.Time.Date,
                    Y1,
                    Y0,
                    Y0.LB,
                    Y0.UB,
                    label.x=F,
                    label.y=F,
                    label.att=F,
                    TP.P){
    
    
    one <- matrix(1/T1,T1,1)
    
    plot( Actual.Time.Date[1:T0], 
          Y0[1:T0], type='l', 
          xlim=range(Actual.Time.Date),
          ylim=YL.True.P, 
          lwd=LWD[1], 
          col=COL[1],
          lty=LTY[1],
          ylab="",xlab="" )
    
    lines(Actual.Time.Date, 
          Y1, type='l', 
          lwd=LWD[2], 
          col=COL[2],
          lty=LTY[2]) 
    
    abline(v=Actual.Time.Date[T0],
           lty=3,
           col=1)
    
    polygon(c(Actual.Time.Date[T0+1:T1],
              Actual.Time.Date[T0+T1:1]),
            c( Y0.LB[1:T1] ,
               Y0.UB[T1:1]),
            col=AREA.COL,
            border=F)
    
    lines(Actual.Time.Date[T0+1:T1], 
          Y0[T0+1:T1], 
          type='l', 
          lwd=LWD[1], 
          col=COL.Overlap,
          lty=LTY[1])
    
    abline(h=0,col=1,lty=3)
    
    if(label.y){
      axis(2,
           at=mean(YL.True.P),
           labels="Outcome",
           tick=F,
           line=1.75,
           cex.axis=1.25)
    }
    
    
    text(Actual.Time.Date[5],
         YL.True.P[1]*(1-TP.P)+YL.True.P[2]*TP.P,
         sprintf("Aver. width\n= %0.3f",
                 abs(mean(Y0.UB-Y0.LB))),
         cex=1.1,
         pos=4)
    
  }
  
  
  load(sprintf("Conformal/Placebo_GP%d_Data_ASC.RData",gp))
  
  T1 <- 36
  T0 <- 181
  
  Actual.Time.Date.Short <- Actual.Time.Date[1:(T1+T0)]
  
  PLOT2(Actual.Time.Date.Short, 
        Y1 = Y1.series, 
        Y0 = Y1.series - ASC.ATT$att$Estimate, 
        Y0.LB = (Y1.series-ASC.ATT$att$upper_bound)[T0+1:T1],
        Y0.UB = (Y1.series-ASC.ATT$att$lower_bound)[T0+1:T1],
        label.y=T,
        TP.P=TP[1])
  
  
  load(sprintf("Conformal/Placebo_GP%d_Data_SCPI.RData",gp))
  
  PLOT2(Actual.Time.Date.Short, 
        Y1 = Y1.series, 
        Y0 = c(SCPI$est.results$Y.pre.fit,
               SCPI$est.results$Y.post.fit), 
        Y0.LB = apply(cbind(SCPI$inference.results$CI.all.gaussian[,1],SCPI$est.results$Y.post.fit),1,min),
        Y0.UB = apply(cbind(SCPI$inference.results$CI.all.gaussian[,2],SCPI$est.results$Y.post.fit),1,max),
        label.y=F,
        TP.P=TP[2])
  
  
  load(sprintf("Result/Placebo_GP%d_Data_SPSC_ATT.RData",gp))
  CPCI <- read.csv(sprintf("Conformal/Placebo_GP%d_Data_SPSC_CPCI.csv",gp))[-167,]
  
  COEF <- GAP <- ATT <- list()
  
  COEF[[1]] <- GMM.gamma.naive
  COEF[[2]] <- GMM.gamma.naive.lambda
  COEF[[3]] <- OLS.gamma.naive
  
  GAP[[1]]  <- Y1.series[1:(T0+T1)] - Wmat.series%*%COEF[[1]]
  GAP[[2]]  <- Y1.series[1:(T0+T1)] - Wmat.series%*%COEF[[2]]
  GAP[[3]]  <- Y1.series[1:(T0+T1)] - Wmat.series%*%COEF[[3]]
  
  ATT[[1]]  <- c(rep(0,T0), rep(GMM.Simple.Coef[1:lengthb],T1))
  ATT[[2]]  <- c(rep(0,T0), rep(GMM.Regular.Coef[1:lengthb],T1))
  ATT[[3]]  <- c(rep(0,T0), rep(OLS.Simple.Coef[1:lengthb],T1))
  
  PLOT2(Actual.Time.Date.Short, 
        Y1 = Y1.series[1:(T0+T1)], 
        Y0 = Wmat.series%*%COEF[[2]], 
        Y0.LB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_LB)*(-1),
        Y0.UB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_UB)*(-1),
        label.y=F,
        label.x=F,
        label.att=T,
        TP.P=TP[3])
   
  XLXL <- c(0,0.15,0.2)
  
  
  par(mar=MAR*c(1,0.25,1,0))
  plot.new()
  YP <- c(0.9,0.6,0.3)
  for(bb in 1:2){
    segments(XLXL[1],YP[bb],
             XLXL[2],YP[bb],
             lwd=LWD[bb], 
             col=COL[bb],
             lty=LTY[bb])
  }
  points( mean(XLXL[1:2]), YP[3], pch=15, cex=2.5, col=AREA.COL)
  
  addImg(caption1, x = 0.6, y = YP[1], width = 0.8)
  addImg(caption2, x = 0.6, y = YP[2], width = 0.8)
  text(XLXL[3], YP[3] ,"95%\nPrediction\nInterval",pos=4, cex=1.1)
  
  
  par(mar=MAR*c(0,1,0,1))
  
  plot.new()
  text(0.5,0.3,"ASC", cex=1.25)
  plot.new()
  text(0.5,0.3,"SCPI", cex=1.25)
  plot.new()
  text(0.5,0.3,expression("SPSC-Ridge, Time-invariant "*g), cex=1.25)
  plot.new()
   
  
  if(head==T){
    mtext(sprintf("Donor Pool = Group %s",gp),
          side=3,
          outer=T,
          line=0,
          at=(1.5)/(3+LastC),
          cex=1.1)
  }
  
  
}



## Figure S6 of the main paper
gp <- 1
png(sprintf("Result_Placebo_GP%0.1d.png",gp),height=2.25,width=10.5,unit="in",res=500)
PPP(gp,c(6.4,7),head=T,c(0.1,0.1,0.1))
dev.off()

gp <- 2
png(sprintf("Result_Placebo_GP%0.1d.png",gp),height=2.25,width=10.5,unit="in",res=500)
PPP(gp,c(6.725,6.95),head=T,c(0.1,0.1,0.1))
dev.off()

gp <- 3
png(sprintf("Result_Placebo_GP%0.1d.png",gp),height=2.25,width=10.5,unit="in",res=500)
PPP(gp,c(6.7,7),head=T,c(0.1,0.1,0.1))
dev.off()

gp <- 4
png(sprintf("Result_Placebo_GP%0.1d.png",gp),height=2.25,width=10.5,unit="in",res=500)
PPP(gp,c(6.5,7),head=T,c(0.1,0.1,0.1))
dev.off()

gp <- 5
png(sprintf("Result_Placebo_GP%0.1d.png",gp),height=2.25,width=10.5,unit="in",res=500)
PPP(gp,c(6.5,7.1),head=T,c(0.1,0.1,0.1))
dev.off()

gp <- 6
png(sprintf("Result_Placebo_GP%0.1d.png",gp),height=2.25,width=10.5,unit="in",res=500)
PPP(gp,c(6.75,6.925),head=T,c(0.1,0.1,0.1))
dev.off()
