rm(list=ls())

library(png)
caption1 <- readPNG("Caption1.png")
caption2 <- readPNG("Caption2.png")
YL.True <- c(5.5,7.0)

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

gp <- 6
ATT.Type <- "constant"
YFYF <- "Result_Unit13"

T1 <- 167
T0 <- 217

PP <- function(gp,head=F){

  load(sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  
  ATT.Type <- "constant"
  LastC <- 0.3
  T1 <- 167
  T0 <- 217
  
  SHADE <- 0.5
  YL.True <- c(5.5,7.2)
  
  LWD <- c(1,2,1,1)
  COL <- c(1,
           rgb(0,0,0,0.5),
           rgb(0,0,0,1),
           1)
  COL.Overlap <- 1
  AREA.COL <- rgb(0,0,0,0.25)
  
  LTY <- c(1,3,1,1) 
   
  
  MAR <- c(1,4,0,0)
  layout(rbind(rbind(5:8,matrix(c(1:4),1,4)),c(9,9,9,10)), 
         widths=c(1,1,1,LastC), heights=c(0.75,8))
  par(mar=MAR)
  
  PLOT2 <- function(Actual.Time.Date,
                    Y1,
                    Y0,
                    Y0.LB,
                    Y0.UB,
                    label.x=F,
                    label.y=F,
                    label.att=F){
    
    
    one <- matrix(1/T1,T1,1)
    
    plot( Actual.Time.Date[1:T0], 
          Y0[1:T0], type='l', 
          xlim=range(Actual.Time.Date),
          ylim=YL.True, 
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
           at=mean(YL.True),
           labels="Outcome",
           tick=F,
           line=1.75,
           cex.axis=1.25)
    }
    
    text(Actual.Time.Date[1],
         5.75,
         sprintf("Average width\n= %0.3f",
                 abs(mean(Y0.UB-Y0.LB))),
         pos=4,
         cex=1.25)
  }
  
  load(sprintf("Conformal/GP%d_Data_ASC.RData",gp))
  
  T1 <- 167
  T0 <- 217
  
  PLOT2(Actual.Time.Date, 
        Y1 = Y1.series, 
        Y0 = Y1.series - ASC.ATT$att$Estimate, 
        Y0.LB = (Y1.series-ASC.ATT$att$upper_bound)[T0+1:T1],
        Y0.UB = (Y1.series-ASC.ATT$att$lower_bound)[T0+1:T1],
        label.y=T)
  
  load(sprintf("Conformal/GP%d_Data_SCPI.RData",gp))
  
  PLOT2(Actual.Time.Date, 
        Y1 = Y1.series, 
        Y0 = c(SCPI$est.results$Y.pre.fit,
               SCPI$est.results$Y.post.fit), 
        Y0.LB = apply(cbind(SCPI$inference.results$CI.all.gaussian[,1],SCPI$est.results$Y.post.fit),1,min),
        Y0.UB = apply(cbind(SCPI$inference.results$CI.all.gaussian[,2],SCPI$est.results$Y.post.fit),1,max),
        label.y=F)

  load(sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  CPCI <- read.csv(sprintf("Conformal/GP%d_Data_SPSC_CPCI.csv",gp))[-167,]
  
  COEF <- GAP <- ATT <- list()
  
  COEF[[1]] <- GMM.gamma.naive
  COEF[[2]] <- GMM.gamma.naive.lambda
  COEF[[3]] <- OLS.gamma.naive
  
  GAP[[1]]  <- Y1.series - Wmat.series%*%COEF[[1]]
  GAP[[2]]  <- Y1.series - Wmat.series%*%COEF[[2]]
  GAP[[3]]  <- Y1.series - Wmat.series%*%COEF[[3]]
  
  ATT[[1]]  <- c(rep(0,T0), rep(GMM.Simple.Coef[1:lengthb],T1))
  ATT[[2]]  <- c(rep(0,T0), rep(GMM.Regular.Coef[1:lengthb],T1))
  ATT[[3]]  <- c(rep(0,T0), rep(OLS.Simple.Coef[1:lengthb],T1))
  
  PLOT2(Actual.Time.Date, 
        Y1 = Y1.series, 
        Y0 = Wmat.series%*%COEF[[2]], 
        Y0.LB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_LB)*(-1),
        Y0.UB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_UB)*(-1),
        label.y=F,
        label.x=F,
        label.att=T)
  
  
  XLXL <- c(0,0.15,0.2)
  
  
  par(mar=MAR*c(1,0.25,1,0))
  plot.new()
  YP <- c(0.9,0.5,0.1)
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
  text(0.5,0.5,"ASC", cex=1.25)
  plot.new()
  text(0.5,0.5,"SCPI", cex=1.25)
  plot.new()
  text(0.5,0.5,"SPSC-Ridge", cex=1.25)
  plot.new()
  
  mtext("Year",
        side=1,
        outer=T,
        line=-1.0,
        at=(1.075)/(2+LastC),
        cex=0.8)
  
  # dev.off()
  
}


PPP <- function(gp,head=F){
  
  
  
  load(sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  
  
  ATT.Type <- "constant"
  LastC <- 0.4
  ATT.Type="constant"
  T1 <- 167
  T0 <- 217
  
  
  SHADE <- 0.5
  YL.True <- c(5.5,7.0)
  
  LWD <- c(1,2,1,1)
  COL <- c(1,
           rgb(0,0,0,0.5),
           rgb(0,0,0,1),
           1)
  COL.Overlap <- 1
  AREA.COL <- rgb(0,0,0,0.25)
  
  LTY <- c(1,3,1,1) 
  
  
  
  # png(sprintf("../plot_paper/%s_PI_GP%s.png",YFYF,gp),width=10,height=4.5,unit="in",res=500)
  
  MAR <- c(3,4,0.5,0)
  layout(rbind(rbind(6:10,matrix(c(1:5),1,5))), 
         widths=c(1,1,1,1,LastC), heights=c(1,8))
  par(mar=MAR,oma=c(0.5,0.5,1.5,0.5))
  
  PLOT2 <- function(Actual.Time.Date,
                    Y1,
                    Y0,
                    Y0.LB,
                    Y0.UB,
                    label.x=F,
                    label.y=F,
                    label.att=F){
    
    
    one <- matrix(1/T1,T1,1)
    
    plot( Actual.Time.Date[1:T0], 
          Y0[1:T0], type='l', 
          xlim=range(Actual.Time.Date),
          ylim=YL.True, 
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
           at=mean(YL.True),
           labels="Outcome",
           tick=F,
           line=1.75,
           cex.axis=1.25)
    }
     
    
    text(Actual.Time.Date[5],
         6.0,
         sprintf("Aver. width\n= %0.3f",
                 abs(mean(Y0.UB-Y0.LB))),
         cex=1.1,
         pos=4)
     
  }
  
  
  load(sprintf("Conformal/GP%d_Data_ASC.RData",gp))
  
  T1 <- 167
  T0 <- 217
  
  PLOT2(Actual.Time.Date, 
        Y1 = Y1.series, 
        Y0 = Y1.series - ASC.ATT$att$Estimate, 
        Y0.LB = (Y1.series-ASC.ATT$att$upper_bound)[T0+1:T1],
        Y0.UB = (Y1.series-ASC.ATT$att$lower_bound)[T0+1:T1],
        label.y=T)
  
  
  load(sprintf("Conformal/GP%d_Data_SCPI.RData",gp))
  
  PLOT2(Actual.Time.Date, 
        Y1 = Y1.series, 
        Y0 = c(SCPI$est.results$Y.pre.fit,
               SCPI$est.results$Y.post.fit), 
        Y0.LB = apply(cbind(SCPI$inference.results$CI.all.gaussian[,1],SCPI$est.results$Y.post.fit),1,min),
        Y0.UB = apply(cbind(SCPI$inference.results$CI.all.gaussian[,2],SCPI$est.results$Y.post.fit),1,max),
        label.y=F)
  
  
  load(sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  CPCI <- read.csv(sprintf("Conformal/GP%d_Data_SPSC_CPCI.csv",gp))[-167,]
  
  COEF <- GAP <- ATT <- list()
  
  COEF[[1]] <- GMM.gamma.naive
  COEF[[2]] <- GMM.gamma.naive.lambda
  COEF[[3]] <- OLS.gamma.naive
  
  GAP[[1]]  <- Y1.series - Wmat.series%*%COEF[[1]]
  GAP[[2]]  <- Y1.series - Wmat.series%*%COEF[[2]]
  GAP[[3]]  <- Y1.series - Wmat.series%*%COEF[[3]]
  
  ATT[[1]]  <- c(rep(0,T0), rep(GMM.Simple.Coef[1:lengthb],T1))
  ATT[[2]]  <- c(rep(0,T0), rep(GMM.Regular.Coef[1:lengthb],T1))
  ATT[[3]]  <- c(rep(0,T0), rep(OLS.Simple.Coef[1:lengthb],T1))
  
  PLOT2(Actual.Time.Date, 
        Y1 = Y1.series, 
        Y0 = Wmat.series%*%COEF[[2]], 
        Y0.LB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_LB)*(-1),
        Y0.UB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_UB)*(-1),
        label.y=F,
        label.x=F,
        label.att=T)
  
  
  
  load(sprintf("Result/Time_GP%d_Data_SPSC_ATT.RData",gp))
  CPCI <- read.csv(sprintf("Conformal/Time_GP%d_Data_SPSC_CPCI.csv",gp))[-167,]
  
  COEF <- GAP <- ATT <- list()
  
  COEF[[1]] <- GMM.gamma.naive
  COEF[[2]] <- GMM.gamma.naive.lambda
  COEF[[3]] <- OLS.gamma.naive
  
  GAP[[1]]  <- Y1.series - Wmat.series%*%COEF[[1]]
  GAP[[2]]  <- Y1.series - Wmat.series%*%COEF[[2]]
  GAP[[3]]  <- Y1.series - Wmat.series%*%COEF[[3]]
  
  ATT[[1]]  <- c(rep(0,T0), rep(GMM.Simple.Coef[1:lengthb],T1))
  ATT[[2]]  <- c(rep(0,T0), rep(GMM.Regular.Coef[1:lengthb],T1))
  ATT[[3]]  <- c(rep(0,T0), rep(OLS.Simple.Coef[1:lengthb],T1))
  
  PLOT2(Actual.Time.Date, 
        Y1 = Y1.series, 
        Y0 = Wmat.series%*%COEF[[2]], 
        Y0.LB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_LB)*(-1),
        Y0.UB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_UB)*(-1),
        label.y=F,
        label.x=F,
        label.att=T)
  
  
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
  text(0.5,0.3,expression("SPSC-Ridge, Time-varying "*g["t"]), cex=1.25)
  plot.new()
  
  # mtext("Year",
  #       side=1,
  #       outer=T,
  #       line=-1,
  #       at=(1.5)/(3+LastC),
  #       cex=0.95)
  
  # dev.off()
  
  if(head==T){
    mtext(sprintf("Donor Pool = Group %s",gp),
          side=3,
          outer=T,
          line=0,
          at=(1.5)/(3+LastC),
          cex=1.1)
  }
  
  
}


## Figure 2 of the main paper
YFYF <- "Result_Unit13"
png(sprintf("Result_GP%0.1d.png",6),height=3,width=11,unit="in",res=500)
PP(6)
dev.off()



## Figure S6 of the main paper
for(gp in 1:6){
  png(sprintf("Result_GP%0.1d_Title.png",gp),height=2.25,width=11.5,unit="in",res=500)
  PPP(gp,head=T)
  dev.off()
}
