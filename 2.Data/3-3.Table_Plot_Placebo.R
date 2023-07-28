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
TABLE <- function(gp){
  
  ATT.Type="constant"
  T1 <- 36
  T0 <- 217 - T1
  
  one <- matrix(1/T1,T1,1)
  
  one <- matrix(1/T1,T1,1)
  
  load(sprintf("Result/Placebo_GP%d_Data_SPSC_ATT.RData",gp))
  
  Result.constant <- rbind( c( as.numeric(t(one)%*%Post.Time.Basis%*%GMM.Simple.Coef[1:lengthb]), 
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
  
  C1 <- c( sprintf("& $%0.3f$",Result.constant[1,3]),
           sprintf("& $%0.3f$",Result.constant[1,1]),
           sprintf("& $%0.3f$",Result.constant[1,2]) )
  
  C2 <- c( sprintf("& $%0.3f$",Result.constant[2,3]),
           sprintf("& $%0.3f$",Result.constant[2,1]),
           sprintf("& $%0.3f$",Result.constant[2,2]) )
  
  C3 <- c( sprintf("& $(%0.3f, %0.3f)$",
                   Result.constant[1,3]-qnorm(0.975)*Result.constant[2,3],
                   Result.constant[1,3]+qnorm(0.975)*Result.constant[2,3]),
           sprintf("& $(%0.3f, %0.3f)$",
                   Result.constant[1,1]-qnorm(0.975)*Result.constant[2,1],
                   Result.constant[1,1]+qnorm(0.975)*Result.constant[2,1]),
           sprintf("& $(%0.3f, %0.3f)$",
                   Result.constant[1,2]-qnorm(0.975)*Result.constant[2,2],
                   Result.constant[1,2]+qnorm(0.975)*Result.constant[2,2]) )
  
  
  
  C11 <- c("Estimate","ASE","95\\% CI")
  cbind(C11,
        rbind(c(C1),c(C2),c(C3)),
        c("\\\\ \\cline{2-5}",
          "\\\\ \\cline{2-5}",
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
        ylim=YL, 
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
         at=mean(YL),
         labels="Outcome",
         tick=F,
         line=1.75,
         cex.axis=1.25)
  } 
  
  text(Actual.Time.Date[5],
       6.98,
       sprintf("Aver. width = %0.4f",
               abs(mean(Y0.UB-Y0.LB))),
       cex=1.1,
       pos=4)
   
}




## Figure S7 of the supplement
ATT.Type <- "constant"
LastC <- 0.4
ATT.Type="constant"
T1 <- 167
T0 <- 217


SHADE <- 0.5
YL <- c(6.7,7.0)

LWD <- c(1,2,1,1)
COL <- c(1,
         rgb(0,0,0,0.5),
         rgb(0,0,0,1),
         1)
COL.Overlap <- 1
AREA.COL <- rgb(0,0,0,0.25)

LTY <- c(1,3,1,1) 


MAR <- c(3,4,2,0)
layout(matrix(c(1,4,2,5,3,6,7,7),2,4), 
       widths=c(1,1,1,0.5), heights=c(8,8))
par(mar=MAR,oma=c(1.5,0.5,0.5,0.5))



for(gp in 1:6){
  
  load(sprintf("Result/Placebo_GP%d_Data_SPSC_ATT.RData",gp))
  CPCI <- read.csv(sprintf("Conformal/Placebo_GP%d_Data_SPSC_CPCI.csv",gp))
  
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
  
  PLOT2(Actual.Time.Date[1:(T0+T1)],
        Y1 = Y1.series[1:(T0+T1)],
        Y0 = Wmat.series%*%COEF[[2]],
        Y0.LB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_LB)*(-1),
        Y0.UB = (-Y1.series[T0+1:T1] + CPCI$SPSC_Reg_UB)*(-1),
        label.y=F,
        label.x=F,
        label.att=T)
  title(main=sprintf("Donor Pool = Group %s",gp), line=1, font.main=1)
  
}


mtext("Year",
      side=1,
      outer=T,
      line=-0.5,
      at=(1.5)/(3+LastC),
      cex=0.95)

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

