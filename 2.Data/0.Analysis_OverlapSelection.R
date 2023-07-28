

Data <- read.csv("Data.csv")
Data$prc_log <- log( Data$mid_itp )
Data$date <- as.numeric( as.Date(Data$date) )+23715

# Firm ID 1 has missing time between date 665 (time 99) - 670 (time 101)
# We impute this value

Data <- rbind(Data[1:99,],
              Data[Data$ID==1 & Data$date==665,],
              Data[100:dim(Data)[1],])
Data$time <- rep(1:411,59)
Data$date <- rep(sort(unique(Data$date)),59)

# CUT
Start <- as.numeric( as.Date("1905-12-31") )+23715
END   <- as.numeric( as.Date("1909-01-02") )+23715
Data  <- Data[Start <= Data$date & Data$date <= END, ]

# 34, 37, 57 ; 37 = Lincoln

Actual.Time <- unique(Data$date)
Actual.Time.Date <- as.Date(Actual.Time-23715,origin="1970-01-01")

Donor.Index <- setdiff(unique(Data$ID),unique( Data$ID[Data$treat_a+Data$treat_c==1] ))



T.Post <- T1  <- sum(Data$treat_c)/3
T.Pre  <- T0  <- dim(Data)[1]/length(unique(Data$ID)) - T1


N <- length(Donor.Index)

Wmat.series <- matrix(0,T0+T1,N)

for(w.iter in 1:length(Donor.Index)){
  d.index <- Donor.Index[w.iter]
  Wmat.series[,w.iter] <- Data[Data$ID==d.index,]$prc_log
}

Wmat.Pre <- Wmat.series[(1:T0),]
Wmat.Post <- Wmat.series[T0+(1:T1),]

overlap <- cbind(1:N,(sapply(1:N,
                             function(bb){
                               RR1 <- range(Wmat.series[1:T0,bb])
                               mean(c(as.numeric(RR1[1] <= Wmat.series[T0+1:T1,bb] &
                                                   Wmat.series[T0+1:T1,bb] <= RR1[2] ) ))
                               
                             })))


layout(matrix(1:2,1,2),widths=c(3,1.2))
MAR <- c(3.5, 3.5, 1, 0.5)

par(mar=MAR)

plot(sort(overlap[,2]),
     xlab="",
     ylab="",
     pch=c(rep(1,12),
           rep(3,13),
           rep(4,12),
           rep(0,12)))
title(xlab="Donor Candidates",
      ylab="Overlap",
      line=2.5)
# abline(h=0.75,
#        col=2,
#        lty=3,
#        lwd=1.5)
# abline(h=0.185,
#        col=2,
#        lty=3,
#        lwd=1.5)
# abline(h=0.0685,
#        col=2,
#        lty=3,
#        lwd=1.5)


par(mar=MAR*c(1,0,1,0))

plot.new()

XL <- c(0.1,0.2)
YL <- c(0.9, 0.7, 0.5, 0.3)

points(rep(XL[1],4),
       YL,
       pch=c(0,4,3,1))
text(rep(XL[2],4),
     YL,
     c("Group 1 (12 donors)",
       "Group 2 (12 donors)",
       "Group 3 (13 donors)",
       "Group 4 (12 donors)"),
     pos=4)
