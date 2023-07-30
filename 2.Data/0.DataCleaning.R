##################################
# Data Cleaning
##################################

Ypos <- c(1,3)
Data <- read.csv("Data.csv")
Data$prc_log <- log( Data$mid_itp )
Data$date <- as.numeric( as.Date(Data$date) )+23715

# Firm ID 1 has missing time between date 665 (time 99) - 670 (time 101)
# We impute this value

(table(Data$ID, Data$date))[1,] 
Data[Data$ID==1 & Data$date==665,]
Data[Data$ID==1 & Data$date==670,]

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

Post.Time.Basis     <- rep(1,T1)

N <- length(Donor.Index)

Wmat.series <- matrix(0,T0+T1,N)

for(w.iter in 1:length(Donor.Index)){
  d.index <- Donor.Index[w.iter]
  Wmat.series[,w.iter] <- Data[Data$ID==d.index,]$prc_log
}

Wmat.series  <- Wmat.series
Ymat.series <- rep(0,T0+T1)
Ymat.series <- cbind(Data$prc_log[Data$ID==34],
                     Data$prc_log[Data$ID==37],
                     Data$prc_log[Data$ID==57])
Y1.series <- apply(Ymat.series[,Ypos],1,mean)
Ymat.Pre  <- Ymat.series[(1:T0),]
Ymat.Post <- Ymat.series[T0+(1:T1),]

## Pre-treatment series 

Wmat.Pre <- Wmat.series[(1:T0),] 
Y1.Pre   <- Y1.series[(1:T0)] 

Wmat.Post <- Wmat.series[T0+(1:T1),] 
Y1.Post   <- Y1.series[T0+(1:T1)] 

Wmat.series <- Wmat.series 
Wmat.Pre    <- Wmat.Pre    
Wmat.Post   <- Wmat.Post   

Wmat.series.work <- Wmat.series

