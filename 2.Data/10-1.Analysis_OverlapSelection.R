rm(list=ls())

##################################
# Functions for the SPSC approach
##################################

source("0.Function_SPSC_Data.R")

ATT.Type <- "constant"
lengthb  <- 1  
Ypos     <- c(1,3)
gY.Bound <- c(6.5,7.1) 

##################################
# Library + Data Cleaning
##################################

library(splines)
library(MASS)
source("0.DataCleaning.R")

##################################
# Overlap Metric
##################################

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

SEQ <- c(12,24,36,48)+0.5
(1:N)[ (50-rank(overlap[,2]) < SEQ[1]) ]
(1:N)[ (SEQ[1] <= 50-rank(overlap[,2]) & 50-rank(overlap[,2]) < SEQ[2]) ]
(1:N)[ (SEQ[2] <= 50-rank(overlap[,2]) & 50-rank(overlap[,2]) < SEQ[3]) ]
(1:N)[ (SEQ[3] <= 50-rank(overlap[,2]) & 50-rank(overlap[,2]) < SEQ[4]) ]
