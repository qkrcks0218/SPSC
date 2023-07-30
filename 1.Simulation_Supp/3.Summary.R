# Merge Result

# LFT <- list.files(sprintf("%s/Result_Raw",getwd()))
# RRR <- read.csv(sprintf("%s/Result_Raw/%s",getwd(),LFT[1]))
# if(length(LFT)>1){
#   for(lft in 2:length(LFT)){
#     RRR <- rbind(RRR,read.csv(sprintf("%s/Result_Raw/%s",getwd(),LFT[lft])))
#   }
# }
# write.csv(RRR,
#           sprintf("Result/Result_Static.csv"),
#           row.names=F) 
# 
# LFT <- list.files(sprintf("%s/Result_Time_Raw",getwd()))
# RRR <- read.csv(sprintf("%s/Result_Time_Raw/%s",getwd(),LFT[1]))
# if(length(LFT)>1){
#   for(lft in 2:length(LFT)){
#     RRR <- rbind(RRR,read.csv(sprintf("%s/Result_Time_Raw/%s",getwd(),LFT[lft])))
#   }
# }
# write.csv(RRR,
#           sprintf("Result/Result_Time.csv"),
#           row.names=F) 


Result.A.S <- read.csv("Result/Result_Static.csv")
Result.A.T <- read.csv("Result/Result_Time.csv")

COVER.Index <- c(c(1,2,3),
                 c(1,2,3)+9,
                 c(1,2,3)+18,
                 c(1,2,3)+27,
                 c(1,2,3)+36,
                 c(c(1,2,3),
                   c(1,2,3)+9,
                   c(1,2,3)+18,
                   c(1,2,3)+27,
                   c(1,2,3)+36)+45,
                 c(c(1,2,3),
                   c(1,2,3)+9,
                   c(1,2,3)+18,
                   c(1,2,3)+27,
                   c(1,2,3)+36)+90)

Length.Index <- COVER.Index + 3
Bias.Index   <- COVER.Index + 6

SPSC.NoR <- (1:15)*3-2
SPSC.RR  <- (1:15)*3-1
SPCI     <- (1:15)*3

############################################
# Table S7 of the supplement
############################################

ROWS <- list()

ROWS[[1]] <- apply(Result.A.S,2,mean)[COVER.Index[SPSC.RR]]
ROWS[[2]] <- apply(Result.A.T,2,mean)[COVER.Index[SPSC.RR]]
ROWS[[3]] <- apply(Result.A.S,2,mean)[COVER.Index[SPCI]]


RRR <- rbind(round(ROWS[[1]],3),
             round(ROWS[[2]],3),
             round(ROWS[[3]],3))


R1 <- "\\multicolumn{2}{|c|}{$\\rho$} & \\multicolumn{5}{c|}{0} & \\multicolumn{5}{c|}{0.5} & \\multicolumn{5}{c|}{1} \\\\ \\hline"
R2 <- "\\multicolumn{2}{|c|}{$\\zeta$} & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 \\\\ \\hline"
C1 <- c("\\multicolumn{1}{|c|}{\\multirow{2}{*}{SPSC-Ridge}}",
        rep("\\multicolumn{1}{|c|}{}",1),
        "\\multicolumn{2}{|c|}{SCPI}")
C2 <- c(" & Time-invariant $\bg$",
        " & Time-varying $\bg_t$","")
CC <- list()
for(tt in 1:14){
  CC[[tt]] <- c(sprintf(" & \\multicolumn{1}{c|}{%0.3f}", RRR[,tt]))
}
tt <- 15
CC[[tt]] <- c(sprintf(" & %0.3f", RRR[,tt]))

CCC <- cbind(CC[[1]],CC[[2]])
for(tt in 3:15){
  CCC <- cbind(CCC,CC[[tt]])
}

C4 <- c(" \\\\ \\cline{2-17}",
        " \\\\ \\hline"," \\\\ \\hline")

RRR.Print <- c(R1,R2,apply(cbind(C1,C2,CCC,C4),1,function(v){paste(v,collapse="")}))

RRR.Print <- as.data.frame(RRR.Print)

print(RRR.Print,row.names=F)  


############################################
# Table S8 of the supplement
############################################


ROWS <- list()

ROWS[[1]] <- apply(Result.A.S,2,mean)[Length.Index[SPSC.RR]]
ROWS[[2]] <- apply(Result.A.T,2,mean)[Length.Index[SPSC.RR]]
ROWS[[3]] <- apply(Result.A.S,2,mean)[Length.Index[SPCI]]


RRR <- rbind(round(ROWS[[1]],3),
             round(ROWS[[2]],3),
             round(ROWS[[3]],3))


R1 <- "\\multicolumn{2}{|c|}{$\\rho$} & \\multicolumn{5}{c|}{0} & \\multicolumn{5}{c|}{0.5} & \\multicolumn{5}{c|}{1} \\\\ \\hline"
R2 <- "\\multicolumn{2}{|c|}{$\\zeta$} & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 \\\\ \\hline"
C1 <- c("\\multicolumn{1}{|c|}{\\multirow{2}{*}{SPSC-Ridge}}",
        rep("\\multicolumn{1}{|c|}{}",1),
        "\\multicolumn{2}{|c|}{SCPI}")
C2 <- c(" & Time-invariant $\bg$",
        " & Time-varying $\bg_t$","")
CC <- list()
for(tt in 1:14){
  CC[[tt]] <- c(sprintf(" & \\multicolumn{1}{c|}{%0.3f}", RRR[,tt]))
}
tt <- 15
CC[[tt]] <- c(sprintf(" & %0.3f", RRR[,tt]))

CCC <- cbind(CC[[1]],CC[[2]])
for(tt in 3:15){
  CCC <- cbind(CCC,CC[[tt]])
}

C4 <- c(" \\\\ \\cline{2-17}",
        " \\\\ \\hline"," \\\\ \\hline")

RRR.Print <- c(R1,R2,apply(cbind(C1,C2,CCC,C4),1,function(v){paste(v,collapse="")}))

RRR.Print <- as.data.frame(RRR.Print)

print(RRR.Print,row.names=F)  


############################################
# Table S9 of the supplement
############################################


ROWS <- list()

ROWS[[1]] <- apply(Result.A.S,2,mean)[Bias.Index[SPSC.RR]]
ROWS[[2]] <- apply(Result.A.T,2,mean)[Bias.Index[SPSC.RR]]
ROWS[[3]] <- apply(Result.A.S,2,mean)[Bias.Index[SPCI]]


RRR <- rbind(round(ROWS[[1]],3),
             round(ROWS[[2]],3),
             round(ROWS[[3]],3))


R1 <- "\\multicolumn{2}{|c|}{$\\rho$} & \\multicolumn{5}{c|}{0} & \\multicolumn{5}{c|}{0.5} & \\multicolumn{5}{c|}{1} \\\\ \\hline"
R2 <- "\\multicolumn{2}{|c|}{$\\zeta$} & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 \\\\ \\hline"
C1 <- c("\\multicolumn{1}{|c|}{\\multirow{2}{*}{SPSC-Ridge}}",
        rep("\\multicolumn{1}{|c|}{}",1),
        "\\multicolumn{2}{|c|}{SCPI}")
C2 <- c(" & Time-invariant $\bg$",
        " & Time-varying $\bg_t$","")
CC <- list()
for(tt in 1:14){
  CC[[tt]] <- c(sprintf(" & \\multicolumn{1}{c|}{%0.3f}", RRR[,tt]))
}
tt <- 15
CC[[tt]] <- c(sprintf(" & %0.3f", RRR[,tt]))

CCC <- cbind(CC[[1]],CC[[2]])
for(tt in 3:15){
  CCC <- cbind(CCC,CC[[tt]])
}

C4 <- c(" \\\\ \\cline{2-17}",
        " \\\\ \\hline"," \\\\ \\hline")

RRR.Print <- c(R1,R2,apply(cbind(C1,C2,CCC,C4),1,function(v){paste(v,collapse="")}))

RRR.Print <- as.data.frame(RRR.Print)

print(RRR.Print,row.names=F)  


############################################
# Table S9 of the supplement
############################################


ROWS <- list()

ROWS[[1]] <- apply(Result.A.S,2,sd)[Bias.Index[SPSC.RR]]
ROWS[[2]] <- apply(Result.A.T,2,sd)[Bias.Index[SPSC.RR]]
ROWS[[3]] <- apply(Result.A.S,2,sd)[Bias.Index[SPCI]]


RRR <- rbind(round(ROWS[[1]],10),
             round(ROWS[[2]],10),
             round(ROWS[[3]],3))


R1 <- "\\multicolumn{2}{|c|}{$\\rho$} & \\multicolumn{5}{c|}{0} & \\multicolumn{5}{c|}{0.5} & \\multicolumn{5}{c|}{1} \\\\ \\hline"
R2 <- "\\multicolumn{2}{|c|}{$\\zeta$} & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 & \\multicolumn{1}{c|}{-1} & \\multicolumn{1}{c|}{-0.5} & \\multicolumn{1}{c|}{0} & \\multicolumn{1}{c|}{0.5} & 1 \\\\ \\hline"
C1 <- c("\\multicolumn{1}{|c|}{\\multirow{2}{*}{SPSC-Ridge}}",
        rep("\\multicolumn{1}{|c|}{}",1),
        "\\multicolumn{2}{|c|}{SCPI}")
C2 <- c(" & Time-invariant $\bg$",
        " & Time-varying $\bg_t$","")
CC <- list()
for(tt in 1:14){
  CC[[tt]] <- c(sprintf(" & \\multicolumn{1}{c|}{%0.3f}", RRR[,tt]))
}
tt <- 15
CC[[tt]] <- c(sprintf(" & %0.3f", RRR[,tt]))

CCC <- cbind(CC[[1]],CC[[2]])
for(tt in 3:15){
  CCC <- cbind(CCC,CC[[tt]])
}

C4 <- c(" \\\\ \\cline{2-17}",
        " \\\\ \\hline"," \\\\ \\hline")

RRR.Print <- c(R1,R2,apply(cbind(C1,C2,CCC,C4),1,function(v){paste(v,collapse="")}))

RRR.Print <- as.data.frame(RRR.Print)

print(RRR.Print,row.names=F)  

