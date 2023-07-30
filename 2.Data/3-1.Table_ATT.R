TABLE <- function(gp){
  
  ATT.Type="constant"
  T1 <- 167
  T0 <- 217
  
  one <- matrix(1/T1,T1,1)
  
  one <- matrix(1/T1,T1,1)
  
  load(sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  
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
  
  
  
  load(sprintf("Result/Time_GP%d_Data_SPSC_ATT.RData",gp))
  
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
  
  C4 <- c( sprintf("& $%0.3f$",Result.constant[1,3]),
           sprintf("& $%0.3f$",Result.constant[1,1]),
           sprintf("& $%0.3f$",Result.constant[1,2]) )
  
  C5 <- c( sprintf("& $%0.3f$",Result.constant[2,3]),
           sprintf("& $%0.3f$",Result.constant[2,1]),
           sprintf("& $%0.3f$",Result.constant[2,2]) )
  
  C6 <- c( sprintf("& $(%0.3f, %0.3f)$",
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
        rbind(c(C1,C4),c(C2,C5),c(C3,C6)),
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

YFYF <- "Result_Unit13"

print(data.frame(cbind(Group, " & ",
                       rbind(TABLE(1),
                             TABLE(2),
                             TABLE(3),
                             TABLE(4),
                             TABLE(5),
                             TABLE(6)))),row.names=F)



Abadie.Estimator <- 
  rbind(c(-0.1197775 ,  0.771310895),
        c(-1.0157451 , -0.049089549),
        c(-0.8568228 ,  0.053757550),
        c(-0.8627475 , -0.008617287),
        c( 0.7041317 ,  1.576279078),
        c(-0.8776551 , -0.012510339) )


TABLE.Paper <- function(gp){
  
  ATT.Type="constant"
  T1 <- 167
  T0 <- 217
  
  one <- matrix(1/T1,T1,1)
  
  one <- matrix(1/T1,T1,1)
  
  load(sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  
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
           sprintf("& $%0.3f$",Abadie.Estimator[gp,1]),
           sprintf("& $%0.3f$",Result.constant[1,1]),
           sprintf("& $%0.3f$",Result.constant[1,2]) )
  
  C2 <- c( sprintf("& $%0.3f$",Result.constant[2,3]),
           sprintf("& -"),
           sprintf("& $%0.3f$",Result.constant[2,1]),
           sprintf("& $%0.3f$",Result.constant[2,2]) )
  
  C3 <- c( sprintf("& $(%0.3f, %0.3f)$",
                   Result.constant[1,3]-qnorm(0.975)*Result.constant[2,3],
                   Result.constant[1,3]+qnorm(0.975)*Result.constant[2,3]),
           sprintf("& -"),
           sprintf("& $(%0.3f, %0.3f)$",
                   Result.constant[1,1]-qnorm(0.975)*Result.constant[2,1],
                   Result.constant[1,1]+qnorm(0.975)*Result.constant[2,1]),
           sprintf("& $(%0.3f, %0.3f)$",
                   Result.constant[1,2]-qnorm(0.975)*Result.constant[2,2],
                   Result.constant[1,2]+qnorm(0.975)*Result.constant[2,2]) )
  
  
  
  C11 <- c("Estimate","ASE","95\\% CI")
  cbind(C11,
        rbind(c(C1),c(C2),c(C3)),
        c("\\\\ \\hline",
          "\\\\ \\hline",
          "\\\\ \\hline"))
  
}

## Table 3 of the main paper
TABLE.Paper(6)
