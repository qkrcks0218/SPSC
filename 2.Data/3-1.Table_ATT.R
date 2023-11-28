Abadie.Estimator <- 
  rbind(c(-0.1197775 ,  0.771310895),
        c(-1.0157451 , -0.049089549),
        c(-0.8568228 ,  0.053757550),
        c(-0.8627475 , -0.008617287),
        c( 0.7041317 ,  1.576279078),
        c(-0.8776551 , -0.012510339) )

TABLE <- function(gp){
  
  ATT.Type="constant"
  T1 <- 167
  T0 <- 217
  
  one <- matrix(1/T1,T1,1)
  
  one <- matrix(1/T1,T1,1)
  
  
  load(sprintf("Conformal/GP%d_Data_ASC.RData",gp))
  load(sprintf("Conformal/GP%d_Data_SCPI.RData",gp))
  
  Co1 <- c(sprintf("& %0.3f",Abadie.Estimator[gp,1]),
           sprintf("& %0.3f",ASC.ATT$average_att$Estimate),
           sprintf("& %0.3f",mean(SCPI$data$Y.post - SCPI$est.results$Y.post.fit)) )
  Co2 <- Co3 <- rep("& -",3)
  
  load(sprintf("Result/GP%d_Data_SPSC_ATT.RData",gp))
  
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
  
  C4 <- c( sprintf("& %0.3f",Result.constant[1,3]),
           sprintf("& %0.3f",Result.constant[1,1]),
           sprintf("& %0.3f",Result.constant[1,2]) )
  
  C5 <- c( sprintf("& %0.3f",Result.constant[2,3]),
           sprintf("& %0.3f",Result.constant[2,1]),
           sprintf("& %0.3f",Result.constant[2,2]) )
  
  C6 <- c( sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,3]-qnorm(0.975)*Result.constant[2,3],
                   Result.constant[1,3]+qnorm(0.975)*Result.constant[2,3]),
           sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,1]-qnorm(0.975)*Result.constant[2,1],
                   Result.constant[1,1]+qnorm(0.975)*Result.constant[2,1]),
           sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,2]-qnorm(0.975)*Result.constant[2,2],
                   Result.constant[1,2]+qnorm(0.975)*Result.constant[2,2]) )
  
  
  C11 <- c("Estimate","ASE","95\\% CI")
  cbind(C11,
        rbind(c(Co1,C1,C4),c(Co2,C2,C5),c(Co3,C3,C6)),
        c("\\\\ \\cline{2-11}",
          "\\\\ \\cline{2-11}",
          "\\\\ \\hline"))
  
}

Group <- c("\\multicolumn{1}{|c|}{\\multirow{3}{*}{{\\begin{tabular}[c]{@{}c@{}}Group 1\\\\ (12)\\end{tabular}}}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{{\\begin{tabular}[c]{@{}c@{}}Group 2\\\\ (12)\\end{tabular}}}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{{\\begin{tabular}[c]{@{}c@{}}Group 3\\\\ (13)\\end{tabular}}}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{{\\begin{tabular}[c]{@{}c@{}}Group 4\\\\ (12)\\end{tabular}}}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{{\\begin{tabular}[c]{@{}c@{}}Group 5\\\\ (5)\\end{tabular}}}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{\\multirow{3}{*}{{\\begin{tabular}[c]{@{}c@{}}Group 6\\\\ (24)\\end{tabular}}}}",
           "\\multicolumn{1}{|c|}{}",
           "\\multicolumn{1}{|c|}{}")

YFYF <- "Result_Unit13"

DDDD <- cbind(Group, " & ",
              rbind(TABLE(1),
                    TABLE(2),
                    TABLE(3),
                    TABLE(4),
                    TABLE(5),
                    TABLE(6)))

print(data.frame(apply(DDDD,1,function(v){paste(v,collapse="")})),row.names=F)


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
  
  load(sprintf("Conformal/GP%d_Data_ASC.RData",gp))
  load(sprintf("Conformal/GP%d_Data_SCPI.RData",gp))
  
  C1 <- c( sprintf("& %0.3f",Result.constant[1,3]),
           sprintf("& %0.3f",Abadie.Estimator[gp,1]),
           sprintf("& %0.3f",ASC.ATT$average_att$Estimate),
           sprintf("& %0.3f",mean(SCPI$data$Y.post - SCPI$est.results$Y.post.fit)),
           sprintf("& %0.3f",Result.constant[1,1]),
           sprintf("& %0.3f",Result.constant[1,2]) )
  
  C2 <- c( sprintf("& %0.3f",Result.constant[2,3]),
           sprintf("& -"),
           sprintf("& -"),
           sprintf("& -"),
           sprintf("& %0.3f",Result.constant[2,1]),
           sprintf("& %0.3f",Result.constant[2,2]) )
  
  C3 <- c( sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,3]-qnorm(0.975)*Result.constant[2,3],
                   Result.constant[1,3]+qnorm(0.975)*Result.constant[2,3]),
           sprintf("& -"),
           sprintf("& -"),
           sprintf("& -"),
           sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,1]-qnorm(0.975)*Result.constant[2,1],
                   Result.constant[1,1]+qnorm(0.975)*Result.constant[2,1]),
           sprintf("& (%0.3f,%0.3f)",
                   Result.constant[1,2]-qnorm(0.975)*Result.constant[2,2],
                   Result.constant[1,2]+qnorm(0.975)*Result.constant[2,2]) )
  
  
  
  C11 <- c("Estimate","ASE","95\\% CI")
  RRR <- cbind(C11,
        rbind(c(C1),c(C2),c(C3)),
        c("\\\\ \\hline",
          "\\\\ \\hline",
          "\\\\ \\hline"))
  
  print(data.frame(apply(RRR,1,function(v){paste(v,collapse="")})),row.names=F)
  
}

## Table 3 of the main paper
TABLE.Paper(6)
