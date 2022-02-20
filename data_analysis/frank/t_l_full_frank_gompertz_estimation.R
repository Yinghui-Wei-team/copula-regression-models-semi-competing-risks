# Real data analysis:NHSBT data
# YW: 15 July 2021 updates:
#     1. Frank Gompertz survival models
#     2. make unified data set for analysis paper2_data.csv
#     3. Output results into a vector
# original script by LS; edited and updated for paper2 by YW
# model: gompertz surival distribution; frank copula; regression on both hazard rates (lambda) and association parameter (theta)
# table 4 gompertz survival models with frank copula
# script name: t - theta; l- lambda; 
#              full - regression models on both theta and lambda
###############################################################
rm(list=ls())
library(copula)
library(mvtnorm)
library(ggplot2)
library(plyr)
library(survival)


########################################################
##################### Load data ########################
########################################################
# YW: need to firstly set working directory to project directory and send through the next two lines
setwd("../../../")
df <- read.csv(file="NHSBT/paper2_data.csv")
dim(df)
attach(df)
names(df)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}


########################################################
############### Frank pseudo likelihood ################
########################################################
start_time = Sys.time()
fpl <- function(para, X, Y, d1, d2, donor, age.grp, gen){
  gamma1 <- para[1]
  p0 <- para[2]
  p1 <- para[3]
  p2 <- para[4]
  p3 <- para[5]
  gamma2 <- para[6]
  q0 <- para[7]
  q1 <- para[8]
  q2 <- para[9]
  q3 <- para[10]
  
  b0 <- para[11]
  b1 <- para[12]
  b2 <- para[13]
  b3 <- para[14]
  
  lambda1 <- exp(p0+p1*age.grp+p2*gen+p3*donor)
  lambda2 <- exp(q0+q1*age.grp+q2*gen+q3*donor)
  
  S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
  S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
  S1[which(S1<0.1^8)] <- 0.1^8
  S2[which(S2<0.1^8)] <- 0.1^8
  
  f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
  f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
  f1[which(f1<0.1^8)] <- 0.1^8
  f2[which(f2<0.1^8)] <- 0.1^8
  
  theta <- exp(b0+b1*age.grp+b2*gen+b3*donor)
  
  C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
  C[which(C < 0.1^8)] <- 0.1^8 
  
  part1 <- d1*d2*(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*S1)-1)-log(exp(theta*S2)-1)+log(f1)+log(f2))
  part2 <- d1*(1-d2)*log(((1-exp(theta*C))*f1)/(1-exp(theta*S1)))
  part3 <- (1-d1)*d2*log(((1-exp(theta*C))*f2)/(1-exp(theta*S2)))
  part4<-((1-d1)*(1-d2))*log(C)
  
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

#fpl(c(0.2,-1,1,1,1, 0.2,-2,2,1,0.01  ,-2,0.1,0.1,0.1),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen)


plfoptim <- optim(c(0.02, -1,-0.01,-0.01,-0.01,  0.02,-1,-0.01,-0.01,-0.01,  2,0.1,0.1,0.1), fpl, method="L-BFGS-B",
                  lower=c(-0.2,-10,-10,-10,-10,  -0.2,-10,-10,-10,-10,  0.01,-1,-1,-1),
                  upper=c(0.2,-1,1,1,1, 0.2,-2,2,1,0.01  ,10,6,3,3), 
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen,
                  control=list(fnscale=-1),hessian=TRUE)

end_time = Sys.time()
run_time = end_time - start_time
run_time


plfoptim$par

########################################################
################## Confidence Intervals ################
########################################################

#Fisher's Information matrix
fisher_info<-solve(-plfoptim$hessian) #inverse -hess
#Standard error = sqrt(var/n)
se<-sqrt(diag(fisher_info)) 

#beta ci#
est_b0 <- plfoptim$par[11]
est_b1 <- plfoptim$par[12]
est_b2 <- plfoptim$par[13]
est_b3 <- plfoptim$par[14]
lwci_b0 <- est_b0 - 1.96*se[11]     
lwci_b1 <- est_b1 - 1.96*se[12]
lwci_b2 <- est_b2 - 1.96*se[13]     
lwci_b3 <- est_b3 - 1.96*se[14]
upci_b0 <- est_b0 + 1.96*se[11]
upci_b1 <- est_b1 + 1.96*se[12] 
upci_b2 <- est_b2 + 1.96*se[13]
upci_b3 <- est_b3 + 1.96*se[14] 

#g ci#
est_g1 <- plfoptim$par[1]
est_g2 <- plfoptim$par[6]
lwci_g1 <- est_g1 - 1.96*se[1]
lwci_g2 <- est_g2 - 1.96*se[6]     
upci_g1 <- est_g1 + 1.96*se[1] 
upci_g2 <- est_g2 + 1.96*se[6]

#p ci#
est_p0 <- plfoptim$par[2]
est_p1 <- plfoptim$par[3]
est_p2 <- plfoptim$par[4]
est_p3 <- plfoptim$par[5]
lwci_p0 <- est_p0 - 1.96*se[2]     
lwci_p1 <- est_p1 - 1.96*se[3]
lwci_p2 <- est_p2 - 1.96*se[4]     
lwci_p3 <- est_p3 - 1.96*se[5]
upci_p0 <- est_p0 + 1.96*se[2]
upci_p1 <- est_p1 + 1.96*se[3] 
upci_p2 <- est_p2 + 1.96*se[4]
upci_p3 <- est_p3 + 1.96*se[5] 

#q ci#
est_q0 <- plfoptim$par[7]
est_q1 <- plfoptim$par[8]
est_q2 <- plfoptim$par[9]
est_q3 <- plfoptim$par[10]
lwci_q0 <- est_q0 - 1.96*se[7]     
lwci_q1 <- est_q1 - 1.96*se[8]
lwci_q2 <- est_q2 - 1.96*se[9]     
lwci_q3 <- est_q3 - 1.96*se[10]
upci_q0 <- est_q0 + 1.96*se[7]
upci_q1 <- est_q1 + 1.96*se[8] 
upci_q2 <- est_q2 + 1.96*se[9]
upci_q3 <- est_q3 + 1.96*se[10] 

#hrs
var_p1 <- fisher_info[3,3]
var_p2 <- fisher_info[4,4]
var_p3 <- fisher_info[5,5]
var_q1 <- fisher_info[8,8]
var_q2 <- fisher_info[9,9]
var_q3 <- fisher_info[10,10]

esthr_l1_age <- exp(est_p1)
esthr_l1_gen <- exp(est_p2)
esthr_l1_donor <- exp(est_p3)

esthr_l2_age <- exp(est_q1)
esthr_l2_gen <- exp(est_q2)
esthr_l2_donor <- exp(est_q3)

var_hr_l1_age <- exp(est_p1)^2 * var_p1
var_hr_l1_gen <- exp(est_p2)^2 * var_p2
var_hr_l1_donor <- exp(est_p3)^2 * var_p3

var_hr_l2_age <- exp(est_q1)^2 * var_q1
var_hr_l2_gen <- exp(est_q2)^2 * var_q2
var_hr_l2_donor <- exp(est_q3)^2 * var_q3


hr_l1_lwci_age <- esthr_l1_age - 1.96*sqrt(var_hr_l1_age)
hr_l1_upci_age <- esthr_l1_age + 1.96*sqrt(var_hr_l1_age)

hr_l1_lwci_gen <- esthr_l1_gen - 1.96*sqrt(var_hr_l1_gen)
hr_l1_upci_gen <- esthr_l1_gen + 1.96*sqrt(var_hr_l1_gen)

hr_l1_lwci_donor <- esthr_l1_donor - 1.96*sqrt(var_hr_l1_donor)
hr_l1_upci_donor <- esthr_l1_donor + 1.96*sqrt(var_hr_l1_donor)

hr_l2_lwci_age <- esthr_l2_age - 1.96*sqrt(var_hr_l2_age)
hr_l2_upci_age <- esthr_l2_age + 1.96*sqrt(var_hr_l2_age)

hr_l2_lwci_gen <- esthr_l2_gen - 1.96*sqrt(var_hr_l2_gen)
hr_l2_upci_gen <- esthr_l2_gen + 1.96*sqrt(var_hr_l2_gen)

hr_l2_lwci_donor <- esthr_l2_donor - 1.96*sqrt(var_hr_l2_donor)
hr_l2_upci_donor <- esthr_l2_donor + 1.96*sqrt(var_hr_l2_donor)


# Results --------------------------------------------------------------------
# YW data needed for paper 2: age
hr_gf_age <-c(esthr_l1_age,hr_l1_lwci_age, hr_l1_upci_age)
hr_gf_age
hr_d_age <-c(esthr_l2_age,hr_l2_lwci_age, hr_l2_upci_age)
hr_d_age

# YW data needed for paper 2: gender
hr_gf_gender <-c(esthr_l1_gen,hr_l1_lwci_gen, hr_l1_upci_gen)
hr_gf_gender
hr_d_gender <-c(esthr_l2_gen,hr_l2_lwci_gen, hr_l2_upci_gen)
hr_d_gender

# YW data needed for paper 2: donor
hr_gf_donor <-c(esthr_l1_donor,hr_l1_lwci_donor, hr_l1_upci_donor)
hr_gf_donor

hr_d_donor <-c(esthr_l2_donor,hr_l2_lwci_donor, hr_l2_upci_donor)
hr_d_donor


# YW data needed for paper 2: regression coefficients on association
association_age <- c(est_b1, lwci_b1, upci_b1)
association_gender<- c(est_b2, lwci_b2, upci_b2)
association_donor<- c(est_b3, lwci_b3, upci_b3)
association_age
association_gender
association_donor

results.age <- c(hr_gf_age, hr_d_age, association_age)
results.gender <- c(hr_gf_gender, hr_d_gender, association_gender)
results.donor <- c(hr_gf_donor, hr_d_donor, association_donor)

results <-rbind(results.age, results.gender, results.donor)
results <- data.frame(results)
names(results) <-c("hr_gf", "l_gf", "u_gf", "hr_d", "l_d", "u_d", "theta", "l_theta", "u_theta")
results <- round(results, 3)
results
# Results --------------------------------------------------------------------

##AIC BIC
loglik <- fpl(plfoptim$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(plfoptim$par)
n<-length(X)
aic<- -2*loglik+2*k
bic<- -2*loglik+log(n)*k
loglik
aic
bic


setwd("R/R code for paper 2/bivariate-copula-models-semi-competing-risks")
results$aic = c(round(aic,1), "NA", "NA")
results$run_time= c(round(run_time,2), "NA", "NA")
row.names(results) <- c("age.gl50", "gender.female","donor.living") 
write.csv(results, "results/real_data_analysis/table4_frank_gompertz.csv")

