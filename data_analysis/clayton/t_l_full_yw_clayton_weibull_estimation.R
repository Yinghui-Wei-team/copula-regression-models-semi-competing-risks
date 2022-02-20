# Real data analysis
# YW: 9 July 2021  Clayton copula exponential survival models
# original script by LS; edited and updated for paper2 by YW
# Table 5: clayton copula weibull survival model
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

########################################################
############### Clayton pseudo likelihood ##############
########################################################
start_time = Sys.time()
cpl <- function(para, X, Y, d1, d2, donor, age.grp, gen){
  alpha1 <- para[1]
  x0 <- para[2]
  x1 <- para[3]
  x2 <- para[4]
  x3 <- para[5]
  alpha2 <- para[6]
  y0 <- para[7]
  y1 <- para[8]
  y2 <- para[9]
  y3 <- para[10]
  
  b0 <- para[11]
  b1 <- para[12]
  b2 <- para[13]
  b3 <- para[14]
  
  beta1 <- exp(x0+x1*age.grp+x2*gen+x3*donor)
  beta2 <- exp(y0+y1*age.grp+y2*gen+y3*donor)
  
  S1 <- exp(-beta1*X^alpha1)
  S2 <- exp(-beta2*Y^alpha2)
  
  f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
  f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
  
  theta <- exp(b0+b1*age.grp+b2*gen+b3*donor)
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  
  part1 <- d1*d2*(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1*S2)^(1+theta)))
  
  part2 <- d1*(1-d2)*(log(C^(1+theta)*f1)-log(S1^(1+theta)))
  
  part3 <- d2*(1-d1)*(log(C^(1+theta)*f2)-log(S2^(1+theta)))
  
  part4 <- ((1-d1)*(1-d2))*log(C)
  
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}


plcoptim <- optim(c(0.02, -1,-0.01,-0.01,-0.01,  0.02,-1,-0.01,-0.01,-0.01,  2,0.1,0.1,0.1), cpl, method="L-BFGS-B",
                  lower=c(-0.2,-10,-10,-10,-10,  -0.2,-10,-10,-10,-10,  0.01,-1,-1,-1),
                  upper=c(1,-1,1,1,1, 1,-2,2,1,0.01  ,10,6,3,3), 
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen,
                  control=list(fnscale=-1),hessian=TRUE)

end_time = Sys.time()
run_time = end_time - start_time
run_time
plcoptim$par

########################################################
################## Confidence Intervals ################
########################################################

#Fisher's Information matrix
fisher_info<-solve(-plcoptim$hessian) #inverse -hess
#Standard error = sqrt(var/n)
se<-sqrt(diag(fisher_info)) 

#beta ci#
est_b0 <- plcoptim$par[11]
est_b1 <- plcoptim$par[12]
est_b2 <- plcoptim$par[13]
est_b3 <- plcoptim$par[14]
lwci_b0 <- est_b0 - 1.96*se[11]     
lwci_b1 <- est_b1 - 1.96*se[12]
lwci_b2 <- est_b2 - 1.96*se[13]     
lwci_b3 <- est_b3 - 1.96*se[14]
upci_b0 <- est_b0 + 1.96*se[11]
upci_b1 <- est_b1 + 1.96*se[12] 
upci_b2 <- est_b2 + 1.96*se[13]
upci_b3 <- est_b3 + 1.96*se[14] 

#a ci#
est_a1 <- plcoptim$par[1]
est_a2 <- plcoptim$par[6]
lwci_a1 <- est_a1 - 1.96*se[1]
lwci_a2 <- est_a2 - 1.96*se[6]     
upci_a1 <- est_a1 + 1.96*se[1] 
upci_a2 <- est_a2 + 1.96*se[6]

#x ci#
est_x0 <- plcoptim$par[2]
est_x1 <- plcoptim$par[3]
est_x2 <- plcoptim$par[4]
est_x3 <- plcoptim$par[5]
lwci_x0 <- est_x0 - 1.96*se[2]     
lwci_x1 <- est_x1 - 1.96*se[3]
lwci_x2 <- est_x2 - 1.96*se[4]     
lwci_x3 <- est_x3 - 1.96*se[5]
upci_x0 <- est_x0 + 1.96*se[2]
upci_x1 <- est_x1 + 1.96*se[3] 
upci_x2 <- est_x2 + 1.96*se[4]
upci_x3 <- est_x3 + 1.96*se[5] 

#x ci#
est_y0 <- plcoptim$par[7]
est_y1 <- plcoptim$par[8]
est_y2 <- plcoptim$par[9]
est_y3 <- plcoptim$par[10]
lwci_y0 <- est_y0 - 1.96*se[7]     
lwci_y1 <- est_y1 - 1.96*se[8]
lwci_y2 <- est_y2 - 1.96*se[9]     
lwci_y3 <- est_y3 - 1.96*se[10]
upci_y0 <- est_y0 + 1.96*se[7]
upci_y1 <- est_y1 + 1.96*se[8] 
upci_y2 <- est_y2 + 1.96*se[9]
upci_y3 <- est_y3 + 1.96*se[10] 

#hrs
var_x1 <- fisher_info[3,3]
var_x2 <- fisher_info[4,4]
var_x3 <- fisher_info[5,5]
var_y1 <- fisher_info[8,8]
var_y2 <- fisher_info[9,9]
var_y3 <- fisher_info[10,10]

esthr_l1_age <- exp(est_x1)
esthr_l1_gen <- exp(est_x2)
esthr_l1_donor <- exp(est_x3)

esthr_l2_age <- exp(est_y1)
esthr_l2_gen <- exp(est_y2)
esthr_l2_donor <- exp(est_y3)

var_hr_l1_age <- exp(est_x1)^2 * var_x1
var_hr_l1_gen <- exp(est_x2)^2 * var_x2
var_hr_l1_donor <- exp(est_x3)^2 * var_x3

var_hr_l2_age <- exp(est_y1)^2 * var_y1
var_hr_l2_gen <- exp(est_y2)^2 * var_y2
var_hr_l2_donor <- exp(est_y3)^2 * var_y3


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
para <- c(est_a1, est_x0,est_x1, est_x2, est_x3, est_a2, est_y0, est_y1, est_y2, est_y3, est_b0, est_b1, est_b2, est_b3)
loglik <- cpl(para,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(para)
n<-length(X)
aic<- -2*loglik+2*k
bic<- -2*loglik+log(n)*k
loglik
aic
bic

results$aic[1] = round(aic,1)
results$run_time[1]= round(run_time,2)
#results$run_time[c(2,3)]=results$aic[c(2,3)]=""
setwd("R/R code for paper 2/bivariate-copula-models-semi-competing-risks")
# to indicate the level in the estimated results
row.names(results) <- c("age.gl50", "gender.female","donor.living") 
write.csv(results, "results/real_data_analysis/table5_clayton_weibull_coef.csv")


