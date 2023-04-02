# Real data analysis:NHSBT data
# YW: 15 July 2021 updates:
#     1. Gumbel Gompertz survival models
#     2. make unified data set for analysis paper2_data.csv
#     3. Output results into a vector
# original script by LS; edited and updated for paper2 by YW
# YW: 2022-02-20: update load data and output results sections
###############################################################

rm(list=ls())
library(copula)
library(mvtnorm)
library(ggplot2)
library(plyr)
library(survival)

copula <- "gumbel"
survival_distribution <- "gompertz"
if(survival_distribution == "exp") {table_ref = "table3"}
if(survival_distribution == "gompertz") {table_ref = "table4"}
if(survival_distribution == "weibull") {table_ref = "table5"}
########################################################
##################### Load data ########################
########################################################
# YW: need to firstly set working directory to project directory and send through the next two lines

# Load data
# set working directory to project directory and send through the next two lines
dir_data <-dir_results <- "../../"
df <- read.csv(file=paste0(dir_data, "NHSBT/paper2_data_v2.csv"))



dim(df)
attach(df)
names(df)


########################################################
############### Gumbel pseudo likelihood ##############
########################################################

gpl <- function(para, X, Y, d1, d2, age.grp, gen, donor){
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
  
  theta <- exp(b0+b1*age.grp+b2*gen+b3*donor)+1
  # print("b2")
  #print(b2)
  
  #print("b3")
  # print(b3)
  #print("b age group")
  #print(exp(b1*age.grp))
  #print("b gen group")
  #print(exp(b2*gen))
  #print("b donor")
  #print(exp(b3*donor))
  
  lambda1 <- exp(p0+p1*age.grp+p2*gen+p3*donor)
  lambda2 <- exp(q0+q1*age.grp+q2*gen+q3*donor)
  
  S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
  S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
  theta[which(theta > 15)]=15
  
  f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
  f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
  
  f1[which(f1 < 0.1^(8))]=0.1^(8)
  f2[which(f2 < 0.1^(8))]=0.1^(8)
  
  #print("f1")
  if(log(f1)< -100) print(log(f1))
  
  #print("f2")
  if(log(f2)< -100) print(log(f2))
  
  S1[which(S1 < 0.1^(8))]=0.1^(8)
  S2[which(S2 < 0.1^(8))]=0.1^(8)
  
  #print("S1")
  if(log(S1)< -100) print(log(S1))
  
  #print("S2")
  if(log(S2)< -100) print(log(S2))
  
  
  C=exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
  
  C[which(C<0.1^(8))]=0.1^(8)
  #print(theta)
  
  #print("C")
  if(log(C)< -100) print(log(C))
  
  part1 <- d1*d2*(log(C)+(theta-1)*log(-log(S1))+(theta-1)*log(-log(S2))+log(theta-1+((-log(S1))^theta+(-log(S2))^theta)^(1/theta))-log(S1)-log(S2)-(2*theta-1)*log(-log(C))+log(f1)+log(f2))
  part2 <- d1*(1-d2)*(log(C)+(theta-1)*log(-log(S1))-log(S1)-(theta-1)*log(-log(C))+log(f1))
  part3 <-((1-d1)*(d2))*(log(C)+(theta-1)*log(-log(S2))-log(S2)-(theta-1)*log(-log(C))+log(f2))
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  #  print(c(part1,part2,part3,part4))
  logpl <- sum(part1+part2+part3+part4) 
  return(logpl)
}

# gpl(c(0.02,   -3,0.08, 0.02,-0.5,  0.05,  -4.5, 1.5,-0.11,-0.6,  -1.8,1,-0.05, 0.05),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,donor=df$donor, age.grp=df$age.grp, gen=df$gen)
# gpl(c(-0.2,   -6,-0.1, -0.2,-1.5,  -0.1,  -7, 0.5,-0.5,-1,      -3,0.1,-0.1, -0.1),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,donor=df$donor, age.grp=df$age.grp, gen=df$gen)
# gpl(c(0.2,    -2,0.4, 0.2,0.5,     0.1,  -3, 3,0.5,1,          -1,2,0.1, 0.1),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,donor=df$donor, age.grp=df$age.grp, gen=df$gen)



# plgoptim <- optim(c(-0.2, -6,-0.1,       -0.2,-1.5, -0.1,  -7, 0.5,-0.5,-1,      -3,0.1,-0.1, -0.1), gpl, method="L-BFGS-B",
#                   lower=c(-0.2,-6,-0.1,  -0.2,-1.5, -0.1,  -7, 0.5,-0.5,-1,      -3,0.1,-0.1, -0.1),
#                   upper=c(0.2, -2, 0.4,   0.2,0.5,  0.1,   -3, 3,   0.5, 1,       -1,2,0.1, 0.1),
#                   X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,donor=df$donor, age.grp=df$age.grp, gen=df$gen,
#                   control=list(fnscale=-1),hessian=TRUE)

start_time= Sys.time()

plgoptim <- optim(c(1, -1,-0.01,-0.01,-0.01,  1, -1,-0.01,-0.01,-0.01,  2,0.1,0.1,0.1), gpl, method="L-BFGS-B",
                  lower=c(-0.2,-6,-0.1,  -0.2,-1.5, -0.1,  -7, 0.5,-0.5,-1,      -3,0.1,-0.1, -0.1),
                  upper=c(0.2, -2, 0.4,   0.2,0.5,  0.1,   -3, 3,   0.5, 1,       -1,2,0.1, 0.1),
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,donor=df$donor, age.grp=df$age.grp, gen=df$gen,
                  control=list(fnscale=-1),hessian=TRUE)

end_time = Sys.time()
run_time = end_time - start_time
run_time

plgoptim$par


########################################################
################## Confidence Intervals ################
########################################################

#Fisher's Information matrix
fisher_info<-solve(-plgoptim$hessian) #inverse -hess
#Standard error = sqrt(var/n)
se<-sqrt(diag(fisher_info)) 

#beta ci#
est_b0 <- plgoptim$par[11]
est_b1 <- plgoptim$par[12]
est_b2 <- plgoptim$par[13]
est_b3 <- plgoptim$par[14]
lwci_b0 <- est_b0 - 1.96*se[11]     
lwci_b1 <- est_b1 - 1.96*se[12]
lwci_b2 <- est_b2 - 1.96*se[13]     
lwci_b3 <- est_b3 - 1.96*se[14]
upci_b0 <- est_b0 + 1.96*se[11]
upci_b1 <- est_b1 + 1.96*se[12] 
upci_b2 <- est_b2 + 1.96*se[13]
upci_b3 <- est_b3 + 1.96*se[14] 

#g ci#
est_g1 <- plgoptim$par[1]
est_g2 <- plgoptim$par[6]
lwci_g1 <- est_g1 - 1.96*se[1]
lwci_g2 <- est_g2 - 1.96*se[6]     
upci_g1 <- est_g1 + 1.96*se[1] 
upci_g2 <- est_g2 + 1.96*se[6]

#p ci#
est_p0 <- plgoptim$par[2]
est_p1 <- plgoptim$par[3]
est_p2 <- plgoptim$par[4]
est_p3 <- plgoptim$par[5]
lwci_p0 <- est_p0 - 1.96*se[2]     
lwci_p1 <- est_p1 - 1.96*se[3]
lwci_p2 <- est_p2 - 1.96*se[4]     
lwci_p3 <- est_p3 - 1.96*se[5]
upci_p0 <- est_p0 + 1.96*se[2]
upci_p1 <- est_p1 + 1.96*se[3] 
upci_p2 <- est_p2 + 1.96*se[4]
upci_p3 <- est_p3 + 1.96*se[5] 

#q ci#
est_q0 <- plgoptim$par[7]
est_q1 <- plgoptim$par[8]
est_q2 <- plgoptim$par[9]
est_q3 <- plgoptim$par[10]
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
loglik <- gpl(plgoptim$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(plgoptim$par)
n<-length(X)
aic<- -2*loglik+2*k
bic<- -2*loglik+log(n)*k
loglik
aic
bic

print(aic)
print(bic)

results$aic = c(round(aic,1), "NA", "NA")
results$run_time= c(round(run_time,2), "NA", "NA")
row.names(results) <- c("age.gl50", "gender.female","donor.living") 
setwd("R/R code for paper 2/bivariate-copula-models-semi-competing-risks")
write.csv(results, paste0("results/real_data_analysis/", table_ref, "_", copula, "_",survival_distribution, ".csv"))
