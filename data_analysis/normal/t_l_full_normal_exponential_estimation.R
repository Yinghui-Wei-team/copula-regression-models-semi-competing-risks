# 11-July-2021
# YW: NHSBT data analysis, Clayton copula exponential survival distribution
# original script by LS; edited and updated for paper2 by YW
# YW: 2022-02-20: update load data and output results sections
###############################################################################

rm(list=ls())
library(copula)
library(mvtnorm)
library(ggplot2)
library(plyr)
library(survival)

copula <- "normal"
survival_distribution <- "exp"
if(survival_distribution == "exp") {table_ref = "table3"}
if(survival_distribution == "gompertz") {table_ref = "table4"}
if(survival_distribution == "weibull") {table_ref = "table5"}

########################################################
##################### Load data ########################
########################################################

# YW: need to firstly set working directory to project directory and send through the next two lines
setwd("../../../")
df <- read.csv(file="NHSBT/paper2_data.csv")
dim(df)
attach(df)
names(df)

# check descriptive statistics
# donor type
dim(df)
table(df$age.grp)
table(df$age.grp)/dim(df)[1]
table(df$gen)
table(df$gen)/dim(df)[1]
table(df$donor)
table(df$donor)/dim(df)[1]


start.time = Sys.time()


#########################################################
############### Normal pseudo likelihood ################
#########################################################
npl<-function(para, X, Y, d1, d2, age.grp, gen, donor){
  a0 <- para[1]
  a1 <- para[2]
  a2 <- para[3]
  a3 <- para[4]
  
  c0 <- para[5]
  c1 <- para[6]
  c2 <- para[7]
  c3 <- para[8]
  
  b0 <- para[9]
  b1 <- para[10]
  b2 <- para[11]
  b3 <- para[12]
  
  lambda1 <- exp(a0+a1*age.grp+a2*gen+a3*donor)
  lambda2 <- exp(c0+c1*age.grp+c2*gen+c3*donor)
  rho <- (exp(2*(b0+b1*age.grp+b2*gen+b3*donor))-1)/(exp(2*(b0+b1*age.grp+b2*gen+b3*donor))+1)
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2;   #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  
  #########################################################
  #################### First Component ####################
  #########################################################
  
  if(sum(df.1)>0){
    
    X.1 <- df[df.1,1]
    Y.1 <- df[df.1,2]
    lambda1.1 <- lambda1[df.1]
    lambda2.1 <- lambda2[df.1]
    S1.1 <- pexp(X.1, rate=lambda1.1) 
    S2.1 <- pexp(Y.1, rate=lambda2.1)
    
    rho.1 <- rho[df.1]
    
    
    part1 <- sum(-0.5*log(1-rho.1^2)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
                                          rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+
                   log(lambda1.1)-lambda1.1*X.1 +log(lambda2.1)-lambda2.1*Y.1)
  } else {
    part1 <- 0;
  }
  
  #########################################################
  ################### Second Component ####################
  #########################################################
  
  if(sum(df.2)>0){
    
    X.2 <- df[df.2,1]
    Y.2 <- df[df.2,2]
    
    lambda1.2 <- lambda1[df.2]
    lambda2.2 <- lambda2[df.2]
    
    S1.2 <- pexp(X.2, rate=lambda1.2) 
    S2.2 <- pexp(Y.2, rate=lambda2.2) 
    
    rho.2 <- rho[df.2]
    
    part2 <- sum(log(pnorm(qnorm(S2.2), mean=rho.2*qnorm(S1.2),
                           sd=sqrt(1-rho.2^2), lower.tail=F)*lambda1.2*exp(-lambda1.2*X.2)))
  } else {
    part2 <- 0;
  }
  
  
  #########################################################
  #################### Third Component ####################
  #########################################################
  
  if(sum(df.3) >0 ){
    
    X.3 <- df[df.3,1]
    Y.3 <- df[df.3,2]
    
    lambda1.3 <- lambda1[df.3]
    lambda2.3 <- lambda2[df.3]
    
    S1.3 <- pexp(X.3, rate=lambda1.3) 
    S2.3 <- pexp(Y.3, rate=lambda2.3) 
    
    rho.3 <- rho[df.3]
    
    part3 <- sum(log(pnorm(qnorm(S1.3), mean=rho.3*qnorm(S2.3),
                           sd=sqrt(1-rho.3^2), lower.tail=F)*lambda2.3*exp(-lambda2.3*Y.3)))
  } else {
    part3 <- 0;
  }
  
  
  #########################################################
  #################### Fourth Component ###################
  #########################################################
  
  if(sum(df.4)>0){
    
    X.4 <- df[df.4,1]
    Y.4 <- df[df.4,2]
    
    lambda1.4 <- lambda1[df.4]
    lambda2.4 <- lambda2[df.4]
    
    S1.4 <- pexp(X.4, rate=lambda1.4) 
    S2.4 <- pexp(Y.4, rate=lambda2.4) 
    
    rho.4 <- rho[df.4]
    
    over.all <- rep(0, length(S1.4))
    for(i in 1:length(over.all)){
      sigma <- matrix(c(1,rho.4[i],rho.4[i],1),nrow=2);
      CDF <- function(V,sigma){
        return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
      }
      over.all[i] <- log(apply(qnorm(cbind(S1.4[i],S2.4[i])),1,CDF,sigma));
    }
    
    part4 <- sum(over.all);
  } else {
    part4 <- 0;
  }
  #print(part4)
  
  #########################################################
  #################### All Components #####################
  ######################################################### 
  
  loglik <- (part1+part2+part3+part4); 
  return(loglik);
}

#tic("n")
# plnoptim <- optim(c(-1,-0.01,-0.01,-0.01,  -1,-0.01,-0.01,-0.01,  0.1,0.1,0.1,0.1), npl, method="L-BFGS-B",
#                   lower=c(-10,-10,-10,-10,  -10,-10,-10,-10,  0.01,-1,-1,-1),
#                   upper=c(-1,1,0.1,0.01, -2,2,0.01,0.01  ,1,1,1,1),
#                   X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen,
#                   control=list(fnscale=-1),hessian=TRUE)
start_time= Sys.time()
plnoptim <- optim(c(-2.74,0.32,0.004,-0.535,  -3.41,1.35,-0.06,-0.61,  0.1,1.03,0.1,0.40), npl, method="L-BFGS-B",
                  lower=c(-10,-10,-10,-10,  -10,-10,-10,-10,  0.01,-1,-1,-1),
                  upper=c(-1,1,0.1,0.01, -2,2,0.01,0.01  ,1,1,1,1),
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen,
                  control=list(fnscale=-1),hessian=TRUE)

end_time = Sys.time()
run_time = end_time - start_time
print(run_time)

plnoptim$par
#toc()

#npl(c(-1,-0.01,-0.01,-0.01,  -1,-0.01,-0.01,-0.01,  2,0.1,0.1,0.1), X=df$X, Y=df$Y, 
 #   d1=df$d1, d2=df$d2,age.grp=df$age.grp,donor=df$donor, gen=df$gen)

########################################################
################## Confidence Intervals ################
########################################################

#Fisher's Information matrix
fisher_info<-solve(-plnoptim$hessian) #inverse -hess
#Standard error = sqrt(var/n)
se<-sqrt(diag(fisher_info)) 

#beta ci#
est_b0 <- plnoptim$par[9]
est_b1 <- plnoptim$par[10]
est_b2 <- plnoptim$par[11]
est_b3 <- plnoptim$par[12]
lwci_b0 <- est_b0 - 1.96*se[9]     
lwci_b1 <- est_b1 - 1.96*se[10]
lwci_b2 <- est_b2 - 1.96*se[11]     
lwci_b3 <- est_b3 - 1.96*se[12]
upci_b0 <- est_b0 + 1.96*se[9]
upci_b1 <- est_b1 + 1.96*se[10] 
upci_b2 <- est_b2 + 1.96*se[11]
upci_b3 <- est_b3 + 1.96*se[12] 

#a ci#
est_a0 <- plnoptim$par[1]
est_a1 <- plnoptim$par[2]
est_a2 <- plnoptim$par[3]
est_a3 <- plnoptim$par[4]
lwci_a0 <- est_a0 - 1.96*se[1]     
lwci_a1 <- est_a1 - 1.96*se[2]
lwci_a2 <- est_a2 - 1.96*se[3]     
lwci_a3 <- est_a3 - 1.96*se[4]
upci_a0 <- est_a0 + 1.96*se[1]
upci_a1 <- est_a1 + 1.96*se[2] 
upci_a2 <- est_a2 + 1.96*se[3]
upci_a3 <- est_a3 + 1.96*se[4] 

#c ci#
est_c0 <- plnoptim$par[5]
est_c1 <- plnoptim$par[6]
est_c2 <- plnoptim$par[7]
est_c3 <- plnoptim$par[8]
lwci_c0 <- est_c0 - 1.96*se[5]     
lwci_c1 <- est_c1 - 1.96*se[6]
lwci_c2 <- est_c2 - 1.96*se[7]     
lwci_c3 <- est_c3 - 1.96*se[8]
upci_c0 <- est_c0 + 1.96*se[5]
upci_c1 <- est_c1 + 1.96*se[6] 
upci_c2 <- est_c2 + 1.96*se[7]
upci_c3 <- est_c3 + 1.96*se[8] 

#hrs
var_a1 <- fisher_info[2,2]
var_a2 <- fisher_info[3,3]
var_a3 <- fisher_info[4,4]
var_c1 <- fisher_info[6,6]
var_c2 <- fisher_info[7,7]
var_c3 <- fisher_info[8,8]

esthr_l1_age <- exp(est_a1)
esthr_l1_gen <- exp(est_a2)
esthr_l1_donor <- exp(est_a3)

esthr_l2_age <- exp(est_c1)
esthr_l2_gen <- exp(est_c2)
esthr_l2_donor <- exp(est_c3)

var_hr_l1_age <- exp(est_a1)^2 * var_a1
var_hr_l1_gen <- exp(est_a2)^2 * var_a2
var_hr_l1_donor <- exp(est_a3)^2 * var_a3

var_hr_l2_age <- exp(est_c1)^2 * var_c1
var_hr_l2_gen <- exp(est_c2)^2 * var_c2
var_hr_l2_donor <- exp(est_c3)^2 * var_c3


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




##AIC BIC
para <- c(est_a0, est_a1, est_a2, est_a3, est_c0, est_c1, est_c2, est_c3, est_b0, est_b1, est_b2, est_b3)
loglik <- npl(para,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(para)
n<-length(X)
aic<- -2*loglik+2*k
bic<- -2*loglik+log(n)*k
loglik
aic
bic

print(est_a0)
print(lwci_a0)
print(upci_a0)

print(est_a1)
print(lwci_a1)
print(upci_a1)

print(est_c0)
print(lwci_c0)
print(upci_c0)

print(est_c1)
print(lwci_c1)
print(upci_c1)

print(est_b0)
print(lwci_b0)
print(upci_b0)

print(est_b1)
print(lwci_b1)
print(upci_b1)

print(esthr_l1_age)
print(hr_l1_lwci_age)
print(hr_l1_upci_age)

print(esthr_l1_gen)
print(hr_l1_lwci_gen)
print(hr_l1_upci_gen)

print(esthr_l1_donor)
print(hr_l1_lwci_donor)
print(hr_l1_upci_donor)

print(esthr_l2_age)
print(hr_l2_lwci_age)
print(hr_l2_upci_age)

print(esthr_l2_gen)
print(hr_l2_lwci_gen)
print(hr_l2_upci_gen)

print(esthr_l2_donor)
print(hr_l2_lwci_donor)
print(hr_l2_upci_donor)

print(aic)
print(bic)

end.time = Sys.time()

run.time = end.time - start.time
run.time

# # gender
# hr_gf_gender <-c(esthr_l1_gen,hr_l1_lwci_gen, hr_l1_upci_gen)
# hr_gf_gender
# hr_d_gender <-c(esthr_l2_gen,hr_l2_lwci_gen, hr_l2_upci_gen)
# hr_d_gender
# 
# 
# # age
# hr_gf_age <-c(esthr_l1_age,hr_l1_lwci_age, hr_l1_upci_age)
# hr_gf_age
# 
# hr_d_age <-c(esthr_l2_age,hr_l2_lwci_age, hr_l2_upci_age)
# hr_d_age
# 
# # donor
# hr_gf_donor <-c(esthr_l1_donor,hr_l1_lwci_donor, hr_l1_upci_donor)
# hr_gf_donor
# 
# hr_d_donor <-c(esthr_l2_donor,hr_l2_lwci_donor, hr_l2_upci_donor)
# hr_d_donor

# Results --------------------------------------------------------------------
#YW data needed for paper 2: age
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



print(aic)


reg_coef <- c(est_a0, lwci_a0, upci_a0, 
              est_a1, lwci_a1, upci_a1, 
              est_a2, lwci_a2, upci_a2, 
              est_a3, lwci_a3, upci_a3,
              est_c0, lwci_c0, upci_c0, 
              est_c1, lwci_c1, upci_c1, 
              est_c2, lwci_c2, upci_c2, 
              est_c3, lwci_c3, upci_c3,
              est_b0, lwci_b0, upci_b0, 
              est_b1, lwci_b1, upci_b1, 
              est_b2, lwci_b2, upci_b2, 
              est_b3, lwci_b3, upci_b3
)

reg_coef <- matrix(reg_coef, ncol=3, byrow=T)
reg_coef <- round(reg_coef, 3)

reg_coef

data.frame(reg_coef, row.names=NULL)

reg_coef

getwd()

results$aic = c(round(aic,1), "NA", "NA")
results$run_time= c(round(run_time,2), "NA", "NA")
row.names(results) <- c("age.gl50", "gender.female","donor.living") 
setwd("R/R code for paper 2/bivariate-copula-models-semi-competing-risks")
write(reg_coef, file=paste0("results/real_data_analysis/parameters_",copula, "_", survival_distribution))
write.csv(results, paste0("results/real_data_analysis/", table_ref, "_", copula, "_",survival_distribution, ".csv"))
