# Real data analysis
# YW: 15 July 2021  Frank copula Weibull survival models
# YW: 20 July 2021: 1. add running time tracker
#                   2. add constraints for terms in likelihood involving log()
rm(list=ls())
library(copula)
library(mvtnorm)
library(ggplot2)
library(plyr)
library(survival)

########################################################
##################### Load data ########################
########################################################
start_time <- Sys.time()

df <- read.csv(file="NHSBT/paper2_data_2021.csv")
# check descriptive statistics
# donor type
dim(df)
attach(df)

names(df)


########################################################
############### Normal pseudo likelihood ###############
########################################################

npl <- function(para, X, Y, d1, d2, donor, age.grp, gen){
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
    beta1.1 <- beta1[df.1]
    beta2.1 <- beta2[df.1]
    
    S1.1 <- 1-exp(-beta1.1*X.1^alpha1)
    S2.1 <- 1-exp(-beta2.1*Y.1^alpha2)
    
    f1.1 <- beta1.1*alpha1*X.1^(alpha1-1)*exp(-beta1.1*X.1^alpha1) 
    f2.1 <- beta2.1*alpha2*Y.1^(alpha2-1)*exp(-beta2.1*Y.1^alpha2) 
    
    
    rho.1 <- rho[df.1]
    
    # YW added 20 July 2021
    S1.1[which(S1.1<0.1^8)] <- 0.1^8
    S1.1[which(S1.1==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    S2.1[which(S2.1<0.1^8)] <- 0.1^8
    S2.1[which(S2.1==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    f1.1[which(f1.1<0.1^8)] <- 0.1^8
    f2.1[which(f2.1<0.1^8)] <- 0.1^8
    temp_var = 1-rho.1^2
    temp_var[which(temp_var< 0.1^8)] <- 0.1^8
    
    # part1 <- sum(-0.5*log(1-rho.1^2)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
    #                                       rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+ 
    #                log(f1.1)+log(f2.1))
    # testing YW
    #a11= sum(temp_var)
   # a22 = sum(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)- rho.1^2*(qnorm(S1.1)^2 
    #                                                      + qnorm(S2.1)^2)))/(2*(1-rho.1^2)))
    #a22.1 = sum(2*rho.1*qnorm(S1.1)*qnorm(S2.1))
   # a22.2 = sum(rho.1^2*(qnorm(S1.1)^2))
    #a22.3 = sum(qnorm(S2.1)^2)
    #a22.4 = 1-rho.1^2
   # a33 = sum(log(f1.1)+log(f2.1))
    
 
    #print(c(a22.1,a22.2,a22.3))
    
    index1 = which(is.na(qnorm(S1.1))==T|is.infinite(qnorm(S1.1))==T)
    #print("qnorm S1.1")
    # if(length(index1)>0) print(c(S1.1[index1],qnorm(S1.1[index1])))
    
    index2 = which(is.na(qnorm(S2.1))==T|is.infinite(qnorm(S2.1))==T)
    #print("s2.1")
    #if(length(index2)>0) print(qnorm(S2.1[index2]))
    
    
    index3 = which(is.na(rho.1==T|is.infinite(rho.1))==T)
   # print("rho.1")
   # if(length(index3)>0) print(rho.1)
    
    
    part1 <- sum(-0.5*log(temp_var)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
                                         rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+ 
                   log(f1.1)+log(f2.1))
  } else {
    part1 <- 0;
  }
  
  
  #########################################################
  ################### Second Component ####################
  #########################################################
  
  if(sum(df.2)>0){
    
    X.2 <- df[df.2,1]
    Y.2 <- df[df.2,2]
    beta1.2 <- beta1[df.2]
    beta2.2 <- beta2[df.2]
    
    S1.2 <- 1-exp(-beta1.2*X.2^alpha1)
    S2.2 <- 1-exp(-beta2.2*Y.2^alpha2)
    
    f1.2 <- beta1.2*alpha1*X.2^(alpha1-1)*exp(-beta1.2*X.2^alpha1) 
    f2.2 <- beta2.2*alpha2*Y.2^(alpha2-1)*exp(-beta2.2*Y.2^alpha2) 
    
    rho.2 <- rho[df.2]
    
    S1.2[which(S1.2<0.1^8)] <- 0.1^8
    S1.2[which(S1.2==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    S2.2[which(S2.2<0.1^8)] <- 0.1^8
    S2.2[which(S2.2==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    f1.2[which(f1.2<0.1^8)] <- 0.1^8
    f2.2[which(f2.2<0.1^8)] <- 0.1^8
    
    part2.1 <- pnorm(qnorm(S2.2), mean=rho.2*qnorm(S1.2),sd=sqrt(1-rho.2^2), lower.tail=F)
    part2.1[which(part2.1<0.1^(10))]=0.1^(10)
    
    part2 <- sum(log(part2.1*f1.2))
  } else {
    part2 <- 0;
  }
  
  
  #########################################################
  #################### Third Component ####################
  #########################################################
  
  if(sum(df.3) >0 ){
    
    X.3 <- df[df.3,1]
    Y.3 <- df[df.3,2]
    
    beta1.3 <- beta1[df.3]
    beta2.3 <- beta2[df.3]
    
    S1.3 <- 1-exp(-beta1.3*X.3^alpha1)
    S2.3 <- 1-exp(-beta2.3*Y.3^alpha2)
    
    f1.3 <- beta1.3*alpha1*X.3^(alpha1-1)*exp(-beta1.3*X.3^alpha1) 
    f2.3 <- beta2.3*alpha2*Y.3^(alpha2-1)*exp(-beta2.3*Y.3^alpha2) 
    
    rho.3 <- rho[df.3]
    
    S1.3[which(S1.3<0.1^8)] <- 0.1^8
    S1.3[which(S1.3==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    S2.3[which(S2.3<0.1^8)] <- 0.1^8
    S2.3[which(S2.3==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    f1.3[which(f1.3<0.1^8)] <- 0.1^8
    f2.3[which(f2.3<0.1^8)] <- 0.1^8
    
    part3.1 <- pnorm(qnorm(S1.3), mean=rho.3*qnorm(S2.3),sd=sqrt(1-rho.3^2), lower.tail=F)
    part3.1[which(part3.1<0.1^(10))]=0.1^(10)
    
    part3 <- sum(log(part3.1*f2.3))
    
  } else {
    part3 <- 0;
  }
  
  
  #########################################################
  #################### Fourth Component ###################
  #########################################################
  
  if(sum(df.4)>0){
    
    X.4 <- df[df.4,1]
    Y.4 <- df[df.4,2]
    
    beta1.4 <- beta1[df.4]
    beta2.4 <- beta2[df.4]
    
    S1.4 <- 1-exp(-beta1.4*X.4^alpha1)
    S2.4 <- 1-exp(-beta2.4*Y.4^alpha2)
    
    S1.4[which(S1.4<0.1^8)] <- 0.1^8
    S1.4[which(S1.4==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    S2.4[which(S2.4<0.1^8)] <- 0.1^8
    S2.4[which(S2.4==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1

    
    rho.4 <- rho[df.4]
    
    over.all <- rep(0, length(S1.4))
    for(i in 1:length(over.all)){
      sigma <- matrix(c(1,rho.4[i],rho.4[i],1),nrow=2);
      CDF <- function(V,sigma){
        return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
      }
      part4.1 <- (apply(qnorm(cbind(S1.4[i],S2.4[i])),1,CDF,sigma))
      part4.1[which(part4.1<0.1^(10))]=0.1^(10)
      over.all[i] <- log(part4.1);    }
    
    part4 <- sum(over.all);
  } else {
    part4 <- 0;
  }
  #print(part4)
  
  #########################################################
  #################### All Components #####################
  ######################################################### 
  
  loglik <- (part1+part2+part3+part4); 
 # print(c("log likelihood=",loglik));
  print(c("part 1=",part1, "part 2=", part2, "part 3=",part3, "part 4=", part4))
  return(loglik);


}

#npl(c(0.7,  -2.5,0.2,0.01,-0.6,    1,   -4,1.3,-0.08,-0.64,  0.35,0.3,0.02,0.03),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen)
#npl(c(0.2,  -10,-1,-1,-2,         0.5,    -10,-2,-2,-2,       0.01,-1,-0.5,-0.5),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen)
#npl(c(1,      -1,1,1,1,           1.5,     -1,2,2,1            ,0.6,0.6,0.5,0.5),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen)


plnoptim <- optim(c(0.7,  -2.5,0.2,0.01,-0.6,    1,   -4,1.3,-0.08,-0.64,  0.35,0.3,0.02,0.03), npl, method="L-BFGS-B",
           lower=c(0.2,  -10,-1,-1,-2,         0.5,    -10,-2,-2,-2,       0.01,-1,-0.5,-0.5),
          upper=c(1,      -1,1,1,1,           1.5,     -1,2,2,1            ,0.6,0.6,0.5,0.5), 
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen,
                  control=list(fnscale=-1),hessian=TRUE)


plnoptim$par

########################################################
################## Confidence Intervals ################
########################################################

#Fisher's Information matrix
fisher_info<-solve(-plnoptim$hessian) #inverse -hess
#Standard error = sqrt(var/n)
se<-sqrt(diag(fisher_info)) 

#beta ci#
est_b0 <- plnoptim$par[11]
est_b1 <- plnoptim$par[12]
est_b2 <- plnoptim$par[13]
est_b3 <- plnoptim$par[14]
lwci_b0 <- est_b0 - 1.96*se[11]     
lwci_b1 <- est_b1 - 1.96*se[12]
lwci_b2 <- est_b2 - 1.96*se[13]     
lwci_b3 <- est_b3 - 1.96*se[14]
upci_b0 <- est_b0 + 1.96*se[11]
upci_b1 <- est_b1 + 1.96*se[12] 
upci_b2 <- est_b2 + 1.96*se[13]
upci_b3 <- est_b3 + 1.96*se[14] 

#a ci#
est_a1 <- plnoptim$par[1]
est_a2 <- plnoptim$par[6]
lwci_a1 <- est_a1 - 1.96*se[1]
lwci_a2 <- est_a2 - 1.96*se[6]     
upci_a1 <- est_a1 + 1.96*se[1] 
upci_a2 <- est_a2 + 1.96*se[6]

#x ci#
est_x0 <- plnoptim$par[2]
est_x1 <- plnoptim$par[3]
est_x2 <- plnoptim$par[4]
est_x3 <- plnoptim$par[5]
lwci_x0 <- est_x0 - 1.96*se[2]     
lwci_x1 <- est_x1 - 1.96*se[3]
lwci_x2 <- est_x2 - 1.96*se[4]     
lwci_x3 <- est_x3 - 1.96*se[5]
upci_x0 <- est_x0 + 1.96*se[2]
upci_x1 <- est_x1 + 1.96*se[3] 
upci_x2 <- est_x2 + 1.96*se[4]
upci_x3 <- est_x3 + 1.96*se[5] 

#x ci#
est_y0 <- plnoptim$par[7]
est_y1 <- plnoptim$par[8]
est_y2 <- plnoptim$par[9]
est_y3 <- plnoptim$par[10]
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

est_hr_l1_age <- exp(est_x1)
est_hr_l1_gen <- exp(est_x2)
est_hr_l1_donor <- exp(est_x3)

est_hr_l2_age <- exp(est_y1)
est_hr_l2_gen <- exp(est_y2)
est_hr_l2_donor <- exp(est_y3)

var_hr_l1_age <- exp(est_x1)^2 * var_x1
var_hr_l1_gen <- exp(est_x2)^2 * var_x2
var_hr_l1_donor <- exp(est_x3)^2 * var_x3

var_hr_l2_age <- exp(est_y1)^2 * var_y1
var_hr_l2_gen <- exp(est_y2)^2 * var_y2
var_hr_l2_donor <- exp(est_y3)^2 * var_y3


hr_l1_lwci_age <- est_hr_l1_age - 1.96*sqrt(var_hr_l1_age)
hr_l1_upci_age <- est_hr_l1_age + 1.96*sqrt(var_hr_l1_age)

hr_l1_lwci_gen <- est_hr_l1_gen - 1.96*sqrt(var_hr_l1_gen)
hr_l1_upci_gen <- est_hr_l1_gen + 1.96*sqrt(var_hr_l1_gen)

hr_l1_lwci_donor <- est_hr_l1_donor - 1.96*sqrt(var_hr_l1_donor)
hr_l1_upci_donor <- est_hr_l1_donor + 1.96*sqrt(var_hr_l1_donor)

hr_l2_lwci_age <- est_hr_l2_age - 1.96*sqrt(var_hr_l2_age)
hr_l2_upci_age <- est_hr_l2_age + 1.96*sqrt(var_hr_l2_age)

hr_l2_lwci_gen <- est_hr_l2_gen - 1.96*sqrt(var_hr_l2_gen)
hr_l2_upci_gen <- est_hr_l2_gen + 1.96*sqrt(var_hr_l2_gen)

hr_l2_lwci_donor <- esthr_l2_donor - 1.96*sqrt(var_hr_l2_donor)
hr_l2_upci_donor <- esthr_l2_donor + 1.96*sqrt(var_hr_l2_donor)



##AIC BIC
para <- c(est_a1, est_x0,est_x1, est_x2, est_x3, est_a2, est_y0, est_y1, est_y2, est_y3, est_b0, est_b1, est_b2, est_b3)
loglik <- npl(para,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(para)
n<-length(X)
aic<- -2*loglik+2*k
bic<- -2*loglik+log(n)*k
loglik
aic
bic



#### print ###
print(est_b0)
print(lwci_b0)
print(upci_b0)

print(est_b1)
print(lwci_b1)
print(upci_b1)

print(est_b2)
print(lwci_b2)
print(upci_b2)

print(est_b3)
print(lwci_b3)
print(upci_b3)

print(est_a1)
print(lwci_a1)
print(upci_a1)

print(est_a2)
print(lwci_a2)
print(upci_a2)

print(est_x0)
print(lwci_x0)
print(upci_x0)

print(est_x1)
print(lwci_x1)
print(upci_x1)

print(est_x2)
print(lwci_x2)
print(upci_x2)

print(est_x3)
print(lwci_x3)
print(upci_x3)

print(est_y0)
print(lwci_y0)
print(upci_y0)

print(est_y1)
print(lwci_y1)
print(upci_y1)

print(est_y2)
print(lwci_y2)
print(upci_y2)

print(est_y3)
print(lwci_y3)
print(upci_y3)


print(est_hr_l1_age)
print(hr_l1_lwci_age)
print(hr_l1_upci_age)

print(est_hr_l1_gen)
print(hr_l1_lwci_gen)
print(hr_l1_upci_gen)

print(est_hr_l1_donor)
print(hr_l1_lwci_donor)
print(hr_l1_upci_donor)


print(est_hr_l2_age)
print(hr_l2_lwci_age)
print(hr_l2_upci_age)

print(est_hr_l2_gen)
print(hr_l2_lwci_gen)
print(hr_l2_upci_gen)

print(est_hr_l2_donor)
print(hr_l2_lwci_donor)
print(hr_l2_upci_donor)

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
print(bic)

end_time <- Sys.time()
run_time = end_time - start_time
run_time


# > names(results) <-c("hr_gf", "l_gf", "u_gf", "hr_d", "l_d", "u_d", "theta", "l_theta", "u_theta")
# > results <- round(results, 3)
# > results
# hr_gf  l_gf  u_gf  hr_d   l_d   u_d theta l_theta u_theta
# results.age    1.005 0.963 1.048 3.694 3.526 3.862 0.284   0.241   0.326
# results.gender 1.021 0.979 1.064 0.904 0.865 0.943 0.028  -0.015   0.071
# results.donor  0.573 0.544 0.602 0.521 0.490 0.552 0.034  -0.026   0.095
# > aic
# [1] 145038.4