# Real data analysis
# YW: 9 July 2021 normal gompertz survival models - regression on both hazards and association parameters
# YW: 14 July 2021: make unified data set for analysis paper2_data.csv
# YW: 17 July 2021: adding constraints for terms in likelihood involving log()
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
survival_distribution <- "gompertz"
if(survival_distribution == "exp") {table_ref = "table3"}
if(survival_distribution == "gompertz") {table_ref = "table4"}
if(survival_distribution == "weibull") {table_ref = "table5"}

start.time = Sys.time()


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




########################################################
############### Normal pseudo likelihood ###############
########################################################

npl <- function(para, X, Y, d1, d2, donor, age.grp, gen){
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
    
  
    # Gompertz survival functions
    S1.1 <- 1-exp(-lambda1.1/gamma1*(exp(gamma1*X.1)-1))
    S2.1 <- 1-exp(-lambda2.1/gamma2*(exp(gamma2*Y.1)-1))

    # Gompertz probability density functions
    f1.1 <- lambda1.1*exp(gamma1*X.1-lambda1.1/gamma1*(exp(gamma1*X.1)-1))
    f2.1 <- lambda2.1*exp(gamma2*Y.1-lambda2.1/gamma2*(exp(gamma2*Y.1)-1))
    
    rho.1 <- rho[df.1]
    
    
    # YW added
    S1.1[which(S1.1<0.1^8)] <- 0.1^8
    S1.1[which(S1.1==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    S2.1[which(S2.1<0.1^8)] <- 0.1^8
    S2.1[which(S2.1==1)] <- 0.999999 # YW added 20 July 2021: avoid qnorm(p)=inf if p=1
    f1.1[which(f1.1<0.1^8)] <- 0.1^8
    f2.1[which(f2.1<0.1^8)] <- 0.1^8
    temp_var = 1-rho.1^2
    temp_var[which(temp_var< 0.1^8)] <- 0.1^8
    
    # part1 <- sum(-0.5*log(1-rho.1^2)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
    #               rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+ 
    #                log(f1.1)+log(f2.1))
    # 
    
    # testing YW
    a11= sum(temp_var)
    a22 = sum(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)- rho.1^2*(qnorm(S1.1)^2 
                  + qnorm(S2.1)^2)))/(2*(1-rho.1^2)))
    a22.1 = sum(2*rho.1*qnorm(S1.1)*qnorm(S2.1))
    a22.2 = sum(rho.1^2*(qnorm(S1.1)^2))
    a22.3 = sum(qnorm(S2.1)^2)
    #a22.4 = 1-rho.1^2
    a33 = sum(log(f1.1)+log(f2.1))
    
    #print(c(a11,a22,a33))
    print(c(a22.1,a22.2,a22.3))
   # print("sum qnorm S1.1")
    #print(sum(qnorm(S1.1)))
   
    #print("sum qnorm S2.1")
   # print(sum(qnorm(S2.1)))
    
   # print("sum rho.1^2")
   # print(sum(rho.1^2))
    
    index1 = which(is.na(qnorm(S1.1))==T|is.infinite(qnorm(S1.1))==T)
   # print(qnorm(S1.1))
   # print(which(is.infinite(qnorm(S1.1))==T))
    #print("index1")
    print("qnorm S1.1")
    if(length(index1)>0) print(c(S1.1[index1],qnorm(S1.1[index1])))
    
    index2 = which(is.na(qnorm(S2.1))==T|is.infinite(qnorm(S2.1))==T)
    #print("s2.1")
    if(length(index2)>0) print(qnorm(S2.1[index2]))
    
  
    index3 = which(is.na(rho.1==T|is.infinite(rho.1))==T)
    print("rho.1")
    if(length(index3)>0) print(rho.1)
    
    
    part1 <- sum(-0.5*log(temp_var)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
                  rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+ 
                  log(f1.1)+log(f2.1))
    
    #print(c(rho.1, f1.1, f2.1))
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
    
    
    S1.2 <- 1-exp(-lambda1.2/gamma1*(exp(gamma1*X.2)-1))
    S2.2 <- 1-exp(-lambda2.2/gamma2*(exp(gamma2*Y.2)-1))
    
    f1.2 <- lambda1.2*exp(gamma1*X.2-lambda1.2/gamma1*(exp(gamma1*X.2)-1))
    f2.2 <- lambda2.2*exp(gamma2*Y.2-lambda2.2/gamma2*(exp(gamma2*Y.2)-1))
    
    rho.2 <- rho[df.2]
    
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
    
    lambda1.3 <- lambda1[df.3]
    lambda2.3 <- lambda2[df.3]
    
    S1.3 <- 1-exp(-lambda1.3/gamma1*(exp(gamma1*X.3)-1))
    S2.3 <- 1-exp(-lambda2.3/gamma2*(exp(gamma2*Y.3)-1))
    
    f1.3 <- lambda1.3*exp(gamma1*X.3-lambda1.3/gamma1*(exp(gamma1*X.3)-1))
    f2.3 <- lambda2.3*exp(gamma2*Y.3-lambda2.3/gamma2*(exp(gamma2*Y.3)-1))
    
    rho.3 <- rho[df.3]
    
    part3.1 <- pnorm(qnorm(S1.3), mean=rho.3*qnorm(S2.3),sd=sqrt(1-rho.3^2), lower.tail=F)
    part3.1[which(part3.1<0.1^(10))]=0.1^(10)
    
    # YW added
    f2.3[which(f2.3<0.1^8)] <- 0.1^8
    
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
    
    lambda1.4 <- lambda1[df.4]
    lambda2.4 <- lambda2[df.4]
    
    S1.4 <- 1-exp(-lambda1.4/gamma1*(exp(gamma1*X.4)-1))
    S2.4 <- 1-exp(-lambda2.4/gamma2*(exp(gamma2*Y.4)-1))
    
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
  #print(c(part1, part2, part3, part4))
  return(loglik);
}

#npl(c(),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen)
#npl(c(),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen)
#npl(c(),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen)


plnoptim <- optim(c(0.01,  -3,0.2,0.01,-0.5,    0.04,   -4,1.3,-0.08,-0.58,  0.35,0.3,0.02,0.03), npl, method="L-BFGS-B",
                  lower=c(-0.2,  -10,-1,-1,-2,         -0.2,    -10,-2,-2,-2,       0.01,-1,-0.5,-0.5),
                  upper=c(0.2,      -1,1,1,1,           0.2,     -1,2,2,1            ,0.6,0.6,0.5,0.5), 
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

#g ci#
est_g1 <- plnoptim$par[1]
est_g2 <- plnoptim$par[6]
lwci_g1 <- est_g1 - 1.96*se[1]
lwci_g2 <- est_g2 - 1.96*se[6]     
upci_g1 <- est_g1 + 1.96*se[1] 
upci_g2 <- est_g2 + 1.96*se[6]

#p ci#
est_p0 <- plnoptim$par[2]
est_p1 <- plnoptim$par[3]
est_p2 <- plnoptim$par[4]
est_p3 <- plnoptim$par[5]
lwci_p0 <- est_p0 - 1.96*se[2]     
lwci_p1 <- est_p1 - 1.96*se[3]
lwci_p2 <- est_p2 - 1.96*se[4]     
lwci_p3 <- est_p3 - 1.96*se[5]
upci_p0 <- est_p0 + 1.96*se[2]
upci_p1 <- est_p1 + 1.96*se[3] 
upci_p2 <- est_p2 + 1.96*se[4]
upci_p3 <- est_p3 + 1.96*se[5] 

#q ci#
est_q0 <- plnoptim$par[7]
est_q1 <- plnoptim$par[8]
est_q2 <- plnoptim$par[9]
est_q3 <- plnoptim$par[10]
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

est_hr_l1_age <- exp(est_p1)
est_hr_l1_gen <- exp(est_p2)
est_hr_l1_donor <- exp(est_p3)

est_hr_l2_age <- exp(est_q1)
est_hr_l2_gen <- exp(est_q2)
est_hr_l2_donor <- exp(est_q3)

var_hr_l1_age <- exp(est_p1)^2 * var_p1
var_hr_l1_gen <- exp(est_p2)^2 * var_p2
var_hr_l1_donor <- exp(est_p3)^2 * var_p3

var_hr_l2_age <- exp(est_q1)^2 * var_q1
var_hr_l2_gen <- exp(est_q2)^2 * var_q2
var_hr_l2_donor <- exp(est_q3)^2 * var_q3


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

hr_l2_lwci_donor <- est_hr_l2_donor - 1.96*sqrt(var_hr_l2_donor)
hr_l2_upci_donor <- est_hr_l2_donor + 1.96*sqrt(var_hr_l2_donor)



##AIC BIC
loglik <- npl(plnoptim$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(plnoptim$par)
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

print(est_g1)
print(lwci_g1)
print(upci_g1)

print(est_g2)
print(lwci_g2)
print(upci_g2)

print(est_p0)
print(lwci_p0)
print(upci_p0)

print(est_p1)
print(lwci_p1)
print(upci_p1)

print(est_p2)
print(lwci_p2)
print(upci_p2)

print(est_p3)
print(lwci_p3)
print(upci_p3)

print(est_q0)
print(lwci_q0)
print(upci_q0)

print(est_q1)
print(lwci_q1)
print(upci_q1)

print(est_q2)
print(lwci_q2)
print(upci_q2)

print(est_q3)
print(lwci_q3)
print(upci_q3)


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

# Results --------------------------------------------------------------------
# YW data needed for paper 2: age
hr_gf_age <-c(est_hr_l1_age,hr_l1_lwci_age, hr_l1_upci_age)
hr_gf_age
hr_d_age <-c(est_hr_l2_age,hr_l2_lwci_age, hr_l2_upci_age)
hr_d_age

# YW data needed for paper 2: gender
hr_gf_gender <-c(est_hr_l1_gen,hr_l1_lwci_gen, hr_l1_upci_gen)
hr_gf_gender
hr_d_gender <-c(est_hr_l2_gen,hr_l2_lwci_gen, hr_l2_upci_gen)
hr_d_gender

# YW data needed for paper 2: donor
hr_gf_donor <-c(est_hr_l1_donor,hr_l1_lwci_donor, hr_l1_upci_donor)
hr_gf_donor

hr_d_donor <-c(est_hr_l2_donor,hr_l2_lwci_donor, hr_l2_upci_donor)
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
print(bic)

end.time = Sys.time()

run.time = end.time - start.time

print("normal copula gompertz survival models")
run.time

results$aic = c(round(aic,1), "NA", "NA")
results$run_time= c(round(run_time,2), "NA", "NA")
row.names(results) <- c("age.gl50", "gender.female","donor.living") 
setwd("R/R code for paper 2/bivariate-copula-models-semi-competing-risks")
write.csv(results, paste0("results/real_data_analysis/", table_ref, "_", copula, "_",survival_distribution, ".csv"))
