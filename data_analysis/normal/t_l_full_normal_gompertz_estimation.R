###############################################################################
# Real data analysis
# YW: 9 July 2021 normal gompertz survival models - regression on both hazards and association parameters
# YW: 14 July 2021: make unified data set for analysis paper2_data.csv
# YW: 17 July 2021: adding constraints for terms in likelihood involving log()
# original script by LS; edited and updated for paper2 by YW
# YW: 2022-02-20: update load data and output results sections
# YW: 2022-12-29: output regression coefficients and add comments
###############################################################################

rm(list=ls())
library(copula);library(mvtnorm);library(plyr)

start.time = Sys.time()

################################################################################
# Set up model specs and load data                                             #
################################################################################
# Model specs
copula <- "normal"
survival_distribution <- "gompertz"
if(survival_distribution == "exp") {table_ref = "table3"}
if(survival_distribution == "gompertz") {table_ref = "table4"}
if(survival_distribution == "weibull") {table_ref = "table5"}

# Load data
# set working directory to project directory and send through the next two lines
dir_data <-dir_results <- "../../"
df <- read.csv(file=paste0(dir_data, "NHSBT/paper2_data.csv"))

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

################################################################################
# Normal pseudo likelihood                                                     #
################################################################################

npl <- function(para, X, Y, d1, d2, donor, age.grp, gen){
  gamma1 <- para[1]  # parameter in Gompertz distribution
  a0 <- para[2]      # regression coefficients for lambda1 (hazard for graft failure)
  a1 <- para[3]
  a2 <- para[4]
  a3 <- para[5]
  gamma2 <- para[6]  # parametr in Gompertz distribution
  c0 <- para[7]      # regression coefficients for lambda2 (hazard for death)
  c1 <- para[8]
  c2 <- para[9]
  c3 <- para[10]
  
  b0 <- para[11]     # regression coefficients for association parameter
  b1 <- para[12]
  b2 <- para[13]
  b3 <- para[14]
  
  lambda1 <- exp(a0+a1*age.grp+a2*gen+a3*donor)  # hazard for graft failure
  lambda2 <- exp(c0+c1*age.grp+c2*gen+c3*donor)  # hazard for death
  
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


plnoptim <- optim(c(0.01,                # gamma1: parameter for Gompertz distribution for graft failure
                    -3,0.2,0.01,-0.5,    # a: regression coefficients for lambda 1 (hazard for graft failure)
                    0.04,                # gamma2: parameter for Gompertz distribution for death
                    -4,1.3,-0.08,-0.58,  # c: regression coefficients for lambda 2 (hazard for death)
                    0.35,0.3,0.02,0.03   # b: regression coefficients for association parameter
                    ), npl, method="L-BFGS-B",
                  lower=c(-0.2,  
                          -10,-1,-1,-2,        
                          -0.2,    
                          -10,-2,-2,-2,       
                          0.01,-1,-0.5,-0.5
                          ),
                  upper=c(0.2,      
                          -1,1,1,1,           
                          0.2,     
                          -1,2,2,1   
                          ,0.6,0.6,0.5,0.5), 
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
est_a0 <- plnoptim$par[2]
est_a1 <- plnoptim$par[3]
est_a2 <- plnoptim$par[4]
est_a3 <- plnoptim$par[5]
lwci_a0 <- est_a0 - 1.96*se[2]     
lwci_a1 <- est_a1 - 1.96*se[3]
lwci_a2 <- est_a2 - 1.96*se[4]     
lwci_a3 <- est_a3 - 1.96*se[5]
upci_a0 <- est_a0 + 1.96*se[2]
upci_a1 <- est_a1 + 1.96*se[3] 
upci_a2 <- est_a2 + 1.96*se[4]
upci_a3 <- est_a3 + 1.96*se[5] 

#q ci#
est_c0 <- plnoptim$par[7]
est_c1 <- plnoptim$par[8]
est_c2 <- plnoptim$par[9]
est_c3 <- plnoptim$par[10]
lwci_c0 <- est_c0 - 1.96*se[7]     
lwci_c1 <- est_c1 - 1.96*se[8]
lwci_c2 <- est_c2 - 1.96*se[9]     
lwci_c3 <- est_c3 - 1.96*se[10]
upci_c0 <- est_c0 + 1.96*se[7]
upci_c1 <- est_c1 + 1.96*se[8] 
upci_c2 <- est_c2 + 1.96*se[9]
upci_c3 <- est_c3 + 1.96*se[10] 

#hrs
var_a1 <- fisher_info[3,3]
var_a2 <- fisher_info[4,4]
var_a3 <- fisher_info[5,5]
var_c1 <- fisher_info[8,8]
var_c2 <- fisher_info[9,9]
var_c3 <- fisher_info[10,10]

est_hr_l1_age <- exp(est_a1)
est_hr_l1_gen <- exp(est_a2)
est_hr_l1_donor <- exp(est_a3)

est_hr_l2_age <- exp(est_c1)
est_hr_l2_gen <- exp(est_c2)
est_hr_l2_donor <- exp(est_c3)

var_hr_l1_age <- exp(est_a1)^2 * var_a1
var_hr_l1_gen <- exp(est_a2)^2 * var_a2
var_hr_l1_donor <- exp(est_a3)^2 * var_a3

var_hr_l2_age <- exp(est_c1)^2 * var_c1
var_hr_l2_gen <- exp(est_c2)^2 * var_c2
var_hr_l2_donor <- exp(est_c3)^2 * var_c3


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

print(est_a0)
print(lwci_a0)
print(upci_a0)

print(est_a1)
print(lwci_a1)
print(upci_a1)

print(est_a2)
print(lwci_a2)
print(upci_a2)

print(est_a3)
print(lwci_a3)
print(upci_a3)

print(est_c0)
print(lwci_c0)
print(upci_c0)

print(est_c1)
print(lwci_c1)
print(upci_c1)

print(est_c2)
print(lwci_c2)
print(upci_c2)

print(est_c3)
print(lwci_c3)
print(upci_c3)


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
run.time

results$aic = c(round(aic,1), "NA", "NA")
results$run_time= c(round(run.time,2), "NA", "NA")
row.names(results) <- c("age.gl50", "gender.female","donor.living") 


print("normal copula gompertz survival models")

################################################################################
# Create a data frame for regression coefficients                              #
################################################################################
# regression coefficients in hazard 1 (lambda1):  est_a0, est_a1, est_a2, est_a3
# regression coefficients in hazard 2 (lambda2):  est_c0, est_c1, est_c2, est_c3
# regression coefficients in association parameter: est_b0, est_b1, est_b2, est_b3, 
# Gompertz parameter: g1 = gamma1, g2 = gamma2, 
reg_coef <- c(est_g1, lwci_g1, upci_g1,
              est_g2, lwci_g2, upci_g2, 
              
              est_a0, lwci_a0, upci_a0, 
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
row.names(reg_coef) <-c("gamma1","gamma2",      # parameters in Gompertz distributions
                        "a0", "a1","a2","a3",   # regression coefficients for hazard 1 (graft failure)
                        "c0","c1","c2","c3",    # regression coefficients for hazard 2 (death)
                        "b0","b1","b2","b3"     # regression coefficients for association parameter
                        )

################################################################################
# Output results                                                               #
################################################################################
write.csv(reg_coef, paste0(dir_results, "parameters_",copula, "_", survival_distribution,".csv"))
write.csv(results, paste0(dir_results, table_ref, "_", copula, "_",survival_distribution, ".csv"))
