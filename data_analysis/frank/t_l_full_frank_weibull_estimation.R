################################################################################
# Real data analysis:NHSBT data
# YW: 15 July 2021  Frank copula Weibull survival models
# 1. use unified data
# 2. output results
# 3. put exp on the transformation below: theta <- exp(b0+b1*age.grp+b2*gen+b3*donor)
# original script by LS; edited and updated for paper2 by YW
# YW: 2022-02-20: update load data and output results sections
################################################################################

rm(list=ls())
library(copula); library(mvtnorm); library(plyr)

################################################################################
# Set up model specs and load data                                             #
################################################################################
# Model specs
copula <- "frank"
survival_distribution <- "weibull"
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
table(d1)
table(d2)
table(d1,d2)

################################################################################
#Frank pseudo likelihood                                                      #
################################################################################
start_time = Sys.time()
fpl <- function(para, X, Y, d1, d2, donor, age.grp, gen){
  alpha1 <- para[1]     # shape parameter in weibull survival distribution for graft failure
  a0 <- para[2]         # regression parameters in lambda1 (hazard for graft failure)
  a1 <- para[3]
  a2 <- para[4]
  a3 <- para[5]
  alpha2 <- para[6]     # shape parameter in weibull survival distribution for death
  c0 <- para[7]         # regression parameters in lambda2 (hazard for death)
  c1 <- para[8]
  c2 <- para[9]
  c3 <- para[10]        
  b0 <- para[11]        # regression parameters for association parameter
  b1 <- para[12]
  b2 <- para[13]
  b3 <- para[14]
  
  beta1 <- exp(a0+a1*age.grp+a2*gen+a3*donor)   # scale parameter in weibull survival distribution for graft failure
  beta2 <- exp(c0+c1*age.grp+c2*gen+c3*donor)   # scale parameter in weibull survival distribution for death
  
  S1 <- exp(-beta1*X^alpha1)          # survival function for graft failure
  S2 <- exp(-beta2*Y^alpha2)          # survival function for death
  S1[which(S1 < 0.1^8)] <- 0.1^8 
  S2[which(S2 < 0.1^8)] <- 0.1^8 
  f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1)   # pdf for graft failure
  f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2)   # pdf for death
  f1[which(f1 < 0.1^8)] <- 0.1^8      
  f2[which(f2 < 0.1^8)] <- 0.1^8 

  theta <- b0+b1*age.grp+b2*gen+b3*donor
  #print(theta)
  #print(log(theta))
  theta[which(theta > 15)]=15

  #print(c(a0,a1,a2))
  
  C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
  C[which(C < 0.1^8)] <- 0.1^8 
  
  #part1 <- d1*d2*(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*S1)-1)-log(exp(theta*S2)-1)+log(f1)+log(f2))
  #rewrite part1 to allow theta to take negative values
  part1 <- d1*d2*(log(theta*(exp(theta*C)-1)/(exp(theta*S1)-1)/(exp(theta*S2)-1))+theta*C+log(f1)+log(f2))
  
  part2 <- d1*(1-d2)*log(((1-exp(theta*C))*f1)/(1-exp(theta*S1)))
  part3 <- (1-d1)*d2*log(((1-exp(theta*C))*f2)/(1-exp(theta*S2)))
  part4<-((1-d1)*(1-d2))*log(C)
  #print(c(part1,part2,part3,part4))
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

plfoptim <- optim(c(0.02,                   # alpha1
                    -1,-0.01,-0.01,-0.01,   # a: regression coefficients on lambda 1 (graft failure)
                    0.02,                   # alpha2
                    -1,-0.01,-0.01,-0.01,   # c: regression coefficients on lambda 2 (death)
                    2,0.1,0.1,0.1           # b: regression coefficients on association parameter
                    ), fpl, method="L-BFGS-B",
                  lower=c(-0.2,             # alpha1
                          -10,-10,-10,-10,  # a
                          -0.2,             # alpha2
                          -10,-10,-10,-10,  # c
                          1,-0.5,-0.4,0     # b
                          ),
                  upper=c(1,                # alpha1
                          -1,1,1,1,         # a
                          1,                # alpha2
                          -2,2,1,0.01,      # c
                          10,6,3,3          # b
                          ), 
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen,
                  control=list(fnscale=-1),hessian=TRUE)

end_time = Sys.time()
run_time = end_time - start_time
run_time

plfoptim$par

################################################################################
# Confidence Intervals                                                         #
################################################################################
#Fisher's Information matrix
fisher_info<-solve(-plfoptim$hessian) #inverse -hess
#Standard error = sqrt(var/n)
se<-sqrt(diag(fisher_info)) 

#beta ci: regression parameters for association parameter#
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

#alpha ci: shape parameters in Weibull distributions#
est_alp1 <- plfoptim$par[1]
est_alp2 <- plfoptim$par[6]
lwci_alp1 <- est_alp1 - 1.96*se[1]
lwci_alp2 <- est_alp2 - 1.96*se[6]     
upci_alp1 <- est_alp1 + 1.96*se[1] 
upci_alp2 <- est_alp2 + 1.96*se[6]

#a ci: regression parameters for hazard 1 (graft failure)#
est_a0 <- plfoptim$par[2]
est_a1 <- plfoptim$par[3]
est_a2 <- plfoptim$par[4]
est_a3 <- plfoptim$par[5]
lwci_a0 <- est_a0 - 1.96*se[2]     
lwci_a1 <- est_a1 - 1.96*se[3]
lwci_a2 <- est_a2 - 1.96*se[4]     
lwci_a3 <- est_a3 - 1.96*se[5]
upci_a0 <- est_a0 + 1.96*se[2]
upci_a1 <- est_a1 + 1.96*se[3] 
upci_a2 <- est_a2 + 1.96*se[4]
upci_a3 <- est_a3 + 1.96*se[5] 

#c ci: regression parameters for hazard 2 (death)#
est_c0 <- plfoptim$par[7]
est_c1 <- plfoptim$par[8]
est_c2 <- plfoptim$par[9]
est_c3 <- plfoptim$par[10]
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

# Results --------------------------------------------------------------------
# age
hr_gf_age <-c(esthr_l1_age,hr_l1_lwci_age, hr_l1_upci_age)
hr_gf_age
hr_d_age <-c(esthr_l2_age,hr_l2_lwci_age, hr_l2_upci_age)
hr_d_age

# gender
hr_gf_gender <-c(esthr_l1_gen,hr_l1_lwci_gen, hr_l1_upci_gen)
hr_gf_gender
hr_d_gender <-c(esthr_l2_gen,hr_l2_lwci_gen, hr_l2_upci_gen)
hr_d_gender

# donor
hr_gf_donor <-c(esthr_l1_donor,hr_l1_lwci_donor, hr_l1_upci_donor)
hr_gf_donor

hr_d_donor <-c(esthr_l2_donor,hr_l2_lwci_donor, hr_l2_upci_donor)
hr_d_donor

# regression coefficients on association
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
para <- c(est_alp1, 
          est_a0,est_a1, est_a2, est_a3, 
          est_alp2, 
          est_c0, est_c1, est_c2, est_c3, 
          est_b0, est_b1, est_b2, est_b3)
loglik <- fpl(para,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(para)
n<-length(df$X)
aic<- -2*loglik+2*k
bic<- -2*loglik+log(n)*k
loglik
aic
bic
# to indicate the level in the estimated results
results$aic = c(round(aic,1), "NA", "NA")
results$run_time= c(round(run_time,2), "NA", "NA")
row.names(results) <- c("age.gl50", "gender.female","donor.living") 

################################################################################
# Create a data frame for regression coefficients                              #
################################################################################
#beta ci#
# regression coefficients in hazard 1:  est_a0, est_a1, est_a2, est_a3
# regression coefficients in hazard 2:  est_c0, est_c1, est_c2, est_c3
# regression coefficients in association parameter: est_b0, est_b1, est_b2, est_b3, 
# Weibull parameter: alp1 = alpha1, alp2 = alpha2, 
# Weibull parameter: beta1 <- exp(a0+a1*age.grp+a2*gen+a3*donor), here p are est_x
# Weibull parameter: beta2 <- exp(c0+c1*age.grp+c2*gen+c3*donor), here q are est_y
reg_coef <- c(est_alp1, lwci_alp1, upci_alp1,
              est_alp2, lwci_alp2, upci_alp2, 
              
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
reg_coef <- round(reg_coef, 2)
reg_coef
data.frame(reg_coef, row.names=NULL)

row.names(reg_coef) <-c("alpha1","alpha2",      # shape parameters in Weibull distributions
                        "a0", "a1","a2","a3",   # regression coefficients for hazard 1 (graft failure)
                        "c0","c1","c2","c3",    # regression coefficients for hazard 2 (death)
                        "b0","b1","b2","b3"     # regression coefficients for association parameter
                        )
################################################################################
# Output results                                                               #
################################################################################
dir_results <- paste0(dir_data, "results/real_data_analysis/revision_1/")
write.csv(reg_coef, paste0(dir_results, "parameters_",copula, "_", survival_distribution,".csv"))
write.csv(results, paste0(dir_results, table_ref, "_", copula, "_",survival_distribution, ".csv"))
print(paste0("saved output for ", copula, " ", survival_distribution, "!"))
