################################################################################
# 9-July-2021
# YW: NHSBT data analysis, Clayton copula exponential survival distribution
# original script by LS; edited and updated for paper2 by YW
# Table 3: exponential survival clayton copula model
# YW: 2022-02-20: update load data and output results sections
# YW: 2022-12-30: output regression coefficients and add comments
################################################################################

rm(list=ls())
library(copula);library(mvtnorm);library(plyr)
start_time = Sys.time()

################################################################################
# Set up model specs and load data                                             #
################################################################################
# Model specs
copula <- "clayton"
survival_distribution <- "exp" 
if(survival_distribution == "exp") {table_ref = "table3"}
if(survival_distribution == "gompertz") {table_ref = "table4"}
if(survival_distribution == "weibull") {table_ref = "table5"}

# Load data
# set working directory to project directory and send through the next two lines
dir_data <-dir_results <- "../../"
df <- read.csv(file=paste0(dir_data, "NHSBT/paper2_data_v2.csv"))

################################################################################
# Clayton pseudo likelihood                                                    #
################################################################################
cpl <- function(para, X, Y, d1, d2, donor, age.grp, gen){
  a0 <- para[1]       # regression parameters in lambda1 (hazard for graft failure)
  a1 <- para[2]
  a2 <- para[3]
  a3 <- para[4]
  
  c0 <- para[5]       # regression parameters in lambda2 (hazard for death)
  c1 <- para[6]
  c2 <- para[7]
  c3 <- para[8]
  
  b0 <- para[9]       # regression parameters for association parameter
  b1 <- para[10]
  b2 <- para[11]
  b3 <- para[12]
  
  lambda1 <- exp(a0+a1*age.grp+a2*gen+a3*donor)    # hazard function for graft failure
  lambda2 <- exp(c0+c1*age.grp+c2*gen+c3*donor)    # hazard function for death
  
  S1 <- exp(-lambda1*X)        # survival function for graft failure
  S2 <- exp(-lambda2*Y)        # survival function for death
  
  theta <- exp(b0+b1*age.grp+b2*gen+b3*donor)      # association parameter
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  part1 <- d1*d2*(log(1+theta)+(1+2*theta)*log(C)-(theta+1)*log(S1)-(theta+1)*log(S2)+log(lambda1)-lambda1*X+log(lambda2)-lambda2*Y)
  part2 <- d1*(1-d2)*((theta+1)*log(C)-(theta+1)*log(S1)+log(lambda1)-lambda1*X)
  part3<-((1-d1)*(d2))*((theta+1)*log(C)-(theta+1)*log(S2)+log(lambda2)-lambda2*Y)
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

plcoptim <- optim(c(-1,-0.01,-0.01,-0.01,  
                    -1,-0.01,-0.01,-0.01,  
                    2,0.1,0.1,0.1
                    ), cpl, method="L-BFGS-B",
                  lower=c(-10,-10,-10,-10,  
                          -10,-10,-10,-10,  
                          0.01,-1,-1,-1
                          ),
                  upper=c(-1,1,0.01,0.01, 
                          -2,2,0.01,0.01,
                          10,6,3,3
                          ), 
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, donor=df$donor, gen=df$gen,
                  control=list(fnscale=-1),hessian=TRUE)

end_time = Sys.time()
run_time = end_time - start_time
run_time
plcoptim$par

################################################################################
# Confidence Intervals                                                         #
################################################################################

#Fisher's Information matrix
fisher_info<-solve(-plcoptim$hessian) #inverse -hess
#Standard error = sqrt(var/n)
se<-sqrt(diag(fisher_info)) 

#association parameter: regression coefficient beta ci#
est_b0 <- plcoptim$par[9]
est_b1 <- plcoptim$par[10]
est_b2 <- plcoptim$par[11]
est_b3 <- plcoptim$par[12]
lwci_b0 <- est_b0 - 1.96*se[9]     
lwci_b1 <- est_b1 - 1.96*se[10]
lwci_b2 <- est_b2 - 1.96*se[11]     
lwci_b3 <- est_b3 - 1.96*se[12]
upci_b0 <- est_b0 + 1.96*se[9]
upci_b1 <- est_b1 + 1.96*se[10] 
upci_b2 <- est_b2 + 1.96*se[11]
upci_b3 <- est_b3 + 1.96*se[12]

# recipient age on association
beta_age <- c(est_b1, lwci_b1, upci_b1)
beta_gender<- c(est_b2, lwci_b2, upci_b2)
beta_donor<- c(est_b3, lwci_b3, upci_b3)
beta_age
beta_gender
beta_donor

#non-terminal event: a ci#
est_a0 <- plcoptim$par[1]
est_a1 <- plcoptim$par[2]
est_a2 <- plcoptim$par[3]
est_a3 <- plcoptim$par[4]
lwci_a0 <- est_a0 - 1.96*se[1]     
lwci_a1 <- est_a1 - 1.96*se[2]
lwci_a2 <- est_a2 - 1.96*se[3]     
lwci_a3 <- est_a3 - 1.96*se[4]
upci_a0 <- est_a0 + 1.96*se[1]
upci_a1 <- est_a1 + 1.96*se[2] 
upci_a2 <- est_a2 + 1.96*se[3]
upci_a3 <- est_a3 + 1.96*se[4] 

#terminal event: c ci#
est_c0 <- plcoptim$par[5]
est_c1 <- plcoptim$par[6]
est_c2 <- plcoptim$par[7]
est_c3 <- plcoptim$par[8]
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

# recipient gender
hr_gf_gender <-c(esthr_l1_gen,hr_l1_lwci_gen, hr_l1_upci_gen)
hr_gf_gender
hr_d_gender <-c(esthr_l2_gen,hr_l2_lwci_gen, hr_l2_upci_gen)
hr_d_gender

# recipient age
hr_gf_age <-c(esthr_l1_age,hr_l1_lwci_age, hr_l1_upci_age)
hr_gf_age

hr_d_age <-c(esthr_l2_age,hr_l2_lwci_age, hr_l2_upci_age)
hr_d_age

# donor type
hr_gf_donor <-c(esthr_l1_donor,hr_l1_lwci_donor, hr_l1_upci_donor)
hr_gf_donor

hr_d_donor <-c(esthr_l2_donor,hr_l2_lwci_donor, hr_l2_upci_donor)
hr_d_donor

# Results --------------------------------------------------------------------
# recipient age
hr_gf_age <-c(esthr_l1_age,hr_l1_lwci_age, hr_l1_upci_age)
hr_gf_age
hr_d_age <-c(esthr_l2_age,hr_l2_lwci_age, hr_l2_upci_age)
hr_d_age

# recipient gender
hr_gf_gender <-c(esthr_l1_gen,hr_l1_lwci_gen, hr_l1_upci_gen)
hr_gf_gender
hr_d_gender <-c(esthr_l2_gen,hr_l2_lwci_gen, hr_l2_upci_gen)
hr_d_gender

# donor type
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
para <- c(est_a0, est_a1, est_a2, est_a3, est_c0, est_c1, est_c2, est_c3, est_b0, est_b1, est_b2, est_b3)
loglik <- cpl(para,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
k<-length(para)
n<-length(df$X)
aic<- -2*loglik+2*k
bic<- -2*loglik+log(n)*k
loglik
aic
bic

results$aic[1] = round(aic,1)
results$run_time[1]= round(run_time,2)

################################################################################
# Create a data frame for regression coefficients                              #
################################################################################

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
reg_coef <- round(reg_coef, 2)
data.frame(reg_coef, row.names=NULL)
reg_coef

################################################################################
# Output results                                                               #
################################################################################
dir_results <- paste0(dir_data, "results/real_data_analysis/revision_1/")
write.csv(reg_coef, paste0(dir_results, "parameters_",copula, "_", survival_distribution,".csv"))
write.csv(results, paste0(dir_results, table_ref, "_", copula, "_",survival_distribution, ".csv"))
print(paste0("saved output for ", copula, " ", survival_distribution, "!"))
