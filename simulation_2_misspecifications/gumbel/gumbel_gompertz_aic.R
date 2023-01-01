#######################################################################################################
# Simulation study: evaluation of misspecification of survival distributions  #
# Content: underlying model - gumbel copula gompertz survival model - evaluation aic
# Original code by LS; reviewed, edited and updated by YW for paper2
# Data are simulated from gumbel copula gompertz distribution
# YW, 24 July 2021: 1. correct bias, mse and re-calculate mse without using loop
#                   2. rename variables and define vectors to save to estimates
#                   3. set up working directory, save output to estimates and summary, debug the code
#                   4. put likelihood to functions outside the loop
#                   5. rewrite some calculations by using vectors to improve efficiency
# YW, 1/1/2023:     1. update output directory and tidy up
#                   2. put starting values out of the loops
#                   3. move likelihood functions from this script to a generic functions folder
#######################################################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(numDeriv)
start_time = Sys.time()

#####################################################################################
#Output directory and output files                                                  #
#####################################################################################
## directory if on own PC
dir_results = "../../results/simulation_results/"

# # directory if on cluster
#dir = "/home/ywei/Simulation/Paper2/Gumbel"
#setwd(dir)

# likelihood function
source("functions/function_sim2.R")

# output files
out_file_summary <- "S2_misspec_underlying_gumbel_gompertz_summary.csv"
out_file_estimates <-"S2-misspec_underlying_gumbel_gompertz_estimates.csv"

#####################################################################################
#################### Gumbel, age, gen from gom chose with aic #######################
#####################################################################################
set.seed(12345) # changed by YW
n <- 3000
runs <- 1

true_b0 <- -2.25
true_b1 <- 1.31

true_g1 <- 0.001 #gomp gamma1
true_g2 <- 0.06  #gomp gamma2
true_p0 <- -3.37 #gomp lambda1
true_p1 <- 0.11 #gomp lambda1
true_q0 <- -4.82 #gomp lambda2
true_q1 <- 1.49 #gomp lambda2

true_theta_d0 <- exp(true_b0)+1
true_theta_d1 <- exp(true_b0+true_b1)+1

t_theta_d0_cop <- gumbelCopula(true_theta_d0)
true_rho_d0 <- rho(t_theta_d0_cop)

t_theta_d1_cop <- gumbelCopula(true_theta_d1)
true_rho_d1 <- rho(t_theta_d1_cop)

#S1 exp, S2 weib
true_hr_l1 <- exp(true_p1)
true_hr_l2 <- exp(true_q1)

true_lambda1 <- rep(0,n)
true_lambda2 <- rep(0,n)
true_t <- rep(0,n)
true_r <- rep(0,n)
U1 <- rep(0,n)
V1 <- rep(0,n)

## YW: added lower and upper bounds of 95%CI
save_hr_l1 <-  hr_l1_lwci <- hr_l1_upci <- rep(0,runs)
save_hr_l2 <- hr_l2_lwci <- hr_l2_upci <- rep(0,runs)
save_rho_d0 <- rho_d0_lwci <- rho_d0_upci<- rep(0,runs)
save_rho_d1 <- rho_d1_lwci <- rho_d1_upci<- rep(0,runs)

# YW: added theta, be consistent with the reporting statistics in the main text
theta_d0 <- theta_d0_lwci <- theta_d0_upci <- rep(0, n)
theta_d1 <- theta_d1_lwci <- theta_d1_upci <- rep(0, n)

# counter for the model selected by aic
counter_exp = counter_wei = counter_gom = 0

hr_1_lw = hr_1_up = hr_1_cross =  hr_2_lw = hr_2_up = hr_2_cross = 0

#################################################################################
# Specification of starting values for optim                                    #
#################################################################################
gumbel_exp_optim_lower = c(-10, -10, -10, -10,  -5,  -5) # lower bound 
gumbel_exp_optim_upper = c(-2.0,  1.5, -2.0,  3.0, -1.0,  4.0)# upper bound 
# YW changed the last two elements starting values from true_b0 and true_b1
gumbel_exp_optim_starting_values = c(-3.5, 0.2, -4.3, 1.4,-2, 1) # starting values 
# checking lower == gumbe_exp_optim_lower

gumbel_wei_optim_lower = c(0.01, -10.00, -10.00,  0.01, -10.00, -10.00,  -6.00,  -5.00) # lower bound 
gumbel_wei_optim_upper = c(1.5, -1.5,  1.0,  1.5, -2.0,  2.5, -1.0,  3.0)# upper bound 
# change last two elements starting values from true_b0 and true_b1
gumbel_wei_optim_starting_values =c(0.6, -3, 0, 0.9, -4, 1.5, -2, 1) # starting values 
# checking lower == gumbel_wei_optim_lower

# gamma1 (g1), p0, true_p1, gamma2 (g2), q0, q1, b0, b1
gumbel_gom_optim_lower = c(-0.25, -6.00, -2.00, -0.20, -6.00, -1.00, -4.00,  0.00) # lower bound
gumbel_gom_optim_upper = c(0.8, -2.0,  2.0,  0.1, -2.0,  3.0, -1.0,  3.0)# upper bound 
#gumbel_gom_optim_starting_values = c(0.01, true_p0, true_p1, true_g2, true_q0, true_q1, true_b0, true_b1) # starting values 
gumbel_gom_optim_starting_values = c(0.01, -3, 0.02, 0.01, -4, 1, -2, 1) # starting values 

# replicate 'runs' times 
for (i in 1:runs){
  ###############################################################
  ######################## generate data ########################
  ###############################################################
  
  #Step 1: generate age categories
  age.grp <- rbinom(n,1,0.40)          #40% are in the older age group in NHSBT data
  
  for(k in 1: n){   #loop to generate U an V from age-varying theta
    m=1                  
    
    #Step 2: generate 1 random variable from Uniform(0,a) distribution 
    u1 <- runif(m,0,1)       
    
    #Step 3: X_true generated from u1 values (T1 from later)
    theta1 <- exp(true_b0 + true_b1 * age.grp[k])+1
    true_lambda1[k] <- exp(true_p0 + true_p1 * age.grp[k])
    true_lambda2[k] <- exp(true_q0 + true_q1 * age.grp[k]) 
    
    #Step 4: Conditional distribution method
    fc<- gumbelCopula(theta1, dim=2) #only allows 1 theta at a time (-> loop)
    uv<- cCopula(cbind(u1, runif(m)), copula = fc, inverse = TRUE) #gives vector (u1,v) - new v
    #this generates v using theta1 and u1 
    u<-uv[,1]  #split u and v from the results of cdm
    v<-uv[,2]
    
    #SAVE:
    U1[k]=u     #add to u and v vectors on the outside
    V1[k]=v
    true_t[k] <- theta1  #save theta for this individual
    
  }
  
  #Step 4: T1 and T2 from inverse survival
  T1 <- 1/true_g1 *log (1-true_g1/true_lambda1 *log(U1))
  T2 <- 1/true_g2 *log (1-true_g2/true_lambda2 *log(V1))
  
  #Step 7: Follow up time C, censoring variable
  C<-runif(n,0,25) 
  
  #Step 8 and 9: Find observed X and Y and indicators
  X <- pmin(T1,T2,C)
  Y <- pmin(T2, C) 
  d1 <- ifelse(T1<=Y,1,0) 
  d2 <- ifelse(T2<=C,1,0) 
  
  #Step 10: Create dataframe, true values of X and Y have association theta=b0+b1*X
  df <- data.frame(X, Y, d1, d2, age.grp)
  df$X[which(df$X==0)] <- 0.1
  df$Y[which(df$Y==0)] <- 0.1
  
  ########################################################
  ################ Gumbel pseudo likelihood ##############
  ##################### Exponential ######################
  ########################################################
  plgoptim_exp <- optim(gumbel_exp_optim_starting_values, gpl_exp, method="L-BFGS-B",
                        lower=gumbel_exp_optim_lower, upper=gumbel_exp_optim_upper, 
                        X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                        control=list(fnscale=-1),hessian=TRUE)
  
  index_lower = which(plgoptim_exp$par == gumbel_exp_optim_lower)
  index_upper = which(plgoptim_exp$par == gumbel_exp_optim_upper)
  
  if(length(index_lower)>0)
  {
    counter_exp_low[index_lower] = counter_exp_low[index_lower]+1
    next
  }
  if(length(index_upper)>0)
  {
    counter_exp_upper[index_upper] = counter_exp_upper[index_upper]+1
    next
  }
  
  ########################################################
  ################ Gumbel pseudo likelihood ##############
  ####################### Weibull ########################
  ########################################################

  plgoptim_wei <- optim(gumbel_wei_optim_starting_values, gpl_wei, method="L-BFGS-B",
                        lower=gumbel_wei_optim_lower, upper=gumbel_wei_optim_upper, 
                        X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                        control=list(fnscale=-1),hessian=TRUE)
  
  
  index_lower = which(plgoptim_exp$par == gumbel_exp_optim_lower)
  index_upper = which(plgoptim_exp$par == gumbel_exp_optim_upper)
  
  if(length(index_lower)>0)
  {
    counter_exp_low[index_lower] = counter_exp_low[index_lower]+1
    next
  }
  if(length(index_upper)>0)
  {
    counter_exp_upper[index_upper] = counter_exp_upper[index_upper]+1
    next
  }
  
  ########################################################
  ############### Gumbel pseudo likelihood ###############
  ####################### Gompertz #######################
  ########################################################
  plgoptim_gom <- optim(gumbel_gom_optim_starting_values, gpl_gom, method="L-BFGS-B",
                        lower=gumbel_gom_optim_lower, upper=gumbel_gom_optim_upper, 
                        X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                        control=list(fnscale=-1),hessian=TRUE)
  
  
  index_lower = which(plgoptim_exp$par == gumbel_exp_optim_lower)
  index_upper = which(plgoptim_exp$par == gumbel_exp_optim_upper)
  
  if(length(index_lower)>0)
  {
    counter_exp_low[index_lower] = counter_exp_low[index_lower]+1
    next
  }
  if(length(index_upper)>0)
  {
    counter_exp_upper[index_upper] = counter_exp_upper[index_upper]+1
    next
  }
  
  ########################################################
  ######################### AICS #########################
  ########################################################
  
  loglik_exp <- gpl_exp(plgoptim_exp$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp)
  k_exp <- length(plgoptim_exp$par)
  aic_exp <- -2*loglik_exp+2*k_exp
  
  loglik_wei <- gpl_wei(plgoptim_wei$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp)
  k_wei <- length(plgoptim_wei$par)
  aic_wei <- -2*loglik_wei+2*k_wei
  
  loglik_gom <- gpl_gom(plgoptim_gom$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp)
  k_gom <- length(plgoptim_gom$par)
  aic_gom <- -2*loglik_gom+2*k_gom
  
  ########################################################
  ####################### RESULTS ########################
  ########################################################
  
  aics <- c(aic_exp, aic_wei, aic_gom)                                             #all AIC values
  cops <- c("Exponential", "Weibull", "Gompertz")                                  #Names of distributions
  index <- which.min(aics)                                                         #gives index distribution with min aic wrt aics     
  print(cops[index])                                                               #print which distribution
  
  # YW: checking
  #index = 1
  if (index==1) {counter_exp=counter_exp+1
  } else if (index==2) {counter_wei=counter_wei+1
  } else {counter_gom=counter_gom+1}
  
  if (index==1){ #chooses exponential
    
    fisher_info <- solve(-plgoptim_exp$hessian) #inverse -hess
    se <- sqrt(diag(fisher_info)) 
    
    #b ci
    est_a0 <- plgoptim_exp$par[1]
    est_a1 <- plgoptim_exp$par[2]
    est_c0 <- plgoptim_exp$par[3]
    est_c1<- plgoptim_exp$par[4]
    est_b0 <- plgoptim_exp$par[5]
    est_b1 <- plgoptim_exp$par[6]
    varb0 <- fisher_info[5,5]
    varb1 <- fisher_info[6,6]
    covb0b1 <- fisher_info[5,6]
    
    #rho for gen=0
    theta_d0[i] <- exp(est_b0)+1
    var_theta_d0 <- exp(2*est_b0)*varb0
    theta_d0_lwci[i] <- theta_d0[i] - 1.96*sqrt(var_theta_d0)
    if (theta_d0_lwci[i]<1) {theta_d0_lwci[i]=1} #if hits lower bound take it as 1
    theta_d0_upci[i] <- theta_d0[i] + 1.96*sqrt(var_theta_d0)
    
    # YW: checking
   # print(c("position 1", theta_d0[i], theta_d0_lwci[i], theta_d0_upci[i]))
    
    rho_d0_cop <- gumbelCopula(theta_d0[i])
    est_rho_d0 <- rho(rho_d0_cop)
    save_rho_d0[i] <- est_rho_d0
    rho_d0_lwci_cop <- gumbelCopula(theta_d0_lwci[i])
    rho_d0_lwci[i] <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- gumbelCopula(theta_d0_upci[i])
    rho_d0_upci[i] <- rho(rho_d0_upci_cop) 
    
    #YW: checking
    #print(c("position 2",rho_d0_lwci_cop, rho_d0_upci_cop))
    #rho for gen=1
    theta_d1[i] <- exp(est_b0+est_b1)+1
    var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
    theta_d1_lwci[i] <- theta_d1[i] - 1.96*sqrt(var_theta_d1)
    theta_d1_upci[i] <- theta_d1[i] + 1.96*sqrt(var_theta_d1)
    if (theta_d1_lwci[i]<1) {theta_d1_lwci[i]=1} #if hits lower bound take it as 1
    
    rho_d1_cop <- gumbelCopula(theta_d1[i])
    est_rho_d1 <- rho(rho_d1_cop)
    save_rho_d1[i] <- est_rho_d1
    rho_d1_lwci_cop <- gumbelCopula(theta_d1_lwci[i])
    rho_d1_lwci[i] <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- gumbelCopula(theta_d1_upci[i])
    rho_d1_upci[i] <- rho(rho_d1_upci_cop) 
    
    #hrs
    vara0 <- fisher_info[1,1]
    vara1 <- fisher_info[2,2] 
    cov_a0a1 <- fisher_info[1,2] 
    varc0 <- fisher_info[3,3] 
    varc1 <- fisher_info[4,4] 
    cov_c0c1 <- fisher_info[3,4] 
    
    est_hr_l1 <- exp(est_a1)
    est_hr_l2 <- exp(est_c1)
    save_hr_l1[i] <- est_hr_l1
    save_hr_l2[i] <- est_hr_l2
    
    var_hr_l1 <- exp(est_a1)^2*vara1
    var_hr_l2 <- exp(est_c1)^2*varc1
    
    hr_l1_lwci[i] <- est_hr_l1 - 1.96*sqrt(var_hr_l1)
    hr_l1_upci[i] <- est_hr_l1 + 1.96*sqrt(var_hr_l1)
    hr_l2_lwci[i] <- est_hr_l2 - 1.96*sqrt(var_hr_l2)
    hr_l2_upci[i] <- est_hr_l2 + 1.96*sqrt(var_hr_l2)
  } else if (index==2){
    
    fisher_info <- solve(-plgoptim_wei$hessian) #inverse -hess
    se <- sqrt(diag(fisher_info)) 
    
    #point and var
    est_a1 <- plgoptim_wei$par[1]
    est_a2 <- plgoptim_wei$par[4]
    est_x0 <- plgoptim_wei$par[2]
    est_x1 <- plgoptim_wei$par[3]
    est_y0 <- plgoptim_wei$par[5]
    est_y1 <- plgoptim_wei$par[6]
    est_b0 <- plgoptim_wei$par[7]
    est_b1 <- plgoptim_wei$par[8]
    varb0 <- fisher_info[7,7]
    varb1 <- fisher_info[8,8]
    covb0b1 <- fisher_info[7,8]
    
    #rho for gen=0
    theta_d0[i] <- exp(est_b0)+1
    var_theta_d0 <- exp(2*est_b0)*varb0
    theta_d0_lwci[i] <- theta_d0[i] - 1.96*sqrt(var_theta_d0)
    if (theta_d0_lwci[i]<1) {theta_d0_lwci[i]=1} #if hits lower bound take it as 1
    theta_d0_upci[i] <- theta_d0[i] + 1.96*sqrt(var_theta_d0)
    
    # YW: checking
    #print(c("position 1", theta_d0[i], theta_d0_lwci[i], theta_d0_upci[i]))
    
    rho_d0_cop <- gumbelCopula(theta_d0[i])
    est_rho_d0 <- rho(rho_d0_cop)
    save_rho_d0[i] <- est_rho_d0
    rho_d0_lwci_cop <- gumbelCopula(theta_d0_lwci[i])
    rho_d0_lwci[i] <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- gumbelCopula(theta_d0_upci[i])
    rho_d0_upci[i] <- rho(rho_d0_upci_cop) 
    
    #YW: checking
    #print(c("position 2",rho_d0_lwci_cop, rho_d0_upci_cop))
    #rho for gen=1
    theta_d1[i] <- exp(est_b0+est_b1)+1
    var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
    theta_d1_lwci[i] <- theta_d1[i] - 1.96*sqrt(var_theta_d1)
    theta_d1_upci[i] <- theta_d1[i] + 1.96*sqrt(var_theta_d1)
    if (theta_d1_lwci[i]<1) {theta_d1_lwci[i]=1} #if hits lower bound take it as 1
    
    rho_d1_cop <- gumbelCopula(theta_d1[i])
    est_rho_d1 <- rho(rho_d1_cop)
    save_rho_d1[i] <- est_rho_d1
    rho_d1_lwci_cop <- gumbelCopula(theta_d1_lwci[i])
    rho_d1_lwci[i] <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- gumbelCopula(theta_d1_upci[i])
    rho_d1_upci[i] <- rho(rho_d1_upci_cop) 
    
    #hrs
    varx1 <- fisher_info[3,3]
    vary1 <- fisher_info[6,6] 
    
    est_hr_l1 <- exp(est_x1) #hr for gf
    est_hr_l2 <- exp(est_y1) #hr for death
    save_hr_l1[i] <- est_hr_l1
    save_hr_l2[i] <- est_hr_l2
    
    var_hr_l1 <- exp(est_x1)^2 * varx1
    var_hr_l2 <- exp(est_y1)^2 * vary1
    
    hr_l1_lwci[i] <- est_hr_l1 - 1.96*sqrt(var_hr_l1)
    hr_l1_upci[i] <- est_hr_l1 + 1.96*sqrt(var_hr_l1)
    hr_l2_lwci[i] <- est_hr_l2 - 1.96*sqrt(var_hr_l2)
    hr_l2_upci[i] <- est_hr_l2 + 1.96*sqrt(var_hr_l2)
  
  } else{ # gompertz is selected
    #hessian
    hessian <- hessian(gpl_gom, plgoptim_gom$par, X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp)
    #fishers info matrix
    fisher_info <- solve(-hessian) 
    #Standard Error
    se <- sqrt(diag(fisher_info)) 
    
    est_b0 <- plgoptim_gom$par[7]
    est_b1 <- plgoptim_gom$par[8]
    est_g1 <- plgoptim_gom$par[1]
    est_g2 <- plgoptim_gom$par[4]
    est_p0 <- plgoptim_gom$par[2]
    est_p1 <- plgoptim_gom$par[3]
    est_q0 <- plgoptim_gom$par[5]
    est_q1 <- plgoptim_gom$par[6]
    varb0 <- fisher_info[7,7]
    varb1 <- fisher_info[8,8]
    covb0b1 <- fisher_info[7,8]
    
    #rho for gen=0
    theta_d0[i] <- exp(est_b0)+1
    var_theta_d0 <- exp(2*est_b0)*varb0
    theta_d0_lwci[i] <- theta_d0[i] - 1.96*sqrt(var_theta_d0)
    if (theta_d0_lwci[i]<1) {theta_d0_lwci[i]=1} #if hits lower bound take it as 1
    theta_d0_upci[i] <- theta_d0[i] + 1.96*sqrt(var_theta_d0)
    
    # YW: checking
    print(c("position 1", theta_d0[i], theta_d0_lwci[i], theta_d0_upci[i]))
    
    rho_d0_cop <- gumbelCopula(theta_d0[i])
    est_rho_d0 <- rho(rho_d0_cop)
    save_rho_d0[i] <- est_rho_d0
    rho_d0_lwci_cop <- gumbelCopula(theta_d0_lwci[i])
    rho_d0_lwci[i] <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- gumbelCopula(theta_d0_upci[i])
    rho_d0_upci[i] <- rho(rho_d0_upci_cop) 
    
    #YW: checking
   # print(c("position 2",rho_d0_lwci_cop, rho_d0_upci_cop))
    #rho for gen=1
    theta_d1[i] <- exp(est_b0+est_b1)+1
    var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
    theta_d1_lwci[i] <- theta_d1[i] - 1.96*sqrt(var_theta_d1)
    theta_d1_upci[i] <- theta_d1[i] + 1.96*sqrt(var_theta_d1)
    if (theta_d1_lwci[i]<1) {theta_d1_lwci[i]=1} #if hits lower bound take it as 1
    
    rho_d1_cop <- gumbelCopula(theta_d1[i])
    est_rho_d1 <- rho(rho_d1_cop)
    save_rho_d1[i] <- est_rho_d1
    rho_d1_lwci_cop <- gumbelCopula(theta_d1_lwci[i])
    rho_d1_lwci[i] <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- gumbelCopula(theta_d1_upci[i])
    rho_d1_upci[i] <- rho(rho_d1_upci_cop) 
    
    #hrs
    varp1 <- fisher_info[3,3]
    varq1 <- fisher_info[6,6] 
    
    est_hr_l1 <- exp(est_p1) #hr for gf
    est_hr_l2 <- exp(est_q1) #hr for death
    save_hr_l1[i] <- est_hr_l1
    save_hr_l2[i] <- est_hr_l2
    
    var_hr_l1 <- exp(est_p1)^2 * varp1
    var_hr_l2 <- exp(est_q1)^2 * varq1
    
    hr_l1_lwci[i] <- est_hr_l1 - 1.96*sqrt(var_hr_l1)
    hr_l1_upci[i]  <- est_hr_l1 + 1.96*sqrt(var_hr_l1)
    hr_l2_lwci[i]  <- est_hr_l2 - 1.96*sqrt(var_hr_l2)
    hr_l2_upci[i]  <- est_hr_l2 + 1.96*sqrt(var_hr_l2)

  }
  
  print(i)
} # end of loop

#bias
hr_l1_bias <- mean(save_hr_l1 -true_hr_l1)
hr_l2_bias <- mean(save_hr_l2 -true_hr_l2)

rho_d0_bias <- mean(save_rho_d0 - true_rho_d0)
rho_d1_bias <- mean(save_rho_d1 - true_rho_d1)

theta_d0_bias <- mean(theta_d0 - true_theta_d0)
theta_d1_bias <- mean(theta_d1 - true_theta_d1)

#coverage
hr_l1_cov <- 100* sum(true_hr_l1 <= hr_l1_upci & true_hr_l1 >= hr_l1_lwci)/runs
hr_l2_cov <- 100* sum(true_hr_l2 <= hr_l2_upci & true_hr_l2 >= hr_l2_lwci)/runs

rho_d0_cov <- 100* sum(true_rho_d0 <= rho_d0_upci & true_rho_d0 >= rho_d0_lwci)/runs
rho_d1_cov <- 100* sum(true_rho_d1 <= rho_d1_upci & true_rho_d1 >= rho_d1_lwci)/runs

theta_d0_cov <- 100*sum(true_theta_d0 <= theta_d0_upci & true_theta_d0 >= theta_d0_lwci)/runs
theta_d1_cov <- 100*sum(true_theta_d1 <= theta_d1_upci & true_theta_d1 >= theta_d1_lwci)/runs

#mse
hr_l1_mse <- mean((save_hr_l1 -true_hr_l1)^2)
hr_l2_mse <- mean((save_hr_l2 -true_hr_l2)^2)

rho_d0_mse <- mean((save_rho_d0-true_rho_d0)^2)
rho_d1_mse <- mean((save_rho_d1-true_rho_d1)^2)

theta_d0_mse <- mean((theta_d0 - true_theta_d0)^2)
theta_d1_mse <- mean((theta_d1 - true_theta_d1)^2)

#counters#
exp_perc <- counter_exp / runs *100
wei_perc <- counter_wei / runs *100
gom_perc <- counter_gom / runs *100

#change meaning of hr
#true hr 1 = 1.5 & true hr 2 = 4 for age
hr_1_perc <- hr_1_up / runs *100
hr_2_perc <- hr_2_up / runs *100

end_time = Sys.time()
run_time = end_time - start_time
run_time

# mean of bias
# hr_l1 represents non-terminal event; hr_l2 represents terminal event
bias <- c(hr_l1_bias, hr_l2_bias, rho_d0_bias, rho_d1_bias, theta_d0_bias, theta_d1_bias)

# Coverage probability
CP <- c(hr_l1_cov,hr_l2_cov, rho_d0_cov, rho_d1_cov, theta_d0_cov, theta_d1_cov)
MSE <- c(hr_l1_mse, hr_l2_mse, rho_d0_mse, rho_d1_mse,theta_d0_mse, theta_d1_mse)

# in the order of exponential, weibull, gompertz.
percentage_chosen = c(exp_perc, wei_perc, gom_perc, "na", "na", "na")

# put results together
items<-c("hr_nt", "hr_t", "rho_reference", "rho_covariates", "theta_reference", "theta_covariates")
Results <- cbind.data.frame(items, bias, CP, MSE, percentage_chosen)
Results[,2:4] <- round(Results[,2:4],3)
Results
rownames(Results)<-NULL
end_time <- Sys.time()
run_time = end_time - start_time
run_time

Estimates = data.frame(hr.l1= save_hr_l1, hr.l1.low= hr_l1_lwci, hr.l1.up = hr_l1_upci,
                       hr.l2= save_hr_l2, hr.l2.low= hr_l2_lwci, hr.l2.up = hr_l2_upci,
                       rho.d0= save_rho_d0, rho.d0.low= rho_d0_lwci, rho.d0.up = rho_d0_upci,
                       rho.d1= save_rho_d1, rho.d1.low= rho_d1_lwci, rho.d1.up = rho_d1_upci)

# output results
write.csv(Results, row.names=F,file=paste0(dir_results, out_file_summary))
write.csv(Estimates, row.names=F,file=paste0(dir_results,out_file_estimates))
print("Simulation 2 for gumbel gompertz model completed successfully!")
# percentage chosen is recorded in the order of exponential, weibull and gompertz. The true model is weibull.