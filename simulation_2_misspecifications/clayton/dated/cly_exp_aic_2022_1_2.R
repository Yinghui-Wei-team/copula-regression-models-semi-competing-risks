#########################################################################################
# Simulation study: evaluation of misspecification of survival distributions   
# Data are simulated from Clayton copula exponential distribution
# Original code by LS; reviewed, edited and updated by YW for paper2
# YW, 24 July 2021: 1. correct bias, mse and re-calculate mse without using loops
#                   2. rename variables and define vectors to save to estimates
#                   3. set up working directory, save output to estimates and summary, 
#                      debug the code
#                   4. put likelihood to functions outside the loop
#                   5. rewrite some calculations by using vectors to improve efficiency
# YW, 1/1/2023:     update output directory and tidy up
#########################################################################################

rm(list=ls())
library(copula); library(mvtnorm); library(numDeriv)
start_time = Sys.time()

#####################################################################################
#Output directory and output files                                                  #
#####################################################################################
## directory if on own PC
dir_results = "../../results/simulation_results/"

# # directory if on cluster
# dir = "/home/ywei/Simulation/Paper2/Clayton"
# setwd(dir)

# likelihood function
source("functions/function_sim2.R")

out_file_summary <- "S2_misspec_underlying_clayton_exp_summary.csv"
out_file_estimates <-"S2_misspec_underlying_clayton_exp_estimates.csv"

#####################################################################################
################## Clayton, age, gen from exp chose with aic ########################
#####################################################################################
set.seed(12345)
n <- 3000
runs <- 3

true_b0 <- 0.62
true_b1 <- 1.04

true_a0 <- -3.42 #exp lambda1
true_a1 <-  0.38 #exp lambda1
true_c0 <- -4.28 #exp lambda2
true_c1 <-  1.41 #exp lambda2

# YW, added
true_values <- c(true_a0, true_a1, true_c0, true_c1, true_b0, true_b1)
n_parameters = length(true_values)

true_theta_d0 <- exp(true_b0)
true_theta_d1 <- exp(true_b0+true_b1)

t_theta_d0_cop <- claytonCopula(true_theta_d0)
true_rho_d0 <- rho(t_theta_d0_cop)

t_theta_d1_cop <- claytonCopula(true_theta_d1)
true_rho_d1 <- rho(t_theta_d1_cop)

#S1 exp, S2 weib
true_hr_l1 <- exp(true_a1)
true_hr_l2 <- exp(true_c1)

true_l1 <- rep(0,n)
true_l2 <- rep(0,n)
true_beta1 <- rep(0,n)
true_beta2 <- rep(0,n)
true_t <- rep(0,n)
true_r <- rep(0,n)
U1 <- rep(0,n)
V1 <- rep(0,n)

hr_1_lw = hr_1_up = hr_1_cross =  hr_2_lw = hr_2_up = hr_2_cross = 0

## YW: added lower and upper bounds of 95%CI
save_hr_l1 <-  hr_l1_lwci <- hr_l1_upci <- rep(0,runs)
save_hr_l2 <- hr_l2_lwci <- hr_l2_upci <- rep(0,runs)
save_rho_d0 <- rho_d0_lwci <- rho_d0_upci<- rep(0,runs)
save_rho_d1 <- rho_d1_lwci <- rho_d1_upci<- rep(0,runs)

# YW: added theta, be consistent with the reporting statistics in the main text
theta_d0 <- theta_d0_lwci <- theta_d0_upci <- rep(0, n)
theta_d1 <- theta_d1_lwci <- theta_d1_upci <- rep(0, n)

# how many times these distributions were selected by the AIC
counter_exp = 0
counter_wei = 0
counter_gom = 0

# count number of runs outside the defined boundaries for each parameter
counter_exp_low <- counter_exp_upper <- rep(0, n_parameters)
counter_wei_low <- counter_wei_upper <- rep(0, n_parameters+2)
counter_gom_low <- counter_gom_upper <- rep(0, n_parameters+2)

#################################################################################
# Specification of starting values for optim                                    #
#################################################################################
# boundaries used for parameters in optim for clayton exponential likelihood 
# Clayton exponential
clayton_exp_optim_lower = c(-10, -10, -10, -10, -1, -5.5)  # lower bound for a0, a1, c0, c1, b0, b1
clayton_exp_optim_upper = c(-1,  1, -1,  2,  3,  3)        # upper bound for a0, a1, c0, c1, b0, b1
clayton_exp_optim_starting_values = c(-3,0.01,-3,0.01,3,0) # starting values for a0, a1, c0, c1, b0, b1

# clayton weibull 
clayton_wei_optim_lower = c(0.01, -10.00, -10.00,   0.01, -10.00, -10.00, -10.00, -15.00)   # lower bound for a0, a1, c0, c1, b0, b1
clayton_wei_optim_upper = c(1.5, -1.0,  1.0,  1.5, -1.0,  3.0,  1.2,  3.0)                  # upper bound for a0, a1, c0, c1, b0, b1
# clayton_wei_optim_starting_values = c(0.67, -2.5, -0.6, 0.94, -3.3, -0.9, true_b0, true_b1) # starting values for a0, a1, c0, c1, b0, b1
clayton_wei_optim_starting_values = c(0.67, -2.5, -0.6, 0.94, -3.3, -0.9, 0.5, 1) # starting values for a0, a1, c0, c1, b0, b1

# boundaries and starting values for g1, p0, p1, g2, q0, q1, b0, b1
clayton_gom_optim_lower = c(-0.2,-5.0, -4.0, -0.2, -6.0, -4.0, -2.0, -2.0) # lower bound 
clayton_gom_optim_upper = c(0.1,-2.0,  2.0,  0.1, -2.0,  3.0,  1.2,  1.7) # upper bound 
clayton_gom_optim_starting_values =  c(-0.1,-4.0, -3.0, -0.1, -5.0, -3.0, -1.0, -1.0) # starting values 

##################################################################################
# likelihood function specification                                              #
##################################################################################

# clayton copula exponential distribution: one covariates
cpl_exp <- function(para, X, Y, d1, d2, age){
  
  a0 <- para[1]
  a1 <- para[2]
  c0 <- para[3]
  c1 <- para[4]
  b0 <- para[5]
  b1 <- para[6]
  
  lambda1 <- exp(a0+a1*age)
  lambda2 <- exp(c0+c1*age)
  S1<-exp(-lambda1*X)
  S2<-exp(-lambda2*Y)
  
  theta <- exp(b0+b1*age)
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  part1 <- d1*d2*(log(1+theta)+(1+2*theta)*log(C)-(theta+1)*log(S1)-(theta+1)*log(S2)+log(lambda1)-lambda1*X+log(lambda2)-lambda2*Y)
  part2 <- d1*(1-d2)*((theta+1)*log(C)-(theta+1)*log(S1)+log(lambda1)-lambda1*X)
  part3<-((1-d1)*(d2))*((theta+1)*log(C)-(theta+1)*log(S2)+log(lambda2)-lambda2*Y)
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  
  return(logpl)
}

cpl_wei <- function(para, X, Y, d1, d2, age){
  alpha1 <- para[1]
  x1 <- para[2]
  x2 <- para[3]
  alpha2 <- para[4]
  y1 <- para[5]
  y2 <- para[6]
  b0 <- para[7]
  b1 <- para[8]
  
  theta <- exp(b0+b1*age)  
  beta1 <- exp(x1+x2*age)
  beta2 <- exp(y1+y2*age)
  
  S1 <- exp(-beta1*X^alpha1)
  S2 <- exp(-beta2*Y^alpha2)
  S1[which(S1<0.1^8)]=0.1^8
  S2[which(S2<0.1^8)]=0.1^8
  
  S1S2 <- S1*S2
  S1S2[which(S1S2<0.1^8)]=0.1^8
  
  f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
  f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
  f1[which(f1<0.1^8)]=0.1^8
  f2[which(f2<0.1^8)]=0.1^8
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  C[which(C<0.1^8)] <- 0.1^8
  part1 <- d1*d2*(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1S2)^(1+theta)))
  
  part2 <- d1*(1-d2)*(log(C^(1+theta)*f1)-log(S1^(1+theta)))
  
  part3 <- d2*(1-d1)*(log(C^(1+theta)*f2)-log(S2^(1+theta)))
  
  part4 <- ((1-d1)*(1-d2))*log(C)
  
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}

cpl_gom <- function(para, X, Y, d1, d2, age){
  gamma1 <- para[1]
  p0 <- para[2]
  p1 <- para[3]
  gamma2 <- para[4]
  q0 <- para[5]
  q1 <- para[6]
  b0 <- para[7]
  b1 <- para[8]
  
  theta <- exp(b0+b1*age)
  lambda1 <- exp(p0+p1*age)
  lambda2 <- exp(q0+q1*age)
  
  S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
  S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
  S1[which(S1<0.1^8)]=0.1^8
  S2[which(S2<0.1^8)]=0.1^8
  
  S1S2 <- S1*S2
  S1S2[which(S1S2<0.1^8)]=0.1^8
  
  f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
  f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
  f1[which(f1<0.1^8)]=0.1^8
  f2[which(f2<0.1^8)]=0.1^8
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  C[which(C<0.1^8)] <- 0.1^8
  
  part1 <- d1*d2*(log((1+theta)*C^(1+2*theta)*f1*f2)-log((S1S2)^(1+theta)))
  
  part2 <- d1*(1-d2)*(log(C^(1+theta)*f1)-log(S1^(1+theta)))
  
  part3 <- d2*(1-d1)*(log(C^(1+theta)*f2)-log(S2^(1+theta)))
  
  part4 <- ((1-d1)*(1-d2))*log(C)
  
  logpl <- sum(part1+part2+part3+part4) 
  return(logpl)
}


################################################################################
# run 'runs' times                                                             #
################################################################################

for (i in 1:runs){
  
  ###############################################################
  ######################## generate data ########################
  ###############################################################
  
  #Step 1: generate age categories
  age <- rbinom(n,1,0.40)          #40% are in the older age group in NHSBT data
  
  for(k in 1:n){   #loop to generate U an V from age-varying theta
    m=1                  
    
    #Step 2: generate 1 random variable from Uniform(0,a) distribution 
    
    u1 <- runif(m,0,1)       
    
    #Step 3: X_true generated from u1 values (T1 from later)
    theta1 <- exp(true_b0 + true_b1 * age[k])
    true_l1[k] <- exp(true_a0 + true_a1 * age[k])
    true_l2[k] <- exp(true_c0 + true_c1 * age[k]) 
    
    #Step 4: Conditional distribution method
    
    fc<- claytonCopula(theta1, dim=2) #only allows 1 theta at a time (-> loop)
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
  T1 <- -log(U1)/true_l1
  T2 <- -log(V1)/true_l2
  
  #Step 7: Follow up time C, censoring variable
  C<-runif(n,0,25) 
  
  #Step 8 and 9: Find observed X and Y and indicators
  X<-pmin(T1,T2,C)
  Y<-pmin(T2, C) 
  d1<-ifelse(T1<=Y,1,0) 
  d2<-ifelse(T2<=C,1,0) 
  
  #Step 10: Create dataframe, true values of X and Y have association theta=b0+b1*X
  df<-data.frame(X, Y, d1, d2, age)
  
  ########################################################
  ############### Clayton pseudo likelihood ##############
  ##################### Exponential ######################
  ########################################################
  
  plcoptim_exp <- optim(clayton_exp_optim_starting_values, cpl_exp, method="L-BFGS-B",
                        lower=clayton_exp_optim_lower,upper=clayton_exp_optim_upper, 
                        X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                        control=list(fnscale=-1),hessian=TRUE)
  
  index_lower = which(plcoptim_exp$par == clayton_exp_optim_lower)
  index_upper = which(plcoptim_exp$par == clayton_exp_optim_upper)
  
  if(length(index_lower)>0)
  {
    counter_exp_low[index_lower] = counter_exp_low[index_lower]+1
    break
  }
  if(length(index_upper)>0)
  {
    counter_exp_upper[index_upper] = counter_exp_upper[index_upper]+1
    break
  }
  #-------------------YW end of re-write for clayton exponential-------------------
  
  ########################################################
  ############### Clayton pseudo likelihood ##############
  ####################### Weibull ########################
  ########################################################
  
  plcoptim_wei <- optim(clayton_wei_optim_starting_values, cpl_wei, method="L-BFGS-B",
                        lower=clayton_wei_optim_lower,
                        upper=clayton_wei_optim_upper, 
                        X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age, control=list(fnscale=-1), hessian=TRUE)
  
  index_lower = which(plcoptim_wei$par == clayton_wei_optim_lower)
  index_upper = which(plcoptim_wei$par == clayton_wei_optim_upper)
  
  
  if(length(index_lower)>0)
  {
    counter_wei_low[index_lower] = counter__wei_low[index_lower]+1
    break
  }
  if(length(index_upper)>0)
  {
    counter_wei_upper[index_upper] = counter_wei_upper[index_upper]+1
    break
  }
  
  
  ########################################################
  ############### Clayton pseudo likelihood ##############
  ####################### Gompertz #######################
  ########################################################
  
  plcoptim_gom <- optim(clayton_gom_optim_starting_values, cpl_gom, method="L-BFGS-B",
                        lower=clayton_gom_optim_lower,
                        upper=clayton_gom_optim_upper, 
                        X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                        control=list(fnscale=-1),hessian=TRUE)
  
  index_lower = which(plcoptim_gom$par == clayton_gom_optim_lower)
  index_upper = which(plcoptim_gom$par == clayton_gom_optim_upper)
  
  if(length(index_lower)>0)
  {
    counter_gom_low[index_lower] = counter_gom_low[index_lower]+1
    break
  }
  if(length(index_upper)>0)
  {
    counter_gom_upper[index_upper] = counter_gom_upper[index_upper]+1
    break
  }
  
  ########################################################
  ######################### AIC  #########################
  ########################################################
  
  loglik_exp <- cpl_exp(plcoptim_exp$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  k_exp <- length(plcoptim_exp$par)
  aic_exp <- -2*loglik_exp+2*k_exp
  
  loglik_wei <- cpl_wei(plcoptim_wei$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  k_wei <- length(plcoptim_wei$par)
  aic_wei <- -2*loglik_wei+2*k_wei

  loglik_gom <- cpl_gom(plcoptim_gom$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  k_gom <- length(plcoptim_gom$par)
  aic_gom <- -2*loglik_gom+2*k_gom
  
  ########################################################
  ####################### RESULTS ########################
  ########################################################
  
  aics <- c(aic_exp, aic_wei, aic_gom)                                             #all AIC values
  cops <- c("Exponential", "Weibull", "Gompertz")                                  #Names of distributions
  index <- which.min(aics)                                                         #gives index distribution with min aic wrt aics     
  print(cops[index])                                                               #print which distribution
  
  if (index==1) {counter_exp=counter_exp+1
  } else if (index==2) {counter_wei=counter_wei+1
    } else {counter_gom=counter_gom+1}
  
 
  if (index==1){ #chooses exponential
    
    fisher_info <- solve(-plcoptim_exp$hessian) #inverse -hess
    se <- sqrt(diag(fisher_info)) 
   
    #b ci
    est_a0 <- plcoptim_exp$par[1]
    est_a1 <- plcoptim_exp$par[2]
    est_c0 <- plcoptim_exp$par[3]
    est_c1<- plcoptim_exp$par[4]
    est_b0 <- plcoptim_exp$par[5]
    est_b1 <- plcoptim_exp$par[6]
    varb0 <- fisher_info[5,5]
    varb1 <- fisher_info[6,6]
    covb0b1 <- fisher_info[5,6]
    
    #rho for age=0
    theta_d0[i] <- exp(est_b0)
    var_theta_d0 <- exp(2*est_b0)*varb0
    theta_d0_lwci[i] <- theta_d0[i] - 1.96*sqrt(var_theta_d0)
    theta_d0_upci[i] <- theta_d0[i] + 1.96*sqrt(var_theta_d0)
    
    rho_d0_cop <- claytonCopula(theta_d0[i])
    est_rho_d0 <- rho(rho_d0_cop)
    rho_d0_lwci_cop <- claytonCopula(theta_d0_lwci[i])
    rho_d0_lwci[i] <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- claytonCopula(theta_d0_upci[i])
    rho_d0_upci[i] <- rho(rho_d0_upci_cop) 
    save_rho_d0[i] <- est_rho_d0
    
    #rho for age=1
    theta_d1[i] <- exp(est_b0+est_b1)
    var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
    theta_d1_lwci[i] <- theta_d1[i] - 1.96*sqrt(var_theta_d1)
    theta_d1_upci[i] <- theta_d1[i] + 1.96*sqrt(var_theta_d1)
    
    rho_d1_cop <- claytonCopula(theta_d1[i])
    est_rho_d1 <- rho(rho_d1_cop)
    rho_d1_lwci_cop <- claytonCopula(theta_d1_lwci[i])
    rho_d1_lwci[i] <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- claytonCopula(theta_d1_upci[i])
    rho_d1_upci[i] <- rho(rho_d1_upci_cop) 
    save_rho_d1[i] <- est_rho_d1
    
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
    
    } else if (index==2){#if Weibull is chosen

      fisher_info <- solve(-plcoptim_wei$hessian) #inverse -hess
      se <- sqrt(diag(fisher_info)) 
      
      #point and var
      est_a1 <- plcoptim_wei$par[1]
      est_a2 <- plcoptim_wei$par[4]
      est_x0 <- plcoptim_wei$par[2]
      est_x1 <- plcoptim_wei$par[3]
      est_y0 <- plcoptim_wei$par[5]
      est_y1 <- plcoptim_wei$par[6]
      est_b0 <- plcoptim_wei$par[7]
      est_b1 <- plcoptim_wei$par[8]
      varb0 <- fisher_info[7,7]
      varb1 <- fisher_info[8,8]
      covb0b1 <- fisher_info[7,8]
      
      #rho for age=0
      theta_d0[i] <- exp(est_b0)
      var_theta_d0 <- exp(2*est_b0)*varb0
      theta_d0_lwci[i] <- theta_d0[i] - 1.96*sqrt(var_theta_d0)
      theta_d0_upci[i] <- theta_d0[i] + 1.96*sqrt(var_theta_d0)
      
      rho_d0_cop <- claytonCopula(theta_d0[i])
      est_rho_d0 <- rho(rho_d0_cop)
      rho_d0_lwci_cop <- claytonCopula(theta_d0_lwci[i])
      rho_d0_lwci[i] <- rho(rho_d0_lwci_cop)
      rho_d0_upci_cop <- claytonCopula(theta_d0_upci[i])
      rho_d0_upci[i] <- rho(rho_d0_upci_cop) 
      save_rho_d0[i] <- est_rho_d0
      
      #rho for age=1
      theta_d1[i] <- exp(est_b0+est_b1)
      var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
      theta_d1_lwci[i] <- theta_d1[i] - 1.96*sqrt(var_theta_d1)
      theta_d1_upci[i] <- theta_d1[i] + 1.96*sqrt(var_theta_d1)
      
      rho_d1_cop <- claytonCopula(theta_d1[i])
      est_rho_d1 <- rho(rho_d1_cop)
      rho_d1_lwci_cop <- claytonCopula(theta_d1_lwci[i])
      rho_d1_lwci[i] <- rho(rho_d1_lwci_cop)
      rho_d1_upci_cop <- claytonCopula(theta_d1_upci[i])
      rho_d1_upci[i] <- rho(rho_d1_upci_cop) 
      save_rho_d1[i] <- est_rho_d1
      
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
    
      } else{# Gompertz is chosen
      #hessian
      hessian <- hessian(cpl_gom, plcoptim_gom$par, X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
      #fishers info matrix
      fisher_info <- solve(-hessian) 
      #Standard Error
      se <- sqrt(diag(fisher_info)) 
      
      est_b0 <- plcoptim_gom$par[7]
      est_b1 <- plcoptim_gom$par[8]
      est_g1 <- plcoptim_gom$par[1]
      est_g2 <- plcoptim_gom$par[4]
      est_p0 <- plcoptim_gom$par[2]
      est_p1 <- plcoptim_gom$par[3]
      est_q0 <- plcoptim_gom$par[5]
      est_q1 <- plcoptim_gom$par[6]
      varb0 <- fisher_info[7,7]
      varb1 <- fisher_info[8,8]
      covb0b1 <- fisher_info[7,8]
      
      #rho for donor=0
      theta_d0[i] <- exp(est_b0)
      var_theta_d0 <- exp(2*est_b0)*varb0
      theta_d0_lwci[i] <- theta_d0[i] - 1.96*sqrt(var_theta_d0)
      theta_d0_upci[i] <- theta_d0[i] + 1.96*sqrt(var_theta_d0)
      
      rho_d0_cop <- claytonCopula(theta_d0[i])
      est_rho_d0 <- rho(rho_d0_cop)
      rho_d0_lwci_cop <- claytonCopula(theta_d0_lwci[i])
      rho_d0_lwci[i] <- rho(rho_d0_lwci_cop)
      rho_d0_upci_cop <- claytonCopula(theta_d0_upci[i])
      rho_d0_upci[i] <- rho(rho_d0_upci_cop) 
      save_rho_d0[i] <- est_rho_d0
      
      #rho for donor=1
      theta_d1[i] <- exp(est_b0+est_b1)
      var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
      theta_d1_lwci[i] <- theta_d1[i] - 1.96*sqrt(var_theta_d1)
      theta_d1_upci[i] <- theta_d1[i] + 1.96*sqrt(var_theta_d1)
      
      rho_d1_cop <- claytonCopula(theta_d1[i])
      est_rho_d1 <- rho(rho_d1_cop)
      rho_d1_lwci_cop <- claytonCopula(theta_d1_lwci[i])
      rho_d1_lwci[i] <- rho(rho_d1_lwci_cop)
      rho_d1_upci_cop <- claytonCopula(theta_d1_upci[i])
      rho_d1_upci[i] <- rho(rho_d1_upci_cop) 
      save_rho_d1[i] <- est_rho_d1
  
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
      hr_l1_upci[i] <- est_hr_l1 + 1.96*sqrt(var_hr_l1)
      hr_l2_lwci[i] <- est_hr_l2 - 1.96*sqrt(var_hr_l2)
      hr_l2_upci[i] <- est_hr_l2 + 1.96*sqrt(var_hr_l2)
      
    }
  
  if (hr_l1_lwci[i] <1 & hr_l1_upci[i]  <1) {hr_1_lw=hr_1_lw+1
  } else if (hr_l1_lwci[i]  >1 & hr_l1_upci[i]  >1) {hr_1_up=hr_1_up+1
  } else {hr_1_cross = hr_1_cross+1}
  
  if (hr_l2_lwci[i]  <1 & hr_l2_upci[i]  <1) {hr_2_lw=hr_2_lw+1
  } else if (hr_l2_lwci[i]  >1 & hr_l2_upci[i]  >1) {hr_2_up=hr_2_up+1
  } else {hr_2_cross = hr_2_cross+1}
  
  print(i)
} # end of loop
  
  #hrs#
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
  
  # put results together and write to CSV file
  # mean of bias
  # hr_l1 represents non-terminal event; hr_l2 represents terminal event
  bias <- c(hr_l1_bias, hr_l2_bias, rho_d0_bias, rho_d1_bias, theta_d0_bias, theta_d1_bias)
  
  # Coverage probability
  CP <- c(hr_l1_cov,hr_l2_cov, rho_d0_cov, rho_d1_cov, theta_d0_cov, theta_d1_cov)
  
  MSE <- c(hr_l1_mse, hr_l2_mse, rho_d0_mse, rho_d1_mse,theta_d0_mse, theta_d1_mse)
  
  # in the order of exponential, weibull, gompertz.
  percentage_chosen = c(exp_perc, wei_perc, gom_perc, "na", "na", "na")
  
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
  print("Simulation 2 for clayton exponential model completed successfully!")
  