################################################################################
# Paper 2: Simulation 1 Model 1
# Model 1 - Clayton copula exponential survival model with covariates on hazard rates
# Original code by LS, reviewed and updated by MW and YW
# YW 25 June 2021 updates: 
# 1.correction for the variance post simulation
# 2.add code to output results into a CSV file
# 3.output simulation time
# MW 14 July updates:
# 1. Correction for the formula for bias.
# 2. Change for the formula for MSE.
# 3. Changed starting values.
# 4. All results saved in a csv file.
# 5. Results for b0 collected
# MW 25 July 2021 updates:
# 1. CI's and coverage calculated outside of the loop.
# 2. Std errors saved for each iteration in csv file together with estimates.
# 3. Run time saved in csv file.
# 4. Renamed b0 to theta for estimation.
# YW 2 January 2023 updates
# 1. Move likelihood function out of the loop
# 2. Move lower and upper bounds out of the loop
# 3. reset results directory, 
#    specify out files, n and runs nearer to the front of the script
###############################################################################

rm(list=ls())
library(copula); library(mvtnorm); library(plyr); library(survival); library(numDeriv)

########################################################
####################### set up #########################
########################################################

# results directory
# directory if on own PC
dir_results <- "../../"
dir = paste0(dir_results, "results/simulation_results")

# directory if working on cluster
# dir = "/home/ywei/Simulation/Paper2/Clayton"
# setwd(dir)

# Number of participants and number of replications of the simulation
n = 3000    # number of participants
runs = 1000    # number of replications

# set outfile name
out_file_estimates <- "S1-Table4-clayton-exponential-covariates-hazards.csv"
out_file_summary <- "S1-Clayton-exponential-covariates-hazards.csv"

# Starting values for parameter estimation for optim():
start_a0 <- 0.5; start_a1 <- 0.5; start_a2 <- 0.5; start_a3 <- 0.5;
start_c0 <- 0.5; start_c1 <- 0.5; start_c2 <- 0.5; start_c3 <- 0.5;
start_t <- 0.5;
starting_values <- c(start_a0, start_a1, start_a2, start_a3, 
                     start_c0, start_c1, start_c2, start_c3, 
                     start_t)
# The last parameter renamed to t (theta), and not b0.

# Lower and upper bounds for regression coefficients
# a0,a1,a2, a3,c0,c1,c2,c3,t
a0_lw <- -10; a0_up <- -2; a1_lw <- -10; a1_up <- 1; a2_lw <- -10; a2_up <- 1; a3_lw <- -10; a3_up <- 1
c0_lw <- -10; c0_up <- -2; c1_lw <- -10; c1_up <- 3;c2_lw <- -10; c2_up <- 1;c3_lw <- -10; c3_up <- 1
t_lw <- 0.01; t_up <- 12
reg_coef_lw <- c(a0_lw,a1_lw,a2_lw, a3_lw,c0_lw,c1_lw,c2_lw, c3_lw, t_lw)
reg_coef_up <- c(a0_up,a1_up,a2_up, a3_up, c0_up,c1_up,c2_up,c3_up, t_up)

# Note: True values are set again inside the function simulation.
true_a0 <- -3.28; true_a1 <- 0.32; true_a2 <- 0; true_a3 <- -0.53; 
true_c0 <- -4.09; true_c1 <- 1.35; true_c2 <- -0.07; true_c3 <- -0.62
true_b0 <- 0.39; true_b1 <- 1.09; true_b2 <- 0.14; true_b3 <- 0.53
true_values <- c(true_a0, true_a1, true_a2, true_a3, 
                 true_c0, true_c1, true_c2, true_c3, true_b0)

##########################################################
############### Clayton pseudo likelihood ################
##########################################################
cpl<-function(para, X, Y, d1, d2, age.grp, gen, donor){
  
  a0 <- para[1]; a1 <- para[2]; a2 <- para[3];  a3 <- para[4]
  
  c0 <- para[5]; c1 <- para[6]; c2 <- para[7]; c3 <- para[8]
  
  theta <- para[9]      
  
  lambda1 <- exp(a0+a1*age.grp+a2*gen+a3*donor)
  lambda2 <- exp(c0+c1*age.grp+c2*gen+c3*donor)
  
  S1 <- exp(-lambda1*X)
  S2 <- exp(-lambda2*Y) #problems when S2 too small -> 0 
  f1 <- lambda1*exp(-lambda1*X)
  f2 <- lambda2*exp(-lambda2*Y)
  
  C=(S1^(-theta)+S2^(-theta)-1)^(-1/theta)
  
  C[which(C<0.1^(8))]=0.1^(8)
  S1[which(S1 < 0.1^(8))]=0.1^(8)
  S2[which(S2 < 0.1^(8))]=0.1^(8)
  
  part1 <- d1*d2*(log(1+theta)+(1+2*theta)*log(C)-(theta+1)*log(S1)-(theta+1)*log(S2)+log(lambda1)-lambda1*X+log(lambda2)-lambda2*Y)
  part2 <- d1*(1-d2)*((theta+1)*log(C)-(theta+1)*log(S1)+log(lambda1)-lambda1*X)
  part3<-((1-d1)*(d2))*((theta+1)*log(C)-(theta+1)*log(S2)+log(lambda2)-lambda2*Y)
  part4<-((1-d1)*(1-d2))*log(C)    
  
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}


simulation <- function(runs, n, starting_values, reg_coef_lw, reg_coef_up, out_file_summary, out_file_estimates){
  
  # set up ---------------------------------------------------------------------
  start_time <- Sys.time()
  set.seed(98452221)
  
  #true values from KTX data --------------------------------------------------
  true_b0 <- 0.39
  true_b1 <- 1.09
  true_b2 <- 0.14
  true_b3 <- 0.53
  
  true_a0 <- -3.28
  true_a1 <- 0.32
  true_a2 <- 0
  true_a3 <- -0.53
  
  true_c0 <- -4.09
  true_c1 <- 1.35
  true_c2 <- -0.07
  true_c3 <- -0.62
  
  true_l1 <- true_l2 <- rep(0,n)
  true_t <- rep(0,n)
  U1 <- V1 <- rep(0,n)
  
  true_hr_l1_age <- exp(true_a1)
  true_hr_l1_gen <- exp(true_a2)
  true_hr_l1_donor <- exp(true_a3)
  
  true_hr_l2_age <- exp(true_c1)
  true_hr_l2_gen <- exp(true_c2)
  true_hr_l2_donor <- exp(true_c3)
  
  ## stuff for later ##
  save_a0 <- save_a1 <- save_a2 <- save_a3 <- rep(0,runs)
  save_se_a0 <- save_se_a1 <- save_se_a2 <- save_se_a3 <- rep(0,runs)
  save_c0 <- save_c1 <- save_c2 <- save_c3 <- rep(0,runs)
  save_se_c0 <- save_se_c1 <- save_se_c2 <- save_se_c3 <- rep(0,runs)
  save_theta <- rep(0, runs)     # This renamed to save_theta from save_b0
  save_se_theta <- rep(0, runs)  
  
  #replicate run 'runs' times -------------------------------------------------
  
  for (i in 1:runs){
    
    ###############################################################
    ######################## generate data ########################
    ###############################################################
    
    #Step 1: generate age categories
    age.grp <- rbinom(n,1,0.40)          #40% are in the older age group in NHSBT data
    donor <- rbinom(n,1,0.30)
    gen <- rbinom(n,1,0.38)
    
    for(k in 1:(n)){   #loop to generate U an V from age-varying theta
      m=1                  
      
      #Step 2: generate 1 random variable from Uniform(0,a) distribution 
      u1 <- runif(m,0,1)       
      
      #Step 3: X_true generated from u1 values (T1 from later)
      theta1 <- exp(true_b0+true_b1*age.grp[k]+true_b2*gen[k]+true_b3*donor[k])
      true_l1s <- exp(true_a0 + true_a1*age.grp[k] + true_a2*gen[k] + true_a3*donor[k]) 
      true_l2s <- exp(true_c0 + true_c1*age.grp[k] + true_c2*gen[k] + true_c3*donor[k])
      
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
      true_l1[k] <- true_l1s
      true_l2[k] <- true_l2s
    }
    
    #Step 4: T1 and T2 from inverse exponential CDF (F^{-1}=-log(1-U)/lambda_1), could use qexp(u1)
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
    df<-data.frame(X, Y, d1, d2, age.grp, gen, donor)
    
    
    plcoptim <- optim(starting_values, cpl, method="L-BFGS-B", 
                      lower=reg_coef_lw,
                      upper=reg_coef_up, 
                      X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen,
                      donor=df$donor, control=list(fnscale=-1), hessian=TRUE)
  
    #############################################################
    ################## Estimates and std errors ################
    ############################################################
    
    hess <- hessian(cpl, plcoptim$par, X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
    
    fisher_info <- solve(-hess)
    #fisher_info <- solve(-plcoptim$hessian) #inverse -hess
    #Standard error = sqrt(var/n)
    
    se <- sqrt(diag(fisher_info)) 
  
    #a
    save_a0[i] <- plcoptim$par[1]
    save_a1[i] <- plcoptim$par[2]
    save_a2[i] <- plcoptim$par[3]
    save_a3[i] <- plcoptim$par[4]
    
    save_se_a0[i] <- se[1]
    save_se_a1[i] <- se[2]
    save_se_a2[i] <- se[3]
    save_se_a3[i] <- se[4]
    
    #c
    save_c0[i] <- plcoptim$par[5]
    save_c1[i] <- plcoptim$par[6]
    save_c2[i] <- plcoptim$par[7]
    save_c3[i] <- plcoptim$par[8]
    
    save_se_c0[i] <- se[5]
    save_se_c1[i] <- se[6]
    save_se_c2[i] <- se[7]
    save_se_c3[i] <- se[8]
    
    # t
    save_theta[i] <- plcoptim$par[9]
    save_se_theta[i] <- se[9]
    
    print(i)
  }
  
  # END OF LOOP

  # HR's #
  save_hr_l1_age <- exp(save_a1)
  save_hr_l1_gen <- exp(save_a2)
  save_hr_l1_donor <- exp(save_a3)
  
  save_hr_l2_age <- exp(save_c1)
  save_hr_l2_gen <- exp(save_c2)
  save_hr_l2_donor <- exp(save_c3)
  
  save_se_hr_l1_age <- save_hr_l1_age * save_se_a1
  save_se_hr_l1_gen <- save_hr_l1_gen * save_se_a2  
  save_se_hr_l1_donor <- save_hr_l1_donor * save_se_a3
  
  save_se_hr_l2_age <- save_hr_l2_age * save_se_c1
  save_se_hr_l2_gen <- save_hr_l2_gen * save_se_c2
  save_se_hr_l2_donor <- save_hr_l2_donor * save_se_c3

  #bias
  hr_l1_bias_age <- mean(save_hr_l1_age) - true_hr_l1_age
  hr_l2_bias_age <- mean(save_hr_l2_age) - true_hr_l2_age
  hr_l1_bias_gen <- mean(save_hr_l1_gen) - true_hr_l1_gen
  hr_l2_bias_gen <- mean(save_hr_l2_gen) - true_hr_l2_gen
  hr_l1_bias_donor <- mean(save_hr_l1_donor) - true_hr_l1_donor
  hr_l2_bias_donor <- mean(save_hr_l2_donor) - true_hr_l2_donor
  
  # CI's
  hr_l1_lwci_age <- save_hr_l1_age - 1.96*save_se_hr_l1_age
  hr_l1_upci_age <- save_hr_l1_age + 1.96*save_se_hr_l1_age
  hr_l1_lwci_gen <- save_hr_l1_gen - 1.96*save_se_hr_l1_gen
  hr_l1_upci_gen <- save_hr_l1_gen + 1.96*save_se_hr_l1_gen
  hr_l1_lwci_donor <- save_hr_l1_donor - 1.96*save_se_hr_l1_donor
  hr_l1_upci_donor <- save_hr_l1_donor + 1.96*save_se_hr_l1_donor
  
  hr_l2_lwci_age <- save_hr_l2_age - 1.96*save_se_hr_l2_age
  hr_l2_upci_age <- save_hr_l2_age + 1.96*save_se_hr_l2_age
  hr_l2_lwci_gen <- save_hr_l2_gen - 1.96*save_se_hr_l2_gen
  hr_l2_upci_gen <- save_hr_l2_gen + 1.96*save_se_hr_l2_gen
  hr_l2_lwci_donor <- save_hr_l2_donor - 1.96*save_se_hr_l2_donor
  hr_l2_upci_donor <- save_hr_l2_donor + 1.96*save_se_hr_l2_donor
  
  #coverage
  counter_hr_l1_age <- sum(true_hr_l1_age <= hr_l1_upci_age & true_hr_l1_age >= hr_l1_lwci_age)
  hr_l1_cov_age <- (counter_hr_l1_age / runs) * 100
  
  counter_hr_l2_age <- sum(true_hr_l2_age <= hr_l2_upci_age & true_hr_l2_age >= hr_l2_lwci_age)
  hr_l2_cov_age <- (counter_hr_l2_age / runs) * 100
  
  counter_hr_l1_gen <- sum(true_hr_l1_gen <= hr_l1_upci_gen & true_hr_l1_gen >= hr_l1_lwci_gen)
  hr_l1_cov_gen <- (counter_hr_l1_gen / runs) * 100
  
  counter_hr_l2_gen <- sum(true_hr_l2_gen <= hr_l2_upci_gen & true_hr_l2_gen >= hr_l2_lwci_gen)
  hr_l2_cov_gen <- (counter_hr_l2_gen / runs) * 100
  
  counter_hr_l1_donor <- sum(true_hr_l1_donor <= hr_l1_upci_donor & true_hr_l1_donor >= hr_l1_lwci_donor)
  hr_l1_cov_donor <- (counter_hr_l1_donor / runs) * 100
  
  counter_hr_l2_donor <- sum(true_hr_l2_donor <= hr_l2_upci_donor & true_hr_l2_donor >= hr_l2_lwci_donor)
  hr_l2_cov_donor <- (counter_hr_l2_donor / runs) * 100
  
  #mse
  hr_l1_mse_age <- mean((save_hr_l1_age - true_hr_l1_age)^2)
  hr_l2_mse_age <- mean((save_hr_l2_age - true_hr_l2_age)^2)
  hr_l1_mse_gen <- mean((save_hr_l1_gen - true_hr_l1_gen)^2)
  hr_l2_mse_gen <- mean((save_hr_l2_gen - true_hr_l2_gen)^2)
  hr_l1_mse_donor <- mean((save_hr_l1_donor - true_hr_l1_donor)^2)
  hr_l2_mse_donor <- mean((save_hr_l2_donor - true_hr_l2_donor)^2)
  
  #a#
  #bias
  a0_bias <- mean(save_a0) - true_a0
  a1_bias <- mean(save_a1) - true_a1
  a2_bias <- mean(save_a2) - true_a2
  a3_bias <- mean(save_a3) - true_a3
  
  # CI's
  uci_a0 <- save_a0 + 1.96*save_se_a0 
  lci_a0 <- save_a0 - 1.96*save_se_a0
  uci_a1 <- save_a1 + 1.96*save_se_a1
  lci_a1 <- save_a1 - 1.96*save_se_a1
  uci_a2 <- save_a2 + 1.96*save_se_a2
  lci_a2 <- save_a2 - 1.96*save_se_a2
  uci_a3 <- save_a3 + 1.96*save_se_a3
  lci_a3 <- save_a3 - 1.96*save_se_a3
  
  #coverage
  counter_a0 <- sum(true_a0 <= uci_a0 & true_a0 >= lci_a0)
  a0_cov <- (counter_a0 / runs) * 100
  
  counter_a1 <- sum(true_a1 <= uci_a1 & true_a1 >= lci_a1)
  a1_cov <- (counter_a1 / runs) * 100
  
  counter_a2 <- sum(true_a2 <= uci_a2 & true_a2 >= lci_a2)
  a2_cov <- (counter_a2 / runs) * 100
  
  counter_a3 <- sum(true_a3 <= uci_a3 & true_a3 >= lci_a3)
  a3_cov <- (counter_a3 / runs) * 100
  
  #mse
  a0_mse <- mean((save_a0 - true_a0)^2)
  a1_mse <- mean((save_a1 - true_a1)^2)
  a2_mse <- mean((save_a2 - true_a2)^2)
  a3_mse <- mean((save_a3 - true_a3)^2)
  
  #b0#
  b0_bias <- mean(log(save_theta)) - true_b0    # ?? This might have to be modified
  # depending on whether we are treating theta as an estimator of exp(true_b0) or as an estimator of true theta.
  
  # CI, using delta method
  uci_b0 <- log(save_theta) + 1.96*save_se_theta*(save_theta^(-1))   
  lci_b0 <- log(save_theta) - 1.96*save_se_theta*(save_theta^(-1))
  
  # coverage
  counter_b0 <- sum(true_b0 <= uci_b0 & true_b0 >= lci_b0)  # ?? This might have to be modified
  # depending on whether we are treating theta as an estimator of exp(true_b0) or as an estimator of true theta.
  b0_cov <- (counter_b0 / runs) * 100
  
  # mse
  b0_mse <- mean((log(save_theta) - true_b0)^2)
  
  #cs#
  #bias
  c0_bias <- mean(save_c0) - true_c0
  c1_bias <- mean(save_c1) - true_c1
  c2_bias <- mean(save_c2) - true_c2
  c3_bias <- mean(save_c3) - true_c3
  
  # CI
  uci_c0 <- save_c0 + 1.96*save_se_c0 
  lci_c0 <- save_c0 - 1.96*save_se_c0
  uci_c1 <- save_c1 + 1.96*save_se_c1
  lci_c1 <- save_c1 - 1.96*save_se_c1
  uci_c2 <- save_c2 + 1.96*save_se_c2
  lci_c2 <- save_c2 - 1.96*save_se_c2
  uci_c3 <- save_c3 + 1.96*save_se_c3
  lci_c3 <- save_c3 - 1.96*save_se_c3
  
  #coverage
  counter_c0 <- sum(true_c0 <= uci_c0 & true_c0 >= lci_c0)
  c0_cov <- (counter_c0 / runs) * 100
  
  counter_c1 <- sum(true_c1 <= uci_c1 & true_c1 >= lci_c1)
  c1_cov <- (counter_c1 / runs) * 100
  
  counter_c2 <- sum(true_c2 <= uci_c2 & true_c2 >= lci_c2)
  c2_cov <- (counter_c2 / runs) * 100
  
  counter_c3 <- sum(true_c3 <= uci_c3 & true_c3 >= lci_c3)
  c3_cov <- (counter_c3 / runs) * 100
  
  #mse
  c0_mse <- mean((save_c0 - true_c0)^2)
  c1_mse <- mean((save_c1 - true_c1)^2)
  c2_mse <- mean((save_c2 - true_c2)^2)
  c3_mse <- mean((save_c3 - true_c3)^2)
  ###############################
  
  counter_a0_low <- sum(save_a0 == a0_lw)
  counter_a0_upper <- sum(save_a0 == a0_up)
  counter_a1_low <- sum(save_a1 == a1_lw)
  counter_a1_upper <- sum(save_a1 == a1_up)
  counter_a2_low <- sum(save_a2 == a2_lw)
  counter_a2_upper <- sum(save_a2 == a2_up)
  counter_a3_low <- sum(save_a3 == a3_lw)
  counter_a3_upper <- sum(save_a3 == a3_up)
  
  counter_c0_low <- sum(save_c0 == c0_lw)
  counter_c0_upper <- sum(save_c0 == c0_up)
  counter_c1_low <- sum(save_c1 == c1_lw)
  counter_c1_upper <- sum(save_c1 == c1_up)
  counter_c2_low <- sum(save_c2 == c2_lw)
  counter_c2_upper <- sum(save_c2 == c2_up)
  counter_c3_low <- sum(save_c3 == c3_lw)
  counter_c3_upper <- sum(save_c3 == c3_up)
  
  counter_t_low <- sum(save_theta == t_lw) 
  counter_t_upper <- sum(save_theta == t_up)
  
  counter_low <- c(counter_a0_low, counter_a1_low, counter_a2_low, counter_a3_low,
                   counter_c0_low, counter_c1_low, counter_c2_low, counter_c3_upper, 
                   counter_t_low)
  
  counter_upper <- c(counter_a0_upper, counter_a1_upper, counter_a2_upper, counter_a3_upper,
                     counter_c0_upper, counter_c1_upper, counter_c2_upper, counter_c3_upper, 
                     counter_t_upper)
  
  counter_hit_boundaries <- data.frame(rbind(counter_low, counter_upper))
  names(counter_hit_boundaries) <- c("a0", "a1", "a2", "a3",
                                     "c0", "c1", "c2", "c3", "t")
  row.names(counter_hit_boundaries) <- c("lower", "upper")
  
  print(counter_hit_boundaries)
  
  # put results together and write to CSV file
  # mean of bias
  # hr_l1 represents non-terminal event; hr_l2 represents terminal event
  mean_bias <- c(a0_bias, a1_bias, a2_bias, a3_bias, 
                 b0_bias,
                 c0_bias, c1_bias, c2_bias, c3_bias,
                 hr_l1_bias_age, hr_l1_bias_donor, hr_l1_bias_gen,
                 hr_l2_bias_age, hr_l2_bias_donor, hr_l2_bias_gen)
  
  # Coverage probability
  CP <- c(a0_cov, a1_cov, a2_cov, a3_cov,
          b0_cov,
          c0_cov, c1_cov, c2_cov, c3_cov,
          hr_l1_cov_age, hr_l1_cov_donor, hr_l1_cov_gen,
          hr_l2_cov_age, hr_l2_cov_donor, hr_l2_cov_gen)
  
  MSE <- c(a0_mse, a1_mse, a2_mse, a3_mse,
           b0_mse,
           c0_mse, c1_mse, c2_mse, c3_mse,
           hr_l1_mse_age, hr_l1_mse_donor, hr_l1_mse_gen,
           hr_l2_mse_age, hr_l2_mse_donor, hr_l2_mse_gen)
  
  # YW: put results together
  items<-c("a0", "a1", "a2", "a3",
           "b0",
           "c0", "c1", "c2", "c3", 
           "NT_age", "NT_donor", "NT_gen",
           "T_age", "T_donor", "T_gen")
  
  end_time <- Sys.time()
  run_time = end_time - start_time

  Results <- cbind.data.frame(items, mean_bias, CP, MSE)
  
  #Results[,2:4] <- round(Results[,2:4],3)
  rownames(Results)<-NULL
  
  print(Results)
  print(paste("Run time", run_time))

  # SAVE IN A FILE THE ENTIRE VECTORS OF ESTIMATES AND SE'S FOR EACH ITERATION:
  Estimates = data.frame(a0.est = save_a0, se.a0.est = save_se_a0, 
                         a1.est = save_a1, se.a1.est = save_se_a1,
                         a2.est = save_a2, se.a2.est = save_se_a2,
                         a3.est = save_a3, se.a3.est = save_se_a3, 
                         c0.est = save_c0, se.c0.est = save_se_c0,
                         c1.est = save_c1, se.c1.est = save_se_c1,
                         c2.est = save_c2, se.c2.est = save_se_c2, 
                         c3.est = save_c3, se.c3.est = save_se_c3, 
                         theta.est = save_theta, se.theta.est = save_se_theta,
                         hr.l1.age.est = save_hr_l1_age, se.hr.l1.age.est = save_se_hr_l1_age,
                         hr.l2.age.est = save_hr_l2_age, se.hr.l2.age.est = save_se_hr_l2_age,
                         hr.l1.gen.est = save_hr_l1_gen, se.hr.l1.gen.est = save_se_hr_l1_gen,
                         hr.l2.gen.est = save_hr_l2_gen, se.hr.l2.gen.est = save_se_hr_l2_gen,
                         hr.l1.donor.est = save_hr_l1_donor, se.hr.l1.donor.est = save_se_hr_l1_donor,
                         hr.l2.donor.est = save_hr_l2_donor, se.hr.l2.donor.est = save_se_hr_l2_donor)
  
  write.csv(Results, row.names = F, file = out_file_summary)
  cat("\n\n", file = out_file_summary, append = TRUE)
  cat("Run time\n", file = out_file_summary, append = TRUE)
  cat(run_time, file = out_file_summary, append = TRUE)
  write.csv(Estimates, row.names = F, file=out_file_estimates)
} # END OF FUNCTION simulation

simulation(runs = runs, n= n, starting_values = starting_values, 
           reg_coef_lw = reg_coef_lw, reg_coef_up = reg_coef_up,
           out_file_summary = out_file_summary,
           out_file_estimates = out_file_estimates)
print("Simulation1 model1 for clayton exponential model completed successfully!")