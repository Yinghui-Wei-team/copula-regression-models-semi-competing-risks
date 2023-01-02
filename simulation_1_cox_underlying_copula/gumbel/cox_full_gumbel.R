################################################################################
# Paper 2: Simulation 1
# Data simulated from: Gumbel copula exponential survival model 
#                      with covariates on hazard rates and association parameter
# Analysis:            Cox model is fitted
# Purpose:             Evaluation for misspecification
################################################################################
# original script by LS; edited and updated for paper2 by YW
# YW: 22 July 2021:  Data from Clayton copula, analysis by Cox model
# YW: 22 July 2021: 1. add running time tracker
#                   2. rename variables, output results
# YW: 30 Dec 2022: Tidy up
# YW: 2 Jan 2023:  change results directory
################################################################################
rm(list=ls())
library(copula);library(mvtnorm); library(plyr);library(survival);library(numDeriv)

# directory if on own PC
dir_results <- "../../results/simulation_results/simulation1/"

# # directory if on cluster
# dir2 = "home/ywei/Simulation/Paper2/Gumbel"
# 
# # set working directory
# setwd(dir2)

out_file_summary <- "S1-Cox model-data from Gumbel copula-summary.csv"
out_file_estimates <- "S1-Cox model-data from Gumbel copula-estimates.csv"
start_time = Sys.time()

########################################################
####################### set up #########################
########################################################
set.seed(541121303)
n <- 3000
runs <- 1000

#true values from KTX data
true_b0 <- -2.30
true_b1 <- 1.35
true_b2 <- 0
true_b3 <- 0

true_a0 <- -3.33
true_a1 <- 0.13
true_a2 <- 0
true_a3 <- -0.51

true_c0 <- -4.16
true_c1 <- 1.30
true_c2 <- -0.11
true_c3 <- -0.64

true_l1 <- rep(0,n)
true_l2 <- rep(0,n)
true_t <- rep(0,n)
true_r <- rep(0,n)
U1 <- rep(0,n)
V1 <- rep(0,n)

true_hr_l1_age <- exp(true_a1)
true_hr_l1_gen <- exp(true_a2)
true_hr_l1_donor <- exp(true_a3)

true_hr_l2_age <- exp(true_c1)
true_hr_l2_gen <- exp(true_c2)
true_hr_l2_donor <- exp(true_c3)

# YW defined
hr_l1_age <- hr_l1_lwci_age <- hr_l1_upci_age <- rep(0,runs)
hr_l1_gen<- hr_l1_lwci_gen <- hr_l1_upci_gen <- rep(0,runs)
hr_l1_donor <- hr_l1_lwci_donor <- hr_l1_upci_donor <- rep(0,runs)

hr_l2_age <- hr_l2_lwci_age <- hr_l2_upci_age <- rep(0,runs)
hr_l2_gen <- hr_l2_lwci_gen <- hr_l2_upci_gen <-  rep(0,runs)
hr_l2_donor <- hr_l2_lwci_donor <- hr_l2_upci_donor <- rep(0,runs)


###############################################################
###################### run 'runs' times #######################
###############################################################

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
    
    theta1 <- exp(true_b0+true_b1*age.grp[k]+true_b2*gen[k]+true_b3*donor[k])+1
    true_l1s <- exp(true_a0 + true_a1*age.grp[k] + true_a2*gen[k] + true_a3*donor[k]) 
    true_l2s <- exp(true_c0 + true_c1*age.grp[k] + true_c2*gen[k] + true_c3*donor[k])
    
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
  
  
  #############################################
  ###############    CPHM   ###################
  #############################################  
  
  ## Non-terminal event ##
  cox_l1 <- coxph(Surv(X, d1) ~ age.grp+gen+donor, data = df)
  sum_l1 <- summary(cox_l1)
  
  # This could have been output to CSV to record sum_l1 and sum_l2, and 
  # post simulation metrics can then be calculated 
  hr_l1_age[i] <- sum_l1$coefficients[4]
  hr_l1_gen[i] <- sum_l1$coefficients[5]
  hr_l1_donor[i] <- sum_l1$coefficients[6]
  
  hr_l1_lwci_age[i] <- sum_l1$conf.int[7]
  hr_l1_lwci_gen[i] <- sum_l1$conf.int[8]
  hr_l1_lwci_donor[i] <- sum_l1$conf.int[9]
  
  hr_l1_upci_age[i] <- sum_l1$conf.int[10]
  hr_l1_upci_gen[i] <- sum_l1$conf.int[11]
  hr_l1_upci_donor[i] <- sum_l1$conf.int[12]
  

  ## Terminal event ##
  cox_l2 <- coxph(Surv(Y, d2) ~ age.grp+gen+donor, data = df)
  sum_l2 <- summary(cox_l2)
  
  hr_l2_age[i] <- sum_l2$coefficients[4]
  hr_l2_gen[i] <- sum_l2$coefficients[5]
  hr_l2_donor[i] <- sum_l2$coefficients[6]
  
  hr_l2_lwci_age[i] <- sum_l2$conf.int[7]
  hr_l2_lwci_gen[i] <- sum_l2$conf.int[8]
  hr_l2_lwci_donor[i] <- sum_l2$conf.int[9]
  
  hr_l2_upci_age[i] <- sum_l2$conf.int[10]
  hr_l2_upci_gen[i] <- sum_l2$conf.int[11]
  hr_l2_upci_donor[i] <- sum_l2$conf.int[12]
  
  
  print(i)
}


#hrs#
#bias: corrected
hr_l1_bias_age <- mean(true_hr_l1_age - hr_l1_age)
hr_l2_bias_age <- mean(true_hr_l2_age - hr_l2_age)
hr_l1_bias_gen <- mean(true_hr_l1_gen - hr_l1_gen)
hr_l2_bias_gen <- mean(true_hr_l2_gen - hr_l2_gen)
hr_l1_bias_donor <- mean(true_hr_l1_donor - hr_l1_donor)
hr_l2_bias_donor <- mean(true_hr_l2_donor - hr_l2_donor)

#coverage: YW revised

hr_l1_cov_age <- 100* sum(true_hr_l1_age <= hr_l1_upci_age & true_hr_l1_age >= hr_l1_lwci_age)/runs
hr_l2_cov_age <- 100* sum(true_hr_l2_age <= hr_l2_upci_age & true_hr_l2_age >= hr_l2_lwci_age)/runs

hr_l1_cov_gen <- 100* sum(true_hr_l1_gen <= hr_l1_upci_gen & true_hr_l1_gen >= hr_l1_lwci_gen)/runs
hr_l2_cov_gen <- 100* sum(true_hr_l2_gen <= hr_l2_upci_gen & true_hr_l2_gen >= hr_l2_lwci_gen)/runs

hr_l1_cov_donor <- 100* sum(true_hr_l1_donor <= hr_l1_upci_donor & true_hr_l1_donor >= hr_l1_lwci_donor)/runs
hr_l2_cov_donor <- 100* sum(true_hr_l2_donor <= hr_l2_upci_donor & true_hr_l2_donor >= hr_l2_lwci_donor)/runs

#mse: corrected
hr_l1_mse_age <- mean((true_hr_l1_age - hr_l1_age)^2)
hr_l2_mse_age <- mean((true_hr_l2_age - hr_l2_age)^2)
hr_l1_mse_gen <- mean((true_hr_l1_gen - hr_l1_gen)^2)
hr_l2_mse_gen <- mean((true_hr_l2_gen - hr_l2_gen)^2)
hr_l1_mse_donor <- mean((true_hr_l1_donor - hr_l1_donor)^2)
hr_l2_mse_donor <- mean((true_hr_l2_donor - hr_l2_donor)^2)


end_time = Sys.time()
run_time = end_time - start_time
run_time

# YW 23 June 2021: put results together and write to CSV file
# mean of bias
bias <- c(hr_l1_bias_age, hr_l1_bias_donor, hr_l1_bias_gen,
          hr_l2_bias_age, hr_l2_bias_donor, hr_l2_bias_gen)

# Coverage probability
CP <- c(hr_l1_cov_age, hr_l1_cov_donor, hr_l1_cov_gen,
        hr_l2_cov_age, hr_l2_cov_donor, hr_l2_cov_gen)

MSE <- c(hr_l1_mse_age, hr_l1_mse_donor, hr_l1_mse_gen,
         hr_l2_mse_age, hr_l2_mse_donor, hr_l2_mse_gen)

# YW: put results together
items<-c("NT_age", "NT_donor", "NT_gen",
         "T_age", "T_donor", "T_gen")
Results <- cbind.data.frame(items, bias, CP, MSE)

Results[,2:4] <- round(Results[,2:4],3)
Results

rownames(Results)<-NULL
end_time <- Sys.time()
run_time = end_time - start_time
run_time

Estimates = data.frame(hr.l1.age.est = hr_l1_age, hr.l1.age.lci = hr_l1_lwci_age, hr.l1.age.uci=hr_l1_upci_age,
                       hr.l2.age.est = hr_l2_age, hr.l2.age.lci = hr_l2_lwci_age, hr.l2.age.uci=hr_l2_upci_age,
                       hr.l1.gen.est = hr_l1_gen, hr.l1.gen.lci = hr_l1_lwci_gen, hr.l1.gen.uci=hr_l1_upci_gen,
                       hr.l2.gen.est = hr_l2_gen, hr.l2.gen.lci = hr_l2_lwci_gen, hr.l2.gen.uci=hr_l2_upci_gen,
                       hr.l1.donor.est = hr_l1_donor, hr.l1.donor.lci = hr_l1_lwci_donor, hr.l1.gen.uci=hr_l1_upci_donor,
                       hr.l2.donor.est = hr_l2_donor, hr.l2.donor.lci = hr_l2_lwci_donor, hr.l2.gen.uci=hr_l2_upci_donor)


# Output estimates from each replication and summary of results -----------------
write.csv(Results, row.names=F,file=paste0(dir_results, out_file_summary))
write.csv(Estimates, row.names=F,file=paste0(dir_results, out_file_estimates))
print(run_time)
print("Simulation1 cox model for gumbel exponential data completed successfully!")

