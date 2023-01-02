################################################################################
# Paper 2: Simulation 1 Model 2
# Data simulated from: Model 2- Clayton copula exponential survival model with 
#                      covariates on hazard rates and the association parameters
# Fitted model:        The underlying clayton model as above
# Purpose:             Evaluating performance when the true model is specified
################################################################################
# Original code by LS, reviewed and updated by YW
# YW 26 June 2021 updates: 
# 1.corrected variance and variance post simulation
# 2.added code to output results into a CSV file
# 3.output simulation time
# YW: 11 July 2021 updates
# 4. Change starting values in optim() to be in the middle of range allowed for the 
#    respective parameters instead of starting from the true values
# 5. correct the calculation for bias 
# 6. Change the formula for MSE
# YW: 2 January 2023 updates:
# 1. Put likelihood, starting values, lower and upper bounds outside the loop
# 2. reset results directory and tidy up
# 3. move calculation of bias, coverage outside loop
################################################################################

rm(list=ls())
library(copula); library(mvtnorm);library(numDeriv)

################################################################################
#                                set up                                        #
################################################################################
# directory if on own PC
dir_results <- "../../"
dir = paste0(dir_results, "results/simulation_results")

# directory if working on cluster
# dir = "/home/ywei/Simulation/Paper2/Clayton"
# setwd(dir)

# set up outfile names
out_file_results <- "S2-Table5-clayton-exponential-covariates-hazards-association.csv"
out_file_estimates <- "S2-Clayton-exponential-covariates-for-hazards-estimates.csv"

start_time <- Sys.time()
set.seed(90089811)
n <- 3000
runs <- 10

# starting values for optim ----------------------------------------------------
# a0, a1, a2, a3, c0, c1, c2, c3, b0, b1, b2, b3
starting_values = c(-6, -5, -5, -5, -6, -5, -5, -5,
                    2, 3, 5, 7)

# lower and upper bound for parameters -----------------------------------------
a0_lw <- -10; a0_up <- -2
a1_lw <- -10; a1_up <- 1
a2_lw <- -10; a2_up <- 1
a3_lw <- -10; a3_up <- 1
c0_lw <- -10; c0_up <- -2
c1_lw <- -10; c1_up <- 3
c2_lw <- -10; c2_up <- 1
c3_lw <- -10; c3_up <- 1

b0_lw <- -5; b0_up <- 10
b1_lw <- -5; b1_up <- 12
b2_lw <- -10; b2_up <- 20
b3_lw <- -5; b3_up <- 20

#true values from KTX data
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

## stuff for later ##
save_a0 <- save_a1 <- save_a2 <- save_a3 <- rep(0,runs)
save_b0 <- save_b1 <- save_b2 <- save_b3 <- rep(0,runs)
save_c0 <- save_c1 <- save_c2 <- save_c3 <- rep(0,runs)
save_hr_l1_age <- save_hr_l2_age <- save_hr_l1_gen <- save_hr_l2_gen <- rep(0,runs)
save_hr_l1_donor <- save_hr_l2_donor <- rep(0,runs)

save_se_a0 <- save_se_a1 <- save_se_a2 <- save_se_a3 <- rep(0,runs)
save_se_c0 <- save_se_c1 <- save_se_c2 <- save_se_c3 <- rep(0,runs)
save_se_b0 <- save_se_b1 <- save_se_b2 <- save_se_b3 <- rep(0,runs)

bias_a0 <- bias_a1 <- bias_a2 <- bias_a3 <- rep(0,runs)
bias_b0 <- bias_b1 <- bias_b2 <- bias_b3 <- rep(0,runs)
bias_c0 <- bias_c1 <- bias_c2 <- bias_c3 <- rep(0,runs)
bias_l1_hr_age <- bias_l1_hr_age <- bias_l1_hr_gen <- bias_l1_hr_gen <- rep(0,runs)
bias_l1_hr_donor <- bias_l1_hr_donor <- rep(0,runs)
bias_l2_hr_age <- bias_l2_hr_age <- rep(0,runs)
bias_l2_hr_gen <- bias_l2_hr_gen <- rep(0,runs)
bias_l2_hr_donor <- bias_l2_hr_donor <- rep(0,runs)

counter_a0 = counter_a1 = counter_a2 = counter_a3 = 0
counter_c0 = counter_c1 = 0
counter_c2 = 0
counter_c3 = 0
counter_b0 = 0
counter_b1 = 0
counter_b2 = 0
counter_b3 = 0
counter_hr_l1_age = 0
counter_hr_l1_gen = 0
counter_hr_l1_donor = 0
counter_hr_l2_age = 0
counter_hr_l2_gen = 0
counter_hr_l2_donor = 0 

counter_a0_low = 0
counter_a1_low = 0
counter_a2_low = 0
counter_a3_low = 0
counter_c0_low = 0
counter_c1_low = 0
counter_c2_low = 0
counter_c3_low = 0
counter_b0_low = 0
counter_b1_low = 0
counter_b2_low = 0
counter_b3_low = 0
counter_a0_upper = 0
counter_a1_upper = 0 
counter_a2_upper = 0
counter_a3_upper = 0 
counter_c0_upper = 0
counter_c1_upper = 0
counter_c2_upper = 0
counter_c3_upper = 0 
counter_b0_upper = 0
counter_b1_upper = 0
counter_b2_upper = 0
counter_b3_upper = 0 

##########################################################
############### Clayton pseudo likelihood ################
##########################################################
cpl<-function(para, X, Y, d1, d2, age.grp, gen, donor){
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
  
  S1 <- exp(-lambda1*X)
  S2 <- exp(-lambda2*Y) #problems when S2 too small -> 0 
  f1 <- lambda1*exp(-lambda1*X)
  f2 <- lambda2*exp(-lambda2*Y)
  
  theta <- exp(b0+b1*age.grp+b2*gen+b3*donor)
  
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
                    lower=c(a0_lw,a1_lw,a2_lw, a3_lw,c0_lw,c1_lw,c2_lw, c3_lw, b0_lw,b1_lw, b2_lw, b3_lw),
                    upper=c(a0_up,a1_up,a2_up, a3_up, c0_up,c1_up,c2_up,c3_up, b0_up,b1_up,b2_up, b3_up), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen,
                    donor=df$donor, control=list(fnscale=-1),hessian=TRUE)
  
  plcoptim$par
  
  cpl(c(a0_up,a1_up,a2_up, a3_up, c0_up,c1_up,c2_up,c3_up, b0_up,b1_up,b2_up, b3_up), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
  
  if(plcoptim$par[1] == a0_lw) {counter_a0_low = counter_a0_low + 1}
  if(plcoptim$par[1] == a0_up) {counter_a0_upper = counter_a0_upper + 1}
  if(plcoptim$par[2] == a1_lw) {counter_a1_low = counter_a1_low + 1}
  if(plcoptim$par[2] == a1_up) {counter_a1_upper = counter_a1_upper + 1}
  if(plcoptim$par[3] == a2_lw) {counter_a2_low = counter_a2_low + 1}
  if(plcoptim$par[3] == a2_up) {counter_a2_upper = counter_a2_upper + 1}
  if(plcoptim$par[4] == a3_lw) {counter_a3_low = counter_a3_low + 1}
  if(plcoptim$par[4] == a3_up) {counter_a3_upper = counter_a3_upper + 1}
  
  if(plcoptim$par[5] == c0_lw) {counter_c0_low = counter_c0_low + 1}
  if(plcoptim$par[5] == c0_up) {counter_c0_upper = counter_c0_upper + 1}
  if(plcoptim$par[6] == c1_lw) {counter_c1_low = counter_c1_low + 1}
  if(plcoptim$par[6] == c1_up) {counter_c1_upper = counter_c1_upper + 1}
  if(plcoptim$par[7] == c2_lw) {counter_c2_low = counter_c2_low + 1}
  if(plcoptim$par[7] == c2_up) {counter_c2_upper = counter_c2_upper + 1}
  if(plcoptim$par[8] == c3_lw) {counter_c3_low = counter_c3_low + 1}
  if(plcoptim$par[8] == c3_up) {counter_c3_upper = counter_c3_upper + 1}
  
  if(plcoptim$par[9] == b0_lw) {counter_b0_low = counter_b0_low + 1}
  if(plcoptim$par[9] == b0_up) {counter_b0_upper = counter_b0_upper + 1}
  if(plcoptim$par[10] == b1_lw) {counter_b1_low = counter_b1_low + 1}
  if(plcoptim$par[10] == b1_up) {counter_b1_upper = counter_b1_upper + 1}
  if(plcoptim$par[11] == b2_lw) {counter_b2_low = counter_b2_low + 1}
  if(plcoptim$par[11] == b2_up) {counter_b2_upper = counter_b2_upper + 1}
  if(plcoptim$par[12] == b3_lw) {counter_b3_low = counter_b3_low + 1}
  if(plcoptim$par[12] == b3_up) {counter_b3_upper = counter_b3_upper + 1}
  
  ########################################################
  ################## Confidence Intervals ################
  ########################################################
  hess <- hessian(cpl, plcoptim$par, X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
  
  fisher_info <- solve(-hess)
  #fisher_info <- solve(-plcoptim$hessian) #inverse -hess
  #Standard error = sqrt(var/n)
  se<-sqrt(diag(fisher_info)) 
  
  #a ci
  a0_est <- plcoptim$par[1]
  a1_est <- plcoptim$par[2]
  a2_est <- plcoptim$par[3]
  a3_est <- plcoptim$par[4]
  save_a0[i] <- a0_est
  save_a1[i] <- a1_est
  save_a2[i] <- a2_est
  save_a3[i] <- a3_est
  save_se_a0[i] <- se[1]
  save_se_a1[i] <- se[2]
  save_se_a2[i] <- se[3]
  save_se_a3[i] <- se[4]
  # uci_a0 <- a0_est + 1.96*se[1] 
  # lci_a0 <- a0_est - 1.96*se[1]
  # uci_a1 <- a1_est + 1.96*se[2]
  # lci_a1 <- a1_est - 1.96*se[2]
  # uci_a2 <- a2_est + 1.96*se[3] 
  # lci_a2 <- a2_est - 1.96*se[3]
  # uci_a3 <- a3_est + 1.96*se[4] 
  # lci_a3 <- a3_est - 1.96*se[4]
  # 
  #c ci
  c0_est <- plcoptim$par[5]
  c1_est <- plcoptim$par[6]
  c2_est <- plcoptim$par[7]
  c3_est <- plcoptim$par[8]
  save_c0[i] <- c0_est
  save_c1[i] <- c1_est
  save_c2[i] <- c2_est
  save_c3[i] <- c3_est
  save_se_c0[i] <- se[5]
  save_se_c1[i] <- se[6]
  save_se_c2[i] <- se[7]
  save_se_c3[i] <- se[8]
  # uci_c0 <- c0_est + 1.96*se[5] 
  # lci_c0 <- c0_est - 1.96*se[5]
  # uci_c1 <- c1_est + 1.96*se[6]
  # lci_c1 <- c1_est - 1.96*se[6]
  # uci_c2 <- c2_est + 1.96*se[7] 
  # lci_c2 <- c2_est - 1.96*se[7]
  # uci_c3 <- c3_est + 1.96*se[8]
  # lci_c3 <- c3_est - 1.96*se[8]
  
  #b ci
  b0_est <- plcoptim$par[9]
  b1_est <- plcoptim$par[10]
  b2_est <- plcoptim$par[11]
  b3_est <- plcoptim$par[12]
  save_b0[i] <- b0_est
  save_b1[i] <- b1_est
  save_b2[i] <- b2_est
  save_b3[i] <- b3_est
  save_se_c0[i] <- se[9]
  save_se_c1[i] <- se[10]
  save_se_c2[i] <- se[11]
  save_se_c3[i] <- se[12]
  # uci_b0 <- b0_est + 1.96*se[9] 
  # lci_b0 <- b0_est - 1.96*se[9]
  # uci_b1 <- b1_est + 1.96*se[10]
  # lci_b1 <- b1_est - 1.96*se[10]
  # uci_b2 <- b2_est + 1.96*se[11]
  # lci_b2 <- b2_est - 1.96*se[11]
  # uci_b3 <- b3_est + 1.96*se[12]
  # lci_b3 <- b3_est - 1.96*se[12]
  
  
  ###HR###
  
  var_a0 <- fisher_info[1,1]
  var_a1 <- fisher_info[2,2]
  var_a2 <- fisher_info[3,3]
  var_a3 <- fisher_info[4,4]
  
  var_c0 <- fisher_info[5,5]
  var_c1 <- fisher_info[6,6]
  var_c2 <- fisher_info[7,7]
  var_c3 <- fisher_info[8,8]
  
  var_hr_l1_age <- exp(a1_est)^2 * var_a1
  var_hr_l1_gen <- exp(a2_est)^2 * var_a2
  var_hr_l1_donor <- exp(a3_est)^2 * var_a3
  var_hr_l2_age <- exp(c1_est)^2 * var_c1
  var_hr_l2_gen <- exp(c2_est)^2 * var_c2
  var_hr_l2_donor <- exp(c3_est)^2 * var_c3
  
  esthr_l1_age <- exp(a1_est)
  esthr_l1_gen <- exp(a2_est)
  esthr_l1_donor <- exp(a3_est)
  
  esthr_l2_age <- exp(c1_est)
  esthr_l2_gen <- exp(c2_est)
  esthr_l2_donor <- exp(c3_est)
  
  save_hr_l1_age[i] <- esthr_l1_age
  save_hr_l1_gen[i] <- esthr_l1_gen
  save_hr_l1_donor[i] <- esthr_l1_donor
  
  save_hr_l2_age[i] <- esthr_l2_age
  save_hr_l2_gen[i] <- esthr_l2_gen
  save_hr_l2_donor[i] <- esthr_l2_donor
  
  # YW commented out 2 January 2023
  # hr_l1_lwci_age <- esthr_l1_age - 1.96*sqrt(var_hr_l1_age)
  # hr_l1_upci_age <- esthr_l1_age + 1.96*sqrt(var_hr_l1_age)
  # 
  # hr_l1_lwci_gen <- esthr_l1_gen - 1.96*sqrt(var_hr_l1_gen)
  # hr_l1_upci_gen <- esthr_l1_gen + 1.96*sqrt(var_hr_l1_gen)
  # 
  # hr_l1_lwci_donor <- esthr_l1_donor - 1.96*sqrt(var_hr_l1_donor)
  # hr_l1_upci_donor <- esthr_l1_donor + 1.96*sqrt(var_hr_l1_donor)
  # 
  # hr_l2_lwci_age <- esthr_l2_age - 1.96*sqrt(var_hr_l2_age)
  # hr_l2_upci_age <- esthr_l2_age + 1.96*sqrt(var_hr_l2_age)
  # 
  # hr_l2_lwci_gen <- esthr_l2_gen - 1.96*sqrt(var_hr_l2_gen)
  # hr_l2_upci_gen <- esthr_l2_gen + 1.96*sqrt(var_hr_l2_gen)
  # 
  # hr_l2_lwci_donor <- esthr_l2_donor - 1.96*sqrt(var_hr_l2_donor)
  # hr_l2_upci_donor <- esthr_l2_donor + 1.96*sqrt(var_hr_l2_donor)
  # 
  
  # #############################################
  # ############### REPORTING ###################
  # #############################################
  # 
  # if(true_a0 <= uci_a0   && true_a0 >= lci_a0)   {counter_a0 = counter_a0+1}
  # if(true_a1 <= uci_a1   && true_a1 >= lci_a1)   {counter_a1 = counter_a1+1}
  # if(true_a2 <= uci_a2   && true_a2 >= lci_a2)   {counter_a2 = counter_a2+1}
  # if(true_a3 <= uci_a3   && true_a3 >= lci_a3)   {counter_a3 = counter_a3+1}
  # 
  # if(true_c0 <= uci_c0   && true_c0 >= lci_c0)   {counter_c0 = counter_c0+1}
  # if(true_c1 <= uci_c1   && true_c1 >= lci_c1)   {counter_c1 = counter_c1+1}
  # if(true_c2 <= uci_c2   && true_c2 >= lci_c2)   {counter_c2 = counter_c2+1}
  # if(true_c3 <= uci_c3   && true_c3 >= lci_c3)   {counter_c3 = counter_c3+1}
  # 
  # if(true_b0 <= uci_b0   && true_b0 >= lci_b0)   {counter_b0 = counter_b0+1}
  # if(true_b1 <= uci_b1   && true_b1 >= lci_b1)   {counter_b1 = counter_b1+1}
  # if(true_b2 <= uci_b2   && true_b2 >= lci_b2)   {counter_b2 = counter_b2+1}
  # if(true_b3 <= uci_b3   && true_b3 >= lci_b3)   {counter_b3 = counter_b3+1}
  # 
  # if(true_hr_l1_age <= hr_l1_upci_age && true_hr_l1_age >= hr_l1_lwci_age) {counter_hr_l1_age=counter_hr_l1_age+1}
  # if(true_hr_l2_age <= hr_l2_upci_age && true_hr_l2_age >= hr_l2_lwci_age) {counter_hr_l2_age=counter_hr_l2_age+1}
  # if(true_hr_l1_gen <= hr_l1_upci_gen && true_hr_l1_gen >= hr_l1_lwci_gen) {counter_hr_l1_gen=counter_hr_l1_gen+1}
  # if(true_hr_l2_gen <= hr_l2_upci_gen && true_hr_l2_gen >= hr_l2_lwci_gen) {counter_hr_l2_gen=counter_hr_l2_gen+1}
  # if(true_hr_l1_donor <= hr_l1_upci_donor && true_hr_l1_donor >= hr_l1_lwci_donor) {counter_hr_l1_donor=counter_hr_l1_donor+1}
  # if(true_hr_l2_donor <= hr_l2_upci_donor && true_hr_l2_donor >= hr_l2_lwci_donor) {counter_hr_l2_donor=counter_hr_l2_donor+1}
  # 
  # bias_a0[i] <- true_a0 - a0_est
  # bias_a1[i] <- true_a1 - a1_est
  # bias_a2[i] <- true_a2 - a2_est
  # bias_a3[i] <- true_a3 - a3_est
  # 
  # bias_c0[i] <- true_c0 - c0_est
  # bias_c1[i] <- true_c1 - c1_est
  # bias_c2[i] <- true_c2 - c2_est
  # bias_c3[i] <- true_c3 - c3_est
  # 
  # bias_b0[i] <- true_b0 - b0_est
  # bias_b1[i] <- true_b1 - b1_est
  # bias_b2[i] <- true_b2 - b2_est
  # bias_b3[i] <- true_b3 - b3_est
  # 
  # bias_l1_hr_age[i] <- true_hr_l1_age - esthr_l1_age
  # bias_l2_hr_age[i] <- true_hr_l2_age - esthr_l2_age
  # bias_l1_hr_gen[i] <- true_hr_l1_gen - esthr_l1_gen
  # bias_l2_hr_gen[i] <- true_hr_l2_gen - esthr_l2_gen
  # bias_l1_hr_donor[i] <- true_hr_l1_donor - esthr_l1_donor
  # bias_l2_hr_donor[i] <- true_hr_l2_donor - esthr_l2_donor
  
  print(i)
}  # END OF LOOP

# YW commented out 2 January 2023
# #hrs#
# #bias
# hr_l1_bias_age <- mean(bias_l1_hr_age)
# hr_l2_bias_age <- mean(bias_l2_hr_age)
# hr_l1_bias_gen <- mean(bias_l1_hr_gen)
# hr_l2_bias_gen <- mean(bias_l2_hr_gen)
# hr_l1_bias_donor <- mean(bias_l1_hr_donor)
# hr_l2_bias_donor <- mean(bias_l2_hr_donor)
# 
# #coverage
# hr_l1_cov_age <- (counter_hr_l1_age / runs) * 100
# hr_l2_cov_age <- (counter_hr_l2_age / runs) * 100
# hr_l1_cov_gen <- (counter_hr_l1_gen / runs) * 100
# hr_l2_cov_gen <- (counter_hr_l2_gen / runs) * 100
# hr_l1_cov_donor <- (counter_hr_l1_donor / runs) * 100
# hr_l2_cov_donor <- (counter_hr_l2_donor / runs) * 100
# 
# #variance
# hr_l1_var_age <- var(save_hr_l1_age)
# hr_l2_var_age <- var(save_hr_l2_age)
# hr_l1_var_gen <- var(save_hr_l1_gen)
# hr_l2_var_gen <- var(save_hr_l2_gen)
# hr_l1_var_donor <- var(save_hr_l1_donor)
# hr_l2_var_donor <- var(save_hr_l2_donor)
# 
# #mse
# hr_l1_mse_age <- mean((save_hr_l1_age - true_hr_l1_age)^2)
# hr_l2_mse_age <- mean((save_hr_l2_age - true_hr_l2_age)^2)
# hr_l1_mse_gen <- mean((save_hr_l1_gen - true_hr_l1_gen)^2)
# hr_l2_mse_gen <- mean((save_hr_l2_gen - true_hr_l2_gen)^2)
# hr_l1_mse_donor <- mean((save_hr_l1_donor - true_hr_l1_donor)^2)
# hr_l2_mse_donor <- mean((save_hr_l2_donor - true_hr_l2_donor)^2)
# 
# #betas#
# #bias
# b0_bias <- mean(bias_b0)
# b1_bias <- mean(bias_b1)
# b2_bias <- mean(bias_b2)
# b3_bias <- mean(bias_b3)
# 
# #coverage
# b0_cov <- (counter_b0 / runs) * 100
# b1_cov <- (counter_b1 / runs) * 100
# b2_cov <- (counter_b2 / runs) * 100
# b3_cov <- (counter_b3 / runs) * 100
# 
# #variance
# b0_var <- var(save_b0)
# b1_var <- var(save_b1)
# b2_var <- var(save_b2)
# b3_var <- var(save_b3)
# 
# #mse
# b0_mse <- mean((save_b0 - true_b0)^2)
# b1_mse <- mean((save_b1 - true_b1)^2)
# b2_mse <- mean((save_b2 - true_b2)^2)
# b3_mse <- mean((save_b3 - true_b3)^2)
# 
# #a#
# #bias
# a0_bias <- mean(bias_a0)
# a1_bias <- mean(bias_a1)
# a2_bias <- mean(bias_a2)
# a3_bias <- mean(bias_a3)
# 
# #coverage
# a0_cov <- (counter_a0 / runs) * 100
# a1_cov <- (counter_a1 / runs) * 100
# a2_cov <- (counter_a2 / runs) * 100
# a3_cov <- (counter_a3 / runs) * 100
# 
# #variance
# a0_var <- var(save_a0)
# a1_var <- var(save_a1)
# a2_var <- var(save_a2)
# a3_var <- var(save_a3)
# 
# #mse
# a0_mse <- mean((save_a0 - true_a0)^2)
# a1_mse <- mean((save_a1 - true_a1)^2)
# a2_mse <- mean((save_a2 - true_a2)^2)
# a3_mse <- mean((save_a3 - true_a3)^2)
# 
# #cs#
# #bias
# c0_bias <- mean(bias_c0)
# c1_bias <- mean(bias_c1)
# c2_bias <- mean(bias_c2)
# c3_bias <- mean(bias_c3)
# 
# #coverage
# c0_cov <- (counter_c0 / runs) * 100
# c1_cov <- (counter_c1 / runs) * 100
# c2_cov <- (counter_c2 / runs) * 100
# c3_cov <- (counter_c3 / runs) * 100
# 
# #variance
# c0_var <- var(save_c0)
# c1_var <- var(save_c1)
# c2_var <- var(save_c2)
# c3_var <- var(save_c3)
# 
# #mse
# c0_mse <- mean((save_c0 - true_c0)^2)
# c1_mse <- mean((save_c1 - true_c1)^2)
# c2_mse <- mean((save_c2 - true_c2)^2)
# c3_mse <- mean((save_c3 - true_c3)^2)

# YW added 2 January 2023 - calculate these outside the loop
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
b0_bias <- mean(save_b0) - true_b0    
b1_bias <- mean(save_b1) - true_b1   
b2_bias <- mean(save_b2) - true_b2 
b3_bias <- mean(save_b3) - true_b3 

# CI's
uci_b0 <- save_b0 + 1.96*save_se_b0 
lci_b0 <- save_b0 - 1.96*save_se_b0
uci_b1 <- save_b1 + 1.96*save_se_b1
lci_b1 <- save_b1 - 1.96*save_se_b1
uci_b2 <- save_b2 + 1.96*save_se_b2
lci_b2 <- save_b2 - 1.96*save_se_b2
uci_b3 <- save_b3 + 1.96*save_se_b3
lci_b3 <- save_b3 - 1.96*save_se_b3

# coverage
#coverage
counter_b0 <- sum(true_b0 <= uci_b0 & true_b0 >= lci_b0)
b0_cov <- (counter_b0 / runs) * 100

counter_b1 <- sum(true_b1 <= uci_b1 & true_b1 >= lci_b1)
b1_cov <- (counter_b1 / runs) * 100

counter_b2 <- sum(true_b2 <= uci_b2 & true_b2 >= lci_b2)
b2_cov <- (counter_b2 / runs) * 100

counter_b3 <- sum(true_b3 <= uci_b3 & true_b3 >= lci_b3)
b3_cov <- (counter_b3 / runs) * 100

#mse
b0_mse <- mean((save_b0 - true_b0)^2)
b1_mse <- mean((save_b1 - true_b1)^2)
b2_mse <- mean((save_b2 - true_b2)^2)
b3_mse <- mean((save_b3 - true_b3)^2)

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

counter_b0_low <- sum(save_b0 == b0_lw) 
counter_b0_upper <- sum(save_b0 == b0_up)
counter_b1_low <- sum(save_b1 == b1_lw)
counter_b1_upper <- sum(save_b1 == b1_up)
counter_b2_low <- sum(save_b2 == b2_lw)
counter_b2_upper <- sum(save_b2 == b2_up)
counter_b3_low <- sum(save_b3 == b3_lw)
counter_b3_upper <- sum(save_b3 == b3_up)

counter_low <- c(counter_a0_low, counter_a1_low, counter_a2_low, counter_a3_low,
                 counter_c0_low, counter_c1_low, counter_c2_low, counter_c3_upper, 
                 counter_b0_low, counter_b1_low, counter_b2_low, counter_b3_upper)

counter_upper <- c(counter_a0_upper, counter_a1_upper, counter_a2_upper, counter_a3_upper,
                   counter_c0_upper, counter_c1_upper, counter_c2_upper, counter_c3_upper, 
                   counter_b0_low, counter_b1_low, counter_b2_low, counter_b3_upper)

counter_hit_boundaries <- data.frame(rbind(counter_low, counter_upper))
names(counter_hit_boundaries) <- c("a0", "a1", "a2", "a3",
                                   "c0", "c1", "c2", "c3",
                                   "b0", "b1", "b2", "b3")
row.names(counter_hit_boundaries) <- c("lower", "upper")

print(counter_hit_boundaries)

# put results together and write to CSV file
# mean of bias
# hr_l1 represents non-terminal event; hr_l2 represents terminal event
mean_bias <- c(a0_bias, a1_bias, a2_bias, a3_bias, 
               c0_bias, c1_bias, c2_bias, c3_bias,
               b0_bias, b1_bias, b2_bias, b3_bias,
               hr_l1_bias_age, hr_l1_bias_donor, hr_l1_bias_gen,
               hr_l2_bias_age, hr_l2_bias_donor, hr_l2_bias_gen)

# Coverage probability
CP <- c(a0_cov, a1_cov, a2_cov, a3_cov,
        c0_cov, c1_cov, c2_cov, c3_cov,
        b0_cov, b1_cov, b2_cov, b3_cov,
        hr_l1_cov_age, hr_l1_cov_donor, hr_l1_cov_gen,
        hr_l2_cov_age, hr_l2_cov_donor, hr_l2_cov_gen)

MSE <- c(a0_mse,a1_mse,a2_mse,a3_mse,
         c0_mse, c1_mse, c2_mse, c3_mse,
         b0_mse, b1_mse, b2_mse, b3_mse,
         hr_l1_mse_age, hr_l1_mse_donor, hr_l1_mse_gen,
         hr_l2_mse_age, hr_l2_mse_donor, hr_l2_mse_gen)

# put results together
items<-c("a0", "a1", "a2", "a3",
         "c0", "c1", "c2", "c3", 
         "b0", "b1", "b2", "b3",
         "NT_age", "NT_donor", "NT_gen",
         "T_age", "T_donor", "T_gen")
Results <- cbind.data.frame(items, mean_bias, CP, MSE)

Results[,2:4] <- round(Results[,2:4],3)
Results

rownames(Results)<-NULL
end_time <- Sys.time()

run_time = end_time - start_time
run_time

Estimates = data.frame(a0.est = save_a0, a1.est = save_a1, a2.est = save_a2, a3.est = save_a3, 
                       c0.est = save_c0, c1.est = save_c1, c2.est = save_c2, c3.est = save_c3, 
                       b0.est = save_b0, 
                       hr.l1.age.est = save_hr_l1_age, hr.l2.age.est = save_hr_l2_age,
                       hr.l1.gen.est = save_hr_l1_gen, hr.l2.gen.est = save_hr_l2_gen,
                       hr.l1.donor.est = save_hr_l1_donor, hr.l2.donor.est = save_hr_l2_donor)

# Output estimates from each replication and summary of results -----------------
write.csv(Results, row.names=F,file=paste0(dir_results, out_file_results))
write.csv(Estimates, row.names=F,file=paste0(dir_results, out_file_estimates))
print(run_time)
print("Simulation1 model2 for clayton exponential model completed successfully!")
