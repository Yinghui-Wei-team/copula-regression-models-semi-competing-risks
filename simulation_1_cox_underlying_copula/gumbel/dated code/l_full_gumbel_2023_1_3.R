##############################################################################################
# Paper 2: Simulation 1
# Gumbel copula exponential survival model with covariates for hazard rates
# YW 25June 2021: added time counter, results to CSV file.
# YW 25 June 2021: corrected variances
# YW 25 June 2021: define save_hr
# YW 31/12/2022 update: change starting values to be between lower and upper bounds
##############################################################################################

rm(list=ls())
library(copula); library(mvtnorm);library(plyr); library(survival); library(numDeriv)

#################################################################################
# set up                                                                        #
#################################################################################
# directory if on PC
dir_results = "../../results/simulation_results/simulation1/"
# dir = "results"
# setwd(dir)
# directory if working on cluster
# dir = "/home/ywei/Simulation/Paper2/Gumbel"
# setwd(dir)

out_file_summary = "S1-gumbel-exponential-covariates-hazards.csv"

start_time <- Sys.time()
#set.seed(5333244)
set.seed(12345)
n <- 3000
runs <- 2
#runs <- 1000

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

## stuff for later ##
save_a0 <- rep(0,runs)
save_a1 <- rep(0,runs)
save_a2 <- rep(0,runs)
save_a3 <- rep(0,runs)
save_c0 <- rep(0,runs)
save_c1 <- rep(0,runs)
save_c2 <- rep(0,runs)
save_c3 <- rep(0,runs)
save_hr_l1_age <- rep(0,runs)
save_hr_l2_age <- rep(0,runs)
save_hr_l1_gen <- rep(0,runs)
save_hr_l2_gen <- rep(0,runs)
save_hr_l1_donor <- rep(0,runs)
save_hr_l2_donor <- rep(0,runs)

bias_a0 <- rep(0,runs)
bias_a1 <- rep(0,runs)
bias_a2 <- rep(0,runs)
bias_a3 <- rep(0,runs)
bias_c0 <- rep(0,runs)
bias_c1 <- rep(0,runs)
bias_c2 <- rep(0,runs)
bias_c3 <- rep(0,runs)
bias_l1_hr_age <- rep(0,runs)
bias_l1_hr_age <- rep(0,runs)
bias_l1_hr_gen <- rep(0,runs)
bias_l1_hr_gen <- rep(0,runs)
bias_l1_hr_donor <- rep(0,runs)
bias_l1_hr_donor <- rep(0,runs)
bias_l2_hr_age <- rep(0,runs)
bias_l2_hr_age <- rep(0,runs)
bias_l2_hr_gen <- rep(0,runs)
bias_l2_hr_gen <- rep(0,runs)
bias_l2_hr_donor <- rep(0,runs)
bias_l2_hr_donor <- rep(0,runs)

counter_a0 = 0
counter_a1 = 0
counter_a2 = 0
counter_a3 = 0
counter_c0 = 0
counter_c1 = 0
counter_c2 = 0
counter_c3 = 0
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
counter_t_low = 0
counter_a0_upper = 0
counter_a1_upper = 0 
counter_a2_upper = 0
counter_a3_upper = 0 
counter_c0_upper = 0
counter_c1_upper = 0
counter_c2_upper = 0
counter_c3_upper = 0 
counter_t_upper = 0


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
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  ########################################################
  ############## Gumbel pseudo likelihood ################
  ########################################################
  gpl<-function(para, X, Y, d1, d2, age.grp, gen, donor){
    a0 <- para[1]
    a1 <- para[2]
    a2 <- para[3]
    a3 <- para[4]
    
    c0 <- para[5]
    c1 <- para[6]
    c2 <- para[7]
    c3 <- para[8]
    
    theta <- para[9]
    
    lambda1 <- exp(a0+a1*age.grp+a2*gen+a3*donor)
    lambda2 <- exp(c0+c1*age.grp+c2*gen+c3*donor)
    
    S1 <- exp(-lambda1*X)
    S2 <- exp(-lambda2*Y) #problems when S2 too small -> 0 
    f1 <- lambda1*exp(-lambda1*X)
    f2 <- lambda2*exp(-lambda2*Y)
    
    S1[which(S1 < 0.1^(8))]=0.1^(8)
    S2[which(S2 < 0.1^(8))]=0.1^(8)
    theta[which(theta > 15)]=15
    #theta[which(theta > 3)]=3
    #S1[which(S1 > 0.9)]=0.9
    #S1[which(S1 > 0.9)]=0.9
    #subLogS1 <- log(S1) > -0.1^(8) ? -0.1^(8) : log(S1)
    subLogS1 <- log(S1) 
    subLogS1[which(subLogS1 > -0.1^(7))]= -0.1^7
    subLogS2<- log(S2) 
    subLogS2[which(subLogS2 > -0.1^(7))]= -0.1^7
    
    C=exp(-((-subLogS1)^(theta)+(-subLogS2)^(theta))^(1/theta))
    
    C[which(C<0.1^(8))]=0.1^(8)
    
    #part1 <- d1*d2*(log(C)+(theta-1)*log(-log(S1))+(theta-1)*log(-log(S2))+log(theta-1+((-log(S1))^theta+(-log(S2))^theta)^(1/theta))-log(S1)-log(S2)-(2*theta-1)*log(-log(C))+log(lambda1)-lambda1*X+log(lambda2)-lambda2*Y)
    part1 <- d1*d2*(log(C*(-subLogS1)^(theta-1)*(-subLogS2)^(theta-1)*(theta-1-log(C))*f1*f2)-log(S1*S2*(-log(C))^(2*theta-1)))
    
    part2 <- d1*(1-d2)*(log(C)+(theta-1)*log(-subLogS1)-subLogS1-(theta-1)*log(-log(C))+log(lambda1)-lambda1*X) 
    
    part3 <-((1-d1)*(d2))*(log(C)+(theta-1)*log(-subLogS2)-subLogS2-(theta-1)*log(-log(C))+log(lambda2)-lambda2*Y)
    
    part4<-((1-d1)*(1-d2))*log(C)
    
    #print(C*(-log(S1))^(theta-1)*(-log(S2))^(theta-1))
    #print((S1))
    #print(C)
    #print((-subLogS1)^(theta-1))
    #print((-log(S2))^(theta-1))
    #print(theta)
    #print((theta-1-log(C)))
    logpl<-sum(part1+part2+part3+part4) 
    return(logpl)
  }
  
  #gpl(c(a0_up,a1_up,a2_up, a3_up, c0_up,c1_up,c2_up,c3_up, b0_up,b1_up,b2_up, b3_up), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
  
  #gpl(c(-0.5,0,-0.5,0,5,2),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp)
  
  a0_lw <- -7
  a0_up <- -2
  a1_lw <- -5
  a1_up <- 1
  a2_lw <- -5
  a2_up <- 1
  a3_lw <- -5
  a3_up <- 1
  
  c0_lw <- -10
  c0_up <- -2
  c1_lw <- -5
  c1_up <- 3
  c2_lw <- -5
  c2_up <- 1
  c3_lw <- -5
  c3_up <- 1
  
  
  t_lw <- 1.01   #-7
  t_up <- 15    #-0.5


  starting_values = c(-5, -2, -2, -2,   # a
                      -5, -1, -2, -2,   # c
                      7                 # b
                      )
  # plgoptim <- optim(c(true_a0, true_a1, true_a2, true_a3, true_c0, true_c1, true_c2, true_c3,
  #                     true_b0), gpl, method="L-BFGS-B", 
  #                   lower=c(a0_lw,a1_lw,a2_lw, a3_lw,c0_lw,c1_lw,c2_lw, c3_lw, t_lw),
  #                   upper=c(a0_up,a1_up,a2_up, a3_up, c0_up,c1_up,c2_up,c3_up, t_up), 
  #                   X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen,
  #                   donor=df$donor, control=list(fnscale=-1),hessian=TRUE)
  
  plgoptim <- optim(starting_values, gpl, method="L-BFGS-B", 
                    lower=c(a0_lw,a1_lw,a2_lw, a3_lw,c0_lw,c1_lw,c2_lw, c3_lw, t_lw),
                    upper=c(a0_up,a1_up,a2_up, a3_up, c0_up,c1_up,c2_up,c3_up, t_up), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen,
                    donor=df$donor, control=list(fnscale=-1),hessian=TRUE)
  
  plgoptim$par
  
  #gpl(c(a0_lw,a1_lw,a2_lw, a3_lw,c0_lw,c1_lw,c2_lw, c3_lw,t_lw), X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
  
  if(plgoptim$par[1] == a0_lw) {counter_a0_low = counter_a0_low + 1}
  if(plgoptim$par[1] == a0_up) {counter_a0_upper = counter_a0_upper + 1}
  if(plgoptim$par[2] == a1_lw) {counter_a1_low = counter_a1_low + 1}
  if(plgoptim$par[2] == a1_up) {counter_a1_upper = counter_a1_upper + 1}
  if(plgoptim$par[3] == a2_lw) {counter_a2_low = counter_a2_low + 1}
  if(plgoptim$par[3] == a2_up) {counter_a2_upper = counter_a2_upper + 1}
  if(plgoptim$par[4] == a3_lw) {counter_a3_low = counter_a3_low + 1}
  if(plgoptim$par[4] == a3_up) {counter_a3_upper = counter_a3_upper + 1}
  
  if(plgoptim$par[5] == c0_lw) {counter_c0_low = counter_c0_low + 1}
  if(plgoptim$par[5] == c0_up) {counter_c0_upper = counter_c0_upper + 1}
  if(plgoptim$par[6] == c1_lw) {counter_c1_low = counter_c1_low + 1}
  if(plgoptim$par[6] == c1_up) {counter_c1_upper = counter_c1_upper + 1}
  if(plgoptim$par[7] == c2_lw) {counter_c2_low = counter_c2_low + 1}
  if(plgoptim$par[7] == c2_up) {counter_c2_upper = counter_c2_upper + 1}
  if(plgoptim$par[8] == c3_lw) {counter_c3_low = counter_c3_low + 1}
  if(plgoptim$par[8] == c3_up) {counter_c3_upper = counter_c3_upper + 1}
  
  if(plgoptim$par[9] == t_lw) {counter_t_low = counter_t_low + 1}
  if(plgoptim$par[9] == t_up) {counter_t_upper = counter_t_upper + 1}

  
  ########################################################
  ################## Confidence Intervals ################
  ########################################################
  hess <- hessian(gpl, plgoptim$par, X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, gen=df$gen, donor=df$donor)
  
  if (sum(is.na(hess))>0) {hess=plgoptim$hessian}    #use hess from hessian not optim (generally better) but if nan then use the one from optim 
  
  fisher_info <- solve(-plgoptim$hess)
  #fisher_info <- solve(-plgoptim$hessian) #inverse -hess
  #Standard error = sqrt(var/n)
  se<-sqrt(diag(fisher_info)) 
  
  
  #a ci
  a0_est <- plgoptim$par[1]
  a1_est <- plgoptim$par[2]
  a2_est <- plgoptim$par[3]
  a3_est <- plgoptim$par[4]
  save_a0[i] <- a0_est
  save_a1[i] <- a1_est
  save_a2[i] <- a2_est
  save_a3[i] <- a3_est
  uci_a0 <- a0_est + 1.96*se[1] 
  lci_a0 <- a0_est - 1.96*se[1]
  uci_a1 <- a1_est + 1.96*se[2]
  lci_a1 <- a1_est - 1.96*se[2]
  uci_a2 <- a2_est + 1.96*se[3] 
  lci_a2 <- a2_est - 1.96*se[3]
  uci_a3 <- a3_est + 1.96*se[4] 
  lci_a3 <- a3_est - 1.96*se[4]
  
  #c ci
  c0_est <- plgoptim$par[5]
  c1_est <- plgoptim$par[6]
  c2_est <- plgoptim$par[7]
  c3_est <- plgoptim$par[8]
  save_c0[i] <- c0_est
  save_c1[i] <- c1_est
  save_c2[i] <- c2_est
  save_c3[i] <- c3_est
  uci_c0 <- c0_est + 1.96*se[5] 
  lci_c0 <- c0_est - 1.96*se[5]
  uci_c1 <- c1_est + 1.96*se[6]
  lci_c1 <- c1_est - 1.96*se[6]
  uci_c2 <- c2_est + 1.96*se[7] 
  lci_c2 <- c2_est - 1.96*se[7]
  uci_c3 <- c3_est + 1.96*se[8]
  lci_c3 <- c3_est - 1.96*se[8]
  
  var_a0 <- fisher_info[1,1]
  var_a1 <- fisher_info[2,2]
  var_a2 <- fisher_info[3,3]
  var_a3 <- fisher_info[4,4]
  
  var_c0 <- fisher_info[5,5]
  var_c1 <- fisher_info[6,6]
  var_c2 <- fisher_info[7,7]
  var_c3 <- fisher_info[8,8]
  
  ###HR###
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
  
  
  #############################################
  ############### REPORTING ###################
  #############################################
  
  if(true_a0 <= uci_a0   && true_a0 >= lci_a0)   {counter_a0 = counter_a0+1}
  if(true_a1 <= uci_a1   && true_a1 >= lci_a1)   {counter_a1 = counter_a1+1}
  if(true_a2 <= uci_a2   && true_a2 >= lci_a2)   {counter_a2 = counter_a2+1}
  if(true_a3 <= uci_a3   && true_a3 >= lci_a3)   {counter_a3 = counter_a3+1}
  
  if(true_c0 <= uci_c0   && true_c0 >= lci_c0)   {counter_c0 = counter_c0+1}
  if(true_c1 <= uci_c1   && true_c1 >= lci_c1)   {counter_c1 = counter_c1+1}
  if(true_c2 <= uci_c2   && true_c2 >= lci_c2)   {counter_c2 = counter_c2+1}
  if(true_c3 <= uci_c3   && true_c3 >= lci_c3)   {counter_c3 = counter_c3+1}
  
  if(true_hr_l1_age <= hr_l1_upci_age && true_hr_l1_age >= hr_l1_lwci_age) {counter_hr_l1_age=counter_hr_l1_age+1}
  if(true_hr_l2_age <= hr_l2_upci_age && true_hr_l2_age >= hr_l2_lwci_age) {counter_hr_l2_age=counter_hr_l2_age+1}
  if(true_hr_l1_gen <= hr_l1_upci_gen && true_hr_l1_gen >= hr_l1_lwci_gen) {counter_hr_l1_gen=counter_hr_l1_gen+1}
  if(true_hr_l2_gen <= hr_l2_upci_gen && true_hr_l2_gen >= hr_l2_lwci_gen) {counter_hr_l2_gen=counter_hr_l2_gen+1}
  if(true_hr_l1_donor <= hr_l1_upci_donor && true_hr_l1_donor >= hr_l1_lwci_donor) {counter_hr_l1_donor=counter_hr_l1_donor+1}
  if(true_hr_l2_donor <= hr_l2_upci_donor && true_hr_l2_donor >= hr_l2_lwci_donor) {counter_hr_l2_donor=counter_hr_l2_donor+1}
  
  bias_a0[i] <- true_a0 - a0_est
  bias_a1[i] <- true_a1 - a1_est
  bias_a2[i] <- true_a2 - a2_est
  bias_a3[i] <- true_a3 - a3_est
  
  bias_c0[i] <- true_c0 - c0_est
  bias_c1[i] <- true_c1 - c1_est
  bias_c2[i] <- true_c2 - c2_est
  bias_c3[i] <- true_c3 - c3_est

  bias_l1_hr_age[i] <- true_hr_l1_age - esthr_l1_age
  bias_l2_hr_age[i] <- true_hr_l2_age - esthr_l2_age
  bias_l1_hr_gen[i] <- true_hr_l1_gen - esthr_l1_gen
  bias_l2_hr_gen[i] <- true_hr_l2_gen - esthr_l2_gen
  bias_l1_hr_donor[i] <- true_hr_l1_donor - esthr_l1_donor
  bias_l2_hr_donor[i] <- true_hr_l2_donor - esthr_l2_donor
  
  print(i)
}


#hrs#
#bias
hr_l1_bias_age <- mean(abs(bias_l1_hr_age))
hr_l2_bias_age <- mean(abs(bias_l2_hr_age))
hr_l1_bias_gen <- mean(abs(bias_l1_hr_gen))
hr_l2_bias_gen <- mean(abs(bias_l2_hr_gen))
hr_l1_bias_donor <- mean(abs(bias_l1_hr_donor))
hr_l2_bias_donor <- mean(abs(bias_l2_hr_donor))
#coverage
hr_l1_cov_age <- (counter_hr_l1_age / runs) * 100
hr_l2_cov_age <- (counter_hr_l2_age / runs) * 100
hr_l1_cov_gen <- (counter_hr_l1_gen / runs) * 100
hr_l2_cov_gen <- (counter_hr_l2_gen / runs) * 100
hr_l1_cov_donor <- (counter_hr_l1_donor / runs) * 100
hr_l2_cov_donor <- (counter_hr_l2_donor / runs) * 100
#variance: corrected by Yw
hr_l1_var_age <- var(save_hr_l1_age)
hr_l2_var_age <- var(save_hr_l2_age)
hr_l1_var_gen <- var(save_hr_l1_gen)
hr_l2_var_gen <- var(save_hr_l2_gen)
hr_l1_var_donor <- var(save_hr_l1_donor)
hr_l2_var_donor <- var(save_hr_l2_donor)
#mse
hr_l1_mse_age <- hr_l1_bias_age^2+hr_l1_var_age
hr_l2_mse_age <- hr_l2_bias_age^2+hr_l2_var_age
hr_l1_mse_gen <- hr_l1_bias_gen^2+hr_l1_var_gen
hr_l2_mse_gen <- hr_l2_bias_gen^2+hr_l2_var_gen
hr_l1_mse_donor <- hr_l1_bias_donor^2+hr_l1_var_donor
hr_l2_mse_donor <- hr_l2_bias_donor^2+hr_l2_var_donor


#a#
#bias
a0_bias <- mean(abs(bias_a0))
a1_bias <- mean(abs(bias_a1))
a2_bias <- mean(abs(bias_a2))
a3_bias <- mean(abs(bias_a3))
#coverage
a0_cov <- (counter_a0 / runs) * 100
a1_cov <- (counter_a1 / runs) * 100
a2_cov <- (counter_a2 / runs) * 100
a3_cov <- (counter_a3 / runs) * 100

#variance: corrected by YW 25june2021
a0_var <- var(save_a0)
a1_var <- var(save_a1)
a2_var <- var(save_a2)
a3_var <- var(save_a3)

#mse
a0_mse <- a0_bias^2+a0_var
a1_mse <- a1_bias^2+a1_var
a2_mse <- a2_bias^2+a2_var
a3_mse <- a3_bias^2+a3_var

#cs#
#bias
c0_bias <- mean(abs(bias_c0))
c1_bias <- mean(abs(bias_c1))
c2_bias <- mean(abs(bias_c2))
c3_bias <- mean(abs(bias_c3))
#coverage
c0_cov <- (counter_c0 / runs) * 100
c1_cov <- (counter_c1 / runs) * 100
c2_cov <- (counter_c2 / runs) * 100
c3_cov <- (counter_c3 / runs) * 100
#variance: corrected by YW
c0_var <- var(save_c0)
c1_var <- var(save_c1)
c2_var <- var(save_c2)
c3_var <- var(save_c3)
#mse
c0_mse <- c0_bias^2+c0_var
c1_mse <- c1_bias^2+c1_var
c2_mse <- c2_bias^2+c2_var
c3_mse <- c3_bias^2+c3_var
###############################

### REPORT ###
print(paste("a0 bias", a0_bias))
print(paste("a1 bias", a1_bias))
print(paste("a2 bias", a2_bias))
print(paste("a3 bias", a3_bias))
print(paste("c0 bias", c0_bias))
print(paste("c1 bias", c1_bias))
print(paste("c2 bias", c2_bias))
print(paste("c3 bias", c3_bias))
print(paste("a0 mse", a0_mse))
print(paste("a1 mse", a1_mse))
print(paste("a2 mse", a2_mse))
print(paste("a3 mse", a3_mse))
print(paste("c0 mse", c0_mse))
print(paste("c1 mse", c1_mse))
print(paste("c2 mse", c2_mse))
print(paste("c3 mse", c3_mse))
print(paste("a0 coverage", a0_cov))
print(paste("a1 coverage", a1_cov))
print(paste("a2 coverage", a2_cov))
print(paste("a3 coverage", a3_cov))
print(paste("c0 coverage", c0_cov))
print(paste("c1 coverage", c1_cov))
print(paste("c2 coverage", c2_cov))
print(paste("c3 coverage", c3_cov))
#################################

print(paste("HR l1 bias age", hr_l1_bias_age))
print(paste("HR l2 bias age", hr_l2_bias_age))
print(paste("HR l1 mse age", hr_l1_mse_age))
print(paste("HR l2 mse age", hr_l2_mse_age))
print(paste("HR l1 cov age", hr_l1_cov_age))
print(paste("HR l2 cov age", hr_l2_cov_age))

print(paste("HR l1 bias gen", hr_l1_bias_gen))
print(paste("HR l2 bias gen", hr_l2_bias_gen))
print(paste("HR l1 mse gen", hr_l1_mse_gen))
print(paste("HR l2 mse gen", hr_l2_mse_gen))
print(paste("HR l1 cov gen", hr_l1_cov_gen))
print(paste("HR l2 cov gen", hr_l2_cov_gen))

print(paste("HR l1 bias donor", hr_l1_bias_donor))
print(paste("HR l2 bias donor", hr_l2_bias_donor))
print(paste("HR l1 mse donor", hr_l1_mse_donor))
print(paste("HR l2 mse donor", hr_l2_mse_donor))
print(paste("HR l1 cov donor", hr_l1_cov_donor))
print(paste("HR l2 cov donor", hr_l2_cov_donor))
#################################

print(paste("counter a0 lower bound", counter_a0_low))
print(paste("counter a1 lower bound", counter_a1_low))
print(paste("counter a2 lower bound", counter_a2_low))
print(paste("counter a3 lower bound", counter_a3_low))
print(paste("counter c0 lower bound", counter_c0_low))
print(paste("counter c1 lower bound", counter_c1_low))
print(paste("counter c2 lower bound", counter_c2_low))
print(paste("counter c3 lower bound", counter_c3_low))
print(paste("counter t lower bound", counter_t_low))

print(paste("counter a0 upper bound", counter_a0_upper))
print(paste("counter a1 upper bound", counter_a1_upper))
print(paste("counter a2 upper bound", counter_a2_upper))
print(paste("counter a3 upper bound", counter_a3_upper))
print(paste("counter c0 upper bound", counter_c0_upper))
print(paste("counter c1 upper bound", counter_c1_upper))
print(paste("counter c2 upper bound", counter_c2_upper))
print(paste("counter c3 upper bound", counter_c3_upper))
print(paste("counter t upper bound", counter_t_upper))

end_time <- Sys.time()

run_time = end_time - start_time

run_time

# YW 23 June 2021: put results together and write to CSV file
# mean of bias
# hr_l1 represents non-terminal event; hr_l2 represents terminal event
mean_bias <- c(a0_bias, a1_bias, a2_bias, a3_bias, 
               c0_bias, c1_bias, c2_bias, c3_bias,
               hr_l1_bias_age, hr_l1_bias_donor, hr_l1_bias_gen,
               hr_l2_bias_age, hr_l2_bias_donor, hr_l2_bias_gen)

# also output variance

variance <-c(a0_var, a1_var,a2_var, a3_var,
             c0_var, c1_var, c2_var, c3_var,
             hr_l1_var_age, hr_l1_var_donor, hr_l1_bias_gen,
             hr_l2_var_age, hr_l2_var_donor, hr_l2_bias_gen)

# Coverage probability
CP <- c(a0_cov, a1_cov, a2_cov, a3_cov,
        c0_cov, c1_cov, c2_cov, c3_cov,
        hr_l1_cov_age, hr_l1_cov_donor, hr_l1_cov_gen,
        hr_l2_cov_age, hr_l2_cov_donor, hr_l2_cov_gen)

MSE <- c(a0_mse,a1_mse,a2_mse,a3_mse,
         c0_mse, c1_mse, c2_mse, c3_mse,
         hr_l1_mse_age, hr_l1_mse_donor, hr_l1_mse_gen,
         hr_l2_mse_age, hr_l2_mse_donor, hr_l2_mse_gen)

# YW: put results together
items<-c("a0", "a1", "a2", "a3",
         "c0", "c1", "c2", "c3", 
         "NT_age", "NT_donor", "NT_gen",
         "T_age", "T_donor", "T_gen")
Results <- cbind.data.frame(items, mean_bias, CP, variance, MSE)

Results[,2:4] <- round(Results[,2:4],4)

rownames(Results)<-NULL

# Output results---------------------------------------------------------------
write.csv(Results, row.names=F,file= paste0(dir_results,out_file_summary))
print("Simulation 1 for gumbel exponential model is completed successfuly! ")