rm(list=ls())
library(copula)
library(mvtnorm)
library(plyr)
library(survival)
library(numDeriv)

#####################################################################################
#################### Gumbel, age, gen from gom chose with aic #######################
#####################################################################################

set.seed(91186965)
n <- 3000
runs <- 1000

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

## stuff for later ##
save_hr_l1 <- rep(0,runs)
save_hr_l2 <- rep(0,runs)
save_rho_d0 <- rep(0,runs)
save_rho_d1 <- rep(0,runs)

bias_l1_hr <- rep(0,runs)
bias_l2_hr <- rep(0,runs)
bias_rho_d0 <- rep(0,runs)
bias_rho_d1 <- rep(0,runs)

counter_hr_l1 = 0
counter_hr_l2 = 0
counter_rho_d0 = 0
counter_rho_d1 = 0

counter_exp = 0
counter_wei = 0
counter_gom = 0

###############################################################
###################### run 'runs' times #######################
###############################################################

for (i in 1:runs){
  
  ###############################################################
  ######################## generate data ########################
  ###############################################################
  
  #Step 1: generate age categories
  age.grp <- rbinom(n,1,0.40)          #40% are in the older age group in NHSBT data
  
  for(k in 1:(n)){   #loop to generate U an V from age-varying theta
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
  
  gpl_exp <- function(para, X, Y, d1, d2, age.grp){
    
    a0 <- para[1]
    a1 <- para[2]
    c0 <- para[3]
    c1 <- para[4]
    b0 <- para[5]
    b1 <- para[6]
    
    lambda1 <- exp(a0+a1*age.grp)
    lambda2 <- exp(c0+c1*age.grp)
    S1<-exp(-lambda1*X)
    S2<-exp(-lambda2*Y)
    f1 <- lambda1*exp(-lambda1*X)
    f2 <- lambda2*exp(-lambda2*Y)
    
    theta <- exp(b0+b1*age.grp)+1
    
    C=exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
    
    part1 <- d1*d2*(log(C)+(theta-1)*log(-log(S1))+(theta-1)*log(-log(S2))+log(theta-1+((-log(S1))^theta+(-log(S2))^theta)^(1/theta))-log(S1)-log(S2)-(2*theta-1)*log(-log(C))+log(lambda1)-lambda1*X+log(lambda2)-lambda2*Y)
    part2 <- ((1-d2)*d1) * (log(C*(-log(S1))^(theta-1)*f1)-log(S1*(-log(C))^(theta-1)))
    part3 <- ((1-d1)*d2) * (log(C*(-log(S2))^(theta-1)*f2)-log(S2*(-log(C))^(theta-1)))
    part4<-((1-d1)*(1-d2))*log(C)
    logpl<-sum(part1+part2+part3+part4) 
    
    #print(sum(is.na(part3)))
    return(logpl)
  }
  
  
  a0_lw <- -10
  a0_up <- -2
  a1_lw <- -10
  a1_up <- 1.5
  
  c0_lw <- -10
  c0_up <- -2
  c1_lw <- -10
  c1_up <- 3
  
  b0_lw <- -5
  b0_up <- -1
  b1_lw <- -5
  b1_up <- 4
  
  plgoptim_exp <- optim(c(-3.5, 0.2, -4.3, 1.4, true_b0, true_b1), gpl_exp, method="L-BFGS-B",
                  lower=c(a0_lw,a1_lw,c0_lw,c1_lw,b0_lw,b1_lw),upper=c(a0_up,a1_up,c0_up,c1_up,b0_up,b1_up), 
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp,
                  control=list(fnscale=-1),hessian=TRUE)
  
  if(plgoptim_exp$par[1] == a0_lw) {break}
  if(plgoptim_exp$par[1] == a0_up) {break}
  if(plgoptim_exp$par[2] == a1_lw) {break}
  if(plgoptim_exp$par[2] == a1_up) {break}
  if(plgoptim_exp$par[3] == c0_lw) {break}
  if(plgoptim_exp$par[3] == c0_up) {break}
  if(plgoptim_exp$par[4] == c1_lw) {break}
  if(plgoptim_exp$par[4] == c1_up) {break}
  if(plgoptim_exp$par[5] == b0_lw) {break}
  if(plgoptim_exp$par[5] == b0_up) {break}
  if(plgoptim_exp$par[6] == b1_lw) {break}
  if(plgoptim_exp$par[6] == b1_up) {break}
  
  ########################################################
  ################ Gumbel pseudo likelihood ##############
  ####################### Weibull ########################
  ########################################################
  
  gpl_wei <- function(para, X, Y, d1, d2, age.grp){
    alpha1 <- para[1]
    x1 <- para[2]
    x2 <- para[3]
    alpha2 <- para[4]
    y1 <- para[5]
    y2 <- para[6]
    b0 <- para[7]
    b1 <- para[8]
    
    theta <- exp(b0+b1*age.grp)+1
    beta1 <- exp(x1+x2*age.grp)
    beta2 <- exp(y1+y2*age.grp)
    
    S1 <- exp(-beta1*X^alpha1)
    S2 <- exp(-beta2*Y^alpha2)
    S1[which(S1<0.1^8)] <- 0.1^8
    S2[which(S2<0.1^8)] <- 0.1^8
    
    
    f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
    f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
    f1[which(f1<0.1^8)] <- 0.1^8
    f2[which(f2<0.1^8)] <- 0.1^8
    
    C=exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
    
    part1 <- d1*d2*(log(C)+(theta-1)*log(-log(S1))+(theta-1)*log(-log(S2))+log(theta-1+((-log(S1))^theta+(-log(S2))^theta)^(1/theta))-log(S1)-log(S2)-(2*theta-1)*log(-log(C))+log(f1)+log(f2))
    part2 <- d1*(1-d2)*(log(C)+(theta-1)*log(-log(S1))-log(S1)-(theta-1)*log(-log(C))+log(f1))
    part3 <-((1-d1)*(d2))*(log(C)+(theta-1)*log(-log(S2))-log(S2)-(theta-1)*log(-log(C))+log(f2))
    part4<-((1-d1)*(1-d2))*log(C)
    logpl<-sum(part1+part2+part3+part4) 
    
    logpl<-sum(part1+part2+part3+part4) 
    return(logpl)
  }

  a1_lw <- 0.01
  a1_up <- 1.5
  a2_lw <- 0.01
  a2_up <- 1.5
  
  x1_lw <- -10
  x1_up <- -1.5
  x2_lw <- -10
  x2_up <- 1
  
  y1_lw <- -10
  y1_up <- -2
  y2_lw <- -10
  y2_up <- 2.5
  
  b0_lw <- -6
  b0_up <- -1
  b1_lw <- -5
  b1_up <- 3
  
  plgoptim_wei <- optim(c(0.6, -3, 0, 0.9, -4, 1.5, true_b0, true_b1), gpl_wei, method="L-BFGS-B",
                    lower=c(a1_lw, x1_lw, x2_lw, a2_lw, y1_lw, y2_lw, b0_lw, b1_lw),
                    upper=c(a1_up, x1_up, x2_up, a2_up, y1_up, y2_up, b0_up, b1_up), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age.grp=df$age.grp, control=list(fnscale=-1), hessian=TRUE)
  
  if(plgoptim_wei$par[1] == a1_lw) {break}
  if(plgoptim_wei$par[1] == a1_up) {break}
  if(plgoptim_wei$par[2] == x1_lw) {break}
  if(plgoptim_wei$par[2] == x1_up) {break}
  if(plgoptim_wei$par[3] == x2_lw) {break}
  if(plgoptim_wei$par[3] == x2_up) {break}
  if(plgoptim_wei$par[4] == a2_lw) {break}
  if(plgoptim_wei$par[4] == a2_up) {break}
  if(plgoptim_wei$par[5] == y1_lw) {break}
  if(plgoptim_wei$par[5] == y1_up) {break}
  if(plgoptim_wei$par[6] == y2_lw) {break}
  if(plgoptim_wei$par[6] == y2_up) {break}
  if(plgoptim_wei$par[7] == b0_lw) {break}
  if(plgoptim_wei$par[7] == b0_up) {break}
  if(plgoptim_wei$par[8] == b1_lw) {break}
  if(plgoptim_wei$par[8] == b1_up) {break}
  
  ########################################################
  ############### Gumbel pseudo likelihood ###############
  ####################### Gompertz #######################
  ########################################################
  
  gpl_gom <- function(para, X, Y, d1, d2, age.grp){
    gamma1 <- para[1]
    p0 <- para[2]
    p1 <- para[3]
    gamma2 <- para[4]
    q0 <- para[5]
    q1 <- para[6]
    b0 <- para[7]
    b1 <- para[8]
    
    theta <- exp(b0+b1*age.grp)+1
    lambda1 <- exp(p0+p1*age.grp)
    lambda2 <- exp(q0+q1*age.grp)
    
    S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
    S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
    S1[which(S1 < 0.1^8)] <- 0.1^8
    S2[which(S2 < 0.1^8)] <- 0.1^8
    
    f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
    f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
    f1[which(f1 < 0.1^8)] <- 0.1^8
    f2[which(f2 < 0.1^8)] <- 0.1^8
    
    C=exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
    
    part1 <- d1*d2*(log(C)+(theta-1)*log(-log(S1))+(theta-1)*log(-log(S2))+log(theta-1+((-log(S1))^theta+(-log(S2))^theta)^(1/theta))-log(S1)-log(S2)-(2*theta-1)*log(-log(C))+log(f1)+log(f2))
    part2 <- d1*(1-d2)*(log(C)+(theta-1)*log(-log(S1))-log(S1)-(theta-1)*log(-log(C))+log(f1))
    part3 <-((1-d1)*(d2))*(log(C)+(theta-1)*log(-log(S2))-log(S2)-(theta-1)*log(-log(C))+log(f2))
    part4<-((1-d1)*(1-d2))*log(C)
    logpl<-sum(part1+part2+part3+part4) 
    
    logpl <- sum(part1+part2+part3+part4) 
    return(logpl)
  }
  

  g1_lw <- -0.25
  g1_up <- 0.8
  
  p0_lw <- -6
  p0_up <- -2
  p1_lw <- -2
  p1_up <- 2
  
  g2_lw <- -0.2
  g2_up <- 0.1
  
  q0_lw <- -6
  q0_up <- -2
  q1_lw <- -1
  q1_up <- 3
  
  b0_lw <- -4
  b0_up <- -1
  b1_lw <- 0
  b1_up <- 3
  
  plgoptim_gom <- optim(c(0.01, true_p0, true_p1, true_g2, true_q0, true_q1, true_b0, true_b1), gpl_gom, method="L-BFGS-B",
                    lower=c(g1_lw,p0_lw,p1_lw,g2_lw,q0_lw, q1_lw, b0_lw,b1_lw),
                    upper=c(g1_up,p0_up,p1_up,g2_up,q0_up, q1_up, b0_up,b1_up), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp,
                    control=list(fnscale=-1),hessian=TRUE)
  
  if(plgoptim_gom$par[1] == g1_lw) {break}
  if(plgoptim_gom$par[1] == g1_up) {break}
  if(plgoptim_gom$par[2] == p0_lw) {break}
  if(plgoptim_gom$par[2] == p0_up) {break}
  if(plgoptim_gom$par[3] == p1_lw) {break}
  if(plgoptim_gom$par[3] == p1_up) {break}
  if(plgoptim_gom$par[4] == g2_lw) {break}
  if(plgoptim_gom$par[4] == g2_up) {break}
  if(plgoptim_gom$par[5] == q0_lw) {break}
  if(plgoptim_gom$par[5] == q0_up) {break}
  if(plgoptim_gom$par[6] == q1_lw) {break}
  if(plgoptim_gom$par[6] == q1_up) {break}
  if(plgoptim_gom$par[7] == b0_lw) {break}
  if(plgoptim_gom$par[7] == b0_up) {break}
  if(plgoptim_gom$par[8] == b1_lw) {break}
  if(plgoptim_gom$par[8] == b1_up) {break}
  
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
    
    #rho for age.grp=0
    theta_d0 <- exp(est_b0)+1
    var_theta_d0 <- exp(2*est_b0)*varb0
    theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
    if (theta_d0_lwci<1) {theta_d0_lwci=1} #if hits lower bound take it as 1
    theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
    
    rho_d0_cop <- gumbelCopula(theta_d0)
    est_rho_d0 <- rho(rho_d0_cop)
    rho_d0_lwci_cop <- gumbelCopula(theta_d0_lwci)
    rho_d0_lwci <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- gumbelCopula(theta_d0_upci)
    rho_d0_upci <- rho(rho_d0_upci_cop) 
    save_rho_d0[i] <- est_rho_d0
    
    #rho for age.grp=1
    theta_d1 <- exp(est_b0+est_b1)+1
    var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
    theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
    if (theta_d1_lwci<1) {theta_d1_lwci=1} #if hits lower bound take it as 1
    theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
    
    rho_d1_cop <- gumbelCopula(theta_d1)
    est_rho_d1 <- rho(rho_d1_cop)
    rho_d1_lwci_cop <- gumbelCopula(theta_d1_lwci)
    rho_d1_lwci <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- gumbelCopula(theta_d1_upci)
    rho_d1_upci <- rho(rho_d1_upci_cop) 
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
    
    hr_l1_lwci <- est_hr_l1 - 1.96*sqrt(var_hr_l1)
    hr_l1_upci <- est_hr_l1 + 1.96*sqrt(var_hr_l1)
    hr_l2_lwci <- est_hr_l2 - 1.96*sqrt(var_hr_l2)
    hr_l2_upci <- est_hr_l2 + 1.96*sqrt(var_hr_l2)
    
    ############### REPORTING ###################
    
    if(true_hr_l1 <= hr_l1_upci && true_hr_l1 >= hr_l1_lwci) {counter_hr_l1=counter_hr_l1+1}
    if(true_hr_l2 <= hr_l2_upci && true_hr_l2 >= hr_l2_lwci) {counter_hr_l2=counter_hr_l2+1}
    
    if(true_rho_d0 <= rho_d0_upci && true_rho_d0 >= rho_d0_lwci) {counter_rho_d0=counter_rho_d0+1}
    if(true_rho_d1 <= rho_d1_upci && true_rho_d1 >= rho_d1_lwci) {counter_rho_d1=counter_rho_d1+1}
    
    bias_l1_hr[i] <- true_hr_l1 - est_hr_l1
    bias_l2_hr[i] <- true_hr_l2 - est_hr_l2
    bias_rho_d0[i] <- true_rho_d0 - est_rho_d0
    bias_rho_d1[i] <- true_rho_d1 - est_rho_d1
    
    
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
    
    #rho for age.grp=0
    theta_d0 <- exp(est_b0)+1
    var_theta_d0 <- exp(2*est_b0)*varb0
    theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
    if (theta_d0_lwci<1) {theta_d0_lwci=1} #if hits lower bound take it as 1
    theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
    
    rho_d0_cop <- gumbelCopula(theta_d0)
    est_rho_d0 <- rho(rho_d0_cop)
    rho_d0_lwci_cop <- gumbelCopula(theta_d0_lwci)
    rho_d0_lwci <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- gumbelCopula(theta_d0_upci)
    rho_d0_upci <- rho(rho_d0_upci_cop) 
    save_rho_d0[i] <- est_rho_d0
    
    #rho for age.grp=1
    theta_d1 <- exp(est_b0+est_b1)+1
    var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
    theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
    if (theta_d1_lwci<1) {theta_d1_lwci=1} #if hits lower bound take it as 1
    theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
    
    rho_d1_cop <- gumbelCopula(theta_d1)
    est_rho_d1 <- rho(rho_d1_cop)
    rho_d1_lwci_cop <- gumbelCopula(theta_d1_lwci)
    rho_d1_lwci <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- gumbelCopula(theta_d1_upci)
    rho_d1_upci <- rho(rho_d1_upci_cop) 
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
    
    hr_l1_lwci <- est_hr_l1 - 1.96*sqrt(var_hr_l1)
    hr_l1_upci <- est_hr_l1 + 1.96*sqrt(var_hr_l1)
    hr_l2_lwci <- est_hr_l2 - 1.96*sqrt(var_hr_l2)
    hr_l2_upci <- est_hr_l2 + 1.96*sqrt(var_hr_l2)
    
    ############### REPORTING ###################
    
    if(true_hr_l1 <= hr_l1_upci && true_hr_l1 >= hr_l1_lwci) {counter_hr_l1=counter_hr_l1+1}
    if(true_hr_l2 <= hr_l2_upci && true_hr_l2 >= hr_l2_lwci) {counter_hr_l2=counter_hr_l2+1}
    
    if(true_rho_d0 <= rho_d0_upci && true_rho_d0 >= rho_d0_lwci) {counter_rho_d0=counter_rho_d0+1}
    if(true_rho_d1 <= rho_d1_upci && true_rho_d1 >= rho_d1_lwci) {counter_rho_d1=counter_rho_d1+1}
    
    bias_l1_hr[i] <- true_hr_l1 - est_hr_l1
    bias_l2_hr[i] <- true_hr_l2 - est_hr_l2
    bias_rho_d0[i] <- true_rho_d0 - est_rho_d0
    bias_rho_d1[i] <- true_rho_d1 - est_rho_d1
    
    
  } else{
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
    
    #rho for age.grp=0
    theta_d0 <- exp(est_b0)+1
    var_theta_d0 <- exp(2*est_b0)*varb0
    theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
    if (theta_d0_lwci<1) {theta_d0_lwci=1} #if hits lower bound take it as 1
    theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
    
    rho_d0_cop <- gumbelCopula(theta_d0)
    est_rho_d0 <- rho(rho_d0_cop)
    rho_d0_lwci_cop <- gumbelCopula(theta_d0_lwci)
    rho_d0_lwci <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- gumbelCopula(theta_d0_upci)
    rho_d0_upci <- rho(rho_d0_upci_cop) 
    save_rho_d0[i] <- est_rho_d0
    
    #rho for age.grp=1
    theta_d1 <- exp(est_b0+est_b1)+1
    var_theta_d1 <- exp(2*(est_b0+est_b1))*(varb0+2*covb0b1+varb1)
    theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
    if (theta_d1_lwci<1) {theta_d1_lwci=1} #if hits lower bound take it as 1
    theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
    
    rho_d1_cop <- gumbelCopula(theta_d1)
    est_rho_d1 <- rho(rho_d1_cop)
    rho_d1_lwci_cop <- gumbelCopula(theta_d1_lwci)
    rho_d1_lwci <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- gumbelCopula(theta_d1_upci)
    rho_d1_upci <- rho(rho_d1_upci_cop) 
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
    
    hr_l1_lwci <- est_hr_l1 - 1.96*sqrt(var_hr_l1)
    hr_l1_upci <- est_hr_l1 + 1.96*sqrt(var_hr_l1)
    hr_l2_lwci <- est_hr_l2 - 1.96*sqrt(var_hr_l2)
    hr_l2_upci <- est_hr_l2 + 1.96*sqrt(var_hr_l2)
    
    if(true_hr_l1 <= hr_l1_upci && true_hr_l1 >= hr_l1_lwci) {counter_hr_l1=counter_hr_l1+1}
    if(true_hr_l2 <= hr_l2_upci && true_hr_l2 >= hr_l2_lwci) {counter_hr_l2=counter_hr_l2+1}
    
    if(true_rho_d0 <= rho_d0_upci && true_rho_d0 >= rho_d0_lwci) {counter_rho_d0=counter_rho_d0+1}
    if(true_rho_d1 <= rho_d1_upci && true_rho_d1 >= rho_d1_lwci) {counter_rho_d1=counter_rho_d1+1}
    
    bias_l1_hr[i] <- true_hr_l1 - est_hr_l1
    bias_l2_hr[i] <- true_hr_l2 - est_hr_l2
    bias_rho_d0[i] <- true_rho_d0 - est_rho_d0
    bias_rho_d1[i] <- true_rho_d1 - est_rho_d1
    
  }
  
  
  print(i)
} # end of loop

#hrs#
#bias
hr_l1_bias <- mean(abs(bias_l1_hr))
hr_l2_bias <- mean(abs(bias_l2_hr))
#coverage
hr_l1_cov <- (counter_hr_l1 / runs) * 100
hr_l2_cov <- (counter_hr_l2 / runs) * 100
#variance
hr_l1_var <- (sd(save_hr_l1))^0.5
hr_l2_var <- (sd(save_hr_l2))^0.5
#mse
hr_l1_mse <- hr_l1_bias^2+hr_l1_var
hr_l2_mse <- hr_l2_bias^2+hr_l2_var

#rho#
rho_d0_bias <- mean(abs(bias_rho_d0))
rho_d1_bias <- mean(abs(bias_rho_d1))
rho_d0_cov <- (counter_rho_d0 / runs) * 100 
rho_d1_cov <- (counter_rho_d1 / runs) * 100
rho_d0_var <- (sd(save_rho_d0))^0.5
rho_d1_var <- (sd(save_rho_d1))^0.5
rho_d0_mse <- rho_d0_bias^2+rho_d0_var
rho_d1_mse <- rho_d1_bias^2+rho_d1_var

#counters#
exp_perc <- counter_exp / runs *100
wei_perc <- counter_wei / runs *100
gom_perc <- counter_gom / runs *100

#################################

print(paste("HR l1 bias", hr_l1_bias))
print(paste("HR l2 bias", hr_l2_bias))
print(paste("HR l1 mse", hr_l1_mse))
print(paste("HR l2 mse", hr_l2_mse))
print(paste("HR l1 cov", hr_l1_cov))
print(paste("HR l2 cov", hr_l2_cov))

#################################

print(paste("rho d0 bias", rho_d0_bias))
print(paste("rho d1 bias", rho_d1_bias))
print(paste("rho d0 mse", rho_d0_mse))
print(paste("rho d1 mse", rho_d1_mse))
print(paste("rho d0 cov", rho_d0_cov))
print(paste("rho d1 cov", rho_d1_cov))

#################################

print(paste("Exponential percentage", exp_perc))
print(paste("Weibull percentage", wei_perc))
print(paste("Gompertz percentage", gom_perc))

