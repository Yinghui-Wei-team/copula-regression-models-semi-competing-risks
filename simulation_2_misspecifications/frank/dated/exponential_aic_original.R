rm(list=ls())
library(copula)
library(mvtnorm)
library(plyr)
library(survival)
library(numDeriv)

#####################################################################################
#################### Frank, age, gen from exp chose with aic ########################
#####################################################################################

set.seed(96662391)
n <- 3000
runs <- 1000

true_b0 <- 3.44
true_b1 <- 5.23

true_a0 <- -3.42 #exp lambda1
true_a1 <-  0.37 #exp lambda1
true_c0 <- -4.27 #exp lambda2
true_c1 <-  1.41 #exp lambda2


true_theta_d0 <- true_b0
true_theta_d1 <- true_b0+true_b1

t_theta_d0_cop <- frankCopula(true_theta_d0)
true_rho_d0 <- rho(t_theta_d0_cop)

t_theta_d1_cop <- frankCopula(true_theta_d1)
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

hr_1_lw = 0
hr_1_up = 0
hr_1_cross = 0
hr_2_lw = 0
hr_2_up = 0
hr_2_cross = 0


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

counter_hr_l1 = 0
counter_hr_l2 = 0
counter_rho_d0 = 0
counter_rho_d1 = 0
###############################################################
###################### run 'runs' times #######################
###############################################################

for (i in 1:runs){
  
  ###############################################################
  ######################## generate data ########################
  ###############################################################
  
  #Step 1: generate age categories
  age <- rbinom(n,1,0.40)          #40% are in the older age group in NHSBT data
  
  for(k in 1:(n)){   #loop to generate U an V from age-varying theta
    m=1                  
    
    #Step 2: generate 1 random variable from Uniform(0,a) distribution 
    
    u1 <- runif(m,0,1)       
    
    #Step 3: X_true generated from u1 values (T1 from later)
    
    theta1 <- true_b0 + true_b1 * age[k]
    true_l1[k] <- exp(true_a0 + true_a1 * age[k])
    true_l2[k] <- exp(true_c0 + true_c1 * age[k]) 
    
    #Step 4: Conditional distribution method
    
    fc<- frankCopula(theta1, dim=2) #only allows 1 theta at a time (-> loop)
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
  df$X[which(df$X==0)]<-0.1
  df$Y[which(df$Y==0)]<-0.1
  
  
  ######################################################
  ############### Frank pseudo likelihood ##############
  #################### Exponential #####################
  ######################################################
  
  fpl_exp <- function(para, X, Y, d1, d2, age){
    
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
    
    theta <- b0+b1*age
    
    C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
    
    #part1 <- d1*d2*(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*S1)-1)-log(exp(theta*S2)-1)+log(lambda1)-lambda1*X+log(lambda2)-lambda2*Y)
    part1 <- d1*d2*(log(theta*exp(theta*C)*(exp(theta*C)-1)*lambda1*exp(-lambda1*X)*lambda2*exp(-lambda2*Y))-log((exp(theta*S1)-1)*(exp(theta*S2)-1)))
    part2 <- d1*(1-d2)*log(((1-exp(theta*C))/(1-exp(theta*S1)))*lambda1*exp(-lambda1*X))
    part3 <- (1-d1)*d2*log(((1-exp(theta*C))/(1-exp(theta*S2)))*lambda2*exp(-lambda2*Y))
    part4<-((1-d1)*(1-d2))*log(C)
    logpl<-sum(part1+part2+part3+part4) 
    
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
  
  b0_lw <- 1
  b0_up <- 10
  b1_lw <- 0
  b1_up <- 10
  
  plfoptim_exp <- optim(c(-3,0.01,-3,0.01,3,0), fpl_exp, method="L-BFGS-B",
                  lower=c(a0_lw,a1_lw,c0_lw,c1_lw,b0_lw,b1_lw),upper=c(a0_up,a1_up,c0_up,c1_up,b0_up,b1_up), 
                  X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                  control=list(fnscale=-1),hessian=TRUE)
  
  if(plfoptim_exp$par[1] == a0_lw) {counter_a0_low = counter_a0_low + 1}
  if(plfoptim_exp$par[1] == a0_up) {counter_a0_upper = counter_a0_upper + 1}
  if(plfoptim_exp$par[2] == a1_lw) {counter_a1_low = counter_a1_low + 1}
  if(plfoptim_exp$par[2] == a1_up) {counter_a1_upper = counter_a1_upper + 1}
  if(plfoptim_exp$par[3] == c0_lw) {counter_c0_low = counter_c0_low + 1}
  if(plfoptim_exp$par[3] == c0_up) {counter_c0_upper = counter_c0_upper + 1}
  if(plfoptim_exp$par[4] == c1_lw) {counter_c1_low = counter_c1_low + 1}
  if(plfoptim_exp$par[4] == c1_up) {counter_c1_upper = counter_c1_upper + 1}
  if(plfoptim_exp$par[5] == b0_lw) {counter_b0_low = counter_b0_low + 1}
  if(plfoptim_exp$par[5] == b0_up) {counter_b0_upper = counter_b0_upper + 1}
  if(plfoptim_exp$par[6] == b1_lw) {counter_b1_low = counter_b1_low + 1}
  if(plfoptim_exp$par[6] == b1_up) {counter_b1_upper = counter_b1_upper + 1}
  
  if(plfoptim_exp$par[1] == a0_lw) {break}
  if(plfoptim_exp$par[1] == a0_up) {break}
  if(plfoptim_exp$par[2] == a1_lw) {break}
  if(plfoptim_exp$par[2] == a1_up) {break}
  if(plfoptim_exp$par[3] == c0_lw) {break}
  if(plfoptim_exp$par[3] == c0_up) {break}
  if(plfoptim_exp$par[4] == c1_lw) {break}
  if(plfoptim_exp$par[4] == c1_up) {break}
  if(plfoptim_exp$par[5] == b0_lw) {break}
  if(plfoptim_exp$par[5] == b0_up) {break}
  if(plfoptim_exp$par[6] == b1_lw) {break}
  if(plfoptim_exp$par[6] == b1_up) {break}
  
  ######################################################
  ############### Frank pseudo likelihood ##############
  ###################### Weibull #######################
  ######################################################
  
  fpl_wei <- function(para, X, Y, d1, d2, age){
    alpha1 <- para[1]
    x1 <- para[2]
    x2 <- para[3]
    alpha2 <- para[4]
    y1 <- para[5]
    y2 <- para[6]
    b0 <- para[7]
    b1 <- para[8]
    
    theta <- b0+b1*age
    beta1 <- exp(x1+x2*age)
    beta2 <- exp(y1+y2*age)
    
    S1 <- exp(-beta1*X^alpha1)
    S2 <- exp(-beta2*Y^alpha2)
    S1[which(S1<0.1^8)]=0.1^8
    S2[which(S2<0.1^8)]=0.1^8
    
    f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
    f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
    f1[which(f1<0.1^8)]=0.1^8
    f2[which(f2<0.1^8)]=0.1^8
    
    C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
    C[which(C<0.1^8)]=0.1^8
    
    #part1 <- d1*d2*(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*S1)-1)-log(exp(theta*S2)-1)+log(f1)+log(f2))
    part1 <- d1*d2*(log(theta*exp(theta*C)*(exp(theta*C)-1)*f1*f2)-log((exp(theta*S1)-1)*(exp(theta*S2)-1)))
    part2 <- d1*(1-d2)*log(((1-exp(theta*C))*f1)/(1-exp(theta*S1)))
    part3 <- (1-d1)*d2*log(((1-exp(theta*C))*f2)/(1-exp(theta*S2)))
    part4<-((1-d1)*(1-d2))*log(C)
    
    #print(theta*exp(theta*C)*(exp(theta*C)-1)*f1*f2)
    logpl<-sum(part1+part2+part3+part4) 
    return(logpl)
  }
  
  #fpl_wei(c(0.7, -3, 0.3, 0.9, -4, 1.4, 3.5, 4.1),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  
  a1_lw <- 0.1
  a1_up <- 1.5
  a2_lw <- 0.1
  a2_up <- 1.5
  
  x1_lw <- -5
  x1_up <- -2
  x2_lw <- -2
  x2_up <- 1
  
  y1_lw <- -8
  y1_up <- -3
  y2_lw <- -2
  y2_up <- 2
  
  b0_lw <- 2
  b0_up <- 8
  b1_lw <- 1
  b1_up <- 8
  
  plfoptim_wei <- optim(c(0.67, -2.5, -0.6, 0.94, -3.3, -0.9, true_b0, true_b1), fpl_wei, method="L-BFGS-B",
                    lower=c(a1_lw, x1_lw, x2_lw, a2_lw, y1_lw, y2_lw, b0_lw, b1_lw),
                    upper=c(a1_up, x1_up, x2_up, a2_up, y1_up, y2_up, b0_up, b1_up), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age, control=list(fnscale=-1), hessian=TRUE)
  
  #fpl_wei(c(0.7, -3, 0.3, 0.9, -4, 1.4, 3.5, 4.1),X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  
  if(plfoptim_wei$par[1] == a1_lw) {counter_a1_low <<- counter_a1_low + 1}
  if(plfoptim_wei$par[1] == a1_up) {counter_a1_upper <<- counter_a1_upper + 1}
  if(plfoptim_wei$par[2] == x1_lw) {counter_x1_low <<- counter_x1_low + 1}
  if(plfoptim_wei$par[2] == x1_up) {counter_x1_upper <<- counter_x1_upper + 1}
  if(plfoptim_wei$par[3] == x2_lw) {counter_x2_low <<- counter_x2_low + 1}
  if(plfoptim_wei$par[3] == x2_up) {counter_x2_upper <<- counter_x2_upper + 1}
  if(plfoptim_wei$par[4] == a2_lw) {counter_a2_low <<- counter_a2_low + 1}
  if(plfoptim_wei$par[4] == a2_up) {counter_a2_upper <<- counter_a2_upper + 1}
  if(plfoptim_wei$par[5] == y1_lw) {counter_y1_low <<- counter_y1_low + 1}
  if(plfoptim_wei$par[5] == y1_up) {counter_y1_upper <<- counter_y1_upper + 1}
  if(plfoptim_wei$par[6] == y2_lw) {counter_y2_low <<- counter_y2_low + 1}
  if(plfoptim_wei$par[6] == y2_up) {counter_y2_upper <<- counter_y2_upper + 1}
  if(plfoptim_wei$par[7] == b0_lw) {counter_b0_low <<- counter_b0_low + 1}
  if(plfoptim_wei$par[7] == b0_up) {counter_b0_upper <<- counter_b0_upper + 1}
  if(plfoptim_wei$par[8] == b1_lw) {counter_b1_low <<- counter_b1_low + 1}
  if(plfoptim_wei$par[8] == b1_up) {counter_b1_upper <<- counter_b1_upper + 1}
  
  if(plfoptim_wei$par[1] == a1_lw) {break}
  if(plfoptim_wei$par[1] == a1_up) {break}
  if(plfoptim_wei$par[2] == x1_lw) {break}
  if(plfoptim_wei$par[2] == x1_up) {break}
  if(plfoptim_wei$par[3] == x2_lw) {break}
  if(plfoptim_wei$par[3] == x2_up) {break}
  if(plfoptim_wei$par[4] == a2_lw) {break}
  if(plfoptim_wei$par[4] == a2_up) {break}
  if(plfoptim_wei$par[5] == y1_lw) {break}
  if(plfoptim_wei$par[5] == y1_up) {break}
  if(plfoptim_wei$par[6] == y2_lw) {break}
  if(plfoptim_wei$par[6] == y2_up) {break}
  if(plfoptim_wei$par[7] == b0_lw) {break}
  if(plfoptim_wei$par[7] == b0_up) {break}
  if(plfoptim_wei$par[8] == b1_lw) {break}
  if(plfoptim_wei$par[8] == b1_up) {break}

  ######################################################
  ############### Frank pseudo likelihood ##############
  ###################### Gompertz ######################
  ######################################################
  
  fpl_gom <- function(para, X, Y, d1, d2, age){
    gamma1 <- para[1]
    p0 <- para[2]
    p1 <- para[3]
    gamma2 <- para[4]
    q0 <- para[5]
    q1 <- para[6]
    
    b0 <- para[7]
    b1 <- para[8]
    
    theta <- b0+b1*age
    lambda1 <- exp(p0+p1*age)
    lambda2 <- exp(q0+q1*age)
    
    S1 <- exp(-lambda1/gamma1*(exp(gamma1*X)-1))
    S2 <- exp(-lambda2/gamma2*(exp(gamma2*Y)-1))
    S1[which(S1<0.1^8)]=0.1^8
    S2[which(S2<0.1^8)]=0.1^8
    
    f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
    f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
    f1[which(f1<0.1^8)]=0.1^8
    f2[which(f2<0.1^8)]=0.1^8
    
    C= -1/theta * log(((1-exp(-theta)-(1-exp(-theta*S1))*(1-exp(-theta*S2))))/(1-exp(-theta)))
    
    part1 <- d1*d2*(log(theta*exp(theta*C)*(exp(theta*C)-1)*f1*f2)-log((exp(theta*S1)-1)*(exp(theta*S2)-1)))
    part2 <- d1*(1-d2)*log(((1-exp(theta*C))/(1-exp(theta*S1)))*f1)
    part3 <- (1-d1)*d2*log(((1-exp(theta*C))/(1-exp(theta*S2)))*f2)
    part4<-((1-d1)*(1-d2))*log(C)
    
    logpl <- sum(part1+part2+part3+part4) 
    return(logpl)
  }
  
  
  g1_lw <- -0.1
  g1_up <- 0.1
  
  p0_lw <- -5
  p0_up <- -1
  p1_lw <- -2
  p1_up <- 1
  
  g2_lw <- -0.1
  g2_up <- 0.1
  
  q0_lw <- -5
  q0_up <- -1
  q1_lw <- -2
  q1_up <- 2
  
  b0_lw <- 1
  b0_up <- 8
  b1_lw <- -2
  b1_up <- 8
  
  plfoptim_gom <- optim(c(-0.01, -3, -0.5, 0.02, -3.5, -0.8, 0.5, 0), fpl_gom, method="L-BFGS-B",
                    lower=c(g1_lw,p0_lw,p1_lw,g2_lw,q0_lw, q1_lw, b0_lw,b1_lw),
                    upper=c(g1_up,p0_up,p1_up,g2_up,q0_up, q1_up, b0_up,b1_up), 
                    X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                    control=list(fnscale=-1),hessian=TRUE)
  
  if(plfoptim_gom$par[1] == g1_lw) {counter_g1_low <<- counter_g1_low + 1}
  if(plfoptim_gom$par[1] == g1_up) {counter_g1_upper <<- counter_g1_upper + 1}
  if(plfoptim_gom$par[2] == p0_lw) {counter_p0_low <<- counter_p0_low + 1}
  if(plfoptim_gom$par[2] == p0_up) {counter_p0_upper <<- counter_p0_upper + 1}
  if(plfoptim_gom$par[3] == p1_lw) {counter_p1_low <<- counter_p1_low + 1}
  if(plfoptim_gom$par[3] == p1_up) {counter_p1_upper <<- counter_p1_upper + 1}
  if(plfoptim_gom$par[4] == g2_lw) {counter_g2_low <<- counter_g2_low + 1}
  if(plfoptim_gom$par[4] == g2_up) {counter_g2_upper <<- counter_g2_upper + 1}
  if(plfoptim_gom$par[5] == q0_lw) {counter_q0_low <<- counter_q0_low + 1}
  if(plfoptim_gom$par[5] == q0_up) {counter_q0_upper <<- counter_q0_upper + 1}
  if(plfoptim_gom$par[6] == q1_lw) {counter_q1_low <<- counter_q1_low + 1}
  if(plfoptim_gom$par[6] == q1_up) {counter_q1_upper <<- counter_q1_upper + 1}
  if(plfoptim_gom$par[7] == b0_lw) {counter_b0_low <<- counter_b0_low + 1}
  if(plfoptim_gom$par[7] == b0_up) {counter_b0_upper <<- counter_b0_upper + 1}
  if(plfoptim_gom$par[8] == b1_lw) {counter_b1_low <<- counter_b1_low + 1}
  if(plfoptim_gom$par[8] == b1_up) {counter_b1_upper <<- counter_b1_upper + 1}

  if(plfoptim_gom$par[1] == g1_lw) {break}
  if(plfoptim_gom$par[1] == g1_up) {break}
  if(plfoptim_gom$par[2] == p0_lw) {break}
  if(plfoptim_gom$par[2] == p0_up) {break}
  if(plfoptim_gom$par[3] == p1_lw) {break}
  if(plfoptim_gom$par[3] == p1_up) {break}
  if(plfoptim_gom$par[4] == g2_lw) {break}
  if(plfoptim_gom$par[4] == g2_up) {break}
  if(plfoptim_gom$par[5] == q0_lw) {break}
  if(plfoptim_gom$par[5] == q0_up) {break}
  if(plfoptim_gom$par[6] == q1_lw) {break}
  if(plfoptim_gom$par[6] == q1_up) {break}
  if(plfoptim_gom$par[7] == b0_lw) {break}
  if(plfoptim_gom$par[7] == b0_up) {break}
  if(plfoptim_gom$par[8] == b1_lw) {break}
  if(plfoptim_gom$par[8] == b1_up) {break}
  
  ########################################################
  ######################### AICS #########################
  ########################################################
  
  loglik_exp <- fpl_exp(plfoptim_exp$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  k_exp <- length(plfoptim_exp$par)
  aic_exp <- -2*loglik_exp+2*k_exp
  
  loglik_wei <- fpl_wei(plfoptim_wei$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  k_wei <- length(plfoptim_wei$par)
  aic_wei <- -2*loglik_wei+2*k_wei
  
  loglik_gom <- fpl_gom(plfoptim_gom$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
  k_gom <- length(plfoptim_gom$par)
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
    
    fisher_info <- solve(-plfoptim_exp$hessian) #inverse -hess
    se <- sqrt(diag(fisher_info)) 
    
    #b ci
    est_a0 <- plfoptim_exp$par[1]
    est_a1 <- plfoptim_exp$par[2]
    est_c0 <- plfoptim_exp$par[3]
    est_c1<- plfoptim_exp$par[4]
    est_b0 <- plfoptim_exp$par[5]
    est_b1 <- plfoptim_exp$par[6]
    varb0 <- fisher_info[5,5]
    varb1 <- fisher_info[6,6]
    covb0b1 <- fisher_info[5,6]
    
    #rho for age.grp=0
    theta_d0 <- est_b0
    var_theta_d0 <- varb0
    theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
    theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
    
    rho_d0_cop <- frankCopula(theta_d0)
    est_rho_d0 <- rho(rho_d0_cop)
    rho_d0_lwci_cop <- frankCopula(theta_d0_lwci)
    rho_d0_lwci <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- frankCopula(theta_d0_upci)
    rho_d0_upci <- rho(rho_d0_upci_cop) 
    save_rho_d0[i] <- est_rho_d0
    
    #rho for age.grp=1
    theta_d1 <- est_b0+est_b1
    var_theta_d1 <- varb0+2*covb0b1+varb1
    theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
    theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
    
    rho_d1_cop <- frankCopula(theta_d1)
    est_rho_d1 <- rho(rho_d1_cop)
    rho_d1_lwci_cop <- frankCopula(theta_d1_lwci)
    rho_d1_lwci <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- frankCopula(theta_d1_upci)
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
    
    fisher_info <- solve(-plfoptim_wei$hessian) #inverse -hess
    se <- sqrt(diag(fisher_info)) 
    
    #point and var
    est_a1 <- plfoptim_wei$par[1]
    est_a2 <- plfoptim_wei$par[4]
    est_x0 <- plfoptim_wei$par[2]
    est_x1 <- plfoptim_wei$par[3]
    est_y0 <- plfoptim_wei$par[5]
    est_y1 <- plfoptim_wei$par[6]
    est_b0 <- plfoptim_wei$par[7]
    est_b1 <- plfoptim_wei$par[8]
    varb0 <- fisher_info[7,7]
    varb1 <- fisher_info[8,8]
    covb0b1 <- fisher_info[7,8]
    
    #rho for age.grp=0
    theta_d0 <- est_b0
    var_theta_d0 <- varb0
    theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
    theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
    
    rho_d0_cop <- frankCopula(theta_d0)
    est_rho_d0 <- rho(rho_d0_cop)
    rho_d0_lwci_cop <- frankCopula(theta_d0_lwci)
    rho_d0_lwci <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- frankCopula(theta_d0_upci)
    rho_d0_upci <- rho(rho_d0_upci_cop) 
    save_rho_d0[i] <- est_rho_d0
    
    #rho for age.grp=1
    theta_d1 <- est_b0+est_b1
    var_theta_d1 <- varb0+2*covb0b1+varb1
    theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
    theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
    
    rho_d1_cop <- frankCopula(theta_d1)
    est_rho_d1 <- rho(rho_d1_cop)
    rho_d1_lwci_cop <- frankCopula(theta_d1_lwci)
    rho_d1_lwci <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- frankCopula(theta_d1_upci)
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
    hessian <- hessian(fpl_gom, plfoptim_gom$par, X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
    #fishers info matrix
    fisher_info <- solve(-hessian) 
    #Standard Error
    se <- sqrt(diag(fisher_info)) 
    
    est_b0 <- plfoptim_gom$par[7]
    est_b1 <- plfoptim_gom$par[8]
    est_g1 <- plfoptim_gom$par[1]
    est_g2 <- plfoptim_gom$par[4]
    est_p0 <- plfoptim_gom$par[2]
    est_p1 <- plfoptim_gom$par[3]
    est_q0 <- plfoptim_gom$par[5]
    est_q1 <- plfoptim_gom$par[6]
    varb0 <- fisher_info[7,7]
    varb1 <- fisher_info[8,8]
    covb0b1 <- fisher_info[7,8]
    
    #rho for age.grp=0
    theta_d0 <- est_b0
    var_theta_d0 <- varb0
    theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
    theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
    
    rho_d0_cop <- frankCopula(theta_d0)
    est_rho_d0 <- rho(rho_d0_cop)
    rho_d0_lwci_cop <- frankCopula(theta_d0_lwci)
    rho_d0_lwci <- rho(rho_d0_lwci_cop)
    rho_d0_upci_cop <- frankCopula(theta_d0_upci)
    rho_d0_upci <- rho(rho_d0_upci_cop) 
    save_rho_d0[i] <- est_rho_d0
    
    #rho for age.grp=1
    theta_d1 <- est_b0+est_b1
    var_theta_d1 <- varb0+2*covb0b1+varb1
    theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
    theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
    
    rho_d1_cop <- frankCopula(theta_d1)
    est_rho_d1 <- rho(rho_d1_cop)
    rho_d1_lwci_cop <- frankCopula(theta_d1_lwci)
    rho_d1_lwci <- rho(rho_d1_lwci_cop)
    rho_d1_upci_cop <- frankCopula(theta_d1_upci)
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
  
  if (hr_l1_lwci <1 & hr_l1_upci <1) {hr_1_lw=hr_1_lw+1
  } else if (hr_l1_lwci >1 & hr_l1_upci >1) {hr_1_up=hr_1_up+1
  } else {hr_1_cross = hr_1_cross+1}
  
  if (hr_l2_lwci <1 & hr_l2_upci <1) {hr_2_lw=hr_2_lw+1
  } else if (hr_l2_lwci >1 & hr_l2_upci >1) {hr_2_up=hr_2_up+1
  } else {hr_2_cross = hr_2_cross+1}
  
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

#change meaning of hr
#true hr 1 = 1.5 & true hr 2 = 4 for age
hr_1_perc <- hr_1_up / runs *100
hr_2_perc <- hr_2_up / runs *100


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

print(paste("HR 1 true direction percentage", hr_1_perc)) #percentage of correct direction
print(paste("HR 2 true direction percentage", hr_2_perc)) #percentage of correct direction

