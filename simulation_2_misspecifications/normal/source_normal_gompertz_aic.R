######################################################################################################
# Simulation study 2: evaluation of misspecification of survival distributions  
# Data are simulated from Normal copula weibull distribution
######################################################################################################
# Original code by LS
# YW, 1 Jan 2023:   1. update output directory and tidy up
#                   2. Put likelihood functions into a generic script under the functions folder
#####################################################################################################
library(copula); library(mvtnorm); library(numDeriv)
start_time = Sys.time()
#####################################################################################################
# output files                                                  
#####################################################################################################
# likelihood function
out_file_estimates <- paste0("S2_aic_normal_gompertz_", part, ".csv")

#####################################################################################
################## Normal, age, gen from gom chose with aic #######################
#####################################################################################
#set.seed(7877320)
set.seed(12345)
n <- 3000
runs <- 1000

true_b0 <- 0.36
true_b1 <- 0.28

true_g1 <- 0.001 #gomp gamma1
true_g2 <- 0.06  #gomp gamma2
true_p0 <- -3.37 #gomp lambda1
true_p1 <- 0.14 #gomp lambda1
true_q0 <- -4.79 #gomp lambda2
true_q1 <- 1.49 #gomp lambda2

true_theta_d0 <- (exp(2*true_b0)-1)/(exp(2*true_b0)+1)
true_theta_d1 <- (exp(2*(true_b0+true_b1))-1)/(exp(2*(true_b0+true_b1))+1)

t_theta_d0_cop <- normalCopula(true_theta_d0)
true_rho_d0 <- rho(t_theta_d0_cop)

t_theta_d1_cop <- normalCopula(true_theta_d1)
true_rho_d1 <- rho(t_theta_d1_cop)

#S1 exp, S2 weib
true_hr_l1 <- exp(true_p1)
true_hr_l2 <- exp(true_q1)

true_lambda1 <- rep(0,n)
true_lambda2 <- rep(0,n)
true_beta1 <- rep(0,n)
true_beta2 <- rep(0,n)
true_t <- rep(0,n)
true_r <- rep(0,n)
U1 <- rep(0,n)
V1 <- rep(0,n)

## stuff for later ##
save_hr_l1 <- rep(0,100)
save_hr_l2 <- rep(0,100)
save_rho_d0 <- rep(0,100)
save_rho_d1 <- rep(0,100)

bias_l1_hr <- rep(0,100)
bias_l2_hr <- rep(0,100)
bias_rho_d0 <- rep(0,100)
bias_rho_d1 <- rep(0,100)

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
  age <- rbinom(n,1,0.40)          #40% are in the older age group in NHSBT data
  
  for(k in 1:(n)){   #loop to generate U an V from age-varying theta
    m=1                  
    
    #Step 2: generate 1 random variable from Uniform(0,a) distribution 
    
    u1 <- runif(m,0,1)       
    
    #Step 3: X_true generated from u1 values (T1 from later)
    
    theta1 <- (exp(2*(true_b0+true_b1*age[k]))-1)/(exp(2*(true_b0+true_b1*age[k]))+1)
    true_lambda1[k] <- exp(true_p0 + true_p1 * age[k])
    true_lambda2[k] <- exp(true_q0 + true_q1 * age[k]) 
    
    #Step 4: Conditional distribution method
    
    fc<- normalCopula(theta1, dim=2) #only allows 1 theta at a time (-> loop)
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
  X<-pmin(T1,T2,C)
  Y<-pmin(T2, C) 
  d1<-ifelse(T1<=Y,1,0) 
  d2<-ifelse(T2<=C,1,0) 
  
  #Step 10: Create dataframe, true values of X and Y have association theta=b0+b1*X
  df<-data.frame(X, Y, d1, d2, age)
  
  if(i<101){ #1-100
    #if(i>100 & i<201){ #101-200
    #if(i>200 & i<301){ #201-300
    #if(i>300 & i<401){ #301-400
    #if(i>400 & i<501){ #401-500
    #if(i>500 & i<601){ #501-600
    #if(i>600 & i<701){ #601-700
    #if(i>700 & i<801){ #701-800
    #if(i>800 & i<901){ #801-900
    #if(i>900) #901-1000
    
    ########################################################
    ############### Normal pseudo likelihood ##############
    ##################### Exponential ######################
    ########################################################
    npl_exp <- function(para, X,Y,d1,d2,age){
      
      a0 <- para[1]
      a1 <- para[2]
      c0 <- para[3]
      c1 <- para[4]
      b0 <- para[5]
      b1<- para[6]
      
      rho <- (exp(2*(b0+b1*age))-1)/(exp(2*(b0+b1*age))+1)
      lambda1 <- exp(a0+a1*age)
      lambda2 <- exp(c0+c1*age) 
      
      df.1 <- d1 & d2     #case part 1  
      df.2 <- d1 & (!d2)  #case part 2
      df.3 <- (!d1)&d2;   #case part 3 
      df.4 <- (!d1)&(!d2) #case part 4
      
      #########################################################
      #################### First Component ####################
      #########################################################
      
      if(sum(df.1)>0){
        
        X.1 <- df[df.1,1]
        Y.1 <- df[df.1,2]
        lambda1.1 <- lambda1[df.1]
        lambda2.1 <- lambda2[df.1]
        S1.1 <- pexp(X.1, rate=lambda1.1) 
        S2.1 <- pexp(Y.1, rate=lambda2.1)
        
        rho.1 <- rho[df.1]
        
        
        #qnorm(S1)=-qnorm(F1) therefore when squaring and multiplying they equal out in part 1
        
        part1 <- sum(-0.5*log(1-rho.1^2)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
                                              rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+
                       log(lambda1.1)-lambda1.1*X.1 +log(lambda2.1)-lambda2.1*Y.1)
      } else {
        part1 <- 0;
      }
      
      
      #########################################################
      ################### Second Component ####################
      #########################################################
      
      if(sum(df.2)>0){
        
        X.2 <- df[df.2,1]
        Y.2 <- df[df.2,2]
        
        lambda1.2 <- lambda1[df.2]
        lambda2.2 <- lambda2[df.2]
        
        S1.2 <- pexp(X.2, rate=lambda1.2) 
        S2.2 <- pexp(Y.2, rate=lambda2.2) 
        
        rho.2 <- rho[df.2]
        
        part2 <- sum(log(pnorm(qnorm(S2.2), mean=rho.2*qnorm(S1.2),
                               sd=sqrt(1-rho.2^2), lower.tail=F)*lambda1.2*exp(-lambda1.2*X.2)))
      } else {
        part2 <- 0;
      }
      
      
      #########################################################
      #################### Third Component ####################
      #########################################################
      
      if(sum(df.3) >0 ){
        
        X.3 <- df[df.3,1]
        Y.3 <- df[df.3,2]
        
        lambda1.3 <- lambda1[df.3]
        lambda2.3 <- lambda2[df.3]
        
        S1.3 <- pexp(X.3, rate=lambda1.3) 
        S2.3 <- pexp(Y.3, rate=lambda2.3) 
        
        rho.3 <- rho[df.3]
        
        part3 <- sum(log(pnorm(qnorm(S1.3), mean=rho.3*qnorm(S2.3),
                               sd=sqrt(1-rho.3^2), lower.tail=F)*lambda2.3*exp(-lambda2.3*Y.3)))
      } else {
        part3 <- 0;
      }
      
      
      #########################################################
      #################### Fourth Component ###################
      #########################################################
      
      if(sum(df.4)>0){
        
        X.4 <- df[df.4,1]
        Y.4 <- df[df.4,2]
        
        lambda1.4 <- lambda1[df.4]
        lambda2.4 <- lambda2[df.4]
        
        S1.4 <- pexp(X.4, rate=lambda1.4) 
        S2.4 <- pexp(Y.4, rate=lambda2.4) 
        
        rho.4 <- rho[df.4]
        
        over.all <- rep(0, length(S1.4))
        for(i in 1:length(over.all)){
          sigma <- matrix(c(1,rho.4[i],rho.4[i],1),nrow=2);
          CDF <- function(V,sigma){
            return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
          }
          over.all[i] <- log(apply(qnorm(cbind(S1.4[i],S2.4[i])),1,CDF,sigma));
        }
        
        part4 <- sum(over.all);
      } else {
        part4 <- 0;
      }
      #print(part4)
      
      #########################################################
      #################### All Components #####################
      ######################################################### 
      
      loglik <- (part1+part2+part3+part4); 
      return(loglik);
    }
    
    a0_lw <- -10
    a0_up <- -2
    a1_lw <- -10
    a1_up <- 2
    c0_lw <- -10
    c0_up <- -2
    c1_lw <- -10
    c1_up <- 2
    
    b0_lw <- 0
    b0_up <- 0.7
    b1_lw <- -0.5
    b1_up <- 0.9
    
    plnoptim_exp <- optim(c(-3,0.01,-3,0.01,0.5,0), npl_exp, method="L-BFGS-B",
                          lower=c(a0_lw,a1_lw,c0_lw,c1_lw,b0_lw,b1_lw),upper=c(a0_up,a1_up,c0_up,c1_up,b0_up,b1_up), 
                          X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                          control=list(fnscale=-1),hessian=TRUE)
    
    
    plnoptim_exp$par
    
    # if(plnoptim_exp$par[1] == a0_lw) {break}
    # if(plnoptim_exp$par[1] == a0_up) {break}
    # if(plnoptim_exp$par[2] == a1_lw) {break}
    # if(plnoptim_exp$par[2] == a1_up) {break}
    # if(plnoptim_exp$par[3] == c0_lw) {break}
    # if(plnoptim_exp$par[3] == c0_up) {break}
    # if(plnoptim_exp$par[4] == c1_lw) {break}
    # if(plnoptim_exp$par[4] == c1_up) {break}
    # if(plnoptim_exp$par[5] == b0_lw) {break}
    # if(plnoptim_exp$par[5] == b0_up) {break}
    # if(plnoptim_exp$par[6] == b1_lw) {break}
    # if(plnoptim_exp$par[6] == b1_up) {break}
    
    ########################################################
    ############### Normal pseudo likelihood ##############
    ####################### Weibull ########################
    ########################################################
    
    npl_wei <- function(para, X,Y,d1,d2,age){
      alpha1 <- para[1]
      x1 <- para[2]
      x2 <- para[3]
      alpha2 <- para[4]
      y1 <- para[5]
      y2 <- para[6]
      b0 <- para[7]
      b1 <- para[8]
      
      rho <- (exp(2*(b0+b1*age))-1)/(exp(2*(b0+b1*age))+1)
      beta1 <- exp(x1+x2*age)
      beta2 <- exp(y1+y2*age)
      
      
      df.1 <- d1 & d2     #case part 1  
      df.2 <- d1 & (!d2)  #case part 2
      df.3 <- (!d1)&d2;   #case part 3 
      df.4 <- (!d1)&(!d2) #case part 4
      
      #########################################################
      #################### First Component ####################
      #########################################################
      
      if(sum(df.1)>0){
        
        
        X.1 <- df[df.1,1]
        Y.1 <- df[df.1,2]
        beta1.1 <- beta1[df.1]
        beta2.1 <- beta2[df.1]
        
        p1.1 <- exp(-beta1.1*X.1^alpha1)
        p2.1 <- exp(-beta2.1*Y.1^alpha2)
        
        p1.1[which(p1.1<0.1^(8))]=0.1^(8)
        p2.1[which(p2.1<0.1^(8))]=0.1^(8)
        
        S1.1 <- 1-p1.1
        S2.1 <- 1-p2.1
        
        
        f1.1 <- beta1.1*alpha1*X.1^(alpha1-1)*exp(-beta1.1*X.1^alpha1) 
        f2.1 <- beta2.1*alpha2*Y.1^(alpha2-1)*exp(-beta2.1*Y.1^alpha2) 
        
        
        rho.1 <- rho[df.1]
        
        
        part1 <- sum(-0.5*log(1-rho.1^2)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
                                              rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+ 
                       log(f1.1)+log(f2.1))
      } else {
        part1 <- 0;
      }
      
      #########################################################
      ################### Second Component ####################
      #########################################################
      
      if(sum(df.2)>0){
        
        X.2 <- df[df.2,1]
        Y.2 <- df[df.2,2]
        beta1.2 <- beta1[df.2]
        beta2.2 <- beta2[df.2]
        
        p1.2 <- exp(-beta1.2*X.2^alpha1)
        p2.2 <- exp(-beta2.2*Y.2^alpha2)
        
        p1.2[which(p1.2<0.1^(8))]=0.1^(8)
        p2.2[which(p2.2<0.1^(8))]=0.1^(8)
        
        S1.2 <- 1-p1.2
        S2.2 <- 1-p2.2
        
        f1.2 <- beta1.2*alpha1*X.2^(alpha1-1)*exp(-beta1.2*X.2^alpha1) 
        f2.2 <- beta2.2*alpha2*Y.2^(alpha2-1)*exp(-beta2.2*Y.2^alpha2) 
        
        rho.2 <- rho[df.2]
        
        part2.1 <- pnorm(qnorm(S2.2), mean=rho.2*qnorm(S1.2),sd=sqrt(1-rho.2^2), lower.tail=F)
        part2.1[which(part2.1<0.1^(10))]=0.1^(10)
        
        part2 <- sum(log(part2.1*f1.2))
      } else {
        part2 <- 0;
      }
      
      #########################################################
      #################### Third Component ####################
      #########################################################
      
      if(sum(df.3) >0 ){
        
        X.3 <- df[df.3,1]
        Y.3 <- df[df.3,2]
        
        beta1.3 <- beta1[df.3]
        beta2.3 <- beta2[df.3]
        
        p1.3 <- exp(-beta1.3*X.3^alpha1)
        p2.3 <- exp(-beta2.3*Y.3^alpha2)
        
        p1.3[which(p1.3<0.1^(8))]=0.1^(8)
        p2.3[which(p2.3<0.1^(8))]=0.1^(8)
        
        S1.3 <- 1-p1.3
        S2.3 <- 1-p2.3
        
        f1.3 <- beta1.3*alpha1*X.3^(alpha1-1)*exp(-beta1.3*X.3^alpha1) 
        f2.3 <- beta2.3*alpha2*Y.3^(alpha2-1)*exp(-beta2.3*Y.3^alpha2) 
        
        rho.3 <- rho[df.3]
        
        part3.1 <- pnorm(qnorm(S1.3), mean=rho.3*qnorm(S2.3),sd=sqrt(1-rho.3^2), lower.tail=F)
        part3.1[which(part3.1<0.1^(10))]=0.1^(10)
        
        part3 <- sum(log(part3.1*f2.3))
        
      } else {
        part3 <- 0;
      }
      
      
      #########################################################
      #################### Fourth Component ###################
      #########################################################
      
      if(sum(df.4)>0){
        
        X.4 <- df[df.4,1]
        Y.4 <- df[df.4,2]
        
        beta1.4 <- beta1[df.4]
        beta2.4 <- beta2[df.4]
        
        S1.4 <- 1-exp(-beta1.4*X.4^alpha1)
        S2.4 <- 1-exp(-beta2.4*Y.4^alpha2)
        
        rho.4 <- rho[df.4]
        
        over.all <- rep(0, length(S1.4))
        for(i in 1:length(over.all)){
          sigma <- matrix(c(1,rho.4[i],rho.4[i],1),nrow=2);
          CDF <- function(V,sigma){
            return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
          }
          part4.1 <- (apply(qnorm(cbind(S1.4[i],S2.4[i])),1,CDF,sigma))
          part4.1[which(part4.1<0.1^(10))]=0.1^(10)
          over.all[i] <- log(part4.1);    }
        
        part4 <- sum(over.all);
      } else {
        part4 <- 0;
      }
      #print(part4)
      
      #########################################################
      #################### All Components #####################
      ######################################################### 
      
      loglik <- (part1+part2+part3+part4); 
      return(loglik);
    }
    
    
    
    a1_lw <- 0.1
    a1_up <- 1.5
    a2_lw <- 0.1
    a2_up <- 1.5
    
    x0_lw <- -5
    x0_up <- -1
    x1_lw <- -2
    x1_up <- 1
    
    y0_lw <- -8
    y0_up <- -3
    y1_lw <- -2.5
    y1_up <- 2.5
    
    b0_lw <- 0.01
    b0_up <- 0.9
    b1_lw <- -0.8
    b1_up <-0.8
    
    plnoptim_wei <- optim(c(0.65, -2.5, -0.6, 0.94, -3.3, -0.9, 0.3, 0.2), npl_wei, method="L-BFGS-B",
                          lower=c(a1_lw, x0_lw, x1_lw, a2_lw, y0_lw, y1_lw, b0_lw, b1_lw),
                          upper=c(a1_up, x0_up, x1_up, a2_up, y0_up, y1_up, b0_up, b1_up), 
                          X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age, control=list(fnscale=-1), hessian=TRUE)
    plnoptim_wei$par
    
    # if(plnoptim_wei$par[1] == a1_lw) {break}
    # if(plnoptim_wei$par[1] == a1_up) {break}
    # if(plnoptim_wei$par[2] == x0_lw) {break}
    # if(plnoptim_wei$par[2] == x0_up) {break}
    # if(plnoptim_wei$par[3] == x1_lw) {break}
    # if(plnoptim_wei$par[3] == x1_up) {break}
    # if(plnoptim_wei$par[4] == a2_lw) {break}
    # if(plnoptim_wei$par[4] == a2_up) {break}
    # if(plnoptim_wei$par[5] == y0_lw) {break}
    # if(plnoptim_wei$par[5] == y0_up) {break}
    # if(plnoptim_wei$par[6] == y1_lw) {break}
    # if(plnoptim_wei$par[6] == y1_up) {break}
    # if(plnoptim_wei$par[7] == b0_lw) {break}
    # if(plnoptim_wei$par[7] == b0_up) {break}
    # if(plnoptim_wei$par[8] == b1_lw) {break}
    # if(plnoptim_wei$par[8] == b1_up) {break}
    # 
    ########################################################
    ############### Normal pseudo likelihood ##############
    ####################### Gompertz #######################
    ########################################################
    
    npl_gom <- function(para, X,Y,d1,d2,age){
      gamma1 <- para[1]
      p0 <- para[2]
      p1 <- para[3]
      gamma2 <- para[4]
      q0 <- para[5]
      q1 <- para[6]
      b0 <- para[7]
      b1 <- para[8]
      
      rho <- (exp(2*(b0+b1*age))-1)/(exp(2*(b0+b1*age))+1)
      lambda1 <- exp(p0+p1*age)
      lambda2 <- exp(q0+q1*age)
      
      
      df.1 <- d1 & d2     #case part 1  
      df.2 <- d1 & (!d2)  #case part 2
      df.3 <- (!d1)&d2;   #case part 3 
      df.4 <- (!d1)&(!d2) #case part 4
      
      #########################################################
      #################### First Component ####################
      #########################################################
      
      if(sum(df.1)>0){
        
        
        X.1 <- df[df.1,1]
        Y.1 <- df[df.1,2]
        lambda1.1 <- lambda1[df.1]
        lambda2.1 <- lambda2[df.1]
        
        q <- exp(-lambda1.1/gamma1*(exp(gamma1*X.1)-1))
        q[which(q<0.1^(8))]=0.1^(8)
        S1.1 <- 1-q
        r <- exp(-lambda2.1/gamma2*(exp(gamma2*Y.1)-1))
        r[which(r<0.1^(8))]=0.1^(8)
        S1.1[which(S1.1<0.1^(8))]=0.1^(8)
        S2.1 <- 1-r
        
        
        f1.1 <- lambda1.1*exp(gamma1*X.1-lambda1.1/gamma1*(exp(gamma1*X.1)-1))
        f2.1 <- lambda2.1*exp(gamma2*Y.1-lambda2.1/gamma2*(exp(gamma2*Y.1)-1))
        f1.1[which(f1.1<0.1^(8))]=0.1^(8)
        f2.1[which(f2.1<0.1^(8))]=0.1^(8)
        
        rho.1 <- rho[df.1]
        
        
        part1 <- sum(-0.5*log(1-rho.1^2)+(((2*rho.1*qnorm(S1.1)*qnorm(S2.1)-
                                              rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2)))/((2*(1-rho.1^2))))+ 
                       log(f1.1)+log(f2.1))
      } else {
        part1 <- 0;
      }
      
      #print(qnorm(S1.1))
      #print(qnorm(S2.1))
      #print(rho.1^2*(qnorm(S1.1)^2 + qnorm(S2.1)^2))
      
      #########################################################
      ################### Second Component ####################
      #########################################################
      
      if(sum(df.2)>0){
        
        X.2 <- df[df.2,1]
        Y.2 <- df[df.2,2]
        lambda1.2 <- lambda1[df.2]
        lambda2.2 <- lambda2[df.2]
        
        q.2 <- exp(-lambda1.2/gamma1*(exp(gamma1*X.2)-1))
        q.2[which(q.2<0.1^(8))]=0.1^(8)
        S1.2 <- 1-q.2
        r.2 <- exp(-lambda2.2/gamma2*(exp(gamma2*Y.2)-1))
        r.2[which(r.2<0.1^(8))]=0.1^(8)
        S2.2 <- 1-r.2
        
        #S1.2 <- 1-exp(-lambda1.2/gamma1*(exp(gamma1*X.2)-1))
        #S2.2 <- 1-exp(-lambda2.2/gamma2*(exp(gamma2*Y.2)-1))
        
        f1.2 <- lambda1.2*exp(gamma1*X.2-lambda1.2/gamma1*(exp(gamma1*X.2)-1))
        f2.2 <- lambda2.2*exp(gamma2*Y.2-lambda2.2/gamma2*(exp(gamma2*Y.2)-1))
        
        rho.2 <- rho[df.2]
        
        part2.1 <- pnorm(qnorm(S2.2), mean=rho.2*qnorm(S1.2),sd=sqrt(1-rho.2^2), lower.tail=F)
        part2.1[which(part2.1<0.1^(10))]=0.1^(10)
        
        part2 <- sum(log(part2.1*f1.2))
      } else {
        part2 <- 0;
      }
      
      
      #########################################################
      #################### Third Component ####################
      #########################################################
      
      if(sum(df.3) >0 ){
        
        X.3 <- df[df.3,1]
        Y.3 <- df[df.3,2]
        
        lambda1.3 <- lambda1[df.3]
        lambda2.3 <- lambda2[df.3]
        
        #S1.3 <- 1-exp(-lambda1.3/gamma1*(exp(gamma1*X.3)-1))
        #S2.3 <- 1-exp(-lambda2.3/gamma2*(exp(gamma2*Y.3)-1))
        
        q.3 <- exp(-lambda1.3/gamma1*(exp(gamma1*X.3)-1))
        q.3[which(q.3<0.1^(8))]=0.1^(8)
        S1.3 <- 1-q.3
        r.3 <- exp(-lambda2.3/gamma2*(exp(gamma2*Y.3)-1))
        r.3[which(r.3<0.1^(8))]=0.1^(8)
        S2.3 <- 1-r.3
        
        f1.3 <- lambda1.3*exp(gamma1*X.3-lambda1.3/gamma1*(exp(gamma1*X.3)-1))
        f2.3 <- lambda2.3*exp(gamma2*Y.3-lambda2.3/gamma2*(exp(gamma2*Y.3)-1))
        
        rho.3 <- rho[df.3]
        
        part3.1 <- pnorm(qnorm(S1.3), mean=rho.3*qnorm(S2.3),sd=sqrt(1-rho.3^2), lower.tail=F)
        part3.1[which(part3.1<0.1^(10))]=0.1^(10)
        
        part3 <- sum(log(part3.1*f2.3))
        
      } else {
        part3 <- 0;
      }
      
      
      #########################################################
      #################### Fourth Component ###################
      #########################################################
      
      if(sum(df.4)>0){
        
        X.4 <- df[df.4,1]
        Y.4 <- df[df.4,2]
        
        lambda1.4 <- lambda1[df.4]
        lambda2.4 <- lambda2[df.4]
        
        S1.4 <- 1-exp(-lambda1.4/gamma1*(exp(gamma1*X.4)-1))
        S2.4 <- 1-exp(-lambda2.4/gamma2*(exp(gamma2*Y.4)-1))
        
        rho.4 <- rho[df.4]
        
        over.all <- rep(0, length(S1.4))
        for(i in 1:length(over.all)){
          sigma <- matrix(c(1,rho.4[i],rho.4[i],1),nrow=2);
          CDF <- function(V,sigma){
            return(pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
          }
          part4.1 <- (apply(qnorm(cbind(S1.4[i],S2.4[i])),1,CDF,sigma))
          part4.1[which(part4.1<0.1^(10))]=0.1^(10)
          over.all[i] <- log(part4.1);    }
        
        part4 <- sum(over.all);
      } else {
        part4 <- 0;
      }
      #print(part4)
      
      #########################################################
      #################### All Components #####################
      ######################################################### 
      #print(part2)
      #print(part3)
      #print(part4)
      loglik <- (part1+part2+part3+part4); 
      return(loglik);
    }
    
    g1_lw <- -0.1
    g1_up <- 0.1
    
    p0_lw <- -5
    p0_up <- -2
    p1_lw <- -2
    p1_up <- 2
    
    g2_lw <- -0.1
    g2_up <- 0.1
    
    q0_lw <- -6
    q0_up <- -3
    q1_lw <- -4
    q1_up <- 2.5
    
    b0_lw <- 0.01
    b0_up <- 0.8
    b1_lw <- -0.8
    b1_up <- 0.8
    
    plnoptim_gom <- optim(c(g1_lw,p0_lw,p1_lw,g2_lw,q0_lw, q1_lw, b0_lw,b1_lw), npl_gom, method="L-BFGS-B",
                          lower=c(g1_lw,p0_lw,p1_lw,g2_lw,q0_lw, q1_lw, b0_lw,b1_lw),
                          upper=c(g1_up,p0_up,p1_up,g2_up,q0_up, q1_up, b0_up,b1_up), 
                          X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age=df$age,
                          control=list(fnscale=-1),hessian=TRUE)
    
    
    plnoptim_gom$par
    
    # if(plnoptim_gom$par[1] == g1_lw) {break}
    # if(plnoptim_gom$par[1] == g1_up) {break}
    # if(plnoptim_gom$par[2] == p0_lw) {break}
    # if(plnoptim_gom$par[2] == p0_up) {break}
    # if(plnoptim_gom$par[3] == p1_lw) {break}
    # if(plnoptim_gom$par[3] == p1_up) {break}
    # if(plnoptim_gom$par[4] == g2_lw) {break}
    # if(plnoptim_gom$par[4] == g2_up) {break}
    # if(plnoptim_gom$par[5] == q0_lw) {break}
    # if(plnoptim_gom$par[5] == q0_up) {break}
    # if(plnoptim_gom$par[6] == q1_lw) {break}
    # if(plnoptim_gom$par[6] == q1_up) {break}
    # if(plnoptim_gom$par[7] == b0_lw) {break}
    # if(plnoptim_gom$par[7] == b0_up) {break}
    # if(plnoptim_gom$par[8] == b1_lw) {break}
    # if(plnoptim_gom$par[8] == b1_up) {break}
    # 
    ########################################################
    ######################### AICS #########################
    ########################################################
    
    loglik_exp <- npl_exp(plnoptim_exp$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
    k_exp <- length(plnoptim_exp$par)
    aic_exp <- -2*loglik_exp+2*k_exp
    
    loglik_wei <- npl_wei(plnoptim_wei$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
    k_wei <- length(plnoptim_wei$par)
    aic_wei <- -2*loglik_wei+2*k_wei
    
    loglik_gom <- npl_gom(plnoptim_gom$par,X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
    k_gom <- length(plnoptim_gom$par)
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
      
      fisher_info <- solve(-plnoptim_exp$hessian) #inverse -hess
      se <- sqrt(diag(fisher_info)) 
      
      #b ci
      est_a0 <- plnoptim_exp$par[1]
      est_a1 <- plnoptim_exp$par[2]
      est_c0 <- plnoptim_exp$par[3]
      est_c1<- plnoptim_exp$par[4]
      est_b0 <- plnoptim_exp$par[5]
      est_b1 <- plnoptim_exp$par[6]
      varb0 <- fisher_info[5,5]
      varb1 <- fisher_info[6,6]
      covb0b1 <- fisher_info[5,6]
      
      #rho for age=0
      theta_d0 <- (exp(2*est_b0)-1)/(exp(2*est_b0)+1)
      var_theta_d0 <- ((16*exp(4*est_b0))/(exp(2*est_b0)+1)^4)*varb0
      theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
      theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
      if (theta_d0_upci>1) {theta_d0_upci=1}
      
      rho_d0_cop <- normalCopula(theta_d0)
      est_rho_d0 <- rho(rho_d0_cop)
      rho_d0_lwci_cop <- normalCopula(theta_d0_lwci)
      rho_d0_lwci <- rho(rho_d0_lwci_cop)
      rho_d0_upci_cop <- normalCopula(theta_d0_upci)
      rho_d0_upci <- rho(rho_d0_upci_cop) 
      save_rho_d0[i] <- est_rho_d0
      
      #rho for age=1
      theta_d1 <- (exp(2*(est_b0+est_b1))-1)/(exp(2*(est_b0+est_b1))+1)
      var_theta_d1 <- ((16*exp(4*(est_b0+est_b1)))/(exp(2*(est_b0+est_b1))+1)^4)*(varb0+2*covb0b1+varb1)
      theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
      theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
      if (theta_d1_upci>1) {theta_d1_upci=1}
      
      rho_d1_cop <- normalCopula(theta_d1)
      est_rho_d1 <- rho(rho_d1_cop)
      rho_d1_lwci_cop <- normalCopula(theta_d1_lwci)
      rho_d1_lwci <- rho(rho_d1_lwci_cop)
      rho_d1_upci_cop <- normalCopula(theta_d1_upci)
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
      
      fisher_info <- solve(-plnoptim_wei$hessian) #inverse -hess
      se <- sqrt(diag(fisher_info)) 
      
      #point and var
      est_a1 <- plnoptim_wei$par[1]
      est_a2 <- plnoptim_wei$par[4]
      est_x0 <- plnoptim_wei$par[2]
      est_x1 <- plnoptim_wei$par[3]
      est_y0 <- plnoptim_wei$par[5]
      est_y1 <- plnoptim_wei$par[6]
      est_b0 <- plnoptim_wei$par[7]
      est_b1 <- plnoptim_wei$par[8]
      varb0 <- fisher_info[7,7]
      varb1 <- fisher_info[8,8]
      covb0b1 <- fisher_info[7,8]
      
      #rho for age=0
      theta_d0 <- (exp(2*est_b0)-1)/(exp(2*est_b0)+1)
      var_theta_d0 <- ((16*exp(4*est_b0))/(exp(2*est_b0)+1)^4)*varb0
      theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
      theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
      if (theta_d0_upci>1) {theta_d0_upci=1}
      
      rho_d0_cop <- normalCopula(theta_d0)
      est_rho_d0 <- rho(rho_d0_cop)
      rho_d0_lwci_cop <- normalCopula(theta_d0_lwci)
      rho_d0_lwci <- rho(rho_d0_lwci_cop)
      rho_d0_upci_cop <- normalCopula(theta_d0_upci)
      rho_d0_upci <- rho(rho_d0_upci_cop) 
      save_rho_d0[i] <- est_rho_d0
      
      #rho for age=1
      theta_d1 <- (exp(2*(est_b0+est_b1))-1)/(exp(2*(est_b0+est_b1))+1)
      var_theta_d1 <- ((16*exp(4*(est_b0+est_b1)))/(exp(2*(est_b0+est_b1))+1)^4)*(varb0+2*covb0b1+varb1)
      theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
      theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
      if (theta_d1_upci>1) {theta_d1_upci=1}
      
      rho_d1_cop <- normalCopula(theta_d1)
      est_rho_d1 <- rho(rho_d1_cop)
      rho_d1_lwci_cop <- normalCopula(theta_d1_lwci)
      rho_d1_lwci <- rho(rho_d1_lwci_cop)
      rho_d1_upci_cop <- normalCopula(theta_d1_upci)
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
      hessian <- hessian(npl_gom, plnoptim_gom$par, X=df$X, Y=df$Y, d1=df$d1, d2=df$d2, age=df$age)
      #fishers info matrix
      fisher_info <- solve(-hessian) 
      #Standard Error
      se <- sqrt(diag(fisher_info)) 
      
      est_b0 <- plnoptim_gom$par[7]
      est_b1 <- plnoptim_gom$par[8]
      est_g1 <- plnoptim_gom$par[1]
      est_g2 <- plnoptim_gom$par[4]
      est_p0 <- plnoptim_gom$par[2]
      est_p1 <- plnoptim_gom$par[3]
      est_q0 <- plnoptim_gom$par[5]
      est_q1 <- plnoptim_gom$par[6]
      varb0 <- fisher_info[7,7]
      varb1 <- fisher_info[8,8]
      covb0b1 <- fisher_info[7,8]
      
      #rho for age=0
      theta_d0 <- (exp(2*est_b0)-1)/(exp(2*est_b0)+1)
      var_theta_d0 <- ((16*exp(4*est_b0))/(exp(2*est_b0)+1)^4)*varb0
      theta_d0_lwci <- theta_d0 - 1.96*sqrt(var_theta_d0)
      theta_d0_upci <- theta_d0 + 1.96*sqrt(var_theta_d0)
      if (theta_d0_upci>1) {theta_d0_upci=1}
      
      rho_d0_cop <- normalCopula(theta_d0)
      est_rho_d0 <- rho(rho_d0_cop)
      rho_d0_lwci_cop <- normalCopula(theta_d0_lwci)
      rho_d0_lwci <- rho(rho_d0_lwci_cop)
      rho_d0_upci_cop <- normalCopula(theta_d0_upci)
      rho_d0_upci <- rho(rho_d0_upci_cop) 
      save_rho_d0[i] <- est_rho_d0
      
      #rho for age=1
      theta_d1 <- (exp(2*(est_b0+est_b1))-1)/(exp(2*(est_b0+est_b1))+1)
      var_theta_d1 <- ((16*exp(4*(est_b0+est_b1)))/(exp(2*(est_b0+est_b1))+1)^4)*(varb0+2*covb0b1+varb1)
      theta_d1_lwci <- theta_d1 - 1.96*sqrt(var_theta_d1)
      theta_d1_upci <- theta_d1 + 1.96*sqrt(var_theta_d1)
      if (theta_d1_upci>1) {theta_d1_upci=1}
      
      rho_d1_cop <- normalCopula(theta_d1)
      est_rho_d1 <- rho(rho_d1_cop)
      rho_d1_lwci_cop <- normalCopula(theta_d1_lwci)
      rho_d1_lwci <- rho(rho_d1_lwci_cop)
      rho_d1_upci_cop <- normalCopula(theta_d1_upci)
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
    
  }
  
  print(i)
} 


print(paste("bias_l1_hr", bias_l1_hr))
print(paste("bias_l2_hr", bias_l2_hr))
print(paste("bias_rho_d0", bias_rho_d0))
print(paste("bias_rho_d1", bias_rho_d1))

print(paste("counter_hr_l1", counter_hr_l1))
print(paste("counter_hr_l2", counter_hr_l2))
print(paste("counter_rho_d0", counter_rho_d0))
print(paste("counter_rho_d1", counter_rho_d1))

print(paste("save_hr_l1", save_hr_l1))
print(paste("save_hr_l2", save_hr_l2))
print(paste("save_rho_d0", save_rho_d0))
print(paste("save_rho_d1", save_rho_d1))

print(paste("counter_exp", counter_exp))
print(paste("counter_wei", counter_wei))
print(paste("counter_gom", counter_gom))

df2 <- data.frame(bias_l1_hr,bias_l2_hr,bias_rho_d0,bias_rho_d1,counter_hr_l1,counter_hr_l2,counter_rho_d0,counter_rho_d1,
                  save_hr_l1,save_hr_l2,save_rho_d0,save_rho_d1,counter_exp,counter_wei,counter_gom)

write.csv(df2, file = paste0(dir_results, out_file_estimates))
end_time = Sys.time()
run_time = end_time - start_time
print(paste0("Simulation2 for normal-gompertz model completed successfully for ", part, "!"))
print(run_time)