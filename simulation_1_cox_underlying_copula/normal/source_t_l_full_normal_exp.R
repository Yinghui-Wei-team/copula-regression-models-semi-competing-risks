###################################################################################################
# Simulation study: Paper 2 Copula model 2, 
# normal copula exponential distribution
# YW 31/12/2022 update: 
#    1. change starting values to be between lower and upper bounds
#    2. create a source file to be called by 4 parts of simulation (1-250, 251-500, 501-750, 751-1000)
###################################################################################################

start_time <- Sys.time()
set.seed(235452333)
n <- 3000
runs <- 1

#true values from KTX data
true_b0 <- 0.35
true_b1 <- 0.28
true_b2 <- 0
true_b3 <- 0
true_a0 <- -3.30
true_a1 <- 0.11
true_a2 <- 0.02
true_a3 <- -0.51
true_c0 <- -4.15
true_c1 <- 1.32
true_c2 <- -0.11
true_c3  <- -0.65

true_l1 <- rep(0,n)
true_l2 <- rep(0,n)
true_t <- rep(0,n)
true_r <- rep(0,n)
U1 <- rep(0,n)
V1 <- rep(0,n)

save_a0 <- rep(0,runs)
save_a1 <- rep(0,runs)
save_a2 <- rep(0,runs)
save_a3 <- rep(0,runs)
save_c0 <- rep(0,runs)
save_c1 <- rep(0,runs)
save_c2 <- rep(0,runs)
save_c3 <- rep(0,runs)
save_b0 <- rep(0,runs)
save_b1 <- rep(0,runs)
save_b2 <- rep(0,runs)
save_b3 <- rep(0,runs)
save_se_a0 <- rep(0,runs)
save_se_a1 <- rep(0,runs)
save_se_a2 <- rep(0,runs)
save_se_a3 <- rep(0,runs)
save_se_c0 <- rep(0,runs)
save_se_c1 <- rep(0,runs)
save_se_c2 <- rep(0,runs)
save_se_c3 <- rep(0,runs)
save_se_b0 <- rep(0,runs)
save_se_b1 <- rep(0,runs)
save_se_b2 <- rep(0,runs)
save_se_b3 <- rep(0,runs)
save_var_a0 <- rep(0,runs)
save_var_a1 <- rep(0,runs)
save_var_a2 <- rep(0,runs)
save_var_a3 <- rep(0,runs)
save_var_c0 <- rep(0,runs)
save_var_c1 <- rep(0,runs)
save_var_c2 <- rep(0,runs)
save_var_c3 <- rep(0,runs)
save_var_b0 <- rep(0,runs)
save_var_b1 <- rep(0,runs)
save_var_b2 <- rep(0,runs)
save_var_b3 <- rep(0,runs)

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

###############################################################
###################### run 'runs' times #######################
###############################################################

for (i in 1:runs){
  
  ###############################################################
  ######################## generate data ########################
  ###############################################################
  #Step 1: generate age categories
  age.grp <- rbinom(n,1,0.40)         
  donor <- rbinom(n,1,0.30)
  gen <- rbinom(n,1,0.38)  
  
  for(k in 1:(n)){   #loop to generate U an V from age-varying theta
    m=1                  
    
    #Step 2: generate 1 random variable from Uniform(0,a) distribution 
    u1 <- runif(m,0,1)       
    
    #Step 3: X_true generated from u1 values (T1 from later)
    theta1 <- (exp(2*(true_b0+true_b1*age.grp[k]+true_b2*gen[k]+true_b3*donor[k]))-1)/(exp(2*(true_b0+true_b1*age.grp[k]+true_b2*gen[k]+true_b3*donor[k]))+1)
    true_l1s <- exp(true_a0 + true_a1*age.grp[k] + true_a2*gen[k] + true_a3*donor[k]) 
    true_l2s <- exp(true_c0 + true_c1*age.grp[k] + true_c2*gen[k] + true_c3*donor[k])
    
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
  
  if (i < rep_start) {print(i)}
  if (i >= rep_start & i <= rep_end){
    #########################################################
    ############### Normal pseudo likelihood ################
    #########################################################
    npl<-function(para, X, Y, d1, d2, age.grp, gen, donor){
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
      rho <- (exp(2*(b0+b1*age.grp+b2*gen+b3*donor))-1)/(exp(2*(b0+b1*age.grp+b2*gen+b3*donor))+1)
      
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
        
        part2.1 <- pnorm(qnorm(S2.2), mean=rho.2*qnorm(S1.2),sd=sqrt(1-rho.2^2), lower.tail=F)
        part2.1[which(part2.1<0.1^(10))]=0.1^(10)
        part2 <- sum(log(part2.1*lambda1.2*exp(-lambda1.2*X.2)))
        
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
        
        part3.1 <- pnorm(qnorm(S1.3), mean=rho.3*qnorm(S2.3),   sd=sqrt(1-rho.3^2), lower.tail=F)
        part3.1[which(part3.1<0.1^(10))]=0.1^(10)
        part3 <- sum(log(part3.1*lambda2.3*exp(-lambda2.3*Y.3)))
        
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
    
    a0_lw <- -10; a0_up <- -1
    a1_lw <- -10; a1_up <- 0.5
    a2_lw <- -10; a2_up <- 0.5
    a3_lw <- -10; a3_up <- 0
    
    c0_lw <- -10; c0_up <- -2.5
    c1_lw <- -10; c1_up <- 2
    c2_lw <- -10; c2_up <- 0.5
    c3_lw <- -10; c3_up <- 0
    
    b0_lw <- 0; b0_up <- 0.9
    b1_lw <- 0; b1_up <- 0.6
    b2_lw <- -0.3; b2_up <- 0.3
    b3_lw <- -0.3;b3_up <- 0.3
    
    # YW 2/Sept/2021, changed starting values for a0, a1, a2, a3, c0, c1,c2, c3, b0, b1, b2, b3
    starting_values = c(-2, -2, -2, -2,
                        -4, -2, -2, -2,
                        0.2, 0.2, 0, 0)
 
    plnoptim <- optim(starting_values, npl, method="L-BFGS-B",
                      lower=c(a0_lw,a1_lw,a2_lw,a3_lw,c0_lw,c1_lw,c2_lw,c3_lw,b0_lw,b1_lw,b2_lw,b3_lw),upper=c(a0_up,a1_up,a2_up,a3_up,c0_up,c1_up,c2_up,c3_up,b0_up,b1_up,b2_up,b3_up), 
                      X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,age.grp=df$age.grp, gen=df$gen,donor=df$donor,
                      control=list(fnscale=-1),hessian=TRUE)
    
    plnoptim$par
    
    if(plnoptim$par[1] == a0_lw) {counter_a0_low = counter_a0_low + 1}
    if(plnoptim$par[1] == a0_up) {counter_a0_upper = counter_a0_upper + 1}
    if(plnoptim$par[2] == a1_lw) {counter_a1_low = counter_a1_low + 1}
    if(plnoptim$par[2] == a1_up) {counter_a1_upper = counter_a1_upper + 1}
    if(plnoptim$par[3] == a2_lw) {counter_a2_low = counter_a2_low + 1}
    if(plnoptim$par[3] == a2_up) {counter_a2_upper = counter_a2_upper + 1}
    if(plnoptim$par[4] == a3_lw) {counter_a3_low = counter_a3_low + 1}
    if(plnoptim$par[4] == a3_up) {counter_a3_upper = counter_a3_upper + 1}
    
    if(plnoptim$par[5] == c0_lw) {counter_c0_low = counter_c0_low + 1}
    if(plnoptim$par[5] == c0_up) {counter_c0_upper = counter_c0_upper + 1}
    if(plnoptim$par[6] == c1_lw) {counter_c1_low = counter_c1_low + 1}
    if(plnoptim$par[6] == c1_up) {counter_c1_upper = counter_c1_upper + 1}
    if(plnoptim$par[7] == c2_lw) {counter_c2_low = counter_c2_low + 1}
    if(plnoptim$par[7] == c2_up) {counter_c2_upper = counter_c2_upper + 1}
    if(plnoptim$par[8] == c3_lw) {counter_c3_low = counter_c3_low + 1}
    if(plnoptim$par[8] == c3_up) {counter_c3_upper = counter_c3_upper + 1}
    
    if(plnoptim$par[9] == b0_lw) {counter_b0_low = counter_b0_low + 1}
    if(plnoptim$par[9] == b0_up) {counter_b0_upper = counter_b0_upper + 1}
    if(plnoptim$par[10] == b1_lw) {counter_b1_low = counter_b1_low + 1}
    if(plnoptim$par[10] == b1_up) {counter_b1_upper = counter_b1_upper + 1}
    if(plnoptim$par[11] == b2_lw) {counter_b2_low = counter_b2_low + 1}
    if(plnoptim$par[11] == b2_up) {counter_b2_upper = counter_b2_upper + 1}
    if(plnoptim$par[12] == b3_lw) {counter_b3_low = counter_b3_low + 1}
    if(plnoptim$par[12] == b3_up) {counter_b3_upper = counter_b3_upper + 1}
    
    ########################################################
    ################## Confidence Intervals ################
    ########################################################
    
    fisher_info <- solve(-plnoptim$hessian) #inverse -hess
    #Standard error = sqrt(var/n)
    se <- sqrt(diag(fisher_info)) 
    save_a0[i] <- plnoptim$par[1]
    save_a1[i] <- plnoptim$par[2]
    save_a2[i] <- plnoptim$par[3]
    save_a3[i] <- plnoptim$par[4]
    save_c0[i] <- plnoptim$par[5]
    save_c1[i] <- plnoptim$par[6]
    save_c2[i] <- plnoptim$par[7]
    save_c3[i] <- plnoptim$par[8]
    save_b0[i] <- plnoptim$par[9]
    save_b1[i] <- plnoptim$par[10]
    save_b2[i] <- plnoptim$par[11]
    save_b3[i] <- plnoptim$par[12]
    save_se_a0[i] <- se[1]
    save_se_a1[i] <- se[2]
    save_se_a2[i] <- se[3]
    save_se_a3[i] <- se[4]
    save_se_c0[i] <- se[5]
    save_se_c1[i] <- se[6]
    save_se_c2[i] <- se[7]
    save_se_c3[i] <- se[8]
    save_se_b0[i] <- se[9]
    save_se_b1[i] <- se[10]
    save_se_b2[i] <- se[11]
    save_se_b3[i] <- se[12]
    save_var_a0[i] <- fisher_info[1,1]
    save_var_a1[i] <- fisher_info[2,2]
    save_var_a2[i] <- fisher_info[3,3]
    save_var_a3[i] <- fisher_info[4,4]
    save_var_c0[i] <- fisher_info[5,5]
    save_var_c1[i] <- fisher_info[6,6]
    save_var_c2[i] <- fisher_info[7,7]
    save_var_c3[i] <- fisher_info[8,8]
    save_var_b0[i] <- fisher_info[9,9]
    save_var_b1[i] <- fisher_info[10,10]
    save_var_b2[i] <- fisher_info[11,11]
    save_var_b3[i] <- fisher_info[12,12]
  }
  print(i)
}
estimates <- data.frame(save_a0, save_se_a0, save_var_a0,
                        save_a1, save_se_a1, save_var_a1,
                        save_a2, save_se_a2, save_var_a2,
                        save_a3, save_se_a3, save_var_a3,
                        save_c0, save_se_c0, save_var_c0,
                        save_c1, save_se_c1, save_var_c1,
                        save_c2, save_se_c2, save_var_c2,
                        save_c3, save_se_c3, save_var_c3,
                        save_b0, save_se_b0, save_var_b0,
                        save_b1, save_se_b1, save_var_b1,
                        save_b2, save_se_b2, save_var_b2,
                        save_b3, save_se_b3, save_var_b3)
end_time <- Sys.time()

run_time = end_time - start_time

run_time

# Output results
out_file_estimates <- paste0("s1_model2_t_l_full_normal_exp (", rep_start,"-", rep_end, ").csv")
# if on own PC
write.csv(estimates, file = paste0(dir_results, out_file_estimates))