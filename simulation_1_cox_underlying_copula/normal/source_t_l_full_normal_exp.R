######################################################################################################
# Simulation study 1: Paper 2 Copula model 2
# data are generated from normal copula - exponential survival model 
# Original code by LS; reviewed, edited and updated by YW for paper2
# YW 31/12/2022 update: 
#    1. change starting values to be between lower and upper bounds
#    2. create a source file to be called by 4 parts of simulation (1-250, 251-500, 501-750, 751-1000)
#    3. put some scalars into vectors and tidy up
#    4. take likelihood function out of the loop
######################################################################################################

start_time <- Sys.time()
set.seed(235452333)
n <- 3000
runs <- 1

#true values from KTX data
true_b0 <- 0.35; true_b1 <- 0.28; true_b2 <- 0; true_b3 <- 0
true_a0 <- -3.30; true_a1 <- 0.11; true_a2 <- 0.02; true_a3 <- -0.51
true_c0 <- -4.15; true_c1 <- 1.32; true_c2 <- -0.11; true_c3  <- -0.65

true_l1 <- true_l2 <- true_t <- true_r <- rep(0,n)
U1 <- V1 <- rep(0,n)

save_a0 <- save_a1 <- save_a2 <- save_a3 <- rep(0,runs)
save_c0 <- save_c1 <- save_c2 <- save_c3 <- rep(0,runs)
save_b0 <- save_b1 <- save_b2 <- save_b3 <- rep(0,runs)
save_se_a0 <- save_se_a1 <- save_se_a2 <- save_se_a3 <- rep(0,runs)
save_se_c0 <- save_se_c1 <- save_se_c2 <- save_se_c3 <- rep(0,runs)
save_se_b0 <- save_se_b1 <- save_se_b2 <- save_se_b3 <- rep(0,runs)
save_var_a0 <- save_var_a1 <- save_var_a2 <- save_var_a3 <- rep(0,runs)
save_var_c0 <- save_var_c1 <- save_var_c2 <- save_var_c3 <- rep(0,runs)
save_var_b0 <- save_var_b1 <- save_var_b2 <- save_var_b3 <- rep(0,runs)

# counting number of times the estimation hit the lower or upper bounds specified for the parameter
# counter_coef_lw <- counter_coef_up <- rep(0,12)

#############################################################################
# Normal pseudo likelihood                                                  #
#############################################################################
likelihood_t_l_full_normal<-function(para, X, Y, d1, d2, age.grp, gen, donor){
  a0 <- para[1]; a1 <- para[2]; a2 <- para[3]; a3 <- para[4]   # regression coefficients for hazard 1 (graft failure)
  c0 <- para[5]; c1 <- para[6]; c2 <- para[7];c3 <- para[8]    # regression coefficients for hazard 2 (death)
  b0 <- para[9]; b1 <- para[10];b2 <- para[11];b3 <- para[12]  # regression coefficients for association parameter
  
  lambda1 <- exp(a0+a1*age.grp+a2*gen+a3*donor)
  lambda2 <- exp(c0+c1*age.grp+c2*gen+c3*donor)
  rho <- (exp(2*(b0+b1*age.grp+b2*gen+b3*donor))-1)/(exp(2*(b0+b1*age.grp+b2*gen+b3*donor))+1)
  
  df.1 <- d1 & d2     #case part 1  
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2;   #case part 3 
  df.4 <- (!d1)&(!d2) #case part 4
  
  ###########################################################################
  # First Component                                                         #
  ###########################################################################
  
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
  
  ##########################################################################
  #Second Component                                                        #
  ##########################################################################
  
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
  
  ###########################################################################
  # Third Component                                                         #
  ###########################################################################
  
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
  
  ###########################################################################
  # Fourth Component                                                        #
  ###########################################################################
  
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
  
  ##########################################################################
  # All Components                                                         #
  ##########################################################################
  
  loglik <- (part1+part2+part3+part4); 
  return(loglik);
}

#################################################################################
# replicate 'runs' times                                                        #
################################################################################

for (i in 1:runs){
  
  ##############################################################################
  #  generate data                                                             #
  ##############################################################################
  #Step 1: generate age categories
  age.grp <- rbinom(n,1,0.40)         
  donor <- rbinom(n,1,0.30)
  gen <- rbinom(n,1,0.38)  
  
  for(k in 1:(n)){#loop to generate U an V from age-varying theta
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
    U1[k]=u              #add to u and v vectors on the outside
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
  
  #Step 10: Create data frame, true values of X and Y have association theta=b0+b1*X
  df<-data.frame(X, Y, d1, d2, age.grp, gen, donor)
  df$X[df$X==0] <- 0.1
  df$Y[df$Y==0] <- 0.1
  
  if (i < rep_start) {print(i)}
  if (i >= rep_start & i <= rep_end){
    # starting values for a0, a1, a2, a3, c0, c1,c2, c3, b0, b1, b2, b3
    starting_values = c(-2, -2, -2, -2,      -4, -2, -2, -2,       0.2,  0.2, 0, 0)
    coef_lw      = c(-10, -10, -10, -10,    -10, -10, -10, -10,    0,    0, -0.3, -0.3)
    coef_up      = c(-1, 0.5, 0.5, 0,       -2.5,  2, 0.5, 0,      0.9, 0.6, 0.3, 0.3)
    plnoptim <- optim(starting_values, likelihood_t_l_full_normal, method="L-BFGS-B",
                      lower = coef_lw,
                      upper = coef_up, 
                      X=df$X, Y=df$Y, d1=df$d1, d2=df$d2,
                      age.grp=df$age.grp, gen=df$gen,donor=df$donor,
                      control=list(fnscale=-1),hessian=TRUE)
    
    plnoptim$par
    
    # temp1 = which(plnoptim$par==counter_coef_lw)
    # counter_coef_lw[temp1] = counter_coef_lw[temp1] + 1
    # 
    # temp2 = which(plnoptim$par==counter_coef_up)
    # counter_coef_up[temp2] = counter_coef_up[temp2] + 1
    
    ############################################################################
    # Confidence Intervals                                                     #
    ############################################################################
    
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
write.csv(estimates, file = paste0(dir_results, out_file_estimates))
print(paste0("Simulation 2 (", rep_start, ",", rep_end, ")", 
             " for normal copula exponential survival model is completed successfully"))
