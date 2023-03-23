##################################################################################
# likelihood function specification  - Clayton Copula                            #
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

# clayton copula Weibull distribution: one covariates
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

# clayton copula Gompertz distribution: one covariates
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


##################################################################################
# likelihood function specification  - Frank Copula                              #
##################################################################################
# Frank Exponential 
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

#Frank Weibull 
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

## Frank Gompertz 
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


##################################################################################
# likelihood function specification  - Gumbel Copula                             #
##################################################################################
#GUMBEL-exponential pseudo likelihood 
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
  
  return(logpl)
}

#Gumbel-Gompertz 
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
  
  f1 <- lambda1*exp(gamma1*X-lambda1/gamma1*(exp(gamma1*X)-1))
  f2 <- lambda2*exp(gamma2*Y-lambda2/gamma2*(exp(gamma2*Y)-1))
  
  C=exp(-((-log(S1))^(theta)+(-log(S2))^(theta))^(1/theta))
  
  part1 <- d1*d2*(log(C)+(theta-1)*log(-log(S1))+(theta-1)*log(-log(S2))+log(theta-1+((-log(S1))^theta+(-log(S2))^theta)^(1/theta))-log(S1)-log(S2)-(2*theta-1)*log(-log(C))+log(f1)+log(f2))
  part2 <- d1*(1-d2)*(log(C)+(theta-1)*log(-log(S1))-log(S1)-(theta-1)*log(-log(C))+log(f1))
  part3 <-((1-d1)*(d2))*(log(C)+(theta-1)*log(-log(S2))-log(S2)-(theta-1)*log(-log(C))+log(f2))
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  
  logpl <- sum(part1+part2+part3+part4) 
  return(logpl)
}

#Gumbel-Weibull 
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
  f1 <- beta1*alpha1*X^(alpha1-1)*exp(-beta1*X^alpha1) 
  f2 <- beta2*alpha2*Y^(alpha2-1)*exp(-beta2*Y^alpha2) 
  subLogS1 <- log(S1) 
  subLogS1[which(subLogS1 > -0.1^(7))]= -0.1^7
  subLogS2<- log(S2) 
  subLogS2[which(subLogS2 > -0.1^(7))]= -0.1^7
  
  C=exp(-((-subLogS1)^(theta)+(-subLogS2)^(theta))^(1/theta))
  
  C[which(C<0.1^(8))]=0.1^(8)
  
  part1 <- d1*d2*(log(C*(-subLogS1)^(theta-1)*(-subLogS2)^(theta-1)*(theta-1-log(C))*f1*f2)-log(S1*S2*(-log(C))^(2*theta-1)))
  part2 <- d1*(1-d2)*(log(C)+(theta-1)*log(-subLogS1)-subLogS1-(theta-1)*log(-log(C))+log(f1)) 
  part3 <-((1-d1)*(d2))*(log(C)+(theta-1)*log(-subLogS2)-subLogS2-(theta-1)*log(-log(C))+log(f2))
  part4<-((1-d1)*(1-d2))*log(C)
  logpl<-sum(part1+part2+part3+part4) 
  return(logpl)
}



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
