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
