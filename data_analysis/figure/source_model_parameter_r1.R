#################################################################
# Part 1. Read in parameters                                    #
################################################################
# Parameters of the Weibull distribution for the non-terminal event
alpha1 <- 0.66
# Regression parameters for beta1: non-terminal event
p0 <- -2.48  #  
p1 <- 0.11    # Age>50
#p2 <- 0.011   # SexF, the CI contains 0
p2 <- 0       # SexF, 0.011 but the CI contains 0
p3 <- -0.59  # DonorL

# Parameters of the Weibull distribution for the terminal event
alpha2 <- 1
# Regression parameters for beta2
q0 <- -4.05
q1 <- 1.42   # Age>50
q2 <- -0.12   # SexF
q3 <- -0.71  # DonorL

# Regression parameters for theta for Frank copula
b0 <- 1
b1 <- 3.23   # Age>50
#b2 <- 0.053   # SexF, the CI contains 0
b2 <- -0.40       # SexF, 0.053 but the CI contains 0
b3 <- 0   # DonorL