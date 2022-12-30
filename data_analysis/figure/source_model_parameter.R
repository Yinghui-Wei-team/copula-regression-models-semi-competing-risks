#################################################################
# Part 1. Read in parameters                                    #
################################################################
# Parameters of the Weibull distribution for the non-terminal event
alpha1 <- 0.696
# Regression parameters for beta1
p0 <- -1.998  #  
p1 <- 0.19    # Age>50
#p2 <- 0.011   # SexF, the CI contains 0
p2 <- 0       # SexF, 0.011 but the CI contains 0
p3 <- -0.577  # DonorL

# Parameters of the Weibull distribution for the terminal event
alpha2 <- 0.986
# Regression parameters for beta2
q0 <- -3.294
q1 <- 1.317   # Age>50
q2 <- -0.08   # SexF
q3 <- -0.651  # DonorL

# Regression parameters for theta for Frank copula
b0 <- 1
b1 <- 0.775   # Age>50
#b2 <- 0.053   # SexF, the CI contains 0
b2 <- 0       # SexF, 0.053 but the CI contains 0
b3 <- 0.141   # DonorL