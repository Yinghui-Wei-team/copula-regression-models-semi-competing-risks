#################################################################
# Part 1. Read in parameters                                    #
################################################################
# Parameters of the Weibull distribution for the non-terminal event
alpha1 <- 0.70
# Regression parameters (a0, a1, a2, a3) for graft failures
p0 <- -2.56  #  
p1 <- 0.18    # Age>50
#p2 <- 0.01   # SexF, the CI contains 0
p2 <- 0      # SexF, 0.011 but the CI contains 0
p3 <- -0.58  # DonorL

# Parameters of the Weibull distribution for the terminal event
alpha2 <- 0.99
# Regression parameters (c0,c1,c2,c3) for death
q0 <- -4.03
q1 <- 1.32   # Age>50
q2 <- -0.08   # SexF
q3 <- -0.65  # DonorL

# Regression parameters for theta for Frank copula
b0 <- 3.20
b1 <- 3.95   # Age>50
#b2 <- 0.053   # SexF, the CI contains 0
b2 <- 0       # SexF, 0.053 but the CI contains 0
b3 <- 0   # DonorL