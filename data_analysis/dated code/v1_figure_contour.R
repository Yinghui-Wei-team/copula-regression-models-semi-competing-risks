# create contour plots for paper2
# script by MW, 25/08/2021; edited by YW, 29/09/2021
library(copula)

#setwd("C:\\Users\\mwojtys\\Desktop\\Lexy Sorrell\\Methodology Paper 2\\Simulations\\Real data analysis")

# setwd("C:/Users/ywei3/University of Plymouth/Lexy Sorrell - Lexy's Work/Methodology Paper 2/Countour plots")
# Estimated parameters for the real data (kidney transplant) analysis

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

#----------------------------------------------------------------------
# 1=(Age>50); 1=genFemale; 1=DonorLiving
# male, old, deceased
#age.grp <- 1; gen<- 0; donor <- 0 # patient with highest risks
#pdf("contour_frankweibull_OldMD.pdf",width=4,height=8,paper="special")

# male, old, living
#age.grp <- 1; gen<- 0; donor <- 1 # 
#pdf("contour_frankweibull_OldML.pdf",width=4,height=8,paper="special")

#age.grp <- 1; gen<- 1; donor <- 0 #
#age.grp <- 1; gen<- 1; donor <- 1 #

#----------------------------------------------------------------------
# female, young, living
#age.grp <- 0; gen<- 1; donor <- 1 # patient with lowest risks
#pdf("contour_frankweibull_YoungFL.pdf",width=4,height=8,paper="special")

#female, young, deceased
age.grp <- 0; gen<- 1; donor <- 0 #
pdf("contour_frankweibull_YoungFD.pdf",width=4,height=8,paper="special")

#age.grp <- 0; gen<- 0; donor <- 1 #
#age.grp <- 0; gen<- 0; donor <- 0 #


beta1 <- exp(p0+p1*age.grp+p2*gen+p3*donor)
beta2 <- exp(q0+q1*age.grp+q2*gen+q3*donor)
theta <- exp(b0+b1*age.grp+b2*gen+b3*donor)


# In R, the scale parameter of Weibull distr. is defined differently from our 
# parameterisation so I transform it:

shape1 <- alpha1    # a = alpha
scale1 <- beta1^(-1/alpha1)   # b = beta^(-1/alpha)
  
shape2 <- alpha2    # a = alpha
scale2 <- beta2^(-1/alpha2)    # b = beta^(-1/alpha)
  
  
# Frank copula

frank_cop <- mvdc(copula = frankCopula(param = theta),
                margins = c("weibull","weibull"),
                paramMargins = list(list(shape = shape1, scale = scale1),
                                    list(shape = shape2, scale = scale2))
                )

#pdf("contour_frankweibull.pdf")


par(mar=c(2,2,1,1)+2)

#contour(frank_cop, dMvdc, xlim=c(-0.02,1),ylim=c(-0.02,2),
#        xlab=expression(t[1]), ylab = expression(t[2])) #dMvdc

contourplot2(frank_cop, dMvdc, xlim=c(0,10),ylim=c(0,10),
        xlab="Time to graft failure in years", 
        ylab = "Time to death in years") 

#contourplot2(frank_cop, dMvdc, xlim=c(0,5),ylim=c(0,5),
#            xlab=expression(t[1]), ylab = expression(t[2]), asp=2, region=F) 

dev.off()


par(mfrow=c(1,1))


####### Hazard function


h = function(t,alpha,beta){
  beta*alpha*t^(alpha-1)
}


age.grp <- 0; gen<- 1; donor <- 1 # patient with lowest risks
beta1.lowrisk <- exp(p0+p1*age.grp+p2*gen+p3*donor)
beta2.lowrisk <- exp(q0+q1*age.grp+q2*gen+q3*donor)

age.grp <- 1; gen<- 0; donor <- 0 # patient with highest risks
beta1.highrisk <- exp(p0+p1*age.grp+p2*gen+p3*donor)
beta2.highrisk <- exp(q0+q1*age.grp+q2*gen+q3*donor)

# For non-terminal event
pdf("results/hazard_frankweibull_NT.pdf",width=4,height=3.5,paper="special")
curve(h(x,alpha=alpha1,beta=beta1.lowrisk), from = 0, to = 10, ylim = c(0,0.25), 
      xlab="Time since transplant in years", ylab = "Hazard function ")
curve(h(x,alpha=alpha1,beta=beta1.highrisk), from = 0, to = 10, add=T, lty = 2)
legend("topright", 
       legend = c("male above 50 with deaceased donor","female below 50 with living donor"), 
       lty = c(2,1))
dev.off()

# For terminal event
pdf("results/hazard_frankweibull_T.pdf",width=4,height=3.5,paper="special")
curve(h(x,alpha=alpha2,beta=beta2.lowrisk), from = 0, to = 10, ylim = c(0,0.25), 
      xlab="Time since transplant in years", 
      ylab = "Hazard function")
curve(h(x,alpha=alpha2,beta=beta2.highrisk), from = 0, to = 10, add= T, lty = 2)
       legend("topright", 
       legend = c("male above 50 with deaceased donor","female below 50 with living donor"), 
       lty = c(2,1))
dev.off()
