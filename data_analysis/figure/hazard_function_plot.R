################################################################################
# Purpose: create contour plots and hazard function plots for paper2
# Programmed by Malgorzata Wojtys and Yinghui Wei
# First created on 25/08/2021
# YW updated 31/12/2022: use latest regression coefficients after re-running th analyses
################################################################################
library(copula)

output_dir <- "../../results/real_data_analysis/figures/"

#source("data_analysis/source_model_parameter.R")

source("data_analysis/figure/source_model_parameter_r1.R")

#Hazard function--------------------------------------------------------
h = function(t,alpha,beta){
  beta*alpha*t^(alpha-1)
}
age.grp <- 0; gen<- 1; donor <- 1 # patient with lowest risks
beta1.lowrisk <- exp(p0+p1*age.grp+p2*gen+p3*donor)
beta2.lowrisk <- exp(q0+q1*age.grp+q2*gen+q3*donor)

age.grp <- 1; gen<- 0; donor <- 0 # patient with highest risks
beta1.highrisk <- exp(p0+p1*age.grp+p2*gen+p3*donor)
beta2.highrisk <- exp(q0+q1*age.grp+q2*gen+q3*donor)
#pdf("results/real_data_analysis/hazard.pdf",width=10,height=7,paper="special")
#par(mfrow=c(1,2))
# For non-terminal event
pdf(paste0(output_dir, "hazard_frankweibull_NT.pdf"),width=5,height=5,paper="special")
curve(h(x,alpha=alpha1,beta=beta1.lowrisk), from = 0, to = 10, ylim = c(0,0.25), lwd = 2,
      col = "maroon3",
      xlab="Time since transplant in years", ylab = "Hazard function "
      #main = "Kidney graft failure"
      )
curve(h(x,alpha=alpha1,beta=beta1.highrisk), from = 0, to = 10, add=T, lty = 2, lwd = 2, 
      col = "royalblue3")
legend("topright", 
       legend = c("male above 50 with deaceased donor","female below 50 with living donor"), 
       lty = c(2,1), lwd=c(2,2), col = c("royalblue3", "maroon3"))
dev.off()

# For terminal event
pdf(paste0(output_dir,"hazard_frankweibull_T.pdf"),width=5,height=5,paper="special")
curve(h(x,alpha=alpha2,beta=beta2.lowrisk), from = 0, to = 10, ylim = c(0,0.25), lwd = 2,
      col = "maroon3",
      xlab="Time since transplant in years", 
      ylab = "Hazard function"
      #main = "Death"
      )
curve(h(x,alpha=alpha2,beta=beta2.highrisk), from = 0, to = 10, add= T, lty = 2, lwd = 2,
      col = "royalblue3")
legend("topright", 
       legend = c("male above 50 with deaceased donor","female below 50 with living donor"), 
       lty = c(2,1), lwd=c(2,2), col = c("royalblue3", "maroon3"))
dev.off()
