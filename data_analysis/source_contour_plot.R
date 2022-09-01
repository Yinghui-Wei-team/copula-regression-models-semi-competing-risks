pdf(outfile,width=4,height=8,paper="special")
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

p <- contourplot2(frank_cop, dMvdc, xlim=c(0,10),ylim=c(0,10),
            col.regions = colorRampPalette(c("royalblue3", "maroon3"), space="Lab"),
             xlab="Time to graft failure in years", 
             ylab = "Time to death in years") 

print(p)

dev.off()