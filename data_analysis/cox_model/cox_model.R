# 9-July-2021
# YW: NHSBT data analysis, Cox models only
# original script by LS; edited and updated for paper2 by YW
# YW: 2022-02-20: update load data
###############################################################################

rm(list=ls())
library(copula)
library(mvtnorm)
library(ggplot2)
library(plyr)
library(survival)

#as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

########################################################
##################### Load data ########################
########################################################
# need to firstly set working directory to project directory and send through the next two lines
setwd("../../../")
df <- read.csv(file="NHSBT/paper2_data.csv")
# check descriptive statistics
# donor type

plot(survfit(Surv(X, d1) ~ 1, data = df), 
     xlab = "Years since transplant", xlim =c(0,23),
     ylab = "Graft survival probability")

plot(survfit(Surv(Y, d2) ~ 1, data = df), 
     xlab = "Years since transplant", xlim =c(0,23),
     ylab = "Overall survival probability")

dim(df)
table(df$age.grp)
table(df$age.grp)/dim(df)[1]
table(df$gen)
table(df$gen)/dim(df)[1]
table(df$donor)
table(df$donor)/dim(df)[1]

attach(df)

# YW: Cox PH models
names(df)

# YW: Graft failure
#cox_gf <- coxph(Surv(X, d1) ~ age.grp + relevel(as.factor(df$gen),ref="3") + donor, data = df)
start_time= Sys.time()
cox_gf <- coxph(Surv(X, d1) ~ age.grp + gen + donor, data = df)
end_time = Sys.time()
run_time = end_time - start_time
run_time
summary(cox_gf)
extractAIC(cox_gf)

# YW: Death
#cox_death <- coxph(Surv(Y, d2) ~ age.grp +relevel(as.factor(df$gen),ref="3") + donor, data = df)
start_time= Sys.time()
cox_death <- coxph(Surv(Y, d2) ~ age.grp +gen + donor, data = df)
end_time = Sys.time()
run_time = end_time - start_time
run_time
summary(cox_death)
extractAIC(cox_death)