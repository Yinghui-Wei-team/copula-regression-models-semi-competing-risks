###############################################################################
# Purpose: NHSBT data analysis, Cox models only
# Date: 2021-July-09
# Programmed by: original script by LS; edited and updated for paper2 by YW
# YW: 2022-02-20: update load data
###############################################################################

rm(list=ls())
library(copula)
library(mvtnorm)
library(ggplot2)
library(plyr)
library(survival)
library(dplyr)
library(ggsurvfit) 

################################################################################
#Part 1. Load data                                                             #
################################################################################
# need to firstly set working directory to project directory and then send through the next two lines
dir_data <- "../../"
df <- read.csv(file=paste0(dir_data, "NHSBT/paper2_data_v2.csv"))

df <- df%>% 
  # age.grp = 0 "<=50 years", 1 ">50 years"
  mutate(age.grp = factor(age.grp)) %>%
  # gen = 1 female, 0 male
  mutate(gen = factor(gen)) %>%
  # donor = 1 living, 0 deceased
  mutate(donor = factor(donor))

################################################################################
#Part 2. KM Curves                                                             #
################################################################################
# graft survival
survfit2(Surv(X, d1) ~ 1, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nDays",
    y = "Graft survival\n"
  ) + xlim(0, 25)

# graft survival: KM curves stratified by age group
survfit2(Surv(X, d1) ~ age.grp, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nDays",
    y = "Graft survival\n"
  ) + xlim(0, 25)

# overall survival
survfit2(Surv(Y, d2) ~ 1, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nDays",
    y = "Overall survival\n"
  ) + xlim(0, 25) + ylim(0, 1)

# Overall survival: KM curves stratified by age group
survfit2(Surv(X, d2) ~ age.grp, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nDays",
    y = "Overall survival\n"
  ) + xlim(0, 25) + ylim(0, 1)

################################################################################
#Part 3. Summary Statistics                                                    #
################################################################################
dim(df)
table(df$age.grp)
table(df$age.grp)/dim(df)[1]
table(df$gen)
table(df$gen)/dim(df)[1]
table(df$donor)
table(df$donor)/dim(df)[1]

################################################################################
#Part 4. Cox PH Models                                                         #
################################################################################
names(df)
# Graft failure
start_time= Sys.time()
cox_gf <- coxph(Surv(X, d1) ~ age.grp + gen + donor, data = df)
end_time = Sys.time()
run_time = end_time - start_time
run_time
summary(cox_gf)
extractAIC(cox_gf)

# Death
start_time= Sys.time()
cox_death <- coxph(Surv(Y, d2) ~ age.grp +gen + donor, data = df)
end_time = Sys.time()
run_time = end_time - start_time
run_time
summary(cox_death)
extractAIC(cox_death)