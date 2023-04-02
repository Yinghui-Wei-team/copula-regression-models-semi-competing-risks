###############################################################################
# Purpose: NHSBT data analysis, Cox models only
# Date: 2021-July-09
# Programmed by: original script by LS; edited and updated for paper2 by YW
# YW: 2022-02-20: update load data
###############################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(ggplot2);library(plyr)
library(survival); library(dplyr);library(ggsurvfit);library(readr);library(survminer)

################################################################################
#Part 1. Load data                                                             #
################################################################################
# need to firstly set working directory to project directory and then send through the next two lines
dir_data <-dir_results <- "../../"

df <- read_rds(file=paste0(dir_data, "NHSBT/paper2_data_v2.rds"))
#df <- read.csv(file=paste0(dir_data, "NHSBT/paper2_data.csv"))
################################################################################
#Part 2. KM Curves                                                             #
################################################################################
# graft survival
km_gs <- survfit2(Surv(X, d1) ~ 1, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nDays",
    y = "Graft survival\n"
  ) + xlim(0, 25)

km_gs
ggsave(file="km_gs.pdf", 
       path = paste0(dir_results, "results/real_data_analysis/figures"),
       plot=km_gs)

# graft survival: KM curves stratified by age group
km_gs_agegrp <- survfit2(Surv(X, d1) ~ age.grp, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nYears since transplant",
    y = "Graft survival\n"
  ) + xlim(0, 25)

km_gs_agegrp
ggsave(file="km_gs_agegrp.pdf", 
       path = paste0(dir_results, "results/real_data_analysis/figures"),
       plot=km_gs_agegrp)

# overall survival
km_os <- survfit2(Surv(Y, d2) ~ 1, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nYears since transplant",
    y = "Overall survival\n"
  ) + xlim(0, 25) + ylim(0, 1)

km_os
ggsave(file="km_os.pdf", 
       path = paste0(dir_results, "results/real_data_analysis/figures"),
       plot=km_os)

# Overall survival: KM curves stratified by age group
km_os_agegrp <- survfit2(Surv(X, d2) ~ age.grp, data = df) %>%
  ggsurvfit() +
  labs(
    x = "\nYears since transplant",
    y = "Overall survival\n"
  ) + xlim(0, 25) + ylim(0, 1)

km_os_agegrp
ggsave(file="km_os_agegrp.pdf", 
       path = paste0(dir_results, "results/real_data_analysis/figures"),
       plot=km_os_agegrp)

df1 <- df[which(df$age.grp  =="<=50 years"),]
index <- which(df1$X>23&df1$d1==0)
length(df1[index,])
df1[index,]

index1 <- which(df1$X>23&df1$d1==1)
length(df1[index1,])
df1[index1,]

df1 <- df[which(df$age.grp  ==">50 years"),]
index <- which(df1$X>23&df1$d1==0)
length(df1[index,])
df1[index,]

index1 <- which(df1$X>23&df1$d1==1)
length(df1[index1,])
df1[index1,]
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
test_ph_gf = cox.zph(cox_gf)
test_ph_gf
ggcoxzph(test_ph_gf)
end_time = Sys.time()
run_time = end_time - start_time
run_time
summary(cox_gf)
extractAIC(cox_gf)

# Death
start_time= Sys.time()
cox_death <- coxph(Surv(Y, d2) ~ age.grp +gen + donor, data = df)
test_ph_death = cox.zph(cox_death)
test_ph_death
ggcoxzph(test_ph_death)
end_time = Sys.time()
run_time = end_time - start_time
run_time
summary(cox_death)
extractAIC(cox_death)