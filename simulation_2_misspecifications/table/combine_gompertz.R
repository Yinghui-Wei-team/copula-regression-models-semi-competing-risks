# To combine results from simulation 2 AIC normal exponential model

rm(list=ls())
library(copula)
library(mvtnorm)
library(plyr)
library(survival)
library(numDeriv)

# directory if on own PC
dir_results <- "../../results/simulation_results/"
out_file_summary <- "simulation2/normal_gomp/S2_aic_normal_gomp_summary.csv"

#####################################################################################
################## Normal, age, gen from gom chose with aic #######################
#####################################################################################
set.seed(7877320)
n <- 3000
runs <- 1000

true_b0 <- 0.36
true_b1 <- 0.28

true_g1 <- 0.001 #gomp gamma1
true_g2 <- 0.06  #gomp gamma2
true_p0 <- -3.37 #gomp lambda1
true_p1 <- 0.14 #gomp lambda1
true_q0 <- -4.79 #gomp lambda2
true_q1 <- 1.49 #gomp lambda2



true_theta_d0 <- (exp(2*true_b0)-1)/(exp(2*true_b0)+1)
true_theta_d1 <- (exp(2*(true_b0+true_b1))-1)/(exp(2*(true_b0+true_b1))+1)

t_theta_d0_cop <- normalCopula(true_theta_d0)
true_rho_d0 <- rho(t_theta_d0_cop)

t_theta_d1_cop <- normalCopula(true_theta_d1)
true_rho_d1 <- rho(t_theta_d1_cop)

#S1 exp, S2 weib
true_hr_l1 <- exp(true_p1)
true_hr_l2 <- exp(true_q1)

# normal gompertz
df1 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part1.csv"))
df2 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part2.csv"))
df3 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part3.csv"))
df4 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part4.csv"))
df5 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part5.csv"))
df6 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part6.csv"))
df7 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part7.csv"))
df8 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part8.csv"))
df9 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part9.csv"))
df10 <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_part10.csv"))

df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)


bias_l1_hr <- df$bias_l1_hr
bias_l2_hr <- df$bias_l2_hr
bias_rho_d0 <- df$bias_rho_d0
bias_rho_d1 <- df$bias_rho_d1

total_counter_hr_l1 <- df$counter_hr_l1
total_counter_hr_l2 <- df$counter_hr_l2
total_counter_rho_d0 <- df$counter_rho_d0
total_counter_rho_d1 <- df$counter_rho_d1

counter_hr_l1 <- sum(total_counter_hr_l1[1]+total_counter_hr_l1[101]+total_counter_hr_l1[201]+
                       total_counter_hr_l1[301]+total_counter_hr_l1[401]+
                       total_counter_hr_l1[501]+total_counter_hr_l1[601]+
                       total_counter_hr_l1[701]+total_counter_hr_l1[801]+total_counter_hr_l1[901])

counter_hr_l2 <- sum(total_counter_hr_l2[1]+total_counter_hr_l2[101]+total_counter_hr_l2[201]+
                       total_counter_hr_l2[301]+total_counter_hr_l2[401]+
                       total_counter_hr_l2[501]+total_counter_hr_l2[601]+
                       total_counter_hr_l2[701]+total_counter_hr_l2[801]+total_counter_hr_l2[901])

counter_rho_d0  <- sum(total_counter_rho_d0[1]+total_counter_rho_d0[101]+total_counter_rho_d0[201]+
                         total_counter_rho_d0[301]+total_counter_rho_d0[401]+
                         total_counter_rho_d0[501]+total_counter_rho_d0[601]+
                         total_counter_rho_d0[701]+total_counter_rho_d0[801]+total_counter_rho_d0[901])

counter_rho_d1  <- sum(total_counter_rho_d1[1]+total_counter_rho_d1[101]+total_counter_rho_d1[201]+
                         total_counter_rho_d1[301]+total_counter_rho_d1[401]+
                         total_counter_rho_d1[501]+total_counter_rho_d1[601]+
                         total_counter_rho_d1[701]+total_counter_rho_d1[801]+total_counter_rho_d1[901])

save_hr_l1 <- df$save_hr_l1
save_hr_l2 <- df$save_hr_l2
save_rho_d0 <- df$save_rho_d0
save_rho_d1 <- df$save_rho_d1

total_counter_exp <- df$counter_exp
total_counter_wei <- df$counter_wei
total_counter_gom <- df$counter_gom


counter_exp  <- sum(total_counter_exp[1]+total_counter_exp[101]+total_counter_exp[201]+
                      total_counter_exp[301]+total_counter_exp[401]+
                      total_counter_exp[501]+total_counter_exp[601]+
                      total_counter_exp[701]+total_counter_exp[801]+total_counter_exp[901])

counter_wei  <- sum(total_counter_wei[1]+total_counter_wei[101]+total_counter_wei[201]+
                      total_counter_wei[301]+total_counter_wei[401]+
                      total_counter_wei[501]+total_counter_wei[601]+
                      total_counter_wei[701]+total_counter_wei[801]+total_counter_wei[901])

counter_gom  <- sum(total_counter_gom[1]+total_counter_gom[101]+total_counter_gom[201]+
                      total_counter_gom[301]+total_counter_gom[401]+
                      total_counter_gom[501]+total_counter_gom[601]+
                      total_counter_gom[701]+total_counter_gom[801]+total_counter_gom[901])

#hrs#
#bias
hr_l1_bias <- mean(abs(bias_l1_hr))
hr_l2_bias <- mean(abs(bias_l2_hr))
#coverage
hr_l1_cov <- (counter_hr_l1 / runs) * 100
hr_l2_cov <- (counter_hr_l2 / runs) * 100
#variance
hr_l1_var <- (sd(save_hr_l1))^0.5
hr_l2_var <- (sd(save_hr_l2))^0.5
#mse
hr_l1_mse <- hr_l1_bias^2+hr_l1_var
hr_l2_mse <- hr_l2_bias^2+hr_l2_var

print(hr_l1_bias)
print(hr_l1_cov)
print(hr_l1_mse)

print(hr_l2_bias)
print(hr_l2_cov)
print(hr_l2_mse)


#rho#
rho_d0_bias <- mean(abs(bias_rho_d0))
rho_d1_bias <- mean(abs(bias_rho_d1))
rho_d0_cov <- (counter_rho_d0 / runs) * 100 
rho_d1_cov <- (counter_rho_d1 / runs) * 100
rho_d0_var <- (sd(save_rho_d0))^0.5
rho_d1_var <- (sd(save_rho_d1))^0.5
rho_d0_mse <- rho_d0_bias^2+rho_d0_var
rho_d1_mse <- rho_d1_bias^2+rho_d1_var


print(rho_d0_bias)
print(rho_d0_cov)
print(rho_d0_mse)

print(rho_d1_bias)
print(rho_d1_cov)
print(rho_d1_mse)

#counters#
exp_perc <- counter_exp / runs *100
wei_perc <- counter_wei / runs *100
gom_perc <- counter_gom / runs *100

print(exp_perc)
print(wei_perc)
print(gom_perc)

# YW 22 March 2023: put results together and write to CSV file
# mean of bias
# hr_l1 represents non-terminal event; hr_l2 represents terminal event
bias <- c(hr_l1_bias, hr_l2_bias, rho_d0_bias, rho_d1_bias)

# Coverage probability
CP <- c(hr_l1_cov,hr_l2_cov, rho_d0_cov, rho_d1_cov)

MSE <- c(hr_l1_mse, hr_l2_mse, rho_d0_mse, rho_d1_mse)

# in the order of exponential, weibull, gompertz.
percentage_chosen = c(exp_perc, wei_perc, gom_perc, "na")

# YW: put results together
items<-c("hr_nt", "hr_t", "rho_reference", "rho_covariates")
Results <- cbind.data.frame(items, bias, CP, MSE, percentage_chosen)

Results[,2:4] <- round(Results[,2:4],3)

Results[,3] <- format(round(Results[,3],1), nsmall=1)
Results[,4] <- format(round(Results[,4],1), nsmall=1)
Results
rownames(Results)<-NULL

# output results
write.csv(Results, row.names=F,file=paste0(dir_results,out_file_summary))