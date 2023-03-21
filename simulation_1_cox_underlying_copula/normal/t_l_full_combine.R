# directory if on own PC
dir_results <- "../../results/simulation_results/"

model2_normal1 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_1_250.csv" ))
model2_normal2 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_251_500.csv" ))
model2_normal3 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_501_750.csv" ))
model2_normal4 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_751_1000.csv" ))

model2_normal <- rbind(model2_normal1, model2_normal2, model2_normal3, model2_normal4)
model2_normal <- model2_normal %>% select(!contains("X", ignore.case = TRUE))
df <- data.frame(model2_normal)

true_b0 <- 0.35; true_b1 <- 0.28; true_b2 <- 0; true_b3 <- 0
true_a0 <- -3.30; true_a1 <- 0.11; true_a2 <- 0.02; true_a3 <- -0.51
true_c0 <- -4.15; true_c1 <- 1.32; true_c2 <- -0.11; true_c3  <- -0.65


save_a0 <- df$save_a0
save_a1 <- df$save_a1
save_a2 <- df$save_a2
save_a3 <- df$save_a3
save_c0 <- df$save_c0
save_c1 <- df$save_c1
save_c2 <- df$save_c2
save_c3 <- df$save_c3
save_b0 <- df$save_b0
save_b1 <- df$save_b1
save_b2 <- df$save_b2
save_b3 <- df$save_b3

save_se_a0 <- df$save_se_a0
save_se_a1 <- df$save_se_a1
save_se_a2 <- df$save_se_a2
save_se_a3 <- df$save_se_a3
save_se_c0 <- df$save_se_c0
save_se_c1 <- df$save_se_c1
save_se_c2 <- df$save_se_c2
save_se_c3 <- df$save_se_c3
save_se_b0 <- df$save_se_b0
save_se_b1 <- df$save_se_b1
save_se_b2 <- df$save_se_b2
save_se_b3 <- df$save_se_b3

save_var_a0 <- df$save_var_a0
save_var_a1 <- df$save_var_a1
save_var_a2 <- df$save_var_a2
save_var_a3 <- df$save_var_a3
save_var_c0 <- df$save_var_c0
save_var_c1 <- df$save_var_c1
save_var_c2 <- df$save_var_c2
save_var_c3 <- df$save_var_c3
save_var_b0 <- df$save_var_b0
save_var_b1 <- df$save_var_b1
save_var_b2 <- df$save_var_b2
save_var_b3 <- df$save_var_b3

#attach(df)

runs = 1000

true_hr_l1_age <- exp(true_a1)
true_hr_l1_gen <- exp(true_a2)
true_hr_l1_donor <- exp(true_a3)
true_hr_l2_age <- exp(true_c1)
true_hr_l2_gen <- exp(true_c2)
true_hr_l2_donor <- exp(true_c3)

counter_a0 = 0
counter_a1 = 0
counter_a2 = 0
counter_a3 = 0
counter_c0 = 0
counter_c1 = 0
counter_c2 = 0
counter_c3 = 0
counter_b0 = 0
counter_b1 = 0
counter_b2 = 0
counter_b3 = 0
counter_hr_l1_age = 0
counter_hr_l1_gen = 0
counter_hr_l1_donor = 0
counter_hr_l2_age = 0
counter_hr_l2_gen = 0
counter_hr_l2_donor = 0

############ ci ####################
uci_a0 <- save_a0 + 1.96*save_se_a0
lci_a0 <- save_a0 - 1.96*save_se_a0
uci_a1 <- save_a1 + 1.96*save_se_a1
lci_a1 <- save_a1 - 1.96*save_se_a1
uci_a2 <- save_a2 + 1.96*save_se_a2
lci_a2 <- save_a2 - 1.96*save_se_a2
uci_a3 <- save_a3 + 1.96*save_se_a3
lci_a3 <- save_a3 - 1.96*save_se_a3

uci_c0 <- save_c0 + 1.96*save_se_c0
lci_c0 <- save_c0 - 1.96*save_se_c0
uci_c1 <- save_c1 + 1.96*save_se_c1
lci_c1 <- save_c1 - 1.96*save_se_c1
uci_c2 <- save_c2 + 1.96*save_se_c2
lci_c2 <- save_c2 - 1.96*save_se_c2
uci_c3 <- save_c3 + 1.96*save_se_c3
lci_c3 <- save_c3 - 1.96*save_se_c3

uci_b0 <- save_b0 + 1.96*save_se_b0
lci_b0 <- save_b0 - 1.96*save_se_b0
uci_b1 <- save_b1 + 1.96*save_se_b1
lci_b1 <- save_b1 - 1.96*save_se_b1
uci_b2 <- save_b2 + 1.96*save_se_b2
lci_b2 <- save_b2 - 1.96*save_se_b2
uci_b3 <- save_b3 + 1.96*save_se_b3
lci_b3 <- save_b3 - 1.96*save_se_b3

################## biases ########################
bias_a0 <- mean(rep(true_a0,runs) - save_a0)
bias_a1 <- mean(rep(true_a1,runs) - save_a1)
bias_a2 <- mean(rep(true_a2,runs) - save_a2)
bias_a3 <- mean(rep(true_a3,runs) - save_a3)
bias_c0 <- mean(rep(true_c0,runs) - save_c0)
bias_c1 <- mean(rep(true_c1,runs) - save_c1)
bias_c2 <- mean(rep(true_c2,runs) - save_c2)
bias_c3 <- mean(rep(true_c3,runs) - save_c3)
bias_b0 <- mean(rep(true_b0,runs) - save_b0)
bias_b1 <- mean(rep(true_b1,runs) - save_b1)
bias_b2 <- mean(rep(true_b2,runs) - save_b2)
bias_b3 <- mean(rep(true_b3,runs) - save_b3)

############### variances ###################
var_a0 <- var(save_a0)
var_a1 <- var(save_a1)
var_a2 <- var(save_a2)
var_a3 <- var(save_a3)
var_c0 <- var(save_c0)
var_c1 <- var(save_c1)
var_c2 <- var(save_c2)
var_c3 <- var(save_c3)
var_b0 <- var(save_b0)
var_b1 <- var(save_b1)
var_b2 <- var(save_b2)
var_b3 <- var(save_b3)

############## mse ##################
mse_a0 <- bias_a0^2+var_a0
mse_a1 <- bias_a1^2+var_a1
mse_a2 <- bias_a2^2+var_a2
mse_a3 <- bias_a3^2+var_a3
mse_c0 <- bias_c0^2+var_c0
mse_c1 <- bias_c1^2+var_c1
mse_c2 <- bias_c2^2+var_c2
mse_c3 <- bias_c3^2+var_c3
mse_b0 <- bias_b0^2+var_b0
mse_b1 <- bias_b1^2+var_b1
mse_b2 <- bias_b2^2+var_b2
mse_b3 <- bias_b3^2+var_b3

######### covarage ###############
for (i in 1:length(save_a0)){
  if(true_a0 <= uci_a0[i]   && true_a0 >= lci_a0[i])   {counter_a0 = counter_a0+1}
  if(true_a1 <= uci_a1[i]   && true_a1 >= lci_a1[i])   {counter_a1 = counter_a1+1}
  if(true_a2 <= uci_a2[i]   && true_a2 >= lci_a2[i])   {counter_a2 = counter_a2+1}
  if(true_a3 <= uci_a3[i]   && true_a3 >= lci_a3[i])   {counter_a3 = counter_a3+1}
  
  if(true_c0 <= uci_c0[i]   && true_c0 >= lci_c0[i])   {counter_c0 = counter_c0+1}
  if(true_c1 <= uci_c1[i]   && true_c1 >= lci_c1[i])   {counter_c1 = counter_c1+1}
  if(true_c2 <= uci_c2[i]   && true_c2 >= lci_c2[i])   {counter_c2 = counter_c2+1}
  if(true_c3 <= uci_c3[i]   && true_c3 >= lci_c3[i])   {counter_c3 = counter_c3+1}
  
  if(true_b0 <= uci_b0[i]   && true_b0 >= lci_b0[i])   {counter_b0 = counter_b0+1}
  if(true_b1 <= uci_b1[i]   && true_b1 >= lci_b1[i])   {counter_b1 = counter_b1+1}
  if(true_b2 <= uci_b2[i]   && true_b2 >= lci_b2[i])   {counter_b2 = counter_b2+1}
  if(true_b3 <= uci_b3[i]   && true_b3 >= lci_b3[i])   {counter_b3 = counter_b3+1}
}

cov_a0 <- (counter_a0 / runs) * 100
cov_a1 <- (counter_a1 / runs) * 100
cov_a2 <- (counter_a2 / runs) * 100
cov_a3 <- (counter_a3 / runs) * 100
cov_c0 <- (counter_c0 / runs) * 100
cov_c1 <- (counter_c1 / runs) * 100
cov_c2 <- (counter_c2 / runs) * 100
cov_c3 <- (counter_c3 / runs) * 100
cov_b0 <- (counter_b0 / runs) * 100
cov_b1 <- (counter_b1 / runs) * 100
cov_b2 <- (counter_b2 / runs) * 100
cov_b3 <- (counter_b3 / runs) * 100


############## HR ##################
esthr_l1_age <- exp(save_a1)
esthr_l1_gen <- exp(save_a2)
esthr_l1_donor <- exp(save_a3)

esthr_l2_age <- exp(save_c1)
esthr_l2_gen <- exp(save_c2)
esthr_l2_donor <- exp(save_c3)

var_hr_l1_age <- exp(save_a1)^2 * save_var_a1
var_hr_l1_gen <- exp(save_a2)^2 * save_var_a2
var_hr_l1_donor <- exp(save_a3)^2 * save_var_a3
var_hr_l2_age <- exp(save_c1)^2 * save_var_c1
var_hr_l2_gen <- exp(save_c2)^2 * save_var_c2
var_hr_l2_donor <- exp(save_c3)^2 * save_var_c3

hr_l1_lwci_age <- esthr_l1_age - 1.96*sqrt(var_hr_l1_age)
hr_l1_upci_age <- esthr_l1_age + 1.96*sqrt(var_hr_l1_age)
hr_l1_lwci_gen <- esthr_l1_gen - 1.96*sqrt(var_hr_l1_gen)
hr_l1_upci_gen <- esthr_l1_gen + 1.96*sqrt(var_hr_l1_gen)
hr_l1_lwci_donor <- esthr_l1_donor - 1.96*sqrt(var_hr_l1_donor)
hr_l1_upci_donor <- esthr_l1_donor + 1.96*sqrt(var_hr_l1_donor)
hr_l2_lwci_age <- esthr_l2_age - 1.96*sqrt(var_hr_l2_age)
hr_l2_upci_age <- esthr_l2_age + 1.96*sqrt(var_hr_l2_age)
hr_l2_lwci_gen <- esthr_l2_gen - 1.96*sqrt(var_hr_l2_gen)
hr_l2_upci_gen <- esthr_l2_gen + 1.96*sqrt(var_hr_l2_gen)
hr_l2_lwci_donor <- esthr_l2_donor - 1.96*sqrt(var_hr_l2_donor)
hr_l2_upci_donor <- esthr_l2_donor + 1.96*sqrt(var_hr_l2_donor)

############ bias ##############
bias_hr_l1_age <- mean(rep(true_hr_l1_age,runs) - esthr_l1_age)
bias_hr_l1_gen <- mean(rep(true_hr_l1_gen,runs) - esthr_l1_gen)
bias_hr_l1_donor <- mean(rep(true_hr_l1_donor,runs) - esthr_l1_donor)
bias_hr_l2_age <- mean(rep(true_hr_l2_age,runs) - esthr_l2_age)
bias_hr_l2_gen <- mean(rep(true_hr_l2_gen,runs) - esthr_l2_gen)
bias_hr_l2_donor <- mean(rep(true_hr_l2_donor,runs) - esthr_l2_donor)

############ var ############
hr_l1_var_age <- var(esthr_l1_age)
hr_l2_var_age <- var(esthr_l2_age)
hr_l1_var_gen <- var(esthr_l1_gen)
hr_l2_var_gen <- var(esthr_l2_gen)
hr_l1_var_donor <- var(esthr_l1_donor)
hr_l2_var_donor <- var(esthr_l2_donor)

######### mse ############
mse_hr_l1_age <- bias_hr_l1_age^2 + hr_l1_var_age
mse_hr_l1_gen <- bias_hr_l1_gen^2 + hr_l1_var_gen
mse_hr_l1_donor <- bias_hr_l1_donor^2 + hr_l1_var_donor
mse_hr_l2_age <- bias_hr_l2_age^2 + hr_l2_var_age
mse_hr_l2_gen <- bias_hr_l2_gen^2 + hr_l2_var_gen
mse_hr_l2_donor <- bias_hr_l2_donor^2 + hr_l2_var_donor

############## cov ################
for (i in 1:length(save_a0)){
  if(true_hr_l1_age <= hr_l1_upci_age[i]   && true_hr_l1_age >= hr_l1_lwci_age[i])   {counter_hr_l1_age = counter_hr_l1_age+1}
  if(true_hr_l1_gen <= hr_l1_upci_gen[i]   && true_hr_l1_gen >= hr_l1_lwci_gen[i])   {counter_hr_l1_gen = counter_hr_l1_gen+1}
  if(true_hr_l1_donor <= hr_l1_upci_donor[i]   && true_hr_l1_donor >= hr_l1_lwci_donor[i])   {counter_hr_l1_donor = counter_hr_l1_donor+1}
  if(true_hr_l2_age <= hr_l2_upci_age[i]   && true_hr_l2_age >= hr_l2_lwci_age[i])   {counter_hr_l2_age = counter_hr_l2_age+1}
  if(true_hr_l2_gen <= hr_l2_upci_gen[i]   && true_hr_l2_gen >= hr_l2_lwci_gen[i])   {counter_hr_l2_gen = counter_hr_l2_gen+1}
  if(true_hr_l2_donor <= hr_l2_upci_donor[i]   && true_hr_l2_donor >= hr_l2_lwci_donor[i])   {counter_hr_l2_donor = counter_hr_l2_donor+1}
}

cov_hr_l1_age <- (counter_hr_l1_age / runs) * 100
cov_hr_l1_gen <- (counter_hr_l1_gen / runs) * 100
cov_hr_l1_donor <- (counter_hr_l1_donor / runs) * 100
cov_hr_l2_age <- (counter_hr_l2_age / runs) * 100
cov_hr_l2_gen <- (counter_hr_l2_gen / runs) * 100
cov_hr_l2_donor <- (counter_hr_l2_donor / runs) * 100

