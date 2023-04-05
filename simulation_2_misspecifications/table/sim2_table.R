###############################################################################################
# Programmed by YW, 22 March 2023
# Purpose: to combine results for simulation 2 for paper 2
# Output: table in a format that can be included to LaTeX file directly
###############################################################################################

library(dplyr); library(xtable); library(stringr)

# directory if on own PC
dir_results <- "../../results/simulation_results/"

# specify a data frame
df_setup <- function(df){
  df <- data.frame(matrix(nrow=3,ncol=16))
  names(df) <- c("model","hr_nt_bias",  "hr_nt_cp", "hr_nt_mse",
                 "hr_t_bias",   "hr_nt_cp", "hr_nt_mse", 
                 "ref_rho_bias", "ref_rho_cp", "ref_rho_mse",
                 "cov_rho_bias", "cov_rho_cp", "cov_rho_mse",
                 "percent_exp", "percent_weibull", "percent_gompertz")
  df$model <- c("Exponential", "Weibull", "Gompertz")
  df
}

# specify the data in the data frame
df_input <- function(df, row_ref, results){
  df[row_ref,2:4] <- results[1,2:4] 
  df[row_ref,5:7] <- results[2,2:4]
  df[row_ref,8:10] <- results[3,2:4]
  df[row_ref,11:13] <- results[4,2:4]
  df[row_ref,14:16] <- results$percentage_chosen[1:3]
  df
}

# convert the data frame into a LaTEX table
function_latex_table <- function(df){
  for(j in 1:ncol(df)){
    if(j!=ncol(df)){
      df[,j] = paste0(df[,j],"&")      # column separator
    }else{
      df[,j] = paste0(df[,j],"\\\\")    # end of the column for a new line
    }
    df[,j] = gsub("-", "$-$", df[,j])  # maths symbol 
  }
  df
}

###############################################################################################
# Simulation 2: miss-specification of survival distributions
###############################################################################################
# Clayton copula 
clayton_exp <- read.csv(file=paste0(dir_results,"simulation2/clayton/s2_aic_clayton_exp_summary.csv"))
clayton_weibull <- read.csv(file=paste0(dir_results,"simulation2/clayton/s2_aic_clayton_weibull_summary.csv"))
clayton_gompertz <- read.csv(file=paste0(dir_results,"simulation2/clayton/s2_aic_clayton_gompertz_summary.csv"))

df <- NULL
df <- df_setup(df)

df <- df_input(df, row_ref=1, results=clayton_exp)
df <- df_input(df, row_ref=2, results=clayton_weibull)
df <- df_input(df, row_ref=3, results=clayton_gompertz)

df <- function_latex_table(df)

write.csv(df, file=paste0(dir_results, "sim2_aic_table_clayton_copula.csv"))

# Frank copula 
frank_exp <- read.csv(file=paste0(dir_results,"simulation2/frank/s2_aic_frank_exp_summary.csv"))
frank_weibull <- read.csv(file=paste0(dir_results,"simulation2/frank/s2_aic_frank_weibull_summary.csv"))
frank_gompertz <- read.csv(file=paste0(dir_results,"simulation2/frank/s2_aic_frank_gompertz_summary.csv"))

df <- NULL
df <- df_setup(df)
# # exponential
# df[1,2:4] <- frank_exp[1,2:4] 
# df[1,5:7] <- frank_exp[2,2:4]
# df[1,8:10] <- frank_exp[3,2:4]
# df[1,11:13] <- frank_exp[4,2:4]
# df[1,14:16] <- frank_exp$percentage_chosen[1:3]

df <- df_input(df, row_ref=1, results=frank_exp)
df <- df_input(df, row_ref=2, results=frank_weibull)
df <- df_input(df, row_ref=3, results=frank_gompertz)

df <- function_latex_table(df)

write.csv(df, file=paste0(dir_results, "sim2_aic_table_frank_copula.csv"))

################################################################################
# Gumbel copula 
################################################################################
gumbel_exp <- read.csv(file=paste0(dir_results,"simulation2/gumbel/s2_aic_gumbel_exp_summary.csv"))
gumbel_weibull <- read.csv(file=paste0(dir_results,"simulation2/gumbel/s2_aic_gumbel_weibull_summary.csv"))
gumbel_gompertz <- read.csv(file=paste0(dir_results,"simulation2/gumbel/s2_aic_gumbel_gompertz_summary.csv"))

df <- NULL
df <- df_setup(df)

df <- df_input(df, row_ref=1, results=gumbel_exp)
df <- df_input(df, row_ref=2, results=gumbel_weibull)
df <- df_input(df, row_ref=3, results=gumbel_gompertz)

df <- function_latex_table(df)

write.csv(df, file=paste0(dir_results, "sim2_aic_table_gumbel_copula.csv"))


################################################################################
# Normal copula 
################################################################################
# normal exponential
normal_exp <- read.csv(file=paste0(dir_results,"simulation2/normal_exp/S2_aic_normal_exp_summary.csv"))

# normal weibull
normal_weibull <- read.csv(file=paste0(dir_results,"simulation2/normal_weib/S2_aic_normal_weib_summary.csv"))

# normal gompertz
normal_gompertz <- read.csv(file=paste0(dir_results,"simulation2/normal_gomp/S2_aic_normal_gomp_summary.csv"))

df <- NULL
df <- df_setup(df)
df <- df_input(df, row_ref=1, results=normal_exp)
df <- df_input(df, row_ref=2, results=normal_weibull)
df <- df_input(df, row_ref=3, results=normal_gompertz)
df[,c(3,6,9,12)] <- format(round(df[,c(3,6,9,12)],1), nsmall=1)

df <- function_latex_table(df)
write.csv(df, file=paste0(dir_results, "sim2_aic_table_normal_copula.csv"))
