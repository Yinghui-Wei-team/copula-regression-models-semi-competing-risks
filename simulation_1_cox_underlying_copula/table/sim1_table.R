###############################################################################################
# Programmed by YW, 3 January 2023
# Purpose: to combine results for simulation 1 for paper 2
# Output: table in a format that can be included to LaTeX file directly
###############################################################################################

library(dplyr); library(xtable)
# directory if on own PC
dir_results <- "../../results/simulation_results/"

function_latex_table <- function(df){
  for(j in 2:ncol(df)){
    if(j!=ncol(df)){
      df[,j] = paste0(df[,j],"&")      # column separator
    }else{
      df[,j] = paste0(df[,j],"&\\")    # end of the column for a new line
    }
    df[,j] = gsub("-", "$-$", df[,j])  # maths symbol 
  }
  df
}
###############################################################################################
# Simulation 1 Cox Model
###############################################################################################
cox_normal <- read.csv(file=paste0(dir_results,"simulation1/s1_cox_model_summary_normal_copula.csv"))
names(cox_normal) <- c("items", "normal_bias", "normal_cp", "normal_mse")

cox_clayton <- read.csv(file=paste0(dir_results,"simulation1/s1_cox_model_summary_clayton_copula.csv"))
names(cox_clayton) <- c("items", "clayton_bias", "clayton_cp", "clayton_mse")
cox_clayton <- cox_clayton %>% select(!contains("item", ignore.case = TRUE))

cox_frank <- read.csv(file=paste0(dir_results,"simulation1/s1_cox_model_summary_frank_copula.csv"))
names(cox_frank) <- c("items", "frank_bias", "frank_cp", "frank_mse")
cox_frank <- cox_frank %>% select(!contains("item", ignore.case = TRUE))

cox_gumbel <- read.csv(file=paste0(dir_results,"simulation1/s1_cox_model_summary_gumbel_copula.csv"))
names(cox_gumbel) <- c("items", "gumbel_bias", "gumbel_cp", "gumbel_mse")
cox_gumbel <- cox_gumbel %>% select(!contains("item", ignore.case = TRUE))

df <- cbind(cox_normal, cox_clayton, cox_frank, cox_gumbel)

df <- function_latex_table(df)

df$items <- gsub("_", ",\\\\mbox{",df$items)
df$items <- paste0("$\\mbox{HR}_", df$items, "}$ &")
df$items
df

row.number <-c(1,3,2,4,6,5)
df <- cbind(row.number,df)
df <- df[order(row.number),]
df <- df %>% dplyr::select(!row.number)
df

write.csv(df, file=paste0(dir_results, "sim1_table_cox.csv"))

###############################################################################################
# Simulation 1 Model 1
###############################################################################################
model1_normal <- read.csv(file=paste0(dir_results,"simulation1/s1_model1_summary_normal_exponential.csv" ))
names(model1_normal) <- c("items", "normal_bias", "normal_cp", "normal_mse", "normal_run_time")
model1_normal <- model1_normal %>% select(!contains("time", ignore.case = TRUE))

model1_clayton <- read.csv(file=paste0(dir_results,"simulation1/s1_model1_summary_clayton_exponential.csv" ))
index <- which(model1_clayton == "b0")
model1_clayton <- model1_clayton[-index,]
names(model1_clayton) <- c("items", "clayton_bias", "clayton_cp", "clayton_mse")
model1_clayton <- model1_clayton %>% select(!contains(c("item","time"), ignore.case = TRUE))

model1_frank <- read.csv(file=paste0(dir_results,"simulation1/s1_model1_summary_frank_exponential.csv" ))
names(model1_frank) <- c("items", "frank_bias", "frank_cp", "frank_mse", "frank_run_time")
model1_frank <- model1_frank %>% select(!contains(c("item","time"), ignore.case = TRUE))

model1_gumbel <- read.csv(file=paste0(dir_results,"simulation1/s1_model1_summary_gumbel_exponential.csv" ))
names(model1_gumbel) <- c("items", "gumbel_bias", "gumbel_cp", "gumbel_mse", "gumbel_run_time")
model1_gumbel <- model1_gumbel %>% select(!contains(c("item","time"), ignore.case = TRUE))

df <- cbind(model1_normal, model1_clayton, model1_frank, model1_gumbel)

df <- df %>% filter(!df$items %in% c("a0", "a1", "a2", "a3", "c0", "c1", "c2", "c3"))

row.number <-c(1,3,2,4,6,5)

df <- cbind(row.number,df)

df <- df[order(row.number),]

df <- df %>% select(!row.number)

df[,2:ncol(df)] <- round(df[,2:ncol(df)],3)
df <- df %>% mutate(clayton_cp = format(round(100*(clayton_cp/100),1), nsmall=1)) %>%
  mutate(normal_cp = format(round(100*(normal_cp/100),1), nsmall=1)) %>%
  mutate(frank_cp = format(round(100*(frank_cp/100),1),nsmall=1)) %>%
  mutate(gumbel_cp = format(round(100*(gumbel_cp/100),1), nsmall=1))

for(j in 2:ncol(df)){
  if(j!=ncol(df)){
    df[,j] = paste0(df[,j],"&")      # column separator
  }else{
    df[,j] = paste0(df[,j],"\\")    # end of the column for a new line
  }
  df[,j] = gsub("-", "$-$", df[,j])  # maths symbol 
}

df$items <- gsub("_", ",\\\\mbox{",df$items)
df$items <- paste0("$\\mbox{HR}_", df$items, "}$ &")
df$items
df

# df <- xtable(df, type = "latex", file = "filename2.tex")

write.csv(df, file=paste0(dir_results, "sim1_table_copula_model1.csv"))


###############################################################################################
# Simulation 1 Model 2
###############################################################################################

# Model 2 normal copula
model2_normal1 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_1_250.csv" ))
model2_normal2 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_251_500.csv" ))
model2_normal3 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_501_750.csv" ))
model2_normal4 <- read.csv(file=paste0(dir_results,"simulation1/s1_model2_t_l_full_normal_exp_751_1000.csv" ))

model2_normal <- rbind(model2_normal1, model2_normal2, model2_normal3, model2_normal4)
model2_normal <- model2_normal %>% select(!contains("X", ignore.case = TRUE))

true_b0 <- 0.35; true_b1 <- 0.28; true_b2 <- 0; true_b3 <- 0
true_a0 <- -3.30; true_a1 <- 0.11; true_a2 <- 0.02; true_a3 <- -0.51
true_c0 <- -4.15; true_c1 <- 1.32; true_c2 <- -0.11; true_c3  <- -0.65


# model2_normal <- read.csv(file=paste0(dir_results,"simulation1/s1_model1_summary_normal_exponential.csv" ))
# names(model1_normal) <- c("items", "normal_bias", "normal_cp", "normal_mse", "normal_run_time")
# model1_normal <- model1_normal %>% select(!contains("time", ignore.case = TRUE))

# Model 2 Clayton copula exponential survival
model2_clayton <- read.csv(file=paste0(dir_results,
                                       "simulation1/s1_model2_summary_clayton_exponential.csv" ))
index <- which(model2_clayton == "a0"|model2_clayton == "a1"|model2_clayton == "a2"|model2_clayton == "a3"|
               model2_clayton == "c0"|model2_clayton == "c1"|model2_clayton == "c2"|model2_clayton == "c3")
model2_clayton <- model2_clayton[-index,]
names(model2_clayton) <- c("items", "clayton_bias", "clayton_cp", "clayton_mse")

# Model 2 frank copula exponential survival
model2_frank <- read.csv(file=paste0(dir_results,
                                       "simulation1/s1_model2_summary_frank_exponential.csv" ))
index <- which(model2_frank == "a0"|model2_frank == "a1"|model2_frank == "a2"|model2_frank == "a3"|
              model2_frank == "c0"|model2_frank == "c1"|model2_frank == "c2"|model2_frank == "c3")
model2_frank <- model2_frank[-index,]
names(model2_frank) <- c("items", "frank_bias", "frank_cp", "frank_mse")
model2_frank <- model2_frank %>% select(!contains(c("item"), ignore.case = TRUE))

# Model 2 gumbel copula exponential survival
model2_gumbel <- read.csv(file=paste0(dir_results,
                                      "simulation1/s1_model2_summary_gumbel_exponential.csv" ))
index <- which(model2_gumbel == "a0"|model2_gumbel == "a1"|model2_gumbel == "a2"|model2_gumbel == "a3"|
                 model2_gumbel == "c0"|model2_gumbel == "c1"|model2_gumbel == "c2"|model2_gumbel == "c3")
model2_gumbel <- model2_gumbel[-index,]
names(model2_gumbel) <- c("items", "gumbel_bias", "gumbel_cp", "gumbel_mse")
model2_gumbel <- model2_gumbel %>% select(!contains(c("item"), ignore.case = TRUE))

#df <- cbind(model1_normal, model1_clayton, model1_frank, model1_gumbel)
df <- cbind(model2_clayton, model2_frank, model2_gumbel)

row.number <-c(7, 8, 9, 10,
               1, 3, 2, 4, 
               6, 5)

df <- cbind(row.number,df)
df <- df[order(row.number),]
df <- df %>% select(!row.number)
