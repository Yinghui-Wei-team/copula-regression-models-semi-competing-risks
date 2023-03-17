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
