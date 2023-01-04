##################################################################################################
# Programmed by YW, 31 December 2022
# Purpose: to combine results for hazard ratios together for paper2
# Output: table in a format that can be included to LaTeX file directly
##################################################################################################
library(dplyr); library(formattable)
dir_results <- "../../results/real_data_analysis/revision_1/"

# Read in results from each copula regression model-------------------------------------------
combined_table_hr <- function(dir_results, table_ref, survival_distribution){
  # normal
  norm <- read.csv(file=paste0(dir_results, "table",table_ref ,"_normal_",survival_distribution, ".csv"))
  
  # put AIC and running time on a separate row - row one
  row.one <- norm[1,]
  row.one$X = paste0("Normal copula ", survival_distribution, " survival model")
  row.one[2:(ncol(row.one)-2)] = ""
  norm <- rbind(row.one, norm)
  norm$aic[2:4] = norm$run_time[2:4] = ""
  norm
  
  # clayton
  cly <- read.csv(file=paste0(dir_results, "table",table_ref,"_clayton_", survival_distribution,".csv"))
  
  # put AIC and running time on a sepeate row - row one
  row.one <- cly[1,]
  row.one$X = paste0("Clayton copula ", survival_distribution, " survival model")
  row.one[2:(ncol(row.one)-2)] = ""
  cly <- rbind(row.one, cly)
  cly$aic[2:4] = cly$run_time[2:4] = ""
  cly
  
  # frank
  frank <- read.csv(file=paste0(dir_results, "table", table_ref, "_frank_", survival_distribution,".csv"))
  frank
  
  # put AIC and running time on a sepeate row - row one
  row.one <- frank[1,]
  row.one$X = paste0("Frank copula ", survival_distribution, " survival model")
  row.one[2:(ncol(row.one)-2)] = ""
  frank <- rbind(row.one, frank)
  frank$aic[2:4] = frank$run_time[2:4] = ""
  frank
  
  # gumbel
  gum <- read.csv(file=paste0(dir_results, "table",table_ref,"_gumbel_",survival_distribution,".csv"))
  gum
  
  # put AIC and running time on a sepeate row - row one
  row.one <- gum[1,]
  row.one$X = paste0("Gumbel copula ", survival_distribution, " survival model")
  row.one[2:(ncol(row.one)-2)] = ""
  gum <- rbind(row.one, gum)
  gum$aic[2:4] = gum$run_time[2:4] = ""
  gum
  
  # combine results into one data frame
  df <- rbind(norm, cly, frank, gum)
  df
  
  # prepare a table for LaTeX file
  df$X[which(df$X=="age.gl50")] = "\\phantom{xxx} Age group: $>50$ &"
  df$X[which(df$X=="gender.female")] = "\\phantom{xxx} Sex: Female &"
  df$X[which(df$X=="donor.living")] = "\\phantom{xxx} Donor type: Living &"
  
  df$X[which(df$X=="results.age")] = "\\phantom{xxx} Age group: $>50$ &"
  df$X[which(df$X=="results.gender")] = "\\phantom{xxx} Sex: Female &"
  df$X[which(df$X=="results.donor")] = "\\phantom{xxx} Donor type: Living &"
  df
  
  df <- df %>% mutate(l_gf = paste0(" (", l_gf,", " )) %>%
    mutate(u_gf = paste0(u_gf, ") &")) %>%
    mutate(l_d = paste0(" (", l_d,", " )) %>%
    mutate(u_d = paste0(u_d, ") &")) %>%
    mutate(l_theta = paste0(" (", l_theta,", " )) %>%
    mutate(u_theta = paste0(u_theta, ") &")) %>%
    mutate(run_time = paste0(run_time, "\\\\")) %>%
    mutate(aic= paste0(aic, "&"))
  
  index <- which(grepl("copula", df[,1], fixed = TRUE)==T)
  
  df[index,2:(ncol(df)-2)]=""
  
  df[index,1] <- paste0("\\hline\\multicolumn{3}{|l}{", df[index,1], "} & &")
  df
  
  # Output table to a CSV file-----------------------------------------------------------------
  write.csv(df, 
            file=paste0(dir_results, "table_hr_", survival_distribution,".csv"), 
            row.names=F)
  return(df)
}

# Hazard ratios from copula exponential models ---------------------------------
survival_distribution = "exp"
table_ref = "3"
combined_table_hr(dir_results, table_ref, survival_distribution)

# Hazard ratios from copula gompertz models ---------------------------------
survival_distribution = "gompertz"
table_ref = "4"
combined_table_hr(dir_results, table_ref, survival_distribution)

# Hazard ratios from copula weibull models ---------------------------------
survival_distribution = "weibull"
table_ref = "5"
combined_table_hr(dir_results, table_ref, survival_distribution)
