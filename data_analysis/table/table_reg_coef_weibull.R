################################################################################
# Programmed by YW, 30 December 2022
# Purpose: to combine results together to produce table S2 for paper2
###############################################################################

dir_data <-dir_results <- "../../results/real_data_analysis/revision_1/"

survival_distribution = "weibull"
# normal weibull
wei_norm <- read.csv(file=paste0(dir_results, "parameters_normal_weibull.csv"))
names(wei_norm) <- c("para", "norm", "norm.lower", "norm.upper")

# clayton weibull
wei_cly <- read.csv(file=paste0(dir_results, "parameters_clayton_weibull.csv"))
names(wei_cly) <- c("para", "clayton", "clayton.lower", "clayton.upper")

# frank weibull
wei_frank <- read.csv(file=paste0(dir_results, "parameters_frank_weibull.csv"))
names(wei_frank) <- c("para", "frank", "frank.lower", "frank.upper")

# gumbel weibull
wei_gum <- read.csv(file=paste0(dir_results, "parameters_gumbel_weibull.csv"))
names(wei_gum) <- c("para", "gum", "gum.lower", "gum.upper")

df <- merge(wei_norm, wei_cly, by.x = "para")

df <- merge(df, wei_frank, by.x = "para")

df <- merge(df, wei_gum, by.x = "para")

tbl <- data.frame(para = as.character(),
                  norm = as.character(),
                  clayton = as.character(),
                  frank = as.character())

tbl <- as.data.frame(matrix(nrow=nrow(df),ncol=5))
names(tbl) <- c("Parameter", "Normal-Weibull", "Clayton-Weibull", "Frank-Weibull", "Gumbel-Weibull")
tbl$Parameter <- df$para
tbl$`Normal-Weibull` <- paste0("&",format(round(df$norm,2),nsmall=2), 
                               " (", format(round(df$norm.lower,2), nsmall=2), ", ", 
                               format(round(df$norm.upper,2), nsmall=2),")")
tbl$`Clayton-Weibull` <- paste0("&", format(round(df$clayton,2),nsmall=2), " (", 
                                format(round(df$clayton.lower,2), nsmall=2), ", ", 
                                format(round(df$clayton.upper,2), nsmall=2),")")
tbl$`Frank-Weibull` <- paste0("&", format(round(df$frank,2),nsmall=2), " (", 
                              format(round(df$frank.lower,2),nsmall=2), ", ", 
                              format(round(df$frank.upper,2),nsmall=2),")")
tbl$`Gumbel-Weibull` <- paste0("&",format(round(df$gum,2),nsmall=2), " (", 
                               format(round(df$gum.lower,2),nsmall=2), ", ", 
                               format(round(df$gum.upper,2),nsmall=2),")", "\\")

tbl$Parameter <- c("$a_0$", "$a_1$", "$a_2$", "$a_3$",     
                   "$alpha_1$", "$alpha_2$",
                   "$b_0$", "$b_1$",  "$b_2$", "$b_3$",
                   "$c_0$", "$c_1$",  "$c_2$", "$c_3$" )
tbl

write.csv(tbl, 
          file=paste0(dir_results, "table_reg_coef_", survival_distribution,".csv"), 
          row.names=F)

print(paste0("Results for regression coefficients for ", survival_distribution, " models have been saved!"))
# tbl2 <- gsub("-","$-$",tbl)
# 
# tbl2
