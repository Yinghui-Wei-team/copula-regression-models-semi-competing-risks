# https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html
#install.packages("devtools")
#library(devtools)
#devtools::install_github("adayim/forestploter")
library(forestploter)

library(forestploter); library(dplyr);library(ggplot2)
output_dir <- "results/simulation_results/"
# df <- read.csv("results/real_data_analysis/table3_combined_hazard_ratios.csv")
# df <- reshape(df, idvar = "covariate", timevar = "model", direction = "wide")
# write.csv(df, file=paste0(output_dir, "df_wide.csv"), row.names=F)

df <- read.csv("results/real_data_analysis/table3_df_wide2.csv")
df <- df %>% rename(Covariate = covariate)
#df$` ` <- paste(rep(" ", 40), collapse = " ")
df$'                                               Hazard Ratios for Graft Failure' <- paste(rep("                  ", nrow(df)), collapse = " ")
df$'                                               Hazard Ratios for Death' <- paste(rep("                  ", nrow(df)), collapse = " ")

#################################################################################
# Set-up theme                                                                 #
#################################################################################
tm <- forest_theme(base_size = 10,
                   refline_lty = "solid",
                   ci_pch = c(15, 18, 16, 17, 19),
                   ci_col = c("#808080", "#00FF00", "royalblue3", "maroon3", "red"),
                   ci_lwd = 2,
                   footnote_col = "blue",
                   legend_name = "Model:   ", legend_position = "bottom",
                   legend_value = c("Cox  ", "Normal  ", "Clayton  ",  "Frank", "Gumbel"),
                   vertline_lty = c("dashed", "dotdash"),
                   vertline_col = c("#d6604d", "#A52A2A"))

##########################################################################################
outcome <- "gf"

p <- forest(df[,c(1, 32)],
            est = list(df$hr_gf.cox,
                       df$hr_gf.normal,
                       df$hr_gf.clayton,
                       df$hr_gf.frank,
                       df$hr_gf.gumbel
            ),
            lower = list(df$low_gf.cox,
                         df$low_gf.normal,
                         df$low_gf.clayton,
                         df$low_gf.frank,
                         df$low_gf.gumbel),
            upper = list(df$high_gf.cox,
                         df$high_gf.normal,
                         df$high_gf.clayton,
                         df$high_gf.frank,
                         df$high_gf.gumbel),
            ci_column = 2,
            ref_line = 1,
            #vert_line = c(0.5, 2),
            x_trans = "log10",
            nudge_y = 0.2,
            xlim = c(0.4,1.5),
            ticks_at = c(0.5, 1,  1.3),
            theme = tm)

p

ggsave(file=paste0("plot_v2_HR_", outcome,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=9, height=6)


##########################################################################################
outcome <- "death"

p <- forest(df[,c(1, 33)],
            est = list(df$hr_d.cox,
                       df$hr_d.normal,
                       df$hr_d.clayton,
                       df$hr_d.frank
            ),
            lower = list(df$low_d.cox,
                         df$low_d.normal,
                         df$low_d.clayton,
                         df$low_d.frank),
            upper = list(df$high_d.cox,
                         df$high_d.normal,
                         df$high_d.clayton,
                         df$high_d.frank),
            ci_column = 2,
            ref_line = 1,
            vert_line = c(0.5, 2),
            x_trans = "log10",
            nudge_y = 0.2,
            xlim = c(0.4,5),
            ticks_at = c(0.5, 1, 2,  5),
            theme = tm)

p

ggsave(file=paste0("plot_HR_", outcome,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=9, height=6)