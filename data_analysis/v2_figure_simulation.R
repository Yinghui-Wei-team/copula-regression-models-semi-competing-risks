# Purpose: produce a 2x2 plot to show bias from simulation study
# Programmed by Yinghui Wei

library(ggplot2); library(dplyr); library(tidyr)
output_dir <- "results/simulation_results/"
df <- read.csv("results/simulation_results/table5_simulation.csv")

df <- df %>% mutate(model = ifelse(model=="copula1", "Normal copula with covariates on hazards",
                                   ifelse(model=="copula2", "Normal copula with covariates on hazards and association",
                                          "Cox model")))
df <- df %>% rename(Model = model)

df_bias_wide <- df %>% select(Estimand, Model, contains("bias"))

df_bias_long <- gather(df_bias_wide, copula, bias, bias_normal:bias_gum, factor_key=TRUE)

df_bias_long <- df_bias_long %>% mutate(copula = gsub("bias_", "", copula)) %>%
  mutate(copula = ifelse(copula == "clay", "Clayton", 
                         ifelse(copula == "normal", "Normal",
                                ifelse(copula == "frank", "Frank", "Gumbel"))))
# Part 1 - Bias
measure = "bias"
p <- ggplot(df_bias_long, aes(x=Estimand, y=bias, colour=Model, shape=Model, size=Model)) + geom_point()
p <- p +  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c("royalblue3", "maroon3", "#808080"))+
  ylab("Bias \n") + xlab("\nOutcomes and covariates") +
  ylim(-1,1) +
  scale_size_manual(values=c(4,4,4)) +  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  geom_hline(yintercept=0, linetype="dashed", color = "#A52A2A")+
  theme(legend.position="bottom", legend.title=element_text(size=13), 
        legend.text=element_text(size=13),
        axis.text = element_text(size = 15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_line(colour = "black")) + 
  facet_wrap(~ copula) + 
  geom_line(aes(group = Model), lwd=1)

p 

ggsave(file=paste0("plot_sim_", measure,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=12, height=6)


# Part 2 - MSE
df_mse_wide <- df %>% select(Estimand, Model, contains("mse"))

df_mse_long <- gather(df_mse_wide, copula, mse, mse_normal:mse_gum, factor_key=TRUE)

df_mse_long <- df_mse_long %>% mutate(copula = gsub("mse_", "", copula)) %>%
  mutate(copula = ifelse(copula == "clay", "Clayton", 
                         ifelse(copula == "normal", "Normal",
                                ifelse(copula == "frank", "Frank", "Gumbel"))))
measure = "mse"
p <- ggplot(df_mse_long, aes(x=Estimand, y=mse, colour=Model, shape=Model, size=Model)) + geom_point()
p <- p +  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c("royalblue3", "maroon3", "#808080"))+
  ylab("Mean square error \n") + xlab("\nOutcome and covariate") +
  ylim(-0.5,0.5) +
  scale_size_manual(values=c(4,4,4)) +  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  geom_hline(yintercept=0, linetype="dashed", color = "#A52A2A")+
  theme(legend.position="bottom", legend.title=element_text(size=13), 
        legend.text=element_text(size=13),
        axis.text = element_text(size = 15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_line(colour = "black")) + 
  facet_wrap(~ copula)  + 
  geom_line(aes(group = Model), lwd=1)

p 

ggsave(file=paste0("plot_sim_", measure,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=12, height=6)


# Part 3 - CP
df_cp_wide <- df %>% select(Estimand, Model, contains("cp"))

df_cp_long <- gather(df_cp_wide, copula, cp, cp_normal:cp_gum, factor_key=TRUE)

df_cp_long <- df_cp_long %>% mutate(copula = gsub("cp_", "", copula)) %>%
  mutate(copula = ifelse(copula == "clay", "Clayton", 
                         ifelse(copula == "normal", "Normal",
                                ifelse(copula == "frank", "Frank", "Gumbel"))))
measure = "cp"
p <- ggplot(df_cp_long, aes(x=Estimand, y=cp, colour=Model, shape=Model, size=Model)) + geom_point()
p <- p +  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c("royalblue3", "maroon3", "#808080"))+
  ylab("Coverage probability \n") + xlab("\nOutcome and covariate") +
  ylim(0,100) +
  scale_size_manual(values=c(4,4,4)) +  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  geom_hline(yintercept=95, linetype="dashed", color = "#A52A2A")+
  theme(legend.position="bottom", legend.title=element_text(size=13), 
        legend.text=element_text(size=13),
        axis.text = element_text(size = 15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_line(colour = "black")) + 
  facet_wrap(~ copula)  + 
  geom_line(aes(group = Model), lwd=1)

p 

ggsave(file=paste0("plot_sim_", measure,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=12, height=6)
