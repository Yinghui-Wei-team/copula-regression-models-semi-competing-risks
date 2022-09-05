# Purpose: produce stacked bar plot to show percentage correctly chosen by AIC for the underlying survival distributions
# Programmed by Yinghui Wei
library(dplyr); library(RColorBrewer);library(ggplot2);library(stringr)

output_dir <- "results/simulation_results/"


dist= rep(c("Exponential", "Weibull", "Gompertz"),3)
model = c(rep("Exponential", 3), rep("Weibull",3), rep("Gompertz", 3))
df <- data.frame(model, dist)
df <- df %>% mutate(chosen=ifelse(model == dist, "Model chosen correctly", "Model chosen incorrectly"))

# Normal copula
df$percent = c(80.9, 8.9, 10.2, 0, 100, 0, 0, 0.2, 99.8)
copula = "normal"
df$copula = "Normal"
df_normal = df

# Clayton copula
df$percent = c(78.8, 11, 10.2, 0, 100, 0, 0, 3.7, 96.3)
df$copula = "Clayton"
df_clayton = df


# Frank copula
df$percent = c(82.5, 14.4, 3.1, 0, 100, 0, 0, 8.7, 91.3)
df$copula = "Frank"
df_frank = df


# Gumbel copula
df$percent = c(79.5, 9.6, 10.9, 0.1, 99.8, 0.1, 0.1, 0.2, 99.8)
df$copula = "Gumbel"
df_gumbel = df

# combine all data into one df
df <- rbind(df_normal, df_clayton, df_frank, df_gumbel)

coul <- brewer.pal(3, "Pastel2")
p <- ggplot(df, aes(x=model,y=percent,  fill=chosen)) + 
  geom_bar(position="fill", stat="identity") + 
  guides(fill=guide_legend(title="")) +
  theme(axis.text = element_text(size =20),
        title =element_text(size=20, face='bold'),
        axis.title = element_text(size =20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom", legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        panel.spacing = unit(2, "lines"),
        strip.text = element_text(size = 20)) + 
  facet_wrap(~ copula, ncol = 2) +
  scale_fill_manual(values = coul) + 
  labs(x = "\nUnderlying survival distributions",
       y = "Percentage of chosen by AIC\n")
p
copula = "all"
ggsave(file=paste0(output_dir,"figures/figure_percent_", copula,".pdf"), 
       plot=p, width=15, height=10)