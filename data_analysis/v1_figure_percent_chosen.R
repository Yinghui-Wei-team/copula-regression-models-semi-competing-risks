library(dplyr); library(RColorBrewer);library(ggplot2);library(stringr)

output_dir <- "results/simulation_results/"
function_bar_plot <- function(df, copula){
  coul <- brewer.pal(3, "Pastel2") 
  p <- ggplot(df, aes(x=model,y=percent,  fill=chosen)) + 
    geom_bar(position="fill", stat="identity") + 
    ggtitle(str_to_title(copula)) + 
    guides(fill=guide_legend(title="")) +
    theme(axis.text = element_text(size =20),
          title =element_text(size=20, face='bold'),
          axis.title = element_text(size =20),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position="bottom", legend.title = element_text(size=15),
          legend.text = element_text(size=15)) + 
    scale_fill_manual(values = coul) +
    labs(title = str_to_title(copula), x = "\nUnderlying survival distributions",
         y = "Percentage of chosen by AIC\n")
  p
  ggsave(file=paste0(output_dir,"figures/figure_percent_", copula,".pdf"), 
         plot=p, width=15, height=10)
}

dist= rep(c("Exponential", "Weibull", "Gompertz"),3)
model = c(rep("Exponential", 3), rep("Weibull",3), rep("Gompertz", 3))
df <- data.frame(model, percent, dist)
df <- df %>% mutate(chosen=ifelse(model == dist, "Model chosen correctly", "Model chosen incorrectly"))

# Normal copula
df$percent = c(80.9, 8.9, 10.2, 0, 100, 0, 0, 0.2, 99.8)
copula = "normal"
function_bar_plot(df, copula)

# Clayton copula
df$percent = c(78.8, 11, 10.2, 0, 100, 0, 0, 3.7, 96.3)
copula = "clayton"
function_bar_plot(df, copula)

# Frank copula
df$percent = c(82.5, 14.4, 3.1, 0, 100, 0, 0, 8.7, 91.3)
copula = "frank"
function_bar_plot(df, copula)

# Gumbel copula
df$percent = c(79.5, 9.6, 10.9, 0.1, 99.8, 0.1, 0.1, 0.2, 99.8)
copula = "gumbel"
function_bar_plot(df, copula)