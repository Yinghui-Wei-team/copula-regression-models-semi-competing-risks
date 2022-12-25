library(ggplot2); library(dplyr)
output_dir <- "results/simulation_results/"
df <- read.csv("results/simulation_results/table5_simulation.csv")

df <- df %>% mutate(model = ifelse(model=="copula1", "Normal copula with covariates on hazards",
                                   ifelse(model=="copula2", "Normal copula with covariates on hazards and association",
                                          "Cox model")))
df <- df %>% rename(Model = model)

# Part 1 - Bias
measure = "bias"
p <- ggplot(df, aes(x=Estimand, y=bias_normal, colour=Model, shape=Model, size=Model)) + geom_point()
p <- p +  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c("royalblue3", "maroon3", "#808080"))+
  ylab("Bias \n") + xlab("\nOutcomes and covariates") +
  ylim(-0.5,0.5) +
  scale_size_manual(values=c(4,4,4)) +  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  geom_hline(yintercept=0, linetype="dashed", color = "#A52A2A")+
  theme(legend.position="bottom", legend.title=element_text(size=13), 
        legend.text=element_text(size=13),
        axis.text = element_text(size = 13),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p

ggsave(file=paste0("plot_sim_", measure,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=12, height=6)


# Part 2 - MSE
measure = "mse"
p <- ggplot(df, aes(x=Estimand, y=mse_normal, colour=Model, shape=Model, size=Model)) + geom_point()
p <- p +  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c("royalblue3", "maroon3", "#808080"))+
  ylab("Mean square error \n") + xlab("\nOutcome and covariate") +
  ylim(-0.5,0.5) +
  scale_size_manual(values=c(4,4,4)) +  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  geom_hline(yintercept=0, linetype="dashed", color = "#A52A2A")+
  theme(legend.position="bottom", legend.title=element_text(size=13), 
        legend.text=element_text(size=13),
        axis.text = element_text(size = 13),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p

ggsave(file=paste0("plot_sim_", measure,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=12, height=6)


# Part 2 - CP
measure = "cp"
p <- ggplot(df, aes(x=Estimand, y=cp_normal, colour=Model, shape=Model, size=Model)) + geom_point()
p <- p +  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c("royalblue3", "maroon3", "#808080"))+
  ylab("Coverage probability \n") + xlab("\nOutcome and covariate") +
  ylim(10,100) +
  scale_size_manual(values=c(4,4,4)) +  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  geom_hline(yintercept=95, linetype="dashed", color = "#A52A2A")+
  theme(legend.position="bottom", legend.title=element_text(size=13), 
        legend.text=element_text(size=13),
        axis.text = element_text(size = 13),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p

ggsave(file=paste0("plot_sim_", measure,".pdf"), path = paste0(output_dir, "figures"),
       plot=p, width=12, height=6)
