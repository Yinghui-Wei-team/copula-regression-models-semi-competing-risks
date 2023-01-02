################################################################################
# Run simulation by part: normal - weibull model                               #
################################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(plyr); library(survival)
rep_start = 701; rep_end = 800; part = round(rep_end/100,0)
part = paste0("part", part, "(", rep_start, ",", rep_end, ")")

# directory if on own PC
dir_results = "../../results/simulation_results/"
source("simulation_2_misspecifications/normal/source_normal_weibull_aic.R")

# directory if on cluster
# dir = "/home/ywei/Simulation/Paper2/Normal"
# setwd(dir)
# source("source_normal_weibull_aic.R")
