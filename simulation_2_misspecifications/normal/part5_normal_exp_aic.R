################################################################################
# Run simulation by part: normal - exponential model                           #
################################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(plyr); library(survival)
part = 5; rep_start = 401; rep_end = 500
part = paste0("part", part, "(", rep_start, ",", rep_end, ")")

# directory if on own PC
dir_results = "../../results/simulation_results/"
source("simulation_2_misspecifications/normal/source_normal_exp_aic.R")

# directory if on cluster
# dir = "/home/ywei/Simulation/Paper2/Normal"
# setwd(dir)
# source("source_normal_exp_aic.R")
