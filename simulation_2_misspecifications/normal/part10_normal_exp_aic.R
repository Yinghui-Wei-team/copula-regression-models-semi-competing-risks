################################################################################
# Run simulation by part: normal - exponential model                           #
################################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(plyr); library(survival)
part = 10; rep_start = 901; rep_end = 1000
part = paste0("part", part, "(", rep_start, ",", rep_end, ")")

## directory if on own PC
#dir_results = "../../results/simulation_results/"
#source("simulation_2_misspecifications/normal/source_normal_exp_aic.R")

## directory if on cluster
dir_results = "/home/ywei/Simulation/Paper2/Normal/simulation2/"
setwd(dir_results)
source("source_normal_exp_aic.R")
