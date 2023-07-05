################################################################################
# Run simulation by part: normal - exponential model                           #
################################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(plyr); library(survival)
part = 9; rep_start = 801; rep_end = 900
part = paste0("part", part, "(", rep_start, ",", rep_end, ")")

## directory if on own PC
#dir_results = "../../results/simulation_results/"
# source("functions/function_sim2.R")
#source("simulation_2_misspecifications/normal/source_normal_exp_aic.R")

## directory if on cluster
dir_results = "/home/ywei/Simulation/Paper2/Normal/simulation2/"
setwd(dir_results)
source("function_sim2.R")
source("source_normal_exp_aic.R")
