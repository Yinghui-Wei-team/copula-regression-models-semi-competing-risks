################################################################################
# Run the second quarter of the replication: normal - exponential model        #
################################################################################
rm(list=ls())
rep_start = 251
rep_end = 500

## directory if on own PC
# dir_results = "../../results/simulation_results/"
# source("simulation_1_cox_underlying_copula/normal/source_t_l_full_normal_exp.R")

## directory if on cluster
dir_results = "/home/ywei/Simulation/Paper2/Normal/"
setwd(dir_results)
source("source_t_l_full_normal_exp.R")
