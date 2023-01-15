################################################################################
# Run the third quarter of the replication: normal - exponential model         #
################################################################################
rm(list=ls())
rep_start = 501
rep_end = 750

## directory if on own PC
# dir_results = "../../results/simulation_results/"
# source("simulation_1_cox_underlying_copula/normal/source_t_l_full_normal_exp.R")

## directory if on cluster
dir_results = "/home/ywei/Simulation/Paper2/Normal/"
setwd(dir_results)
source("source_t_l_full_normal_exp.R")
