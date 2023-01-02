################################################################################
# Run the fourth quarter of the replication: normal - exponential model        #
################################################################################
rm(list=ls())
rep_start = 751
rep_end = 1000

# directory if on own PC
dir_results = "../../results/simulation_results/"
source("simulation_1_cox_underlying_copula/normal/source_t_l_full_normal_exp.R")

# directory if on cluster
# dir = "/home/ywei/Simulation/Paper2/Normal"
# setwd(dir)
# source("source_t_l_full_normal_exp.R")
