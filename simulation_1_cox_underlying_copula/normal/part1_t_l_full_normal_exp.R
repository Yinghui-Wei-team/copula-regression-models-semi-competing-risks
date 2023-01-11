################################################################################
# Run the first quarter of the replication: normal - exponential model         #
################################################################################
rm(list=ls())
rep_start = 1
rep_end = 250

# directory if on own PC
# dir_results = "../../results/simulation_results/"
# source("simulation_1_cox_underlying_copula/normal/source_t_l_full_normal_exp.R")

## directory if on cluster
dir_results = "/home/ywei/Simulation/Paper2/Normal"
setwd(dir_results)
##if source file is under the above directory (dir)
source("source_t_l_full_normal_exp.R")



