################################################################################
# Run the first quarter of the replication: normal - exponential model         #
################################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(plyr); library(survival)
rep_start = 1
rep_end = 250

# directory if on own PC
dir_results = "../../results/simulation_results/"
source("simulation_1_cox_underlying_copula/normal/source_t_l_full_normal_exp.R")

# directory if on cluster
# dir = "/home/ywei/Simulation/Paper2/Normal"
# setwd(dir)
# source("source_t_l_full_normal_exp.R")
