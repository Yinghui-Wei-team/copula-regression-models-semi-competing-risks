################################################################################
# Run simulation by part: normal - gompertz model                             #
################################################################################
rm(list=ls())
library(copula); library(mvtnorm); library(plyr); library(survival)
rep_start = 101; rep_end = 200; part = round(rep_end/100,0)
part = paste0("part", part, "(", rep_start, ",", rep_end, ")")

## directory if on own PC
#dir_results = "../../results/simulation_results/"
#source("simulation_2_misspecifications/normal/source_normal_gompertz_aic.R")

## directory if on cluster
dir_results = "/home/ywei/Simulation/Paper2/Normal/simulation2/"
setwd(dir_results)
source("source_normal_gompertz_aic.R")
