# Purpose: create contour plots and hazard function plots for paper2
# script by MW, 25/08/2021; edited by YW, 29/09/2021
# YW updated 31/08/2022
library(copula)

output_dir <- "results/real_data_analysis/"

source("data_analysis/source_model_parameter.R")

#################################################################
# Part 2. Contour Plots                                         #
################################################################
#----------------------------------------------------------------------
# 1=(Age>50); 1=genFemale; 1=DonorLiving
# Reference: male, old, deceased
age.grp <- 1; gen<- 0; donor <- 0 # patient with highest risks
outfile = paste0(output_dir, "contour_frankweibull_OldMD.pdf")
source("data_analysis/source_contour_plot.R")

# male, old, living
age.grp <- 1; gen<- 0; donor <- 1 # 
outfile = paste0(output_dir, "contour_frankweibull_OldML.pdf")
source("data_analysis/source_contour_plot.R")

#----------------------------------------------------------------------
# female, young, living
age.grp <- 0; gen<- 1; donor <- 1 # patient with lowest risks
outfile = paste0(output_dir, "contour_frankweibull_YoungFL.pdf")
source("data_analysis/source_contour_plot.R")

#female, young, deceased
age.grp <- 0; gen<- 1; donor <- 0 #
outfile = paste0(output_dir, "contour_frankweibull_YoungFD.pdf")
source("data_analysis/source_contour_plot.R")

