#Biomass Script for Apul_Plasticity_tp0
#recreating in a regular script instead of RMD file

#Library Installs____________________________________
# install packages if you dont already have them
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("plotrix")) install.packages("plotrix")

# load packages
library(tidyverse)
library(plotrix)

#Load Data__________________________________________
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity/data/timepoint_0")
Data <- read.csv("0_biomass/Biomass_Data.csv")
Data <- na.omit(Data)


# calculated mass per ml
#different volumes for sym (4ml) and host (5ml)
#In timepoint 0 there is no sym fraction calculated so do not need to identify values
sym <- 5
host <- 4

#Load tissue homogenate volume
homog_vol <- read.csv("0_homogenate_vols/0_homogenate_vols.csv", header=TRUE)

# Load Surface area data
sa <- read.csv("output/1_surface_area.csv")