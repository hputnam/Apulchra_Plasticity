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
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")
Data <- read.csv("data/timepoint_0/0_biomass/Biomass_Data.csv")
Data <- na.omit(Data)


# calculated mass per ml
#different volumes for sym (4ml) and host (5ml)
#In timepoint 0 there is no sym fraction calculated so do not need to identify values
sym <- 5
host <- 4

#Load tissue homogenate volume
homog_vol <- read.csv("data/timepoint_0/0_homogenate_vols/0_homogenate_vols.csv", header=TRUE)

# Load Surface area data
sa <- read.csv("output/0_surface_area.csv")

# Coral sample metadata
metadata <- read_csv("metadata/coral_metadata.csv") %>% select(1:3)

#colony_id column messed up in naming so go back and name it again
homog_vol <- homog_vol %>%
  mutate(colony_id = ï..colony_id)

# Join homogenate volumes and surface area with sample metadata
metadata <- full_join(metadata, homog_vol) %>%
  full_join(sa)
