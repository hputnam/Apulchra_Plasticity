#PI Curve Rate Script for Apul Plasticity TP0 
#Scripts adapted from urol_e5 timepoint_1 Scripts

## install packages if you dont already have them in your library
if (!require("devtools")) install.packages("devtools")
if (!require("furrr")) install.packages("furrr")
if (!require("future")) install.packages("future")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("lubridate")) install.packages("lubridate")
if (!require("cowplot")) install.packages("cowplot")
if (!require("LoLinR")) install_github('colin-olito/LoLinR') 

## load libraries
library(devtools)
library(LoLinR)
library(tidyverse)
library(lubridate)
library(cowplot)

## libraries for parallel processing
library(future)
library(furrr)

#Set Working Directory
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")

#edit down metadata file and output just for ACR
pi_metadata <- read.csv("data/timepoint_0/0_pi_curves/All_PI_Curve_Sample_Info.csv")

acr_only_data <- pi_metadata %>%
  mutate(Date = ï..Date) %>%
  filter(Date == 20191024 | Date == 20191025 | Date == 20191029) %>%
  filter(Species == "Acropora" | Species == "Blank") %>%
  write_csv(path = "data/timepoint_0/0_pi_curves/acr_only_pi_curve_metadata.csv")
  

path.p <- "data/timepoint_0/0_pi_curves/" #the location of all your respirometry files 

# List data files
file.names <- list.files(path = path.p, pattern = "csv$")  # list all csv file names in the folder
file.names <- file.names[!grepl("metadata", file.names)]   # omit metadata from files to be read in as data


# Load PI curve sample metadata (i.e., which corals were in which runs)
sample.info <- read_csv(file = "data/timepoint_0/0_pi_curves/acr_only_pi_curve_metadata.csv")

# Load PI curve run metadata (i.e., light levels and interval times for each run)
run.info <- read_csv(file = "data/1_pi_curves/1_pi_curves_run_metadata.csv")
