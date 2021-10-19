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


path.p <- "data/timepoint_0/0_pi_curves/" #the location of all your respirometry files 

# List data files
file.names <- list.files(path = path.p, pattern = "csv$")  # list all csv file names in the folder
file.names <- file.names[!grepl("metadata", file.names)]   # omit metadata from files to be read in as data
