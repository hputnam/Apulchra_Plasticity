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
#setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")

path.p <- "data/timepoint_0/0_pi_curves" #the location of all your respirometry files 

# List data files
file.names <- list.files(path = path.p, pattern = "csv$")  # list all csv file names in the folder
file.names <- file.names[!grepl("metadata", file.names)]   # omit metadata from files to be read in as data

# Load PI curve sample metadata (i.e., which corals were in which runs)
sample.info <- read_csv(file = "data/timepoint_0/0_pi_curves/0_pi_curves_sample_metadata.csv")


# Load PI curve run metadata (i.e., light levels and interval times for each run)
run.info <- read_csv(file = "data/timepoint_0/0_pi_curves/0_pi_curves_run_metadata.csv")


# Join all coral and run metadata
metadata <- full_join(sample.info, run.info) %>%
  mutate(Date = as_date(as.character(Date), format = "%Y%m%d", tz = "Tahiti"))

# Select only certain columnns
metadata <- metadata %>%
  select(colony_id, Run, Chamber.Vol.L, Date, Start.time, Stop.time, Light_Value)


# Read in all data files - 
#fixed this line of code by deleting column name artifacts in data sheets
df <- tibble(file.name = file.names) %>%
  mutate(colony_id = gsub("_.*", "", file.name),                              # Get colony_id from filename
         info = map(colony_id, ~filter(metadata, colony_id == .)),           # Get associated sample info
         data0 = map(file.name, ~read_csv(file.path(path.p, .), skip = 1)))   # Get associated O2 data

# Select only Time, Value, and Temp columns from O2 data
df <- df %>%
  mutate(data0 = map(data0, ~select(., Time, Value, Temp))) 

####Before running lines below we need to make sure that the start and stop time formats match the raw data time format


#Use the time breaks in the sample info to link O2 data with light levels
df <- df %>%
  mutate(intervals = map2(data0, info, function(.x, .y) {
    split(.x, f = cut(as.numeric(.x$Time), breaks = as.numeric(c(.y$Start.time, last(.y$Stop.time))),
                      labels = as.character(.y$Light_Value)))})) %>%
  mutate(data = map(intervals, ~ unnest(tibble(.), .id = "Light_Value")))

