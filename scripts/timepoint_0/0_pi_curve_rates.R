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
if (!require("timechange")) install.packages('timechange') 

## load libraries
library(devtools)
library(LoLinR)
library(tidyverse)
library(lubridate)
library(cowplot)
library(timechange)

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
  mutate(Date = lubridate::as_date(as.character(Date), format = "%Y%m%d"))

# Select only certain columnns
metadata <- metadata %>%
  select(colony_id, Run, Chamber.Vol.L, Date, Start.time, Stop.time, Light_Value)

#we have a 7 hour time difference between Moorea time zone (notebook and metadata start and stop times) and the computer pi curve data (computer was likely in eastern time zone). Here we need to correct for the time difference by adding 7 hours to the metadata file  

metadata <- metadata %>%
  mutate(Start.time = strptime(Start.time, format="%H:%M:%S", tz="Pacific/Tahiti"))%>%
  mutate(New.start.time = timechange::time_force_tz(Start.time, tz = "Pacific/Tahiti", tzout="Atlantic/Bermuda"))%>%
  mutate(New.start.time = format(as.POSIXlt(New.start.time), format = "%H:%M:%S"))%>%
  mutate(Stop.time = strptime(Stop.time, format="%H:%M:%S", tz="Pacific/Tahiti"))%>%
  mutate(New.stop.time = timechange::time_force_tz(Stop.time, tz = "Pacific/Tahiti", tzout="Atlantic/Bermuda"))%>%
  mutate(New.stop.time = format(as.POSIXlt(New.stop.time), format = "%H:%M:%S"))

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
#DONE ON 20211028

#BREAK DOWN IS HAPPENING HERE WHERE DATA has no info
#Use the time breaks in the sample info to link O2 data with light levels

#Fixed time issue, now we need to figure out why there is an "unexpected bracket"
df <- df %>%
  mutate(intervals = map2(data0, info, function(.x, .y) {
    split(.x, f = cut(as.numeric(.x$Time), breaks = as.numeric(c(.y$New.start.time, last(.y$New.stop.time))),
                      labels = as.character(.y$Light_Value)))})) %>%
  mutate(data = map(intervals, ~ unnest(tibble(.), .id = "Light_Value")))


### Thin data
# Set thinning parameter
thin_par <- 20

# Thin data for all samples
df <- df %>%
  mutate(thin_data = map(data, ~ slice(., seq(1, nrow(.), thin_par))))

# Create plots for full dataset and thinned data
df <- df %>%
  mutate(data_plot = map2(data, colony_id, ~ ggplot(.x, aes(x = Time, y = Value)) + 
                            facet_wrap(~ as.numeric(Light_Value), scales = "free") +
                            geom_point() +
                            labs(title = .y)),
         thin_data_plot = map2(thin_data, colony_id, ~ ggplot(.x, aes(x = Time, y = Value)) + 
                                 facet_wrap(~ as.numeric(Light_Value), scales = "free") +
                                 geom_point() +
                                 labs(title = .y)))

# Example of plots
cowplot::plot_grid(df$data_plot[[1]], df$thin_data_plot[[1]], nrow = 2,
                   labels = c("Example plot: all data", "Example plot: thinned data"))
