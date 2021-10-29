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

# Read in all data files - 
#fixed this line of code by deleting column name artifacts in data sheets
df <- tibble(file.name = file.names) %>%
  mutate(colony_id = gsub("_.*", "", file.name),                              # Get colony_id from filename
         info = map(colony_id, ~filter(metadata, colony_id == .)),           # Get associated sample info
         data0 = map(file.name, ~read_csv(file.path(path.p, .), skip = 1)))   # Get associated O2 data

# Select only Time, Value, and Temp columns from O2 data
df <- df %>%
  mutate(data0 = map(data0, ~select(., Time, Value, Temp))) 

df <- df %>%
  mutate(intervals = map2(data0, info, function(.x, .y) {
    split(.x, f = cut(as.numeric(.x$Time), breaks = as.numeric(c(.y$Start.time, last(.y$Stop.time))),
                      labels = as.character(.y$Light_Value)))})) %>%
  mutate(data = map(intervals, ~ unnest(tibble(.), .id = "Light_Value")))


### Thin data
# Set thinning parameter
thin_par <- 10

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

# Fit regressions to each interval for each sample

# Define function for fitting LoLinR regressions to be applied to all intervals for all samples
fit_reg <- function(df) {
  rankLocReg(xall = as.numeric(df$Time), yall = df$Value, 
             alpha = 0.2, method = "pc", verbose = FALSE)
}

# Setup for parallel processing
future::plan(multiprocess)

# Map LoLinR function onto all intervals of each sample's thinned dataset
df <- df %>%
  mutate(regs = furrr::future_map(thin_data, function(.) {       # future_map executes function in parallel
    group_by(., Light_Value) %>%
      do(rankLcRg = fit_reg(.))
  }))

## Now 'regs' contains the fitted local regressions for each interval of each sample's thinned dataset

# Define function to pull out and plot regression diagnostics
plot_rankLcRg <- function(colony_id, interval_number) {
  df %>%
    filter(colony_id == colony_id) %>%
    pluck("regs", 1, "rankLcRg", interval_number) %>%
    plot()
}

### Extract slope of best regression for each interval for each sample

df.out <- df %>% 
  unnest(regs) %>%
  mutate(micromol.L.s = map_dbl(rankLcRg, ~ pluck(., "allRegs", "b1", 1)))

# Adjust by chamber volume and normalize to surface area

### Merge rates with sample info
pr <- left_join(
  select(df.out, colony_id, Light_Value, micromol.L.s),
  distinct(metadata, colony_id, Run, Chamber.Vol.L)
)

### Correct for chamber volume and blanks
pr <- pr %>%
  mutate(micromol.s = micromol.L.s * Chamber.Vol.L)

# Get blank values -- average for each run and light value in case multiple blanks
blanks <- pr %>%
  filter(grepl("BLK", colony_id)) %>%
  group_by(Light_Value) %>%
  summarise(avgmicromol.s.blank = mean(micromol.s))

# Plot rates vs. irradiance for each sample

ggplot(blanks, aes(x = as.numeric(Light_Value), y = avgmicromol.s.blank)) +
  geom_point() 


#### Average all blanks and use that to subtract from all
# Join blank values with rest of data and subtract values from samples for same run and light value
pr <- left_join(pr, blanks) %>%
  mutate(micromol.s.adj = micromol.s - avgmicromol.s.blank) %>%
  # After correcting for blank values, remove blanks from data
  filter(!grepl("BLK", colony_id))

# Import surface area data
sa <- read.csv("output/0_surface_area.csv")
sa <- sa %>%
  mutate_at("colony_id", str_replace, "_", "-")

# Join surface area with rest of data
pr <- left_join(pr, select(sa, colony_id, surface.area.cm2))

# Normalize rates by surface area
pr <- pr %>%
  mutate(micromol.cm2.s = micromol.s.adj / surface.area.cm2,
         micromol.cm2.h = micromol.cm2.s * 3600)


# Plot rates vs. irradiance for each sample

ggplot(pr, aes(x = as.numeric(Light_Value), y = micromol.cm2.h)) +
  geom_point() +
  facet_wrap(~colony_id, scale = "free_y", ncol = 4)


# Select variables to write to file
pr.out <- pr %>% 
  select(colony_id, Light_Value, Run, micromol.cm2.s, micromol.cm2.h) %>%
  mutate(timepoint="timepoint0")

# Write to output file
write.csv(pr.out, "output/0_pi_curve_rates.csv")











