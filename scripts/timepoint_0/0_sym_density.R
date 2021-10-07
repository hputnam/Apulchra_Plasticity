#Sym Density Script for Time Point 0 in Apul_Plasticity

## install packages if you dont already have them
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("plotrix")) install.packages("plotrix")

# load packages
library(plotrix)
library(tidyverse)

#Set working directory and import data
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")

# Cell count data
sym_counts <- read_csv("data/timepoint_0/0_sym_counts/0_sym_counts_data.csv")

# Surface area data
sa <- read.csv("output/0_surface_area.csv")

# Tissue homogenate volume data
homog_vols <- read_csv("data/timepoint_0/0_homogenate_vols/0_homogenate_vols.csv") %>% select(1:2)

# Coral sample metadata
metadata <- read_csv("metadata/coral_metadata.csv") %>% select(1:3)

# Join homogenate volumes and surface area with sample metadata
metadata <- full_join(metadata, homog_vols) %>%
  full_join(sa)

# Calculate mean counts for each sample
sym_counts1 <- sym_counts %>%
  select(colony_id, Squares.Counted, matches("Count[0-9]")) %>%
  gather("rep", "count", -colony_id, -Squares.Counted) %>%
  group_by(colony_id, Squares.Counted) %>%
  summarise(mean_count = mean(count, na.rm = TRUE))

# Join mean counts with sample metadata
sym_counts2 <- full_join(sym_counts1, metadata)

# Normalize counts by homogenat volume and surface area
sym_counts3 <- sym_counts2 %>%
  mutate(cells.mL = mean_count * 10000 / Squares.Counted,
         cells = cells.mL * homog_vol_ml,
         cells.cm2 = cells / surface.area.cm2)

# Calculate mean counts for each sample - SO FAR NOT WORKING
sym_counts4 <- sym_counts3 %>%
  select(colony_id, Squares.Counted, matches("Count[0-9]")) %>%
  gather("rep", "count", -colony_id, -Squares.Counted) %>%
  group_by(colony_id, Squares.Counted) %>%
  summarise(mean_count = mean(count, na.rm = TRUE))

# Join mean counts with sample metadata
sym_counts <- full_join(sym_counts, metadata)

# Normalize counts by homogenat volume and surface area
sym_counts <- sym_counts %>%
  mutate(cells.mL = mean_count * 10000 / Squares.Counted,
         cells = cells.mL * homog_vol_ml,
         cells.cm2 = cells / surface.area.cm2)