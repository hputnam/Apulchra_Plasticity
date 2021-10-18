#Chlorophyll Concentration Code (merging E5_ROL with urol_e5 scripts for Oct 2019 data)
#Script edited from urol_e5 repository to fit TP_0 (Ariana Huffmeyer and Dennis Conetta)

rm(list=ls()) #clears workspace 

## install packages if you dont already have them
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("plotrix")) install.packages("plotrix")

# load packages
library(tidyverse)
library(plotrix)
library(broom)
library(ggplot2)

#set working directory and save tibble

setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")

chl_data <- read.csv("data/timepoint_0/0_chlorophyll/0_20191031_Chl.csv")

# Define function to read in chl data

read_chl <- function(file) {
  chl_data <- read_csv(file, skip = 24, n_max = 24) %>%
    select(-1) %>%
    magrittr::set_colnames(c("row", 1:12, "wavelength")) %>%
    fill(row) %>%
    gather("col", "absorbance", -wavelength, -row) %>%
    unite("well", c(row, col), sep = "")
}

# List chlorophyll data files
chl_path <- "data/timepoint_0/0_chlorophyll/"                                        # Path to chlorophyll data directory
all_chl_files <- list.files(path = chl_path, pattern = "*.csv")          # List all files in directory
chl_platemaps <- list.files(path = chl_path, pattern = "platemap")       # List platemap files
chl_data_files <- setdiff(all_chl_files, chl_platemaps)                  # List absorbance data files

# Read in all files into tibble
df <- tibble(file = chl_data_files) %>%
  mutate(platemap = map(file, ~ read_csv(paste0(chl_path, tools::file_path_sans_ext(.), "_platemap.csv"))),
         chl_data = map(file, ~ read_chl(paste0(chl_path, .))))

# Merge platemap and data for each plate
df <- df %>%
  mutate(merged = map2(platemap, chl_data, ~ right_join(.x, .y)))

# average all technical replicates for each plate/sample/wavelength, including all acetone blanks together (per plate)
df <- df %>%
  unnest(merged) %>%
  filter(!is.na(colony_id)) %>%                         # remove empty wells (colony_id is NA)
  group_by(file, colony_id, wavelength) %>%
  summarise(n = n(), mean_abs = mean(absorbance)) %>%
  spread(wavelength, mean_abs)


# get the acetone blank 750 absorbace for each file (i.e., plate), and subtract from 630 and 663 values for each sample
df <- df %>%
  group_by(file) %>%
  mutate(blank750 = `750`[colony_id == "BK"]) %>%
  ungroup() %>%
  mutate(adj630 = `630` - blank750,
         adj663 = `663` - blank750)

# calculate chla and chlc2 values based on equations from Jeffrey and Humphrey 1975
# units ?g/ml
df <- df %>%
  mutate(chla.ug.ml = 11.43 * adj663 - 0.64 * adj630,
         chlc2.ug.ml = 27.09 * adj630 - 3.63 * adj663)


# Load homogenate volume

homog.vol <- read_csv("data/timepoint_0/0_homogenate_vols/0_homogenate_vols.csv") %>%
  select(colony_id, homog_vol_ml)

#merge df and homog.vol to normalize to total volume of homogenate
chl_data <- full_join(df, homog.vol)


# Load surface area
sa <- read_csv("output/0_surface_area.csv")
chl_data <- full_join(chl_data, sa)


# Multiply chlorophyll by the homogenate volume and divide by surface area
chl_data <- chl_data %>%
  mutate(chla.ug.cm2 = chla.ug.ml * homog_vol_ml / surface.area.cm2,
         chlc2.ug.cm2 = chlc2.ug.ml * homog_vol_ml / surface.area.cm2)

# remove blanks and NAs
chl_data <- filter(chl_data, !colony_id %in% c("NA", "BK", "Water"))

#read in coral_metadata and merge with species
metadata <- read.csv("metadata/coral_metadata.csv")

#rename sample_id column because it is having issues
metadata <- metadata %>%
  mutate(colony_id = ï..colony_id) %>%
  select(colony_id, species, site, genotype, purpose)
  

chl_data <- full_join(chl_data, metadata)

#filter out for only Acropora species
chl_data <- chl_data %>%
  filter(species == "Acropora")

#graphing the Chlorophyll data for each site

#Plot of the CHL_a Data for each Site
Fig.1 <- chl_data %>%
  ggplot(aes(x = site, y = chla.ug.cm2, color = site)) +
  coord_cartesian(ylim = c(0, 5))+
  labs(x = "Site", y = "Chl A Concentration (ug/cm2)", color = "Site") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

#Quick Stats on Site vs chla [] for all the corals
model2 <- aov(log10(chla.ug.cm2) ~ site, data = chl_data)
anova(model2)
TukeyHSD(model2)
par(mfrow=c(2,2))
boxplot(model2$residuals)
hist(model2$residuals)
plot(model2$fitted.values, model2$residuals)

#Output Stats to Referencein Concatenation (Tukey and ANOVA in same output file)
chla_res <-anova(model2)
chla_tkyhsd_res <- TukeyHSD(model2)

# add ANOVA in CHL A data
cat("A) ANOVA results of Chl a (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chla_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in CHL A data
cat("A) TukeyHSD results of Chl a (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chla_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)



#Plot of the CHL_a Data for each Site
Fig.2 <- chl_data %>%
  ggplot(aes(x = site, y = chlc2.ug.cm2, color = site)) +
  coord_cartesian(ylim = c(0, 8)) +
  labs(x = "Site", y = "Chl C Concentration (ug/cm2)", color = "Site") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

#Quick Stats on Site vs chlc [] for all the corals
model3 <- aov(log10(chlc2.ug.cm2) ~ site, data = chl_data)
anova(model3)
TukeyHSD(model3)
par(mfrow=c(2,2))
boxplot(model3$residuals)
hist(model3$residuals)
plot(model3$fitted.values, model3$residuals)

#Output Stats to Reference Later in Concatenation (both TUkey and ANOVA in same file)
chlc2_res <- anova(model3)
chlc2_tkyhsd_res <- TukeyHSD(model3)

#add ANOVA in CHLC data
cat("B) ANOVA results of Chl c (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chlc2_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#add TukeyHSD in CHLC data
cat("B) TukeyHSD results of Chl c (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chlc2_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


# write chlorophyll data to output csv file
chl_data %>%
  select(colony_id, site, chla.ug.cm2) %>%
  mutate(timepoint="timepoint0")%>%
  write_csv(path = "output/0_chlorophylla.csv")

# write chlorophyll data to output csv file
chl_data %>%
  select(colony_id, site,  chlc2.ug.cm2) %>%
  mutate(timepoint="timepoint0")%>%
  write_csv(path = "output/0_chlorophyllc.csv")







