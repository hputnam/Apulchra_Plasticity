#Protein Analysis Script fitted for TP0
#Modified script from urol_e5 to fit Apul Oct 2019 data (Ariana Huffmeyer and Dennis Conetta)
#HOLOBIONT PROTEIN 

rm(list=ls()) #clears workspace 

### install packages if you dont already have them
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("broom")) install.packages("broom")

# load packages
library(tidyverse)
library(broom)

#setworking directory 
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")

#read data file
holo_prot_data <- read.csv("data/timepoint_0/0_protein/Holobiont/0_20191031_Plate2_BCA_statistics.csv")

# List protein data files
prot_path = "data/timepoint_0/0_protein/Holobiont/"                                              # Path to prot data directory
all_prot_files <- list.files(path = prot_path, pattern = "*.csv")          # List all files in directory
prot_platemaps <- list.files(path = prot_path, pattern = "platemap")       # List platemap files
prot_data_files <- setdiff(all_prot_files, prot_platemaps)                 # List data files

# Read in all files into tibble ; in this case Plate 1 is the Host ONLY and Plate 2 is the HOLOBIONT for these same corals

df <- tibble(file = prot_data_files) %>%
  separate(file, into = c("trip", "date", "plate"), remove = FALSE) %>%
  unite(plate, trip, date, plate) %>%
  mutate(platemap = map(plate, ~read_csv(paste0(prot_path, ., "_BCA_platemap.csv"))),
         prot_data = map(file, ~read_csv(paste0(prot_path, .)) %>% rename(well = Well)))

# Merge platemap and data for each plate
df <- df %>%
  mutate(merged = map2(platemap, prot_data, ~ right_join(.x, .y)))


##PLOT STANDARD CURVE
# Create standard curve following kit instructions
standards <- tribble(
  ~std, ~BSA_ug.mL,
  "A",        2000,
  "B",        1500,
  "C",        1000,
  "D",         750,
  "E",         500,
  "F",         250,
  "G",         125,
  "H",          25,
  "I",           0
)

std_curve <- df %>%
  unnest(merged) %>%
  filter(grepl("Standard", colony_id)) %>%
  select(plate, well, colony_id, abs562 = `A562`) %>%
  rename(std = colony_id) %>%
  mutate(std = str_sub(std, 9, 9)) %>%
  #group_by(std) %>%
  #summarise(abs562 = mean(abs562)) %>%                       # calculate mean of standard duplicates
  #mutate(abs562.adj = abs562 - abs562[std == "I"]) %>%       # subtract blank absorbace value from all
  left_join(standards)


## Fit linear model for standard curve
mod <- lm(BSA_ug.mL ~ abs562, data = std_curve)
coef(mod)
fitted <- mod %>% broom::augment()

## Fit nonlinear model for standard curve - DOESNT WORK FOR ME (DC - 20211013)
#mod <- nls(formula = BSA_ug.mL ~ (z + a * exp(b * abs562)), start = list(z = 0, a = 1, b = 1), data = std_curve)
#fitted <- mod %>% broom::augment()

# Plot standard curve
std_curve_plot <- std_curve %>%
  ggplot(aes(x = abs562, y = BSA_ug.mL)) +
  geom_point(color = "red", size = 3)


#STD CURVE LINE TO PLOT 
std_curve_plot + 
  geom_line(data = fitted, aes(x = abs562, y = .fitted)) +
  labs(title = "Standard curve")


# Calculate protein concentration for all samples using standard curve
holo_prot <- df %>%
  unnest(merged) %>%
  filter(!grepl("Standard", colony_id)) %>%                     # Get just samples (not standards)
  select(plate, well, colony_id, abs562 = `A562`) %>%        # Select only needed columns
  filter(!is.na(colony_id)) %>%                                 # Filter out empty wells
  mutate(prot_ug.mL = map_dbl(abs562, ~ predict(mod, newdata = data.frame(abs562 = .))))    # Use standard curve to convert absorbance to protein

std_curve_plot + 
  geom_point(data = holo_prot, aes(x = abs562, y = prot_ug.mL), pch = "X", cex = 5, alpha = 0.3) +
  geom_line(data = fitted, aes(x = abs562, y = .fitted)) +
  labs(title = "All samples projected on standard curve")


# Surface area data
sa <- read.csv("output/0_surface_area.csv")
# Tissue homogenate volume data
homog_vols <- read_csv("data/timepoint_0/0_homogenate_vols/0_homogenate_vols.csv") %>% select(1:2)

# Coral sample metadata
metadata <- read_csv("metadata/coral_metadata.csv") 

# Join homogenate volumes and surface area with sample metadata
metadata <- full_join(metadata, homog_vols) %>%
  full_join(sa)

# Join prot data with metadata
holo_prot <- left_join(holo_prot, metadata) %>%
  mutate(prot_ug = prot_ug.mL * homog_vol_ml,
         prot_ug.cm2 = prot_ug / surface.area.cm2,
         prot_mg.cm2 = prot_ug.cm2 / 1000)

holo_prot <- holo_prot %>% 
  filter(species == "Acropora") %>%
  filter(genotype != "NA")
  


holo_prot <- holo_prot %>% 
  group_by(colony_id, site, genotype) %>%
  summarise(prot_ug = mean(prot_ug, na.rm = T),
   prot_ug.cm2 = mean(prot_ug.cm2, na.rm = T)) 

#Plot of the Data for each Site
Fig.4 <- holo_prot %>%
  ggplot(aes(x = site, y = prot_ug.cm2, color = site)) +
  coord_cartesian(ylim = c(0, 600))+
  labs(x = "Site", y = "Total Holobiont protein (µg/cm2)", color = "Site") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

#Quick Stats on Site vs Protein Content for all the corals
model4 <- aov(log10(prot_ug.cm2) ~ site, data = holo_prot)
anova(model4)
TukeyHSD(model4)
par(mfrow=c(2,2))
boxplot(model4$residuals)
hist(model4$residuals)
plot(model4$fitted.values, model4$residuals)

# All sites
holo_prot %>%
  group_by(colony_id, site) %>%
  summarise(prot_ug = mean(prot_ug, na.rm = T),
            prot_ug.cm2 = mean(prot_ug.cm2, na.rm = T)) %>%
  select(colony_id, site, prot_ug, prot_ug.cm2) %>%
  mutate(timepoint="timepoint0")%>%
  write_csv(., path = "output/0_holobiont_protein.csv")

# Nursery 4 genotypes
holo_prot_4geno <- holo_prot %>%
  filter(genotype == "Genotype15"| genotype == "Genotype4"| genotype == "Genotype6"|genotype == "Genotype8") %>%
  write_csv(., path = "output/0_holobiont_protein_4geno.csv")
