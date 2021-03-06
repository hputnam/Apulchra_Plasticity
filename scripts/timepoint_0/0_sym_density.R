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
sym_counts2 <- full_join(sym_counts1, metadata, by="colony_id")
sym_counts2 <- sym_counts2%>% filter(species=="Acropora")

# Normalize counts by homogenate volume and surface area
sym_counts3 <- sym_counts2 %>%
  mutate(cells.mL = mean_count * 10000 / Squares.Counted,
         cells = cells.mL * homog_vol_ml,
         cells.cm2 = cells / surface.area.cm2)

#Plot of the Data for each Site
Fig.5 <- sym_counts3 %>%
  ggplot(aes(x = site, y = cells.cm2, color = site)) +
  coord_cartesian(ylim = c(0, 2.5e6))+
  labs(x = "Site", y = "Symbiont Density (cells/cm2)", color = "Site") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

#Quick Stats on Site vs Protein Content for all the corals
model5 <- aov(log10(cells.cm2) ~ site, data = sym_counts3)
anova(model5)
TukeyHSD(model5)
par(mfrow=c(2,2))
boxplot(model5$residuals)
hist(model5$residuals)
plot(model5$fitted.values, model5$residuals)

#save for concatenation (both Tukey and ANOVA to same output file)
sym.counts_res <- anova(model5)
sym.counts_tkyhsd_res <- TukeyHSD(model5)

# add SYM COUNTS
cat("C) ANOVA results of Sym.Counts (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym.counts_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


# add TukeyHSD in SYM COUNTS
cat("C) TukeyHSD results of Sym.Counts (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym.counts_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#Output of Summarized Symbiont cells by cm2
sym_counts3 %>%
  select(colony_id, site, cells.cm2) %>%
  mutate(timepoint="timepoint0")%>%
  write_csv(., path = "output/0_sym_counts.csv")


#re-join with metadata to select genotype information
# Coral sample metadata
sym_counts4 <- sym_counts3 %>%
  select(colony_id, site, cells.cm2) %>%
  mutate(timepoint="timepoint0")

metadata1 <- read_csv("metadata/coral_metadata.csv") %>% 
  select(1:3,6) %>%
  filter(species == "Acropora")

sym_counts_4geno <- full_join(sym_counts4, metadata1)

# Nursery 4 genotypes
sym_counts_4geno <- sym_counts_4geno %>%
  filter(genotype == "Genotype15"| genotype == "Genotype4"| genotype == "Genotype6"|genotype == "Genotype8") %>%
  select(colony_id, site, genotype, cells.cm2, timepoint) %>%
  write_csv(., path = "output/0_sym_counts_4geno.csv")

