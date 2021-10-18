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
Data <- read.csv("data/timepoint_0/0_biomass/0_biomass_data.csv")

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
Data <- full_join(Data, homog_vol) %>%
  full_join(sa) %>% full_join(metadata)

#Filter out just the Acropora species
Data <- Data%>% filter(species=="Acropora")

Data <- Data %>%
  mutate(dry.pan.mass.g.vol.corr = dry.pan.mass.g/vol.added.ml* homog_vol_ml,
         burnt.pan.mass.g.vol.corr = burnt.pan.mass.g/vol.added.ml *homog_vol_ml)


# Calculate Dry Biomass
Data <- Data %>%
  mutate(dry.biomass.g = (dry.pan.mass.g.vol.corr - initial.mass.g),
         DW.mg.cm2 = ((dry.biomass.g)*1000)/ surface.area.cm2)

Data <- Data %>%
  select(colony_id, site, species, timepoint, dry.pan.mass.g.vol.corr, burnt.pan.mass.g.vol.corr, surface.area.cm2, DW.mg.cm2) %>%
  mutate(burnt.biomass.g = (dry.pan.mass.g.vol.corr - burnt.pan.mass.g.vol.corr),
         AFDW.mg.cm2 = ((burnt.biomass.g)*1000)/ surface.area.cm2)

#Plot of the Data for each Site
Fig.6 <- Data %>%
  ggplot(aes(x = site, y = AFDW.mg.cm2, color = site)) +
  coord_cartesian(ylim = c(0, 3))+
  labs(x = "Site", y = "Ash Free Dry Weight (mg/cm2)", color = "Site") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

#Quick Stats on Site vs Protein Content for all the corals
model <- aov(log10(AFDW.mg.cm2) ~ site, data = Data)
anova(model)
TukeyHSD(model)
par(mfrow=c(2,2))
boxplot(model$residuals)
hist(model$residuals)
plot(model$fitted.values, model$residuals)

#save for concatenation to same file (both Tukey and ANOVA in same file)
afdw_res <- anova(model)
afdw_tkyhsd_res <- TukeyHSD(model)

# add ANOVA in AFDW
cat("F) ANOVA results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in AFDW
cat("F) TukeyHSD results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#condensing down output file to just the columns we want
Data <- Data %>%
  select(colony_id, site, AFDW.mg.cm2, timepoint) %>%
  write_csv(path = "output/0_biomass_output.csv")


