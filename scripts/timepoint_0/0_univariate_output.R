#Summary of TP0 Univariate Response Variables (Except AQY, Am, and R)
library(tidyverse)
library(ggpubr)
library(showtext)
library(dplyr)

#Set working directory
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")


#All saved Figures composed into one figure
chla <- read.csv("output/0_chlorophylla.csv")
chla$variable <- colnames(chla[3])
colnames(chla)[3] <- "value"


chlc2 <- read.csv("output/0_chlorophyllc.csv")
chlc2$variable <- colnames(chlc2[3])
colnames(chlc2)[3] <- "value"
#CHL C stats
model3 <- aov(log10(chlc2.ug.cm2) ~ site, data = chl_data)
anova(model3)
TukeyHSD(model3)

prot.holo <- read.csv("output/0_holobiont_protein.csv")
prot.holo$variable <- colnames(prot.holo[3])
colnames(prot.holo)[3] <- "value"

prot.host <- read.csv("output/0_host_protein.csv")
prot.host$variable <- colnames(prot.host[3])
colnames(prot.host)[3] <- "value"

afdw <- read.csv("output/0_biomass_output.csv")
afdw$variable <- colnames(afdw[3])
colnames(afdw)[3] <- "value"

cell.dens <- read.csv("output/0_sym_counts.csv")
cell.dens$variable <- colnames(cell.dens[3])
colnames(cell.dens)[3] <- "value"

# Coral sample metadata
metadata <- read_csv("metadata/coral_metadata.csv") 

Data <- rbind(prot.holo, prot.host, afdw, cell.dens,chla,chlc2)

Uni.Fig <- Data %>%
  ggplot(aes(x = site, y = value, color = site)) +
  labs(x = "Site", color = "Site") +
  facet_wrap(vars(variable), scales = "free_y") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Uni.Fig

data_new <- Data                              # Replicate data
data_new$group <- factor(data_new$variable,      # Reordering group factor levels
                         levels = c("AFDW.mg.cm2", "prot_ug.cm2", "avg_prot_ug.cm2",
                                    "cells.cm2","chla.ug.cm2","chlc2.ug.cm2"))

unique(Data$variable)

Uni.Fig <- data_new %>%
  ggplot(aes(x = site, y = value, color = site)) +
  labs(x = "Site", color = "Site") +
  facet_wrap(vars(group), scales = "free_y") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Uni.Fig


ggsave("output/TP0_Univariate_Figs.pdf", Uni.Fig, width=10, height=8)


#_________________________________________________________________________________________________________
## TO concatonate the output of all the ANOVA with corresponding TukeyHSD data analyses for each univariate metric
#Table of UNIVARIATE RESPONSEs ANOVA and HSD Results
#ORDER to Cat Scripts in: 
  #1) Chl a
  #2) chl c 
  #3) Sym Counts
  #4) Host Protein
  #5) Holobiont Protein
  #6) AFDW

# Creating the Title of the New File
cat("Table_TP0_Univariates.vs.Site_ANOVA_HSD\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt")

# add ANOVA in CHL A data
cat("A) ANOVA results of Chl a (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chla_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in CHL A data
cat("A) TukeyHSD results of Chl a (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chla_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#add ANOVA in CHLC data
cat("B) ANOVA results of Chl c (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chlc2_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#add TukeyHSD in CHLC data
cat("B) TukeyHSD results of Chl c (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chlc2_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


# add SYM COUNTS
cat("C) ANOVA results of Sym.Counts (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym.counts_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


# add TukeyHSD in SYM COUNTS
cat("C) TukeyHSD results of Sym.Counts (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym.counts_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add ANOVA HOST PROTEIN
cat("D) ANOVA results of Host Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in HOST PROTEIN
cat("D) TukeyHSD results of Host Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add HOLOBIONT PROTEIN
cat("E) ANOVA results of Holobiont Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(holo.prot_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in HOLOBIONT PROTEIN
cat("E) TukeyHSD results of Holobiont Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(holo.prot_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add ANOVA in AFDW
cat("F) ANOVA results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in AFDW
cat("F) TukeyHSD results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#_________________________________________________________________________________________________________

#Table of UNIVARIATE RESPONSEs (TUKEYHSD Results)

# Creating the Title of the New File
cat("Table_TP0_Univariates.vs.Site_TUKEYHSD\n\n", file = "output/Table_TP0_Univariates.vs.Site_TUKEYHSD.txt")

# add TukeyHSD in CHL A data
cat("A) TukeyHSD results of Chl a (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chla_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#add TukeyHSD in CHLC data
cat("B) TukeyHSD results of Chl c (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(chlc2_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in SYM COUNTS
cat("C) TukeyHSD results of Sym.Counts (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym.counts_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


# add TukeyHSD in HOST PROTEIN
cat("D) TukeyHSD results of Host Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in HOLOBIONT PROTEIN
cat("E) TukeyHSD results of Holobiont Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(holo.prot_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# add TukeyHSD in AFDW
cat("F) TukeyHSD results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host.prot_tkyhsd_res, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#______________________________________________________________________________________________
## METADATA FILE CREATION
#Creating Metadata of all metrics from TP0 for the 4 genotypes going to be compared with the rest of the timeseries

#read in all files to full join
path.p <- "output" #the location of all your respirometry files 
file.names <- list.files(path = path.p, pattern = "*4geno.csv")  # list all csv file names in the folder

#Couldn't get the more efficient coding to work well...so I manually did it
#df <- path.p %>%
  #list.files(path = ., pattern = "*4geno.csv") %>%
  #lapply(read_csv) %>% 
  #bind_rows 

#df <- tibble(file.name = file.names) %>%
  #mutate(metric = gsub("_.*", "", file.name),
         #info = map(metric, ~filter(metadata, colony_id == .)),           # Get associated sample info
         #data0 = map(file.name, ~read_csv(file.path(path.p, .), skip = 1))) %>%
  #lapply(read_csv) %>% 
 #bind_rows 

#manual way to merge it
biomass <- read.csv("output/0_biomass_4geno.csv")
chla <- read.csv("output/0_chlorophylla_4geno.csv")
chlc <- read.csv("output/0_chlorophyllc_4geno.csv")
holo_prot <- read.csv("output/0_holobiont_protein_4geno.csv")
host_prot <- read.csv("output/0_host_protein_4geno.csv")
sym_counts <- read.csv("output/0_sym_counts_4geno.csv")

#All the individual merges
geno <- full_join(biomass, chla)

geno <- full_join(geno, chlc)

geno <- full_join(geno, holo_prot)

geno <- left_join(geno, host_prot) #absorbing the values from just holo_prot since they have the same column name

geno_metadata <- full_join(geno, sym_counts)


#Converting Apul to ACR to match the timeseries
geno_metadata [1, 1] <- "ACR-237"
geno_metadata [2, 1] <- "ACR-234"
geno_metadata [3, 1] <- "ACR-243"
geno_metadata [4, 1] <- "ACR-244"

#Edit the current metadata sheet to match the current timeseries metadata
geno_metadata <- geno_metadata %>%
  mutate(Genotype = genotype) %>%
  mutate(month = "19-Oct") %>%
  mutate(nutrient = "NA") %>%
  mutate(Am = "NA") %>%
  mutate(AQY = "NA") %>%
  mutate(Rd = "NA") %>%
  mutate(site_code = "Nursery") %>%
  mutate(chla.ug.cell = geno_metadata$chla.ug.cm2 / geno_metadata$cells.cm2) %>%
  mutate(chlc2.ug.cell = geno_metadata$chlc2.ug.cm2/ geno_metadata$cells.cm2) %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, AFDW.mg.cm2, chla.ug.cm2, chlc2.ug.cm2, Am, AQY, Rd, cells.cm2, chla.ug.cell, chlc2.ug.cell)



###Compiled data from Jan/Nov has host and sym AFDW separated while we have it as one big thing, it is also missing protein.
##Need to create new column in data that is Holobiont AFDW
# colony_id, Genotype, timepoint, month, nutrient, site_code
timeseries_data <- read.csv("data/data_jan_nov_SA.csv")

timeseries_data <- timeseries_data %>%
  mutate(AFDW.mg.cm2 = timeseries_data$Host_AFDW.mg.cm2 + timeseries_data$Sym_AFDW.mg.cm2) %>%
  mutate(colony_id = ï..colony_id) %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, AFDW.mg.cm2, chla.ug.cm2, chlc2.ug.cm2, Am, AQY, Rd, cells.cm2, chla.ug.cell, chlc2.ug.cell)

Apul_Plast_Metadata <- Apul_Plast_Metadata %>%
  rbind(geno_metadata, timeseries_data) %>%
  distinct() %>% #removing duplicate data
  write_csv(., path = "data/complete_timeseries_data.csv")


####UNIVARIATE Figures and Analysis ---- NOT WORKING
library(RColorBrewer)
#Plot Biomass

#Subset full data set for just biomass data
timeseries_biomass <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, AFDW.mg.cm2)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_biomass$site_code <- factor(timeseries_biomass$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) #%>%
timeseries_biomass$Genotype <- factor(timeseries_biomass$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

#Make Boxplot - not all of this works so far
ggplot(timeseries_biomass, aes(x = site_code, y = AFDW.mg.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(aes(shape = Genotype), size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Month") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Holobiont Biomass (mg/cm2)"))))+
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=14),
    legend.text=element_text(size=14),
    axis.title=element_text(face="bold", size=14),
    axis.text=element_text(size=10, color="black"), 
    strip.text.x=element_text(face="italic", size=14)
  )

ggplot(Apul_Plast_Metadata, aes(x = site_code, y = AFDW.mg.cm2, fill = timepoint)) +
  geom_boxplot() +
  geom_point(aes(shape = Genotype), size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"), limits = c("Nursery", "Mahana \nLow", "Hilton \nMedium", "Manava \High"))
