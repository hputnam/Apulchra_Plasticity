#Summary of TP0 Univariate Response Variables (Except AQY, Am, and R)
library(tidyverse)
library(ggpubr)
library(showtext)
library(dplyr)

#Set working directory
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")

####NURSERY VS WILD SITES TP0 ____________________________________________________________
#reading in all data files for TP0 (not filtered for genotypes)
biomass_TP0 <- read.csv("output/0_biomass.csv")
chla_TP0 <- read.csv("output/0_chlorophylla.csv")
chlc_TP0 <- read.csv("output/0_chlorophyllc.csv")
holo_prot_TP0 <- read.csv("output/0_holobiont_protein.csv")
host_prot_TP0 <- read.csv("output/0_host_protein.csv")
sym_counts_TP0 <- read.csv("output/0_sym_counts.csv")

#need to do for AQY, Am and Rd (final of all metrics)

#All the individual merges of original TP0 files
TP0 <- full_join(biomass_TP0, chla_TP0)

TP0 <- full_join(TP0, chlc_TP0)

TP0 <- full_join(TP0, holo_prot_TP0)

TP0 <- left_join(TP0, host_prot_TP0) #absorbing the values from just holo_prot since they have the same column name

TP0_metadata <- full_join(TP0, sym_counts_TP0)

#Creating Chl A and C columns per cell as well as Total CHL columns
TP0_metadata <- TP0_metadata %>%
  mutate(chla.ug.cell = TP0_metadata$chla.ug.cm2 / TP0_metadata$cells.cm2)  #chla per symbiont cell

TP0_metadata <- TP0_metadata %>%
  mutate(chlc2.ug.cell = TP0_metadata$chlc2.ug.cm2/ TP0_metadata$cells.cm2) #chlc per symbiont cell

TP0_metadata <- TP0_metadata %>%
  mutate(total_chl.ug.cm2 = TP0_metadata$chla.ug.cm2 + TP0_metadata$chlc2.ug.cm2) #total chlorophyll per cm2 of frag

TP0_metadata <- TP0_metadata %>%
  mutate(total_chl.ug.cell = TP0_metadata$chla.ug.cell + TP0_metadata$chlc2.ug.cell)  #total chlorophyll per sym cell 

TP0_metadata <- TP0_metadata %>%
  mutate(holo_prot_ug.cm2 = prot_ug.cm2) #renaming to holobiont prot column

TP0_metadata <- TP0_metadata %>%
  mutate(host_prot_ug.cm2 = avg_prot_ug.cm2)  #renaming to host prot column

TP0_metadata <- TP0_metadata %>%
  mutate(sym_prot_ug.cm2 = TP0_metadata$holo_prot_ug.cm2 - TP0_metadata$host_prot_ug.cm2) #calculating sym prot

TP0_metadata <- TP0_metadata %>%
  select(colony_id, site, timepoint, AFDW.mg.cm2, host_prot_ug.cm2, sym_prot_ug.cm2, total_chl.ug.cm2, total_chl.ug.cell, cells.cm2) #%>% #selecting the six metrics wanting to output
  
TP0_metadata <- TP0_metadata %>% #writing output of metadata for TP0 (all nursery - n=10 and wild sites n=3 per site)
  write.csv("output/TP0_metadata")
  
#Individually create plots and save them before compiling into one large output plot
#biomass fig
AFDW_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = AFDW.mg.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Ash Free Dry Weight (mg/cm2)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
AFDW_TP0

#Host Protein fig
Host_Prot_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = host_prot_ug.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Host Protein (ug/cm2)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Host_Prot_TP0  

#Sym Protein fig
Sym_Prot_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = sym_prot_ug.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Symbiont Protein (ug/cm2)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Sym_Prot_TP0

#Total CHL per cm2 fig
Tot_CHL.cm2_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = total_chl.ug.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Total Chlorophyll per Fragment (ug/cm2)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Tot_CHL.cm2_TP0

#Total CHL per Symbiont fig
Tot_CHL.cell_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = total_chl.ug.cell, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Total Chlorophyll per Symbiont (ug/cell)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Tot_CHL.cell_TP0

#Sym Density fig
sym_dens_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = cells.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle("Symbiont Density (cells/cm2)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
sym_dens_TP0


#Configuring all the saved figures for TP0 into final stage of itself
TP0Fig <- ggarrange(AFDW_TP0, Host_Prot_TP0,Sym_Prot_TP0, sym_dens_TP0,Tot_CHL.cell_TP0,Tot_CHL.cm2_TP0, ncol = 3, nrow = 2)

ggsave("output/TP0_Univariate_Figs.pdf", TP0Fig, width=12, height=12)

#STATS on all
model1 <- aov(AFDW.mg.cm2 ~ site, TP0_metadata) #save model for biomass
model2 <- aov(host_prot_ug.cm2 ~ site, TP0_metadata) #save model for Host Prot
model3 <- aov(sym_prot_ug.cm2 ~ site, TP0_metadata) #save model for Sym Prot
model4 <- aov(cells.cm2 ~ site, TP0_metadata) #save model for Sym Density
model5 <- aov(total_chl.ug.cm2 ~ site, TP0_metadata) #save model for Total Chl/cell
model6 <- aov(total_chl.ug.cell ~ site, TP0_metadata) #save model for Total Chl/cm2

#ANOVA on all of them with TukeyHSD Posthoc tests

#______________________________________________________________________________________________
#CREATING TIMESERIES METADATA AND FIGS
## METADATA FILE CREATION
#Creating Metadata of all metrics from TP0 for the 4 genotypes going to be compared with the rest of the timeseries


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
  mutate(AFDW.mg.cm2 = Host_AFDW.mg.cm2 + Sym_AFDW.mg.cm2) %>%
  #mutate(colony_id = ?..colony_id) %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, AFDW.mg.cm2, chla.ug.cm2, chlc2.ug.cm2, Am, AQY, Rd, cells.cm2, chla.ug.cell, chlc2.ug.cell)

Apul_Plast_Metadata <- rbind(geno_metadata, timeseries_data) #STILL NEED TO ADD Host and Holobiont Protein to it

Apul_Plast_Metadata$group <- paste0(Apul_Plast_Metadata$Genotype,Apul_Plast_Metadata$timepoint, Apul_Plast_Metadata$site_code)

Apul_Plast_Metadata <- Apul_Plast_Metadata %>%
  distinct(group, .keep_all = TRUE) %>% #removing duplicate data
  write_csv(path = "data/complete_timeseries_data.csv")


####UNIVARIATE Figures and Analysis ---- BIOMASS
library(RColorBrewer)
##Plot Biomass##

#Subset full data set for just biomass data
timeseries_biomass <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, AFDW.mg.cm2)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_biomass$site_code <- factor(timeseries_biomass$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) #%>%
timeseries_biomass$Genotype <- factor(timeseries_biomass$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

#Make Boxplot - BIOMASS
biomass_fig <- ggplot(timeseries_biomass, aes(x = site_code, y = AFDW.mg.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Holobiont Biomass (mg/cm2)"))))+
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )



####UNIVARIATE Figures and Analysis ---- CHL A (ug/cm2)
##Plot CHL A (ug/cm2)##

#Subset full data set for just biomass data
timeseries_chla_cm2 <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, chla.ug.cm2)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_chla_cm2$site_code <- factor(timeseries_chla_cm2$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_chla_cm2$Genotype <- factor(timeseries_chla_cm2$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

chla_ug.cm2_fig <- ggplot(timeseries_chla_cm2, aes(x = site_code, y = chla.ug.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Chlorophyll A per Frag (ug/cm2)")))) +
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  ylim(0, 4) +
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )

####UNIVARIATE Figures and Analysis ---- CHL C (ug/cm2)
##Plot CHL C (ug/cm2)##

#Subset full data set for just biomass data
timeseries_chlc2_cm2 <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, chlc2.ug.cm2)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_chlc2_cm2$site_code <- factor(timeseries_chlc2_cm2$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_chlc2_cm2$Genotype <- factor(timeseries_chlc2_cm2$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

chlc2_ug.cm2_fig <- ggplot(timeseries_chlc2_cm2, aes(x = site_code, y = chlc2.ug.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Chlorophyll C per Frag (ug/cm2)"))))+
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  ylim(0, 4) +
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )

####UNIVARIATE Figures and Analysis ---- CHL A (ug/cells)
##Plot CHL A (ug/cells)##

#Subset full data set for just biomass data
timeseries_chla_cells <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, chla.ug.cell)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_chla_cells$site_code <- factor(timeseries_chla_cells$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_chla_cells$Genotype <- factor(timeseries_chla_cells$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

library(scales)

chla_cells.cm2_fig <- ggplot(timeseries_chla_cells, aes(x = site_code, y = chla.ug.cell, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Chlorophyll A per Sym (ug/cells)")))) +
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  ylim(0.000000e-06, 6.500000e-06) +
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )


####UNIVARIATE Figures and Analysis ---- CHL C (ug/cells)
##Plot CHL C (ug/cells)##

#Subset full data set for just biomass data
timeseries_chlc2_cells <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, chlc2.ug.cell)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_chlc2_cells$site_code <- factor(timeseries_chlc2_cells$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_chlc2_cells$Genotype <- factor(timeseries_chlc2_cells$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

chlc2_cells.cm2_fig <- ggplot(timeseries_chlc2_cells, aes(x = site_code, y = chlc2.ug.cell, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Chlorophyll C per Sym (ug/cells)"))))+
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  ylim(0.000000e-06, 6.500000e-06) +
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )


####UNIVARIATE Figures and Analysis ---- SYM Density (cells/cm2)
##Plot CHL C (ug/cells)##

#Subset full data set for just biomass data
timeseries_sym_counts <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, cells.cm2)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_sym_counts$site_code <- factor(timeseries_sym_counts$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_sym_counts$Genotype <- factor(timeseries_sym_counts$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

sym_counts_fig <- ggplot(timeseries_sym_counts, aes(x = site_code, y = cells.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  scale_y_continuous(labels = scales::scientific) +
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Symbiont Density (cells/cm2)"))))+
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )


###Configure all univariate functions into one major figure

Fig <- ggarrange(biomass_fig, chla_ug.cm2_fig,chlc2_ug.cm2_fig, sym_counts_fig,chla_cells.cm2_fig,chlc2_cells.cm2_fig, ncol = 3, nrow = 2)
ggsave("Output/Univariate_Figs.pdf", Fig, width=16, height=8)

Fig

colnames(Apul_Plast_Metadata)

fixed_2way <- aov(AFDW.mg.cm2 ~site_code*timepoint, Apul_Plast_Metadata)
summary(fixed_2way)



library(lme4)
mixed_2way <- lmer(AFDW.mg.cm2 ~site_code*timepoint + (1 | Genotype), Apul_Plast_Metadata)
summary(mixed_2way)

