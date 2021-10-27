#Summary of TP0 Univariate Response Variables (Except AQY, Am, and R)
library(tidyverse)
library(ggpubr)
library(showtext)
library(dplyr)

#Set working directory
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")


#All saved Figures composed into one figure

#chl a per frag
chla <- read.csv("output/0_chlorophylla.csv")
chla$variable <- colnames(chla[3])
colnames(chla)[3] <- "value"

#chl a per cell - open and merge all the output data to get CHL a and c2 per cell
ChlorA <- read.csv("output/0_chlorophylla.csv")
ChlorC <- read.csv("output/0_chlorophyllc.csv")
Sym_Density <- read.csv("output/0_sym_counts.csv")


TP0_final_data <- full_join(ChlorA, ChlorC)
TP0_final_data <- full_join(TP0_final_data, Sym_Density)

Chl_per_cell <- TP0_final_data %>%
  mutate(chla.ug.cell = TP0_final_data$chla.ug.cm2 / TP0_final_data$cells.cm2) %>%
  mutate(chlc2.ug.cell = TP0_final_data$chlc2.ug.cm2/ TP0_final_data$cells.cm2) %>%
  write.csv("output/chl_per_cell_TP0.csv")

#pulling out only the values needed to get chl a per cell values
chla_ug.cell <- read.csv("output/chl_per_cell_TP0.csv")
chla_ug.cell$variable <- colnames(chla_ug.cell[8])
colnames(chla_ug.cell)[8] <- "value"  

chla_ug.cell <- select(chla_ug.cell, colony_id, site, value, timepoint, variable)

#chl c2 per cell values
chlc2_ug.cell <- read.csv("output/chl_per_cell_TP0.csv")
chlc2_ug.cell$variable <- colnames(chlc2_ug.cell[9])
colnames(chlc2_ug.cell)[9] <- "value"  

chlc2_ug.cell <- select(chlc2_ug.cell, colony_id, site, value, timepoint, variable)

#Chl C per fragment
chlc2 <- read.csv("output/0_chlorophyllc.csv")
chlc2$variable <- colnames(chlc2[3])
colnames(chlc2)[3] <- "value"

#Holobiont Protein
prot.holo <- read.csv("output/0_holobiont_protein.csv")
prot.holo$variable <- colnames(prot.holo[3])
colnames(prot.holo)[3] <- "value"

#Host Protein
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

Data <- rbind(prot.holo, prot.host, afdw, cell.dens, chla_ug.cell, chlc2_ug.cell, chla, chlc2) #add 


#Creating labels/titles for the facets
variable_titles <- c(
  'prot_ug.cm2' = "Holobiont Protein (ug/cm2)",
  'avg_prot_ug.cm2' = "Host Protein (ug/cm2)",
  'AFDW.mg.cm2' = "Ash Free Dry Weight (mg/cm2)",
  'cells.cm2' = "Symbiont Density (cells/cm2)",
  'chla.ug.cell' = "Chl A per Symb. (ug/cells)",
  'chlc2.ug.cell' = "Chl C per Symb. (ug/cm2)",
  'chla.ug.cm2' = "Chl A per Frag (ug/cm2)",
  'chlc2.ug.cm2' = "Chl C per Frag (ug/cm2)"
)
names(variable_titles)

#Try running in Fig facet_wrap() at the end to try to get to work
labeller = labeller(variable = variable_titles) #to insert into figure provided it actually works



#Faceted figure for all univariate responses at TP0
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
                         levels = c("AFDW.mg.cm2", "prot_ug.cm2", "avg_prot_ug.cm2","cells.cm2", 
                                    "chla.ug.cell", "chlc2.ug.cell", "chla.ug.cm2","chlc2.ug.cm2"))  

unique(Data$variable)

Uni.Fig2 <- data_new %>%
  ggplot(aes(x = site, y = value, color = site)) +
  labs(x = "Site", color = "Site") +
  facet_wrap(vars(group), scales = "free_y") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Uni.Fig2


ggsave("output/TP0_Univariate_Figs.pdf", Uni.Fig2, width=12, height=12)

#______________________________________________________________________________________________
## METADATA FILE CREATION
#Creating Metadata of all metrics from TP0 for the 4 genotypes going to be compared with the rest of the timeseries

#read in all files to full join
#path.p <- "output" #the location of all your respirometry files 
#file.names <- list.files(path = path.p, pattern = "*4geno.csv")  # list all csv file names in the folder

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

