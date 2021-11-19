#Summary of TP0 Univariate Response Variables (Except AQY, Am, and R)
library(tidyverse)
library(ggpubr)
library(showtext)
library(dplyr)
library(multcompView)
library(RColorBrewer)


# Analysis of Timepoint 0 samples from Nursery, naturally site 1 and site 2
####NURSERY VS WILD SITES TP0 ____________________________________________________________
#reading in all data files for TP0 (not filtered for genotypes)
biomass_TP0 <- read.csv("output/0_biomass.csv")
biomass_TP0 <- biomass_TP0 %>% 
  mutate_at("colony_id", str_replace, "_", "-") %>%
  mutate_at("colony_id", str_replace, "Apul", "ACR")

chla_TP0 <- read.csv("output/0_chlorophylla.csv")
chla_TP0 <- chla_TP0 %>% 
  mutate_at("colony_id", str_replace, "_", "-") %>%
  mutate_at("colony_id", str_replace, "Apul", "ACR")

chlc_TP0 <- read.csv("output/0_chlorophyllc.csv")
chlc_TP0 <- chlc_TP0 %>% 
  mutate_at("colony_id", str_replace, "_", "-") %>%
  mutate_at("colony_id", str_replace, "Apul", "ACR")

holo_prot_TP0 <- read.csv("output/0_holobiont_protein.csv")
holo_prot_TP0 <- holo_prot_TP0 %>% 
  mutate_at("colony_id", str_replace, "_", "-") %>%
  mutate_at("colony_id", str_replace, "Apul", "ACR")

host_prot_TP0 <- read.csv("output/0_host_protein.csv")
host_prot_TP0 <- host_prot_TP0 %>% 
  mutate_at("colony_id", str_replace, "_", "-") %>%
  mutate_at("colony_id", str_replace, "Apul", "ACR")

sym_counts_TP0 <- read.csv("output/0_sym_counts.csv")
sym_counts_TP0 <- sym_counts_TP0 %>% 
  mutate_at("colony_id", str_replace, "_", "-") %>%
  mutate_at("colony_id", str_replace, "Apul", "ACR")

pi_curve_TP0 <- read.csv("output/0_pi_curve_pars_nls.csv")
pi_curve_TP0 <- pi_curve_TP0 %>%
  mutate_at("colony_id", str_replace, "Apul", "ACR")

#need to do for AQY, Am and Rd (final of all metrics)

#All the individual merges of original TP0 files
TP0 <- full_join(biomass_TP0, chla_TP0)

TP0 <- full_join(TP0, chlc_TP0)

TP0 <- full_join(TP0, holo_prot_TP0)

TP0 <- full_join(TP0, host_prot_TP0) #absorbing the values from just holo_prot since they have the same column name

TP0 <- full_join(TP0, pi_curve_TP0)

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
  select(colony_id, site, timepoint, AFDW.mg.cm2, host_prot_ug.cm2, total_chl.ug.cm2, total_chl.ug.cell, Am, AQY, Rd, cells.cm2) #%>% #selecting the six metrics wanting to output

TP0_metadata %>% #writing output of metadata for TP0 (all nursery - n=10 and wild sites n=3 per site)
  write.csv("output/TP0_metadata")

#STATS on all

# Tukey test to study each pair of treatment :
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

#Biomass__________________________________________
model.AFDW <- aov(log10(AFDW.mg.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.AFDW$residuals)
qqnorm(model.AFDW$residuals)
qqline(model.AFDW$residuals)
plot(model.AFDW$fitted.values, model.AFDW$residuals)
biomass_anova<-summary(model.AFDW) #anova for biomass
biomass_tkyhsd<-TukeyHSD(model.AFDW) #posthoc tests within anova for biomass
AFDW.HSD.Letters<-generate_label_df(biomass_tkyhsd , "site")#generate Tukey HSD letters using function
AFDW.HSD.Letters$name <- "AFDW.mg.cm2" 
AFDW.HSD.Letters$site <- rownames(AFDW.HSD.Letters)

#Host Protein ________________________________
model.HostProt <- aov(log10(host_prot_ug.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.HostProt$residuals)
qqnorm(model.HostProt$residuals)
qqline(model.HostProt$residuals)
plot(model.HostProt$fitted.values, model.HostProt$residuals)
host_prot_anova<-anova(model.HostProt) #anova for biomass
host_prot_tkyhsd<-TukeyHSD(model.HostProt) #posthoc tests within anova for biomass
host_prot.HSD.Letters<-generate_label_df(host_prot_tkyhsd , "site")#generate Tukey HSD letters using function
host_prot.HSD.Letters$name <- "host_prot_ug.cm2"
host_prot.HSD.Letters$site <- rownames(host_prot.HSD.Letters)

#Sym Density______________________________________
model.CellDens <- aov(log10(cells.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.CellDens$residuals)
qqnorm(model.CellDens$residuals)
qqline(model.CellDens$residuals)
plot(model.CellDens$fitted.values, model.CellDens$residuals)
sym_dens_anova<-anova(model.CellDens) #anova for Sym Density
sym_dens_tkyhsd<-TukeyHSD(model.CellDens) #posthoc tests within anova for Sym Density
sym_dens.HSD.Letters<-generate_label_df(sym_dens_tkyhsd , "site")#generate Tukey HSD letters using function
sym_dens.HSD.Letters$name <- "cells.cm2"
sym_dens.HSD.Letters$site <- rownames(sym_dens.HSD.Letters)

#Total Chlorophyll per Frag (ug/cm2)_______________
model.Chlcm <- aov(log10(total_chl.ug.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.Chlcm$residuals)
qqnorm(model.Chlcm$residuals)
qqline(model.Chlcm$residuals)
plot(model.Chlcm$fitted.values, model.Chlcm$residuals)
tot_chl_anova<-anova(model.Chlcm) #anova for Total Chl
tot_chl_tkyhsd<-TukeyHSD(model.Chlcm) #posthoc tests within anova for Total Chl
tot_chl.HSD.Letters<-generate_label_df(tot_chl_tkyhsd , "site")#generate Tukey HSD letters using function
tot_chl.HSD.Letters$name <- "total_chl.ug.cm2" 
tot_chl.HSD.Letters$site <- rownames(tot_chl.HSD.Letters)

#Total Chlorophyll per cell (ug/cells)____________
model.Chl.cell <- aov(log10(total_chl.ug.cell) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.Chl.cell$residuals)
qqnorm(model.Chl.cell$residuals)
qqline(model.Chl.cell$residuals)
plot(model.Chl.cell$fitted.values, model.Chl.cell$residuals)
tot_chl_per.cell_anova<-anova(model.Chl.cell) #anova for Total Chl per cell
tot_chl_per.cell_tkyhsd<-TukeyHSD(model.Chl.cell) #posthoc tests within anova for Total Chl per cell
tot_chl_per.cell.HSD.Letters<-generate_label_df(tot_chl_per.cell_tkyhsd , "site")#generate Tukey HSD letters using function
tot_chl_per.cell.HSD.Letters$name <- "total_chl.ug.cell" 
tot_chl_per.cell.HSD.Letters$site <- rownames(tot_chl_per.cell.HSD.Letters)

#Am (Max Photosynthesis)____________
model.Am <- aov(log10(Am) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.Am$residuals)
qqnorm(model.Am$residuals)
qqline(model.Am$residuals)
plot(model.Am$fitted.values, model.Chl.cell$residuals)
Am_anova<-anova(model.Am) #anova for Am
Am_tkyhsd<-TukeyHSD(model.Am) #posthoc tests within anova for Am
Am.HSD.Letters<-generate_label_df(Am_tkyhsd , "site")#generate Tukey HSD letters using function
Am.HSD.Letters$name <-"Am" 
Am.HSD.Letters$site <- rownames(Am.HSD.Letters)

#AQY (Max Photosynthetic Rate)____________
model.AQY <- aov(log10(AQY) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.AQY$residuals)
qqnorm(model.AQY$residuals)
qqline(model.AQY$residuals)
plot(model.AQY$fitted.values, model.Chl.cell$residuals)
AQY_anova<-anova(model.AQY) #anova for AQY
AQY_tkyhsd<-TukeyHSD(model.AQY) #posthoc tests within anova for AQY
AQY.HSD.Letters<-generate_label_df(AQY_tkyhsd , "site")#generate Tukey HSD letters using function
AQY.HSD.Letters$name <- "AQY" 
AQY.HSD.Letters$site <- rownames(AQY.HSD.Letters)

#Rd (Respiration Rate)____________
model.Rd <- aov(log10(Rd) ~ site, data = TP0_metadata) #save model
par(mfrow=c(1,3))
hist(model.Rd$residuals)
qqnorm(model.Rd$residuals)
qqline(model.Rd$residuals)
plot(model.Rd$fitted.values, model.Chl.cell$residuals)
Rd_anova<-anova(model.Rd) #anova for Rd
Rd_tkyhsd<-TukeyHSD(model.Rd) #posthoc tests within anova for Rd
Rd.HSD.Letters<-generate_label_df(Rd_tkyhsd , "site")#generate Tukey HSD letters using function
Rd.HSD.Letters$name <- "Rd"
Rd.HSD.Letters$site <- rownames(Rd.HSD.Letters)


Letters <- rbind(AFDW.HSD.Letters,host_prot.HSD.Letters,sym_dens.HSD.Letters,tot_chl.HSD.Letters,
                 tot_chl_per.cell.HSD.Letters,Am.HSD.Letters, AQY.HSD.Letters,Rd.HSD.Letters)


Letters$group <- paste0(Letters$site, Letters$name)



TP0_metadata_long <- pivot_longer(TP0_metadata, cols=AFDW.mg.cm2:cells.cm2)
TP0_metadata_long$group <- paste0(TP0_metadata_long$site, TP0_metadata_long$name)

yvalue<- TP0_metadata_long %>% group_by(name,site) %>%
  summarise(mean = mean(value))# obtain letter position for y axis using means
yvalue$group <-paste0(yvalue$site,yvalue$name)

  
final<-left_join(yvalue,Letters, by="group", keep=F) #merge dataframes
final<-left_join(final,TP0_metadata_long, by="group") #merge dataframes

unique(final$name)

labs <- as_labeller(c("AFDW.mg.cm2", "Am","AQY","cells.cm2","host_prot_ug.cm2","Rd","total_chl.ug.cell","total_chl.ug.cm2") )

TP0Fig <- final %>%
  ggplot(aes(x = site.x, y = value)) +
  geom_blank() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(vars(name), scales = "free_y", nrow = 2,ncol = 4)+
  geom_jitter(width = 0.1, size=0.5) + # Plot all points
  theme(plot.title = element_text(hjust = 0.5, face='bold'))+
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.1) +
  stat_summary(fun = mean, geom = "point", color = "black") +           # Plot mean
  geom_text(aes(x = site, y = mean, label = Letters),vjust=0,hjust=-1.5) +
  theme(plot.title = element_text(vjust=-0.6))

ggsave("output/TP0_Univariate_Figs.jpg", TP0Fig, width=12, height=8)



#ANOVA on all variables and TukeyHSD Posthoc tests

#Biomass__________________________________________
# Concatenate ANOVA results in txt file (AFDW)
cat("A) ANOVA results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(biomass_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# Concatenate TukeyHSD results in txt file (AFDW)
cat("A) TukeyHSD results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(biomass_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Host Protein Stats________________________________
## Concatenate ANOVA results in txt file (Host Protein)
cat("B) ANOVA results of Host Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host_prot_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Host Protein)
cat("B) TukeyHSD results of Host Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host_prot_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Sym Density______________________________________
## Concatenate ANOVA results in txt file (Sym Density)
cat("C) ANOVA results of Symbiont Density (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym_dens_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Sym Density)
cat("C) TukeyHSD results of Symbiont Density (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym_dens_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


## Concatenate ANOVA results in txt file (Total CHL per Frag)
cat("D) ANOVA results of Total Chlorophyll per Fragment (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Total CHL per Frag)
cat("D) TukeyHSD results of Total Chlorophyll per Fragment (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Total Chlorophyll per cell (ug/cells)____________
## Concatenate ANOVA results in txt file (Total CHL per cell)
cat("E) ANOVA results of Total Chlorophyll per Symbiont (ug/cell) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_per.cell_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Total CHL per cell)
cat("E) TukeyHSD results of Total Chlorophyll per Symbiont (ug/cell) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_per.cell_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Am (Max Photosynthesis)____________
## Concatenate ANOVA results in txt file (Max Photosynthesis)
cat("F) ANOVA results of Max Photosynthesis Levels at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(Am_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Max Photosynthesis)
cat("F) TukeyHSD results of Max Photosynthesis Levels at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(Am_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#AQY (Max Photosynthetic Rate)____________
## Concatenate ANOVA results in txt file (Max Photosynthetic Rate)
cat("G) ANOVA results of Max Photosynthetic Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(AQY_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Max Photosynthetic Rate)
cat("G) TukeyHSD results of Max Photosynthetic Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(AQY_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Rd (Respiration Rate)____________
## Concatenate ANOVA results in txt file (Respiration Rate)
cat("H) ANOVA results of Respiration Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(Rd_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Respiration Rate)
cat("H) TukeyHSD results of Respiration Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(Rd_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)





