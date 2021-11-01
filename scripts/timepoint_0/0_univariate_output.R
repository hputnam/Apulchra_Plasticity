#Summary of TP0 Univariate Response Variables (Except AQY, Am, and R)
library(tidyverse)
library(ggpubr)
library(showtext)
library(dplyr)

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
  select(colony_id, site, timepoint, AFDW.mg.cm2, host_prot_ug.cm2, sym_prot_ug.cm2, total_chl.ug.cm2, total_chl.ug.cell, Am, AQY, Rd, cells.cm2) #%>% #selecting the six metrics wanting to output
  
TP0_metadata <- TP0_metadata %>% #writing output of metadata for TP0 (all nursery - n=10 and wild sites n=3 per site)
  write.csv("output/TP0_metadata")
  
#Individually create plots and save them before compiling into one large output plot
#biomass fig

#Biomass Fig
AFDW_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = AFDW.mg.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ylim(0, 3.5) +
  ggtitle("Ash Free Dry Weight (mg/cm2)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) + # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
AFDW_TP0

#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation <- data.frame(
  x = c(1, 2, 3),
  y = c(1.75, 3, 2.25),
  label = c("B", "AB", "B")
)

#add text to figure
AFDW_TP0 <- AFDW_TP0 + geom_text(data = annotation, aes (x=x, y=y, label=label),
                     color = "black",
                     size = 5, fontface= "bold")
AFDW_TP0

#Host Protein fig
Host_Prot_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = host_prot_ug.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Host Protein (ug/cm2)") +
  ylim(0, 400) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) + # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Host_Prot_TP0  

#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation1 <- data.frame(
  x = c(1, 2, 3),
  y = c(200, 350, 300),
  label = c("A", "B", "B")
)

#add text to figure
Host_Prot_TP0 <- Host_Prot_TP0 + geom_text(data = annotation1, aes (x=x, y=y, label=label),
                                 color = "black",
                                 size = 5, fontface= "bold")
Host_Prot_TP0


#Sym Protein fig - NOT USING IN THE END
#Sym_Prot_TP0 <- TP0_metadata %>%
  #ggplot(aes(x = site, y = sym_prot_ug.cm2, color = site)) +
  #labs(x = "Site", y = "", color = "Site") +
  #ggtitle("Symbiont Protein (ug/cm2)") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #geom_jitter(width = 0.1) +                                            # Plot all points
  #stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               #geom = "errorbar", color = "black", width = 0.5) +
  #stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
#Sym_Prot_TP0


#Total CHL per cm2 fig
Tot_CHL.cm2_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = total_chl.ug.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Total Chlorophyll per Fragment (ug/cm2)") +
  ylim(0, 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Tot_CHL.cm2_TP0

#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation2 <- data.frame(
  x = c(1, 2, 3),
  y = c(5, 11, 5),
  label = c("A", "B", "B")
)

#add text to figure
Tot_CHL.cm2_TP0 <- Tot_CHL.cm2_TP0 + geom_text(data = annotation2, aes (x=x, y=y, label=label),
                                           color = "black",
                                           size = 5, fontface= "bold")
Tot_CHL.cm2_TP0

#Total CHL per Symbiont fig
Tot_CHL.cell_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = total_chl.ug.cell, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ylim(0, 2.5e-05) +
  ggtitle("Total Chlorophyll per Symbiont (ug/cell)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean
Tot_CHL.cell_TP0

#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation3 <- data.frame(
  x = c(1, 2, 3),
  y = c(1.0e-05, 2.0e-05, 1.0e-05),
  label = c("B", "AB", "B")
)

#add text to figure
Tot_CHL.cell_TP0 <- Tot_CHL.cell_TP0 + geom_text(data = annotation3, aes (x=x, y=y, label=label),
                                               color = "black",
                                               size = 5, fontface= "bold")
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
#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation4 <- data.frame(
  x = c(1, 2, 3),
  y = c(7.50e+05, 1.25e+06, 9.00e+05),
  label = c("B", "AB", "B")
)

#add text to figure
sym_dens_TP0 <- sym_dens_TP0 + geom_text(data = annotation4, aes (x=x, y=y, label=label),
                                                 color = "black",
                                                 size = 5, fontface= "bold")
sym_dens_TP0

#Am (Max Photosynthesis) fig
Am_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = Am, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Max Photosynthesis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

Am_TP0
#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation5 <- data.frame(
  x = c(1, 2, 3),
  y = c(1, 2.1, 2),
  label = c("A", "B", "B")
)

#add text to figure
Am_TP0 <- Am_TP0 + geom_text(data = annotation5, aes (x=x, y=y, label=label),
                                         color = "black",
                                         size = 5, fontface= "bold")
Am_TP0

#AQY (Max Photosynthetic Rate) fig
AQY_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = AQY, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Max Photosynthetic Rate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

AQY_TP0
#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation6 <- data.frame(
  x = c(1, 2, 3),
  y = c(0.004, 0.004, 0.004),
  label = c("A", "A", "A")
)

#add text to figure
AQY_TP0 <- AQY_TP0 + geom_text(data = annotation6, aes (x=x, y=y, label=label),
                             color = "black",
                             size = 5, fontface= "bold")
AQY_TP0


#Rd (Respiration Rate) Fig
Rd_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = Rd, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ggtitle("Respiration Rate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black")           # Plot mean

Rd_TP0

#adding significance to indicators manually to plot based on TukeyHSD Posthoc
annotation7 <- data.frame(
  x = c(1, 2, 3),
  y = c(0.4, 0.7, 0.6),
  label = c("A", "B", "B")
)

#add text to figure
Rd_TP0 <- Rd_TP0 + geom_text(data = annotation7, aes (x=x, y=y, label=label),
                               color = "black",
                               size = 5, fontface= "bold")
Rd_TP0


#Configuring all the saved figures for TP0 into final stage of itself
TP0Fig <- ggarrange(AFDW_TP0, Host_Prot_TP0, Rd_TP0, sym_dens_TP0, Am_TP0, AQY_TP0,Tot_CHL.cm2_TP0, Tot_CHL.cell_TP0, ncol = 4, nrow = 2)

ggsave("output/TP0_Univariate_Figs.pdf", TP0Fig, width=16, height=12)

#STATS on all
model1 <- aov(log10(AFDW.mg.cm2) ~ site, data = TP0_metadata) #save model for biomass
model2 <- aov(log10(host_prot_ug.cm2) ~ site, data = TP0_metadata) #save model for Host Prot
model3 <- aov(log10(cells.cm2) ~ site, data = TP0_metadata) #save model for Sym Density
model4 <- aov(log10(total_chl.ug.cm2) ~ site, data = TP0_metadata) #save model for Total Chl/cell
model5 <- aov(log10(total_chl.ug.cell) ~ site, data = TP0_metadata) #save model for Total Chl/cm2
model6 <- aov(log10(Am) ~ site, data = TP0_metadata) #save model for Total Chl/cm2
model7 <- aov(log10(AQY) ~ site, data = TP0_metadata) #save model for Total Chl/cm2
model8 <- aov(log10(Rd) ~ site, data = TP0_metadata) #save model for Total Chl/cm2


#ANOVA on all of them with TukeyHSD Posthoc tests

#Biomass__________________________________________
biomass_anova<-anova(model1) #anova for biomass
biomass_tkyhsd<-TukeyHSD(model1) #posthoc tests within anova for biomass

# Concatenate ANOVA results in txt file (AFDW)
cat("A) ANOVA results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(biomass_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

# Concatenate TukeyHSD results in txt file (AFDW)
cat("A) TukeyHSD results of AFDW (mg/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(biomass_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Host Protein Stats________________________________
host_prot_anova<-anova(model2) #anova for biomass
host_prot_tkyhsd<-TukeyHSD(model2) #posthoc tests within anova for biomass

## Concatenate ANOVA results in txt file (Host Protein)
cat("B) ANOVA results of Host Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host_prot_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Host Protein)
cat("B) TukeyHSD results of Host Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(host_prot_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Sym Protein______________________________________
#Did not end up using since only host prot run on timeseries data down the line
#sym_prot_anova<-anova(model3) #anova for biomass
#sym_prot_tkyhsd<-TukeyHSD(model3) #posthoc tests within anova for biomass

## Concatenate ANOVA results in txt file (Sym Protein)
#cat("C) ANOVA results of Symbiont Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
#capture.output(sym_prot_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
#cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Sym Protein)
#cat("C) TukeyHSD results of Symbiont Protein (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
#capture.output(sym_prot_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
#cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#Sym Density______________________________________
sym_dens_anova<-anova(model4) #anova for Sym Density
sym_dens_tkyhsd<-TukeyHSD(model1) #posthoc tests within anova for Sym Density

## Concatenate ANOVA results in txt file (Sym Density)
cat("C) ANOVA results of Symbiont Density (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym_dens_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Sym Density)
cat("C) TukeyHSD results of Symbiont Density (cells/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(sym_dens_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#Total Chlorophyll per Frag (ug/cm2)_______________
tot_chl_anova<-anova(model5) #anova for Total Chl
tot_chl_tkyhsd<-TukeyHSD(model5) #posthoc tests within anova for Total Chl

## Concatenate ANOVA results in txt file (Total CHL per Frag)
cat("D) ANOVA results of Total Chlorophyll per Fragment (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Total CHL per Frag)
cat("D) TukeyHSD results of Total Chlorophyll per Fragment (ug/cm2) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#Total Chlorophyll per cell (ug/cells)____________
tot_chl_per.cell_anova<-anova(model6) #anova for Total Chl per cell
tot_chl_per.cell_tkyhsd<-TukeyHSD(model6) #posthoc tests within anova for Total Chl per cell

## Concatenate ANOVA results in txt file (Total CHL per cell)
cat("E) ANOVA results of Total Chlorophyll per Symbiont (ug/cell) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_per.cell_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Total CHL per cell)
cat("E) TukeyHSD results of Total Chlorophyll per Symbiont (ug/cell) at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(tot_chl_per.cell_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#Am (Max Photosynthesis)____________
Am_anova<-anova(model6) #anova for Am
Am_tkyhsd<-TukeyHSD(model6) #posthoc tests within anova for Am

## Concatenate ANOVA results in txt file (Max Photosynthesis)
cat("F) ANOVA results of Max Photosynthesis Levels at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(Am_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Max Photosynthesis)
cat("F) TukeyHSD results of Max Photosynthesis Levels at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(Am_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)


#AQY (Max Photosynthetic Rate)____________
AQY_anova<-anova(model7) #anova for AQY
AQY_tkyhsd<-TukeyHSD(model7) #posthoc tests within anova for AQY

## Concatenate ANOVA results in txt file (Max Photosynthetic Rate)
cat("G) ANOVA results of Max Photosynthetic Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(AQY_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Max Photosynthetic Rate)
cat("G) TukeyHSD results of Max Photosynthetic Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(AQY_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

#Rd (Respiration Rate)____________
Rd_anova<-anova(model8) #anova for Rd
Rd_tkyhsd<-TukeyHSD(model8) #posthoc tests within anova for Rd

## Concatenate ANOVA results in txt file (Respiration Rate)
cat("H) ANOVA results of Respiration Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(AQY_anova, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)

## Concatenate TukeyHSD results in txt file (Respiration Rate)
cat("H) TukeyHSD results of Respiration Rate at TP0 sites\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
capture.output(AQY_tkyhsd, file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)
cat("\n\n", file = "output/Table_TP0_Univariates.vs.Site_ANOVA_HSD.txt", append = TRUE)



#______________________________________________________________________________________________
#CREATING TIMESERIES METADATA AND FIGS
## METADATA FILE CREATION
#Creating Metadata of all metrics from TP0 for the 4 genotypes going to be compared with the rest of the timeseries


#manual way to merge it
biomass <- read.csv("output/0_biomass_4geno.csv")
biomass <- biomass %>% 
  mutate_at("colony_id", str_replace, "_", "-")
chla <- read.csv("output/0_chlorophylla_4geno.csv")
chla <- chla %>% 
  mutate_at("colony_id", str_replace, "_", "-")
chlc <- read.csv("output/0_chlorophyllc_4geno.csv")
chlc <- chlc %>% 
  mutate_at("colony_id", str_replace, "_", "-")
holo_prot <- read.csv("output/0_holobiont_protein_4geno.csv")
holo_prot <- holo_prot %>% 
  mutate_at("colony_id", str_replace, "_", "-")
host_prot <- read.csv("output/0_host_protein_4geno.csv")
host_prot <- host_prot %>% 
  mutate_at("colony_id", str_replace, "_", "-")
sym_counts <- read.csv("output/0_sym_counts_4geno.csv")
sym_counts <- sym_counts %>% 
  mutate_at("colony_id", str_replace, "_", "-")
alpha <- read.csv("output/0_alpha_4geno.csv")
alpha <- alpha %>% 
  mutate_at("colony_id", str_replace, "_", "-")
Pmax <- read.csv("output/0_Pmax_4geno.csv")
Pmax <- Pmax %>% 
  mutate_at("colony_id", str_replace, "_", "-")
Rd <- read.csv("output/0_Rd_4geno.csv")
Rd <- Rd %>% 
  mutate_at("colony_id", str_replace, "_", "-")

#All the individual merges
geno <- full_join(biomass, chla)
geno <- full_join(geno, chlc)
geno <- full_join(geno, holo_prot)
geno <- left_join(geno, host_prot) #absorbing the values from just holo_prot since they have the same column name
geno <- full_join(geno, sym_counts)
geno <- full_join(geno, alpha)
geno <- full_join(geno, Pmax)
geno <- full_join(geno, Rd)


geno_metadata <- geno
geno_metadata <- geno_metadata %>% 
  mutate_at("colony_id", str_replace, "Apul", "ACR")


######
#START WORKIGN HERE
#####

#Edit the current metadata sheet to match the current timeseries metadata
geno_metadata <- geno_metadata %>%
  mutate(Genotype = genotype) %>%
  mutate(month = "19-Oct") %>%
  mutate(nutrient = "NA") %>%
  mutate(site_code = "Nursery") %>%
  mutate(chla.ug.cell = geno_metadata$chla.ug.cm2 / geno_metadata$cells.cm2) %>%
  mutate(chlc2.ug.cell = geno_metadata$chlc2.ug.cm2/ geno_metadata$cells.cm2) %>%
  select(colony_id, Genotype, timepoint, month, nutrient, AFDW.mg.cm2, host_prot_ug.cm2, chla.ug.cm2, chlc2.ug.cm2, Pmax, alpha, Rd, cells.cm2, chla.ug.cell, chlc2.ug.cell)



###Compiled data from Jan/Nov has host and sym AFDW separated while we have it as one big thing, it is also missing protein.
##Need to create new column in data that is Holobiont AFDW
# colony_id, Genotype, timepoint, month, nutrient, site_code
timeseries_data <- read.csv("data/data_jan_nov_SA.csv")


timeseries_data <- timeseries_data %>%
  mutate(AFDW.mg.cm2 = Host_AFDW.mg.cm2 + Sym_AFDW.mg.cm2) %>%
  mutate(colony_id = ï..colony_id) %>%
  mutate(host_prot_ug.cm2 = prot_ug.cm2) %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, AFDW.mg.cm2, host_prot_ug.cm2, chla.ug.cm2, chlc2.ug.cm2, Am, AQY, Rd, cells.cm2, chla.ug.cell, chlc2.ug.cell)

Apul_Plast_Metadata <- rbind(geno_metadata, timeseries_data) #STILL NEED TO ADD Host and Holobiont Protein to it

Apul_Plast_Metadata$group <- paste0(Apul_Plast_Metadata$Genotype,Apul_Plast_Metadata$timepoint, Apul_Plast_Metadata$site_code)

Apul_Plast_Metadata <- Apul_Plast_Metadata %>%
  mutate(total_chl.ug.cm2 = Apul_Plast_Metadata$chla.ug.cm2 + Apul_Plast_Metadata$chlc2.ug.cm2) %>%
  mutate(total_chl.ug.cell = Apul_Plast_Metadata$chla.ug.cell + Apul_Plast_Metadata$chlc2.ug.cell) %>%
  distinct(group, .keep_all = TRUE) %>% #removing duplicate data
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, AFDW.mg.cm2, host_prot_ug.cm2, total_chl.ug.cm2, total_chl.ug.cell, Am, AQY, Rd, cells.cm2, group) %>%
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



####UNIVARIATE Figures and Analysis ---- TOTAL CHL (ug/cm2)
##Plot TOTAL CHL (ug/cm2)##

#Subset full data set for just biomass data
timeseries_tot_chl_cm2 <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, total_chl.ug.cm2)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_tot_chl_cm2$site_code <- factor(timeseries_tot_chl_cm2$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_tot_chl_cm2$Genotype <- factor(timeseries_tot_chl_cm2$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

tot_chl_cm2_fig <- ggplot(timeseries_tot_chl_cm2, aes(x = site_code, y = total_chl.ug.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Total Chlorophyll per Frag (ug/cm2)")))) +
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  ylim(0, 8) +
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )

####UNIVARIATE Figures and Analysis ---- Total CHL per Symbiont (ug/cell)
##Plot CHL per Symbiont (ug/cell)##

#Subset full data set for just biomass data
timeseries_tot_chl_cell <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, total_chl.ug.cell)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_tot_chl_cell$site_code <- factor(timeseries_tot_chl_cell$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_tot_chl_cell$Genotype <- factor(timeseries_tot_chl_cell$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

tot_chl_cell_fig <- ggplot(timeseries_tot_chl_cell, aes(x = site_code, y = total_chl.ug.cell, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Total Chlorophyll per Symbiont (ug/cell)"))))+
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  #ylim(0, 4) +
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )


####UNIVARIATE Figures and Analysis ---- SYM Density (cells/cm2)
##Plot Sym Density (cells/cm2)##

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


####UNIVARIATE Figures and Analysis ---- Host Prot (ug/cm2)
##Plot Host Prot (ug/cm2)##

#Subset full data set for just biomass data
timeseries_host_prot <- Apul_Plast_Metadata %>%
  select(colony_id, Genotype, timepoint, month, nutrient, site_code, host_prot_ug.cm2)

#Reorder data from beginning to try to get them all in a row at least:
timeseries_host_prot$site_code <- factor(timeseries_host_prot$site_code, levels = c("Nursery", "Mahana Low", "Hilton Medium", "Manava High")) 
timeseries_host_prot$Genotype <- factor(timeseries_host_prot$Genotype, levels = c("Genotype 4", "Genotype 6", "Genotype 8", "Genotype 15")) 

host_prot_fig <- ggplot(timeseries_host_prot, aes(x = site_code, y = host_prot_ug.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch =21, size=2, position = position_jitterdodge(0.2)) +
  scale_x_discrete(labels=c("Nursery" = "Nursery", "Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                            "Manava High" = "Manava \nHigh"))+
  scale_y_continuous(labels = scales::scientific) +
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("#EDF8B1", "#7FCDBB", "#2C7FB8"), limits = c("timepoint0", "timepoint1", "timepoint4"), labels = c("October 2019", "January 2020", "November 2020")) +
  xlab("Site") + 
  ylab(expression(bold(paste("Host Total Protein (ug/cm2)"))))+
  geom_vline(xintercept = 1.5, linetype = "longdash")+
  theme_classic() + 
  theme(
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(size=10),
    axis.title=element_text(face="bold", size=10),
    axis.text=element_text(size=8, color="black"), 
    strip.text.x=element_text(face="italic", size=10)
  )

host_prot_fig
###Configure all univariate functions into one major figure

Fig <- ggarrange(biomass_fig, host_prot_fig, sym_counts_fig, tot_chl_cell_fig, tot_chl_cm2_fig, ncol = 3, nrow = 2)
ggsave("Output/Univariate_Figs.pdf", Fig, width=16, height=8)

Fig

colnames(Apul_Plast_Metadata)

fixed_2way <- aov(AFDW.mg.cm2 ~site_code*timepoint, Apul_Plast_Metadata)
summary(fixed_2way)



library(lme4)
mixed_2way <- lmer(AFDW.mg.cm2 ~site_code*timepoint + (1 | Genotype), Apul_Plast_Metadata)
summary(mixed_2way)

