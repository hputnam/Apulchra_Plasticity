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
  select(colony_id, site, timepoint, AFDW.mg.cm2, host_prot_ug.cm2, sym_prot_ug.cm2, total_chl.ug.cm2, total_chl.ug.cell, Am, AQY, Rd, cells.cm2) #%>% #selecting the six metrics wanting to output

TP0_metadata %>% #writing output of metadata for TP0 (all nursery - n=10 and wild sites n=3 per site)
  write.csv("output/TP0_metadata")

#Individually create plots and save them before compiling into one large output plot
#biomass fig

#Biomass Fig
AFDW_TP0 <- TP0_metadata %>%
  ggplot(aes(x = site, y = AFDW.mg.cm2, color = site)) +
  labs(x = "Site", y = "", color = "Site") +
  ylim(0, 3.5) +
  ylab("mg/cm2") +
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
  ylab("µg/cm2") +
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
  ylab("µg/cm2") +
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
  ylab("µg/cell") +
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
  ylab("cells/cm2") +
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
  ggtitle("Max Photosynthesis (Am)") +
  ylab("µmol/cm2/h") +
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
  ggtitle("alpha (AQY)") +
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
  ggtitle("Respiration Rate (Rd)") +
  ylab("µmol/cm2/h") +
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
#Biomass__________________________________________
model.AFDW <- aov(log10(AFDW.mg.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.AFDW$residuals)
qqnorm(model.AFDW$residuals)
qqline(model.AFDW$residuals)
plot(model.AFDW$fitted.values, model.AFDW$residuals)
biomass_anova<-summary(model.AFDW) #anova for biomass
biomass_tkyhsd<-TukeyHSD(model.AFDW) #posthoc tests within anova for biomass

#Host Protein ________________________________
model.HostProt <- aov(log10(host_prot_ug.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.HostProt$residuals)
qqnorm(model.HostProt$residuals)
qqline(model.HostProt$residuals)
plot(model.HostProt$fitted.values, model.HostProt$residuals)
host_prot_anova<-anova(model.HostProt) #anova for biomass
host_prot_tkyhsd<-TukeyHSD(model.HostProt) #posthoc tests within anova for biomass

#Sym Density______________________________________
model.CellDens <- aov(log10(cells.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.CellDens$residuals)
qqnorm(model.CellDens$residuals)
qqline(model.CellDens$residuals)
plot(model.CellDens$fitted.values, model.CellDens$residuals)
sym_dens_anova<-anova(model.CellDens) #anova for Sym Density
sym_dens_tkyhsd<-TukeyHSD(model.CellDens) #posthoc tests within anova for Sym Density


#Total Chlorophyll per Frag (ug/cm2)_______________
model.Chlcm <- aov(log10(total_chl.ug.cm2) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.Chlcm$residuals)
qqnorm(model.Chlcm$residuals)
qqline(model.Chlcm$residuals)
plot(model.Chlcm$fitted.values, model.Chlcm$residuals)
tot_chl_anova<-anova(model.Chlcm) #anova for Total Chl
tot_chl_tkyhsd<-TukeyHSD(model.Chlcm) #posthoc tests within anova for Total Chl


#Total Chlorophyll per cell (ug/cells)____________
model.Chl.cell <- aov(log10(total_chl.ug.cell) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.Chl.cell$residuals)
qqnorm(model.Chl.cell$residuals)
qqline(model.Chl.cell$residuals)
plot(model.Chl.cell$fitted.values, model.Chl.cell$residuals)
tot_chl_per.cell_anova<-anova(model.Chl.cell) #anova for Total Chl per cell
tot_chl_per.cell_tkyhsd<-TukeyHSD(model.Chl.cell) #posthoc tests within anova for Total Chl per cell

#Am (Max Photosynthesis)____________
model.Am <- aov(log10(Am) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.Am$residuals)
qqnorm(model.Am$residuals)
qqline(model.Am$residuals)
plot(model.Am$fitted.values, model.Chl.cell$residuals)
Am_anova<-anova(model.Am) #anova for Am
Am_tkyhsd<-TukeyHSD(model.Am) #posthoc tests within anova for Am

#AQY (Max Photosynthetic Rate)____________
model.AQY <- aov(log10(AQY) ~ site, data = TP0_metadata) #save model 
par(mfrow=c(1,3))
hist(model.AQY$residuals)
qqnorm(model.AQY$residuals)
qqline(model.AQY$residuals)
plot(model.AQY$fitted.values, model.Chl.cell$residuals)
AQY_anova<-anova(model.AQY) #anova for AQY
AQY_tkyhsd<-TukeyHSD(model.AQY) #posthoc tests within anova for AQY

#Rd (Respiration Rate)____________
model.Rd <- aov(log10(Rd) ~ site, data = TP0_metadata) #save model
par(mfrow=c(1,3))
hist(model.Rd$residuals)
qqnorm(model.Rd$residuals)
qqline(model.Rd$residuals)
plot(model.Rd$fitted.values, model.Chl.cell$residuals)
Rd_anova<-anova(model.Rd) #anova for Rd
Rd_tkyhsd<-TukeyHSD(model.Rd) #posthoc tests within anova for Rd


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





