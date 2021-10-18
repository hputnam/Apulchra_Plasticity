#Summary of TP0 Univariate Response Variables (Except AQY, Am, and R)
library(tidyverse)
library(ggpubr)
library(showtext)

#All saved Figures composed into one figure
chla <- read.csv("output/0_chlorophylla.csv")
chla$variable <- colnames(chla[3])
colnames(chla)[3] <- "value"

chlc2 <- read.csv("output/0_chlorophyllc.csv")
chlc2$variable <- colnames(chlc2[3])
colnames(chlc2)[3] <- "value"

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
