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

# add AFDW
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

