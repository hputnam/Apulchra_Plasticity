#Summary of TP0 Univariate Response Variables (Except AQY, Am, and R)
library(tidyverse)
library(ggpubr)
library(showtext)

#All saved Figures composed into one figure
Univariate_Summary_Fig <- ggarrange(Fig.1,Fig.2,Fig.3,Fig.4,Fig.5,Fig.6, ncol = 3, nrow = 2)

#Make Summary Table of all ANOVA's and TukeyHSD
#All data was log transformed to meet the assumptions of normality


ggsave("output/TP0_Univariate_Figs.pdf", Univariate_Summary_Fig, width=10, height=8)
