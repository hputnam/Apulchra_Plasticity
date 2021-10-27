# set working directory and load necessary packages
library(vegan)
library(tidyverse)
library(reshape2)
library(ggbiplot)
library(ggpubr)

# set seed
set.seed(54321)

# import data
data<- read.csv("data/complete_timeseries_data_test.csv") #read in file

#check values for outliers

#Scale and center datasets
data_scaled <- scale(data[5:13], center = T, scale = T) # scaled variables

#Identify Factors 
fac <- data[1:4]

# PCA of all data
pca.out <- prcomp(data_scaled, center=FALSE, scale=FALSE) #run PCA
summary(pca.out) #view results
biplot(pca.out) #plot results
dev.off()

#combine PC1 and PC2 into dataframe for plotting and calculation of PCA distance for plasticity section below
PCs <- as.data.frame(cbind(pca.out$x[,1], pca.out$x[,2]))
PCs.meta <- cbind(fac, PCs)
PCs.meta <- PCs.meta %>%
  select(-colony_id)

PCs.meta.wide <- pivot_wider(PCs.meta, values_from = c(V1, V2), names_from =timepoint)

TP1 <- ggplot(PCs.meta.wide,aes(x=V1_timepoint0,y=V2_timepoint0))+
  geom_point(aes(), size=2, color="black")+ 
  geom_segment(aes(x=V1_timepoint0, y=V2_timepoint0, xend=V1_timepoint1, yend=V2_timepoint1, color=site))+
  facet_wrap(~Genotype)+
  xlim(-4,4)+
  ylim(-4,4)+
  xlab(label = "PC1")+
  ylab(label = "PC2")+
  ggtitle("A) January 2020") +
  theme_classic() 

TP4 <- ggplot(PCs.meta.wide,aes(x=V1_timepoint0,y=V2_timepoint0))+
  geom_point(aes(), size=2, color="black")+ 
  geom_segment(aes(x=V1_timepoint0, y=V2_timepoint0, xend=V1_timepoint4, yend=V2_timepoint4, color=site))+
  facet_wrap(~Genotype)+
  xlim(-4,4)+
  ylim(-4,4)+
  xlab(label = "PC1")+
  ylab(label = "PC2")+
  ggtitle("B) November 2020") +
  theme_classic() 

Plasticity <- ggarrange(TP1, TP4,  common.legend = TRUE)
ggsave("output/Plasticity_pca.pdf", width = 8, height = 4, units = "in")


#Chaculate distance between origin and transplant points
#plasticity (amount of change in X,Y space) between the genotypes for each time and site relative to the origin at TP0
PCA.dist <- as.data.frame(as.matrix(dist(PCs)))

#join with metadata
PCA.dist <- cbind(fac, PCA.dist)


#filter and groupby to get only the data of interest




#then need to extract distance for matched set that we want 
# e.g., Genotype 15 T0 - Genotype 15 all times and all sites
#Genotype 15 T0 - Genotype 15 T1 S1 
#Genotype 15 T0 - Genotype 15 T1 S2 
#Genotype 15 T0 - Genotype 15 T1 S3 
#Genotype 15 T0 - Genotype 15 T2 S1 
#Genotype 15 T0 - Genotype 15 T2 S2 
#Genotype 15 T0 - Genotype 15 T2 S3 

#view summary stats
range(T3.PCA.dist.df$value)
boxplot(T3.PCA.dist.df$value)

# Test for differences in plasticity
T3.PCA.plast.mod <-aov(value ~Species*Origin, data=T3.PCA.dist.df)
anova(aov(value ~Species*Origin, data=T3.PCA.dist.df))

