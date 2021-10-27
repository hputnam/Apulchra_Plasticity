# set working directory and load necessary packages
library(tidyverse)
library(ggpubr)
library(GGally)
library(lemon)


# import data
data<- read.csv("data/complete_timeseries_data_test.csv") #read in file

#check values for outliers
ggpairs(data[,5:13])

#Scale and center datasets
data_scaled <- scale(data[5:13], center = T, scale = T) # scaled variables

#Identify Factors 
fac <- data[1:4]

# PCA of all variables
pca.out <- prcomp(data_scaled, center=FALSE, scale=FALSE) #run PCA
summary(pca.out) #view results

#combine PC1 and PC2 into dataframe for plotting and calculation of PCA distance for plasticity section below
PCs <- as.data.frame(cbind(pca.out$x[,1], pca.out$x[,2]))
PCs.meta <- cbind(fac, PCs)
PCs.meta <- PCs.meta %>%
  select(-colony_id)

#convert to wide format to enable initial and subsequent comparisons
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


#Calculate distance between origin and transplant points
#plasticity (amount of change in X,Y space) between the genotypes for each time and site relative to the origin at TP0
PCA.dist <- as.data.frame(as.matrix(dist(PCs)))

#join with metadata
PCA.dist <- cbind(fac, PCA.dist)

#filter and select to get only the data of interest
G15 <- PCA.dist %>%
  filter(Genotype=="Genotype15")
G15 <- G15[4:nrow(G15),1:5]
colnames(G15)[5] <- "Distance"

G8 <- PCA.dist %>%
  filter(Genotype=="Genotype8")
G8 <- G8[4:nrow(G8),c(1:4,14)]
colnames(G8)[5] <- "Distance"

G4 <- PCA.dist %>%
  filter(Genotype=="Genotype4")
G4 <- G4[4:nrow(G4),c(1:4,11)]
colnames(G4)[5] <- "Distance"

G6 <- PCA.dist %>%
  filter(Genotype=="Genotype6")
G6 <- G6[4:nrow(G6),c(1:4,8)]
colnames(G6)[5] <- "Distance"

Plast.Data <- rbind(G15,G4,G6,G8)
Plast.Data$group <- paste0(Plast.Data$Genotype, Plast.Data$site)
Plast.Data$groupTS <- paste0(Plast.Data$timepoint, Plast.Data$site)

# Test for differences in plasticity
plast.mod <-aov(Distance ~site*timepoint, data=Plast.Data)
summary(plast.mod)
hist(plast.mod$residuals)
qqnorm(plast.mod$residuals)
qqline(plast.mod$residuals)


# Plast.Data %>%
#   ggplot(aes(timepoint,Distance)) +
#   geom_boxplot(aes(fill=site)) +
#   geom_point(aes(group=groupTS, shape=Genotype), color="black", position=position_dodge(width=0.75))+
#   geom_line(aes(group=group, color=site), alpha = .4) +
#   theme_classic() 


pj <- position_jitterdodge(jitter.width=0.05, seed=9,
                           jitter.height = 0,
                           dodge.width = 0.75)
RNorms <- Plast.Data %>%
ggplot(aes(x=timepoint, y=Distance, fill=site)) +
  geom_boxplot(outlier.colour = NA, width =0.75, alpha=0.2) +
  lemon::geom_pointpath(aes(colour=site, shape=Genotype,group=interaction(Genotype,site)),
                        position = pj, alpha=0.4)+
  theme_classic() 

RNorm<-Plasticity <- ggarrange(RNorms,  common.legend = F)
ggsave("output/ReactionNorms.pdf", width = 4, height = 4, units = "in")
