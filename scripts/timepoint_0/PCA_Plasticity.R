# set working directory and load necessary packages
library(tidyverse)
library(ggpubr)
library(GGally)
library(lemon)
library(vegan)
library(ggfortify)
library(ggplot2)

##### IMPORT DATA ##### 

# import data
data<- read.csv("data/complete_timeseries_data.csv") #read in file

#check values for outliers
data%>%
  ggpairs(., columns=4:11, ggplot2::aes(colour=site)) 
data%>%
  ggpairs(., columns=4:11, ggplot2::aes(colour=timepoint)) 
data%>%
  ggpairs(., columns=4:11) 

##### MAKING PCAs #####  

#preparing dataset for making a PCA and performing a PERMANOVA
pca_data<-data[complete.cases(data), ]

scaled_pca_data<-prcomp(pca_data[c(4:11)], scale=TRUE, center=TRUE) 

#make PCA
pca_plot<-ggplot2::autoplot(scaled_pca_data, data=pca_data, frame.colour="Genotype", loadings=FALSE,  colour="Genotype", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  theme_classic()+
  theme(legend.text = element_text(size=18), 
        legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pca_plot

pca_plot<-ggplot2::autoplot(scaled_pca_data, data=pca_data, frame.colour="site", loadings=FALSE,  colour="site", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  theme_classic()+
  theme(legend.text = element_text(size=18), 
        legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pca_plot

pca_plot<-ggplot2::autoplot(scaled_pca_data, data=pca_data, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  theme_classic()+
  theme(legend.text = element_text(size=18), 
        legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pca_plot

#test with PERMANOVA
# scale data
vegan_data <- scale(pca_data[ ,4:11])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis(vegan_data ~ site*timepoint, data = pca_data, method='eu')
adonis(vegan_data ~ Genotype, data = pca_data, method='eu')




##### PLASTICITY #####
#generate t0 data for each colony to make statistical comparisons

#HERE WE NEED TO CORRECT THE FORMAT TO COPY THE VALUES FOR NURSERY TO EACH SITE
data<- read.csv("data/complete_timeseries_data.csv") #read in file

#alter data file here




#Scale and center datasets
data_scaled <- scale(data[4:11], center = T, scale = T) # scaled variables

#Identify Factors 
fac <- data[1:3]

# PCA of all variables
pca.out <- prcomp(data_scaled, center=FALSE, scale=FALSE) #run PCA
summary(pca.out) #view results

#combine PC1 and PC2 into dataframe for plotting and calculation of PCA distance for plasticity section below
PCs <- as.data.frame(cbind(pca.out$x[,1], pca.out$x[,2]))
PCs.meta <- cbind(fac, PCs)
#PCs.meta <- PCs.meta %>%
  #select(-colony_id)

#convert to wide format to enable initial and subsequent comparisons
PCs.meta.wide <- pivot_wider(PCs.meta, values_from = c(V1, V2), names_from =timepoint)



#TEMPORARY FIX FOR NURSERY FORMAT - WILL DELETE AFTER ABOVE ISSUE IS FIXED FOR COPYING NURSERY VALUES
PCs.meta.wide%>%write_csv("data/pcs_meta_wide.csv")

#manually copied nursery values to timepoint 0 for all genotypes
PCs.meta.wide<-read.csv("data/pcs_meta_wide.csv")





#Creating PCA Plot from October to January Comparison
TP1 <- ggplot(PCs.meta.wide,aes(x=V1_timepoint0,y=V2_timepoint0))+
  geom_point(aes(), size=2, color="black")+ 
  geom_segment(aes(x=V1_timepoint0, y=V2_timepoint0, xend=V1_timepoint1, yend=V2_timepoint1, color=site))+
  facet_wrap(~Genotype)+
  xlim(-4,5)+
  ylim(-4,5)+
  xlab(label = "PC1")+
  ylab(label = "PC2")+
  ggtitle("A) January 2020") +
  theme_classic() 

#Creating PCA Plot from October to November Comparison
TP4 <- ggplot(PCs.meta.wide,aes(x=V1_timepoint0,y=V2_timepoint0))+
  geom_point(aes(), size=2, color="black")+ 
  geom_segment(aes(x=V1_timepoint0, y=V2_timepoint0, xend=V1_timepoint4, yend=V2_timepoint4, color=site))+
  facet_wrap(~Genotype)+
  xlim(-4,5)+
  ylim(-4,5)+
  xlab(label = "PC1")+
  ylab(label = "PC2")+
  ggtitle("B) November 2020") +
  theme_classic() 

Plasticity <- ggarrange(TP1, TP4,  common.legend = TRUE) 
ggsave("output/Plasticity_pca.pdf", width = 8, height = 4, units = "in")


#Calculate distance between origin and transplant points
#plasticity (amount of change in X,Y space) between the genotypes for each time and site relative to the origin at TP0

PCA.dist <- as.data.frame(as.matrix(dist(PCs)))



#convert PCs.meta.wide back to long format (this is to fix the formatting issues where we want to replicate nursery values)  - AH come back to this after the nursery values can be copied in the code

#join with metadata
PCA.dist <- cbind(fac, PCA.dist)

#filter and select to get only the data of interest 
###NOT SURE IF THIS IS CORRECT FROM HERE ON SINCE CHANGING THE INPUT DATA###
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

#try kruskal wallis tests for low sample size
kruskal.test(Distance ~ site, data=Plast.Data)
kruskal.test(Distance ~ timepoint, data=Plast.Data)


# Plast.Data %>%
#   ggplot(aes(timepoint,Distance)) +
#   geom_boxplot(aes(fill=site)) +
#   geom_point(aes(group=groupTS, shape=Genotype), color="black", position=position_dodge(width=0.75))+
#   geom_line(aes(group=group, color=site), alpha = .4) +
#   theme_classic() 


pj <- position_jitterdodge(jitter.width=0.05, seed=9,
                           jitter.height = 0,
                           dodge.width = 0.75)

#create group mean dataset 
gd <- Plast.Data %>%
  group_by(site, timepoint)%>%
  summarise(mean=mean(Distance))

#install.packages('Rcpp')
#library(Rcpp)

#Create final plot for Reaction Norms
###ERROR IN CODE - say my split_indices function is not provided by Rcpp function##
RNorms <- Plast.Data %>%
ggplot(aes(x=timepoint, y=Distance, fill=site)) +
  geom_boxplot(outlier.colour = NA, width =0.75, alpha=0.2) +
  lemon::geom_pointpath(aes(colour=site, shape=Genotype,group=interaction(Genotype,site)),
                        position = pj, alpha=0.4)+
  geom_point(aes(colour=site, x=timepoint, y=mean), data=gd)+
  geom_line(aes(colour=site, group=site, x=timepoint, y=mean), data=gd, size=1, alpha=1)+
  theme_classic(); RNorms 

Plasticity <- ggarrange(RNorms,  common.legend = F)
ggsave("output/ReactionNorms.pdf", width = 4, height = 4, units = "in")
