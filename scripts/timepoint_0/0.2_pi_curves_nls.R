## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("phytotools" %in% rownames(installed.packages()) == 'FALSE') install.packages('phytotools') 
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse') 
if ("broom" %in% rownames(installed.packages()) == 'FALSE') install.packages('broom') 

#Read in required libraries

library("devtools")
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('phytotools')
library("tidyverse")
library("broom")


# Import data
Data <- read.csv(file = 'output/0_pi_curve_rates.csv')
Data <- Data[,-1]

#Define data
Data$PAR <- as.numeric(Data$Light_Value)
Data$Pc <- as.numeric(Data$micromol.cm2.h)
Data<-Data%>%
  filter(!Pc=="NA")

# Define PI curve function as a nonlinear Least Squares regression of a quadratic fit, test nls fit
#Aquatic Photosynthesis, Falkowski   
#Pmax = max photosynthesis (AKA Am)  
#alpha = quantum yeild (AKA AQY quantum yield)  
#I/E = irradiance (AKA PAR)  
#Rd = dark respiration   

nls_data <- Data %>% 
  group_by(colony_id) %>%
  nest(-colony_id) %>%
  mutate(model1 = map(data, ~ 
                        nls(Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd), data=., start=list(Am=(max(.$Pc)-min(.$Pc)),  AQY=0.005, Rd=-min(.$Pc))) %>%
                        #nls(Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd), data=., start=list(Am=1,  AQY=0.01, Rd=0.2)) %>%
                        tidy %>%
                        dplyr::select(term, estimate) %>% 
                        spread(term, estimate))) %>%
  unnest(model1) %>%
  unnest(data) %>%
  group_by(colony_id)%>%
  summarise(Am=mean(Am), AQY=mean(AQY), Rd=mean(Rd))%>%
  mutate(timepoint="timepoint0")

write_csv(nls_data, "output/0_pi_curve_pars_nls.csv")
  
  
md <- read_csv("metadata/coral_metadata.csv")
md <- md %>% 
  mutate_at("colony_id", str_replace, "_", "-")
  
df <- left_join(nls_data, md)



# Pmax Nursery 4 genotypes
Pmax_4geno <- df %>%
  filter(genotype == "Genotype15"| genotype == "Genotype4"| genotype == "Genotype6"|genotype == "Genotype8") %>%
  select(colony_id, site, genotype, Am, timepoint) %>%
  write_csv(., path = "output/0_Pmax_4geno.csv")
  

# alpha Nursery 4 genotypes
alpha_4geno <- df %>%
  filter(genotype == "Genotype15"| genotype == "Genotype4"| genotype == "Genotype6"|genotype == "Genotype8") %>%
  select(colony_id, site, genotype, AQY, timepoint) %>%
  write_csv(., path = "output/0_alpha_4geno.csv")


# Pmax Nursery 4 genotypes
Pd_4geno <- df %>%
  filter(genotype == "Genotype15"| genotype == "Genotype4"| genotype == "Genotype6"|genotype == "Genotype8") %>%
  select(colony_id, site, genotype, Rd, timepoint) %>%
  write_csv(., path = "output/0_Rd_4geno.csv")
  


