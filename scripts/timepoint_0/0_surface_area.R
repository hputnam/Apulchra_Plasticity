#Must run surface area scripts first in order to run biomass

#Library Loading__________________________________________
## install packages if you dont already have them
if (!require("tidyverse")) install.packages("tidyverse")

# load packages
library(tidyverse)

#load wax data
setwd("C:/Users/Dennis/Documents/Github/Apulchra_Plasticity")

wax.data <- read.csv("data/timepoint_0/0_surface_area/0_surface_area_data.csv", header=TRUE)

wax.data$delta.mass.g <- wax.data$weight2.g-wax.data$weight1.g
stnds <- subset(wax.data, Sample=="Standard")
#stnds <- stnds[-1,] # the largest data point was an artifact outlier. It was removed, as it had bubbles of air escape the wooden sphere, leaving gaps in the wax

#calculate the surface area of the spherical standards from the diameter
stnds$rad <- stnds$Diameter/2
stnds$surface.area.cm2 <- 4*pi*(stnds$rad)^2

# calculate the curve coefficients for slope and intercept to apply as the standard
stnd.curve <- lm(surface.area.cm2~delta.mass.g, data=stnds)
plot(surface.area.cm2~delta.mass.g, data=stnds)
stnd.curve$coefficients
summary(stnd.curve)$r.squared #for tp0 the r2 value is .918, so pretty good

#Calculate surface area using the standard curve
smpls <- subset(wax.data, Sample=="Coral")
smpls$surface.area.cm2 <- stnd.curve$coefficients[2] * smpls$delta.mass.g + stnd.curve$coefficients[1]

#select the samples only
smpls <- smpls %>%
  select(-Sample, -Diameter)

#rename column colony_id because weird things happened to it.
smpls <- smpls %>%
  mutate(colony_id = ï..colony_id)


smpls %>%
  count(colony_id) %>% arrange(n)

#check the range to make sure your samples fall within the range of the standards
range(smpls$surface.area.cm2)
range(stnds$surface.area.cm2)

#Save the output for use in normilzation for phys assays
smpls%>%
  mutate(timepoint="timepoint0")%>%
  select(colony_id, surface.area.cm2, timepoint)%>%
  write_csv("output/0_surface_area.csv")
