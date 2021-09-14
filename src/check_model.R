##
## check model
## Author: Elizabeth Goult
## Comments:

# clear workspace 
rm(list=ls())
graphics.off()

#libraries
require(pomp)
require(tidyverse)
source("src/CreateSCxSIRmod.R")

# create pomp object ---------------------------
po<-CreateSCxSIRmod()

# deterministic simulation ----------------------
tj0<-trajectory(po, format= "data.frame")[-1,]
plt.tj.c<-ggplot(tj0, aes(x = time , y = round(Y_bac, digits = 15)))+
    geom_line()+
    theme_classic()
plt.tj.i <- ggplot(tj0, aes(x = time , y = Y_vir ))+
    geom_line()+
    theme_classic()
plt.tj.ci<-ggplot(tj0, aes(x = time , y = Y_bac_vir ))+
    geom_line()+
    theme_classic()
plt.comb<-ggplot(tj0 %>%
                select(Y_bac, Y_vir, Y_bac_vir, time) %>% 
                pivot_longer(cols = starts_with("Y_")), 
                aes(x = time , y = value, color = name ))+
    geom_line()+
    theme_classic()

print(plt.tj.i)
dev.new()
print(plt.tj.c)
dev.new()
print(plt.tj.ci)