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
source("src/CreateBacteriaxVirusmod.R")

# create pomp object ---------------------------
po<-CreateBacteriaxVirusmod()

# deterministic simulation ----------------------
tj0<-trajectory(po, format= "data.frame")
plt.tj.c<-ggplot(tj0, aes(x = time , y = Y_C ))+
    geom_line()+
    theme_classic()
plt.tj.i <- ggplot(tj0, aes(x = time , y = Y_I ))+
    geom_line()+
    theme_classic()
plt.tj.ci<-ggplot(tj0, aes(x = time , y = Y_CI ))+
    geom_line()+
    theme_classic()
plt.comb<-ggplot(tj0 %>%
                select(Y_C, Y_I, Y_CI, time) %>% 
                pivot_longer(cols = starts_with("Y_")), 
                aes(x = time , y = value, color = name ))+
    geom_line()+
    theme_classic()

print(plt.tj.i)
dev.new()
print(plt.tj.c)
dev.new()
print(plt.tj.ci)