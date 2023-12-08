#whole fly bacteria loads from 2%, 16% from conventionally reared flies 

#package download + library upload
#install.packages("ggplot2")
#install.packages("tidyverse")
#install.packages("plyr")
#install.packages("lme4")
#install.packages("car")
#install.packages("multcomp")
#install.packages("ggbeeswarm")
#install.packages("lattice")
#install.packages("ggpubr")
#install.packages("emmeans")
library(emmeans)
library(ggpubr)
#adds error bars 
install.packages("Hmisc")
library(lattice)
library(Hmisc)
library(ggbeeswarm)
library(multcomp)
library(lme4)
library(ggplot2)
#install.packages('dplyr')
library(dplyr)
#install.packages('tidyr')
library(tidyr)
#install.packages('reshape')
library(reshape)
#install.packages('tibble')
library(tibble)
library(ggpmisc)
library(lme4)



#csv file read in
load <- read.csv("bacteria_loads.csv", header = T)


#reorder Diet levels and rename 2% and 16% samples
load$Diet <- factor(load$Diet, levels=c("2%", "16%"))

#colors for graph
two_diet = c( "#3399FF", "#FF9966" )

#box plot of bacteria loads

plot<- ggplot(data=load, aes(x=Plate_Type, y = log2(CFU_fly), fill = Diet)) +
  geom_boxplot(aes(fill = Diet), outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = Diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "#3399FF", "16%" = "#FF9966"))+
  scale_shape_manual(values = c(15, 16))+
  labs(y="CFU per Fly (Log2)", x = "Media", fill = "Diet")+
  theme_classic()

plot+
  theme(
    #changes font size and face of the title of plot
    plot.title=element_text(hjust = 0.5,size=14,face="bold"),
    #edits the line thickness and color of x and y axis
    axis.line = element_line(size = 1.5),
    #edits the thickness of tick marks on x and y axis
    axis.ticks = element_line(size =1.5, color ="black"),
    #moves tick marks from outside to inside using - number for unit
    axis.ticks.length = unit(-.25, "cm"),
    #changes siz and bold face of x and y axis units
    axis.text = element_text(size = 18, face= 1.5),
    #changes size and bold face of title of axes
    axis.title = element_text(size = 20, face = 2),
    #changes the size and face of the legend text
    legend.text = element_text(size=20),
    #changes tile of legend size and face
    legend.title = element_text(size=22, face = 2),
    #changes margins of axis title and axis text
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 00, l = 0)),
    panel.spacing = unit(3, "cm"))



######Stats

#linear model of 2%, 16%, and LF diets
load_lmer <- lmer(log2(CFU_fly) ~ Diet+Plate_Type *Diet + (1|Block), data = load)

load_emm <- emmeans(load_lmer, pairwise ~ Diet| Plate_Type)
load_emm

