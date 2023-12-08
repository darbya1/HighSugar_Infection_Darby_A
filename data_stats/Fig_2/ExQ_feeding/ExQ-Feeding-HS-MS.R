#set working directory
setwd("C:/Users/amd439/Box/Darby_Drea/Projects/Sucrose-diet-variation-2021/Ex-Q/")

#package download + library upload

library(lattice)
library(Hmisc)
library(ggbeeswarm)
library(multcomp)
library(lme4)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(emmeans)


#All 3 blocks with six diets 0%-24%
exq <- read.csv("Ex-Q-Sucrose-All-Blocks.csv", header= T)

#read Day as factor not a numeric number
exq$Day <- as.character(exq$Day)

#sort increasing order of diets
exq$Diet <- factor(exq$Diet,
                   levels = c("0%", "2%", "4%", "8%", "16%", "24%"))
#assign colors for diets
six_color = c( "#99CCFF","#3399FF", "#003399",  "#FFCC66", "#FF9966", "#CC6633")




#barplot of all six diets
exqbar<-ggplot(exq, aes(x=Day, y= Mean, fill=Diet))+
  geom_bar(stat="summary", fun.y = "mean", position = "dodge",colour = "black", linewidth = 1) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position= position_dodge(.90), width = .4, size = .9)+
  labs(y =  "Mean Dye consumed (ug) per fly", legend = "Sucrose w/v") + scale_fill_manual(values = six_color)+
  theme_classic()#+annotate("text", x = 1.22, y = , label = "***", size = 14)#+
  

exqbar<-exqbar + #ggtitle(charttitle)+
  theme( #changes font size and face of the title of plot
    #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
    axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
                axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
                axis.text.x = element_text(size = 18, family = "sans"), 
                axis.text.y = element_text(size = 18, family = "sans"),
               legend.text = element_text(size=20),
               #changes tile of legend size and face
               legend.title = element_text(size=25, face = 2))
  
exqbar


#two-way analysis of  means dye
two.wayA <-aov(Mean ~ Diet + Day, data = exq)
summary(two.wayA)

#pairwise comparisons between diet within day 
emm <- emmeans(two.wayA, pairwise ~ Diet|Day)
emm