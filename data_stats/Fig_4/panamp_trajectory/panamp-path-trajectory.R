#pathogen trajectory assay done with Providencia retgerri and Serratia marcescens 2698B fed on 2% or 16% diets using AMP deletion lines

#library upload
library(ggpubr)
library(lattice)
library(Hmisc)
library(ggbeeswarm)
library(multcomp)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape)
library(tibble)
library(ggpmisc)
library(nlme)
library(emmeans)

#set working directory
setwd("C:/Users/amd439/Box/Darby_Drea/Projects/Sucrose-diet-variation-2021/pathogen-trajectory/")


#PanAMP trajectory

#PanAmp_trajectory: Serratia 2698B
burden3 <- read.csv("PanAMP-Path-trajectory-16-vs-2-sucrose/panamp-path-trajectory-2698B.csv", header = T)

#PanAMP trajectory: Providencia
burden4 <- read.csv("PanAMP-Path-trajectory-16-vs-2-sucrose/panamp-path-trajectory-providencia.csv", header = T)



######################################################################
##Generalized least squares (GLS) requires nlme package also
#install.packages("nlme")
library(nlme)
library(emmeans)

###S. marcescens in AMP deletion lines
############################################################################################################################################

#make diet and hour post infection characters, so that barplot can be made appropriately
burden3$diet <- as.factor(burden3$diet)
burden3$dilution <- as.factor(burden3$dilution)

#reorder diet levels
burden3$diet <- factor(burden3$diet, levels=c("2%", "16%"), labels=c("2%", "16%"))

#adjust order of time series
burden3$time_f<- factor(burden3$time_point,
                        levels = c("0", "2", "4", "6","8"))
burden3$log2CFU <- log2(burden3$CFU_fly)
##code for boxplot + geom_point position with jitter in geom_point

AMP_SM<- ggplot(data=burden3, aes(x=time_f, y = log2CFU, fill = diet)) +
  geom_boxplot(aes(fill = diet), outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  #scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  labs(y="Log2 CFU per Fly", x = "Time (hour)", fill = "Sucrose (w/v)")+
  #scale_y_continuous(limits = c(10, 30))+ 
  theme_classic()
AMP_SM


#removes outlier points geom_boxplot(outlier.shape = NA)
AMP_SM + 
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

######################################################################


#GlS allows to fit a 2-way-ANPVA but allows for different variance across components 

SM_m3 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_point), data = burden3)

#weighs time and diet interaction at each time point
SM_m4 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_f*diet), data = burden3)

#comparing models 
anova(SM_m3, SM_m4)

#contrasts of each time point between each diet
SM_emm<-emmeans(SM_m4, pairwise ~ diet | time_f)
SM_emm

###P. rettgeri  in AMP deletion lines
############################################################################################################################################

#make diet and hour post infection characters, so that barplot can be made appropriately
burden4$diet <- as.factor(burden4$diet)
burden4$dilution <- as.factor(burden4$dilution)

#reorder diet levels
burden4$diet <- factor(burden4$diet, levels=c("2%", "16%"), labels=c("2%", "16%"))

#adjust order of time series
burden4$time_f<- factor(burden4$time_point,
                        levels = c("0", "2", "4", "6","8", "10", "12"))
burden4$log2CFU <- log2(burden4$CFU_fly)

##code for boxplot + geom_point position with jitter in geom_point

AMP_PR<- ggplot(data=burden4, aes(x=time_f, y = log2CFU, fill = diet)) +
  geom_boxplot(aes(fill = diet), outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  #scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  labs(y="Log2 CFU per Fly", x = "Time (hour)", fill = "Sucrose (w/v)")+
  #scale_y_continuous(limits = c(10, 30))+ 
  theme_classic()
AMP_PR


#removes outlier points geom_boxplot(outlier.shape = NA)
AMP_PR + 
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

######################################################################

#GlS allows to fit a 2-way-ANPVA but allows for different variance across components 

PR_m1 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_point), data = burden4)

#weighs time and diet interaction at each time point
PR_m2 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_f*diet), data = burden4)

#comparing models 
anova(PR_m1,PR_m2)

#contrasts of each time point between each diet
PR_emm<- emmeans(PR_m2, pairwise ~ diet | time_f)

