#pathogen trajectory assay done with Providencia retgerri and Serratia marcescens 2698B fed on 2% or 16% diets

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
library(reshape)
library(tibble)
library(ggpmisc)

#set working directory
setwd("C:/Users/amd439/Box/Darby_Drea/Projects/Sucrose-diet-variation-2021/pathogen-trajectory/")


#CS_trajectory: Providencia
burden <- read.csv("CS-Path-trajectory-16-vs-2-sucrose/CS-trajectory-providencia.csv", header = T)

#CS_trajectory: Serratia marcescens 2698B
burden2 <- read.csv("CS-Path-trajectory-16-vs-2-sucrose/CS-trajectory-serratia-2698B.csv", header = T)


###Providencia rettgeri in CS wildtype
#make diet and hour post infection characters, so that barplot can be made appropriately
burden$diet <- as.factor(burden$diet)
burden$dilution <- as.factor(burden$dilution)

#reorder diet levels
burden$diet <- factor(burden$diet, levels=c("2%", "16%"), labels=c("2%", "16%"))

#adjust order of time series
burden$time_f<- factor(burden$time_point,
                    levels = c("0", "2", "4", "6","8", "10", "12","14", "16", "24", "36", "48"))
burden$log2CFU <- log2(burden$CFU_fly)
##code for boxplot + geom_point position with jitter in geom_point

CS_PR<- ggplot(data=burden, aes(x=time_f, y = log2CFU, fill = diet)) +
  geom_boxplot(aes(fill = diet), outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  #scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  labs(y="Log2 CFU per Fly", x = "Time (hour)", fill = "Sucrose (w/v)")+
  scale_y_continuous(limits = c(10, 20))+ 
  theme_classic()
CS_PR


#removes outlier points geom_boxplot(outlier.shape = NA)
CS_PR + 
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
##Generalized least squares (GLS) requires nlme package also
#install.packages("nlme")
library(nlme)
library(emmeans)

table(burden$diet, burden$time_f)

#GlS allows to fit a 2-way-ANPVA but allows for different variance across components 

m3 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_point), data = burden)

#weighs time and diet interaction at each time point
m4 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_f*diet), data = burden)

summary(m4)

#comparing models 
anova(m3, m4)

#contrasts of each time point between each diet
emmeans(m4, pairwise ~ diet | time_f)


###Serratia marcescens in CS wildtype
#make diet and hour post infection characters, so that barplot can be made appropriately
burden2$diet <- as.factor(burden2$diet)
burden2$dilution <- as.factor(burden2$dilution)

#reorder diet levels
burden2$diet <- factor(burden2$diet, levels=c("2%", "16%"), labels=c("2%", "16%"))

#adjust order of time series
burden2$time_f<- factor(burden2$time_point,
                       levels = c("0", "2", "4", "6","8", "10", "12","14", "16", "24", "36", "48"))
burden2$log2CFU <- log2(burden2$CFU_fly)
##code for boxplot + geom_point position with jitter in geom_point

CS_SM<- ggplot(data=burden2, aes(x=time_f, y = log2CFU, fill = diet)) +
  geom_boxplot(aes(fill = diet), outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  #scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  labs(y="Log2 CFU per Fly", x = "Time (hour)", fill = "Sucrose (w/v)")+
  scale_y_continuous(limits = c(10, 20))+ 
  theme_classic()
CS_SM


#removes outlier points geom_boxplot(outlier.shape = NA)
CS_SM + 
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
##Generalized least squares (GLS) requires nlme package also

table(burden2$diet, burden2$time_f)

#GlS allows to fit a 2-way-ANPVA but allows for different variance across components 
m3 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_point), data = burden2)

#weighs time and diet interaction at each time point
m4 = gls(log2(CFU_fly +1) ~ time_f* diet, weights = varIdent(form = ~1|time_f*diet), data = burden2)

summary(m4)

#comparing models 
anova(m3, m4)

#contrasts of each time point between each diet
emmeans(m4, pairwise ~ diet | time_f)
