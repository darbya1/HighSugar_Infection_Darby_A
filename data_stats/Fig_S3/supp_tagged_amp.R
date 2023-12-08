library(lattice)
library(Hmisc)
library(ggbeeswarm)
library(multcomp)
library(car)
library(lme4)
#likely one of the only packages needed
library(ggplot2)
library(tidyverse)
library(dplyr)
library(emmeans)
library(lme4)
library(dplyr)
#set working directory


amp<- read.csv("supp_tagged_AMP.csv", header = T)

amp$infection <-  factor(amp$infection,
                         levels = c("uninfected", "PBS", "E. cloacae"))
dro<- subset(amp, amp$AMP == "Dro")
cec<- subset(amp, amp$AMP != "Dro")


##### linear model testing difference in tagged drosocin 
lm_dro<- lm(Mr_peptide ~ infection, data = dro)


#difference between infection treatments
emm_dro <- emmeans(lm_dro, pairwise~infection)
emm_dro


dro_box<-ggplot(dro, aes(x=infection, y= Mr_peptide, fill = infection))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(pch = 21, position = position_jitterdodge())+
  scale_fill_manual(values = c("uninfected" = "white", "PBS" = "lightgrey", "E. cloacae"= "darkgrey"))+
  scale_shape_manual(values = c(15, 16))+
  labs(y="Tagged Drosocin (Mr/fly)", x = "Infection", fill = "Infection")+
  theme_classic()


dro_box + theme( #changes font size and face of the title of plot
  #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
  axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
  axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
  axis.text.x = element_text(size = 18, family = "sans"), 
  axis.text.y = element_text(size = 18, family = "sans"),
  legend.text = element_text(size=20),
  #changes tile of legend size and face
  legend.title = element_text(size=25, face = 2))

##### linear model testing difference in tagged cecropinA1 
lm_cec<- lm(Mr_peptide ~ infection, data = cec)


#difference between infection treatments
emm_cec <- emmeans(lm_cec, pairwise~infection)
emm_cec


cec_box<-ggplot(cec, aes(x=infection, y= Mr_peptide, fill = infection))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(pch = 21, position = position_jitterdodge())+
  scale_fill_manual(values = c("uninfected" = "white", "PBS" = "lightgrey", "E. cloacae"= "darkgrey"))+
  scale_shape_manual(values = c(15, 16))+
  labs(y="Tagged Cecropin (Mr/fly)", x = "Infection", fill = "Infection")+
  theme_classic()


cec_box + theme( #changes font size and face of the title of plot
  #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
  axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
  axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
  axis.text.x = element_text(size = 18, family = "sans"), 
  axis.text.y = element_text(size = 18, family = "sans"),
  legend.text = element_text(size=20),
  #changes tile of legend size and face
  legend.title = element_text(size=25, face = 2))
