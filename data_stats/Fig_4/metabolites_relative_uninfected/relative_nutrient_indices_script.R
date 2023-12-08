#relative metabolite levels to uninfected flies

setwd("C:/Users/amd439/Box/Darby_Drea/")

#load in libraries
library(ggplot2)
library(lattice)
library(Hmisc)
library(ggbeeswarm)
library(multcomp)
library(lme4)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(emmeans)
library(lme4)


#load in csv file
df=read.csv("differences_by_block.csv", header= T)

#assign levels for diet
df$Diet <-as.factor(df$Diet)
df$Diet <- factor(df$Diet, levels = c("2%", "16%"))

#colors for plot; 2% = blue; 16% = orange
plot_color = c( "#3399FF", "#FF9966" )

#barplots


#relative glucose 
glucose_bar <- ggplot(data = df, aes(x=Infection, y= glucose_diff, fill = Diet ))+
  geom_bar(stat="summary", fun.y = "mean", position = "dodge",colour = "black", linewidth = 1) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position= position_dodge(.90), width = .4, size = .9)+
  #geom_point(aes(x = Diet, color = "black"), position = "jitter")+
  scale_fill_manual(values = plot_color)+
  scale_color_manual(values = "black")+
  labs(y = "Glucose relative to uninfected", x = "Infection") +
  theme_classic()

#relative trehalose 
trehalose_bar <- ggplot(data = df, aes(x=Infection, y= trehalose_diff, fill = Diet ))+
  geom_bar(stat="summary", fun.y = "mean", position = "dodge",colour = "black", size = 1) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position= position_dodge(.90), width = .4, size = .9)+
 # geom_point(aes(x = Diet, color = "black"), position = "jitter")+
  scale_fill_manual(values = plot_color)+
  scale_color_manual(values = "black")+
  labs(y = "Trehalose relative to uninfected",x = "Infection)") +
  theme_classic()


#relative glycogen 
glycogen_bar <- ggplot(data = df, aes(x= Infection, y= glycogen_diff, fill = Diet ))+
  geom_bar(stat="summary", fun.y = "mean", position = "dodge",colour = "black", size = 1) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position= position_dodge(.90), width = .4, size = .9)+
 # geom_point(aes(x = Diet, color = "black"), position = "jitter")+
  scale_fill_manual(values = plot_color)+
  labs(y = "Glycogen relative to uninfected", x = "Infection") +
  scale_color_manual(values = "black")+
  theme_classic()

#######relative change bar_plots....................................................................................

glucose_bar + theme( #changes font size and face of the title of plot
  #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
  axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
  axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
  axis.text.x = element_text(size = 18, family = "sans"), 
  axis.text.y = element_text(size = 18, family = "sans"),
  legend.text = element_text(size=20),
  #changes tile of legend size and face
  legend.title = element_text(size=25, face = 2))

trehalose_bar + theme( #changes font size and face of the title of plot
  #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
  axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
  axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
  axis.text.x = element_text(size = 18, family = "sans"), 
  axis.text.y = element_text(size = 18, family = "sans"),
  legend.text = element_text(size=20),
  #changes tile of legend size and face
  legend.title = element_text(size=25, face = 2))

glycogen_bar + theme( #changes font size and face of the title of plot
  #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
  axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
  axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
  axis.text.x = element_text(size = 18, family = "sans"), 
  axis.text.y = element_text(size = 18, family = "sans"),
  legend.text = element_text(size=20),
  #changes tile of legend size and face
  legend.title = element_text(size=25, face = 2))

#######Stats########

## linear mixed effects model lmer(nutrient~ Diet * Infection + Diet + (1|Block), data = df)

##glucose   #############################   

glu_lmer1<- lmer(glucose_diff ~ Diet * Infection + Diet + (1|Block), data = df)

#pairwise comparisons across all treatments
glumer_emm<-emmeans(glu_lmer1, pairwise ~ Diet: Infection)
glumer_emm

#pairwise comparison between diet within infection
glumer_emm1<-emmeans(glu_lmer1, pairwise ~ Diet|Infection)
glumer_emm1

#pairwise comparisons between infection within diet
glumer_emm2 <- emmeans(glu_lmer1, pairwise ~ Infection|Diet)
glumer_emm2

#######trehalose  

tre.lmer  <- lmer(trehalose_diff~  Diet * Infection + Diet + (1|Block), data = df)

#gives contrasts for each condition
tre_emm1<- emmeans(tre.lmer, pairwise ~ Infection : Diet)
tre_emm1

#pairwise contrast of infection within diet
tre_emm2<- emmeans(tre.lmer, pairwise ~ Infection | Diet)
tre_emm2

#pairwise contrast of diet within infection
tre_emm3 <- emmeans(tre.lmer, pairwise ~ Diet | Infection)
tre_emm3


###########glycogen 

gly.lmer  <- lmer(glycogen.standard~ Diet * Infection + Diet + (1|Block)  , data = df)

#gives contrasts for each condition
gly_emm1<- emmeans(gly.lm, pairwise ~ Infection : Diet)

gly_emm1

#pairwise contrast of infection within diet
gly_emm2<- emmeans(gly.lm, pairwise ~ Infection | Diet)

gly_emm2

#pairwise contrast of diet within infection
gly_emm3 <- emmeans(gly.lm, pairwise ~ Diet | Infection)

gly_emm3