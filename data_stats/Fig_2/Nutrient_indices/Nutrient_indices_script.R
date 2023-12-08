##barplots + stats for nutrient indices in uninfected and providencia/serratia infected flies 24-hpi

##Set working directory
setwd("C:/Users/amd439/Box/Darby_Drea/Projects/Sucrose-diet-variation-2021/nutrient_stores/Thesis Data")

#load in required libraries
library(ggplot2)
library(emmeans)
library(lsmeans)
library(multcompView)
#CLD function
library(multcomp)
library(emmeans)
library(ggplot2)
library(scales)
library(ggrepel)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)


#load in csv files
df = read.csv("pooled_blocks.csv", header = T)

#relevel diet 
df$Diet <-as.factor(df$Diet)
df$Diet <- factor(df$Diet, levels = c("2%", "16%"))

#assign color for plots
plot_color = c( "#3399FF", "#FF9966" )


#glucose barplot
glucose_bar <- ggplot(data = df, aes(x=Infection, y= glucose_std_2, fill = Diet ))+
  geom_bar(stat="summary", fun.y = "mean", position = "dodge",colour = "black", size = 1) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position= position_dodge(.90), width = .4, size = .9)+
  #geom_point(aes(x = Diet, color = "black"), position = "jitter")+
  scale_fill_manual(values = plot_color)+
  scale_color_manual(values = "black")+
  labs(y = "Glucose/Protein ratio per fly", x = "Infection") +
  theme_classic()

#trehalose barplot
trehalose_bar <- ggplot(data = df, aes(x=Infection, y= trehalose.standard, fill = Diet ))+
  geom_bar(stat="summary", fun.y = "mean", position = "dodge",colour = "black", size = 1) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position= position_dodge(.90), width = .4, size = .9)+
 # geom_point(aes(x = Diet, color = "black"), position = "jitter")+
  scale_fill_manual(values = plot_color)+
  scale_color_manual(values = "black")+
  labs(y = "Trehalose/Protein ratio per fly",x = "Infection)") +
  theme_classic()


#glycogen barplot
glycogen_bar <- ggplot(data = df, aes(x= Infection, y= glycogen.standard, fill = Diet ))+
  geom_bar(stat="summary", fun.y = "mean", position = "dodge",colour = "black", size = 1) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position= position_dodge(.90), width = .4, size = .9)+
 # geom_point(aes(x = Diet, color = "black"), position = "jitter")+
  scale_fill_manual(values = plot_color)+
  labs(y = "Glycogen/protein ratio per fly", x = "Infection") +
  scale_color_manual(values = "black")+
  theme_classic()

###barplot outputs...................................................................................

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




####### stats

##general model lmer(nutrient~ Diet * Infection + Diet + (1|Block), data = df)

## glucose   #############################   

glu_lmer<- lmer(glucose_std_2~ Diet * Infection + Diet + (1|Block), data = df)

#posthoc pairwise comparisons across all infection and diet treatments
glumer_emm<-emmeans(glu_lmer, pairwise ~ Diet: Infection)

#posthoc pairwise comparisons between diet within infection
glumer_emm1<-emmeans(glu_lmer, pairwise ~ Diet|Infection)
glumer_emm1

#posthoc pairwise comparisons between infection within diet
glumer_emm2 <- emmeans(glu_lmer, pairwise ~ Infection|Diet)
glumer_emm2

####### trehalose stats  #############################################


#linear mixed effects regression
tre.lmer  <- lmer(trehalose.standard~  Diet * Infection + Diet + (1|Block)  , data = df)
tre.lm <-lm(trehalose)
#gives contrasts for each condition
tre_emm1<- emmeans(tre.lmer, pairwise ~ Infection : Diet)
tre_emm1


#pairwise contrast of infection within diet
tre_emm2<- emmeans(tre.lmer, pairwise ~ Infection | Diet)

tre_emm2

#pairwise contrast of diet within infection
tre_emm3 <- emmeans(tre.lmer, pairwise ~ Diet | Infection)
tre_emm3


########### glycogen  stats 

gly.lmer  <- lmer(glycogen.standard~ Diet * Infection + Diet + (1|Block)  , data = df)

#gives contrasts for each condition
gly_emm1<- emmeans(gly.lmer, pairwise ~ Infection : Diet)
gly_emm1

#pairwise contrast of infection within diet
gly_emm2<- emmeans(gly.lmer, pairwise ~ Infection | Diet)
gly_emm2

#pairwise contrast of diet within infection
gly_emm3 <- emmeans(gly.lmer, pairwise ~ Diet | Infection)
gly_emm3
