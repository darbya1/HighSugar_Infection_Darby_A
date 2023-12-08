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


#Drosocin pooled blocks without background abs
dro= read.csv("pooled_drosocin.csv", header = T)

#make the numeric # of diet as a factor in the data frame
dro$diet = as.factor(dro$diet)

#assign increasing order of diets
dro$diet <- factor(dro$diet,
                  levels = c("2%", "16%"))

##########################################################################################33

#split Serratia and Providencia infections from each other due to flies being separate experimental blocks
dro_PR <- subset(dro, dro$Block == "A"| dro$Block == "B"| dro$Block == "C")
dro_SM <- subset(dro, dro$Block == "E"| dro$Block == "F"| dro$Block == "G")

##GLM of binomial data to test presence absence of flies with tagged peptide
droPR.glb <- glm(LR ~ diet + infection*diet, data = dro_PR, family = binomial)
summary(droPR.glb)

droPR_emm <- emmeans(droPR.glb, pairwise ~ diet|infection)
droPR_emm

############################################################33
##figure of probability of drosocin appearing in infected P. rettgeri flies

df <- with(dro_PR,
           expand.grid(diet = unique(diet), infection = unique(infection)))

fits <- predict(droPR.glb, newdata = df, se.fit = TRUE)

df# Get odds from modrl
df$probability <- exp(fits$fit) 
df$plus <- exp(fits$fit + 1.96 * fits$se.fit)
df$minus <- exp(fits$fit - 1.96 * fits$se.fit)

# Convert odds to probabilities
df$probability <- df$probability / (1 + df$probability)
df$plus <- df$plus / (1 + df$plus)
df$minus <- df$minus / (1 + df$minus)


# Plot probabilities
droPR_plot<- ggplot(df, aes(infection, probability)) +
  geom_errorbar(aes(ymin = minus, ymax = plus, colour = diet), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(aes(fill = diet), shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  labs(x = "Infection", y = "Probability of AMP Present")+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  theme_light(base_size = 16) +
  scale_y_continuous(name = "Probability of peptide present", limits = c(0, 1),
                     labels = scales::percent)


droPR_plot +
  theme( #changes font size and face of the title of plot
    #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
    axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
    axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
    axis.text.x = element_text(size = 18, family = "sans"), 
    axis.text.y = element_text(size = 18, family = "sans"),
    legend.text = element_text(size=20),
    #changes tile of legend size and face
    legend.title = element_text(size=25, face = 2))



##### linear mixed effects model looking at difference of Dro peptide between infection and diet
lmer_DPR<- lmer(molar_mass_dro ~ diet * infection + diet +(1|Block), data = dro_PR)


#difference between diet within infection
emm_DPR <- emmeans(lmer_DPR, pairwise~diet|infection)
emm_DPR

#difference between diet within infection
emm_DPR1 <- emmeans(lmer_DPR, pairwise~infection|diet)
emm_DPR1

emm_DPR2 <- emmeans(lmer_DPR,~infection:diet)
cld(emm_DPR2)


dro_box<-ggplot(dro_PR, aes(x=infection, y= Mr_dro, fill= diet))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = Diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  #scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  labs(y="Tagged Drosocin (Mr/fly)", x = "Infection", fill = "Sucrose (w/v)")+
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



################ Serratia-dro           #############################################################3

##GLM of binomial data to test Presence absence of flies with tagged peptide
droSM.glm1 <- glm(LR ~ infection + diet, data = dro_SM, family = binomial )
summary(droSM.glm1)


##figure of probabilities 
df1 <- with(dro_SM,
           expand.grid(diet = unique(diet), infection = unique(infection)))

fits <- predict(droSM.glm1, newdata = df1, se.fit = TRUE)

# Get odds from model
df1$probability <- exp(fits$fit) 
df1$plus <- exp(fits$fit + 1.96 * fits$se.fit)
df1$minus <- exp(fits$fit - 1.96 * fits$se.fit)

# Convert odds to probabilities
df1$probability <- df1$probability / (1 + df1$probability)
df1$plus <- df1$plus / (1 + df1$plus)
df1$minus <- df1$minus / (1 + df1$minus)


# Plot Probabilities
droSM_plot<- ggplot(df1, aes(infection, probability)) +
  geom_errorbar(aes(ymin = minus, ymax = plus, colour = diet), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(aes(fill = diet), shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  labs(x = "Infection", y = "Probability of AMP Present")+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  theme_light(base_size = 16) +
  scale_y_continuous(name = "Probability of peptide present", limits = c(0, 1),
                     labels = scales::percent)


droSM_plot +
  theme( #changes font size and face of the title of plot
    #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
    axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
    axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
    axis.text.x = element_text(size = 18, family = "sans"), 
    axis.text.y = element_text(size = 18, family = "sans"),
    legend.text = element_text(size=20),
    #changes tile of legend size and face
    legend.title = element_text(size=25, face = 2))

#####linear mixed effects model to test difference of drosocin

lmer_DSM <- lmer(Mr_dro ~ diet *infection + diet + (1|Block), data = dro_SM)

summary(lmer_DSM)

#pairwise comparisons between diet within infection
emm_DSM <- emmeans(lmer_DSM, pairwise~ diet|infection)
emm_DSM


#pairwise comparisons within infection between diet
emm_DSM1 <- emmeans(lmer_DSM, pairwise ~ infection|diet)
emm_DSM1


SMdro_box<-ggplot(dro_SM, aes(x=infection, y= Mr_dro, fill= diet))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = Diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  #scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  labs(y="Tagged Drosocin (Mr/fly)", x = "Infection", fill = "Sucrose (w/v)")+
  #scale_y_continuous(limits = c(0, 2500))+
  theme_classic()


SMdro_box + theme( #changes font size and face of the title of plot
  #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
  axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
  axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
  axis.text.x = element_text(size = 18, family = "sans"), 
  axis.text.y = element_text(size = 18, family = "sans"),
  legend.text = element_text(size=20),
  #changes tile of legend size and face
  legend.title = element_text(size=25, face = 2))



#########################################################################################################3
#cecropin pooled sampels
cec =read.csv("Cecropin_pooled_blocks.csv", header = T)

#make the numeric # of diet as a factor in the data frame
cec$diet = as.factor(cec$diet)

#assign increasing order of diets
cec$diet <- factor(cec$diet,
                   levels = c("2%", "16%"))

##sort infected samples into their own data frame
cec_infect <- subset(cec, cec$infection != "naive"   )

## Binomial GLM testing probabilities of flies that have tagged cecropin A1 peptide present 
cec.glm <- glm(LR~diet*infection+diet, data = cec_infect, family = binomial)
cec.gemm <-emmeans(cec.glm, pairwise ~ diet|infection)
cec.gemm




##figure of probabilities 
df2 <- with(cec_infect,
            expand.grid(diet = unique(diet), infection = unique(infection)))

fits2 <- predict(cec.glm, newdata = df2, se.fit = TRUE)

# Get odds from model
df2$probability <- exp(fits2$fit) 
df2$plus <- exp(fits2$fit + 1.96 * fits2$se.fit)
df2$minus <- exp(fits2$fit - 1.96 * fits2$se.fit)

# Convert odds to probabilities
df2$probability <- df2$probability / (1 + df2$probability)
df2$plus <- df2$plus / (1 + df2$plus)
df2$minus <- df2$minus / (1 + df2$minus)


# Plot Probabilities
cec_plot<- ggplot(df2, aes(infection, probability)) +
  geom_errorbar(aes(ymin = minus, ymax = plus, colour = diet), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(aes(fill = diet), shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  labs(x = "Infection", y = "Probability of AMP Present")+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  theme_classic() +
  scale_y_continuous(name = "Probability of peptide present", limits = c(0, 1),
                     labels = scales::percent)


cec_plot +
  theme( #changes font size and face of the title of plot
    #plot.title=element_text(hjust = 0.5,size=14,face="bold"),
    axis.title.x = element_text(size = 20, family= "sans", face = "bold", margin = margin(t = 20, r = 20, b = 00, l = 0)), 
    axis.title.y = element_text(size = 18, family = "sans", face = "bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), 
    axis.text.x = element_text(size = 18, family = "sans"), 
    axis.text.y = element_text(size = 18, family = "sans"),
    legend.text = element_text(size=20),
    #changes tile of legend size and face
    legend.title = element_text(size=25, face = 2))

##glm testing difference in CecA1 levels in samples with detecable CecA1

#sort samples with detectable CecA1
cec_LR <- subset(cec_infect, cec_infect$LR == 1)

cec.in.glm1 <- glm(Mr_cec ~ diet * infection +diet, data = cec_LR)
summary(cec.in.glm1)

#pairwise comparisons between diet within infections
cec.in.emm1 <- emmeans(cec.in.glm1, pairwise ~ diet|infection)
cec.in.emm1


cec_box<-ggplot(cec_LR, aes(x=infection, y= Mr_cec, fill= diet))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(pch = 21, aes(shape = Diet), position = position_jitterdodge())+
  scale_fill_manual(values = c("2%" = "dodgerblue2", "16%" = "orange2"))+
  scale_shape_manual(values = c(15, 16))+
  #scale_color_manual(values = c("2%" = "blue2", "16%" = "orange2"))+
  labs(y="Tagged Cecropin (Mr/fly)", x = "Infection", fill = "Sucrose (w/v)")+
  theme_classic()


cec_box + #ggtitle(Drosocin_title)+
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



