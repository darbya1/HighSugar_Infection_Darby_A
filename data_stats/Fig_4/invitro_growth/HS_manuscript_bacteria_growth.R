#measure effect on addition of casamino acids and sugar (glucose or sucrose) on final OD600 at 18 hours of Providencia rettgeri and Serratia marcescens

#set working directory
setwd("C:/Users/amd439/Box/Darby_Drea/Projects/Sucrose-diet-variation-2021/Bacteria_growth_curves")

library(emmeans)
library(multcomp)
library(lme4)

##load csv file
OD_18 <- read.csv("PR_SM_18hr_OD_all_blocks.csv", header = T)

#create data frame of only Providencia and only Serratia
PR <- subset(OD_18, OD_18$bacteria =="Providencia")
SM <- subset(OD_18, OD_18$bacteria !="Providencia")

##################Providencia

# lmer(OD600 at 18hrs ~ sugar * CAS + sugar +(1|block))

PR.lmer <- lmer(OD600_18hrs ~ sugar* CAS + sugar + (1|block), data = PR)

#pairwise comparisons between all media
PR.emm <-emmeans(PR.lmer, pairwise ~ sugar:CAS)
PR.emm
#compact letter display of difference in p-values
cld(PR.emm)

#pairwise comparison between sugar within media that have CAS present or absent
PR1.emm1 <- emmeans(PR.lmer, pairwise~ sugar|CAS)
PR1.emm1


##################Serratia 

SM.lmer <- lmer(OD600_18hrs ~ sugar*CAS+sugar +(1|block), data = SM)

#pairwise comparisons between all media
SM.emm <-emmeans(SM.lmer,pairwise ~ sugar:CAS)
cld(SM.emm)

#pairwise comaprisons between sugar within CAS
SM.emm1<- emmeans(SM.lmer, pairwise ~ sugar|CAS)
SM.emm1 
