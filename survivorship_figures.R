#set working directory
setwd("C:/Users/amd439/Box/Darby_Drea/Projects/Sucrose-diet-variation-2021")

#load in required libraries
#install.packages("survival")
library(survival)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("survminer")
library(survminer)
#install.packages("coxme")
library(coxme)
#install.packages("multcomp")
library(multcomp)
#install.packages("GGally")
library(GGally)
library(ggpubr)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("dplyr")
library(dplyr)
#install.packages("emmeans")
library(emmeans)

#Providencia infection
df= read.csv("providencia120.csv", header = T)

#Lactococcus lactis infection
df1 = read.csv("lactis120ABCD.csv", header = T)

#Enterococcus faecalis  infection
df2 = read.csv("enterococcus120abcd.csv")

#Serratia marcescens 2698B infection
df3 = read.csv("serratia-abcde.csv")

##organize PBS vs infected treatments into separate data frames

##Providencia----------------------------------------------------
df$diet = as.factor(df$diet)
#split PBS from infected samples
PRinterinfection <- split(df, f = df$treatment)

#PBS injected- sterile control for all Providencia blocks
PBS_PR<- PRinterinfection$`PBS`
#Providencia infections for all blocks
infection_PR<- PRinterinfection$`Providencia`

##Lactococcus----------------------------------------------------
df1$diet = as.factor(df1$diet)

#split PBS from infected samples
LLinterinfection <- split(df1, f = df1$treatment)

#PBS injected- sterile control for all Lactoccocus blocks
PBS_LL<- LLinterinfection$`PBS`
#Providencia infections for all blocks
infection_LL<- LLinterinfection$`L. lactis`

#Enterococcus----------------------------------------------------------------------------
df2$diet = as.factor(df2$diet)
#split PBS from infected samples
EFinterinfection <- split(df2, f = df2$treatment)

#PBS injected- sterile control for all Enterococus blocks
PBS_EF<- EFinterinfection$`PBS`
#Providencia infections for all blocks
infection_EF<- EFinterinfection$`E. faecalis`

#Serratia -------------------------------------------------------------------------------
df3$diet = as.factor(df3$diet)
#split PBS from infected samples
SMinterinfection <- split(df3, f = df3$treatment)

#PBS injected- sterile control for all Serratia blocks
PBS_SM<- SMinterinfection$`PBS`
#Serratia infection
infection_SM<- SMinterinfection$`2698B`

##sort increasing order of diets---------------------------------------------------------

##Providencia
PBS_PR$diet <- factor(PBS_PR$diet,
                            levels = c("0%", "2%", "4%", "8%", "16%", "24%"))
infection_PR$diet <- factor(infection_PR$diet,
                            levels = c("0%", "2%", "4%", "8%", "16%", "24%"))

##Lactococcus
PBS_LL$diet <- factor(PBS_LL$diet,
                      levels = c("0%", "2%", "4%", "8%", "16%", "24%"))
infection_LL$diet <- factor(infection_LL$diet,
                            levels = c("0%", "2%", "4%", "8%", "16%", "24%"))

##Enterococcus
PBS_EF$diet <- factor(PBS_EF$diet,
                      levels = c("0%", "2%", "4%", "8%", "16%", "24%"))
infection_EF$diet <- factor(infection_EF$diet,
                            levels = c("0%", "2%", "4%", "8%", "16%", "24%"))

##Serratia
PBS_SM$diet <- factor(PBS_SM$diet,
                      levels = c("0%", "2%", "4%", "8%", "16%", "24%"))
infection_SM$diet <- factor(infection_SM$diet,
                            levels = c("0%", "2%", "4%", "8%", "16%", "24%"))


#Pool all PBS samples
PBS_ALL <- rbind(PBS_SM, PBS_EF, PBS_PR, PBS_LL)


#assigned color pallete
six_color = c( "#99CCFF","#3399FF", "#003399",  "#FFCC66", "#FF9966", "#CC6633")


#Fit and plots ----------------------------------------------------------------------------------------------

##PBS from all infection Blocks
fit_PBS<- survfit(Surv(hours_til_death,censor==1)~ diet, data = PBS_ALL)

plot_PBS<- ggsurvplot(fit_PBS,
                  break.time.by=24,
                  pval=F,
                  legend= "right",
                  legend.title="Sucrose Content (w/v)",
                  main="Diets",
                  legend.labs=
                    c("0% Sucrose", "2% Sucrose", "4% Sucrose", 
                      "8% Sucrose", "16% Sucrose", "24% Sucrose"),
                  ylab="Proportion Surviving",
                  xlab="Time (hrs)",
                  size = 4,
                  palette = six_color,
                  font.y = c(25, face = "bold"),
                  font.x = c(25, face = "bold"),
                  font.tickslab = c(20, "plain"),
                  font.legend = c(18),
                  risk.table = F,
                  risk.table.col = "strata",
                  risk.table.height = .5,
                  sur.plot.height = 1,
                  #linetype= "strata"
)

ggsave(
  filename,
  plot = plot_PBSOk(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  ...
)


##Providencia--------------------------------------------------------------------------------------------
fit_PR<- survfit(Surv(hours_til_death,censor==1)~ diet, data = infection_PR)
plot_PR<- ggsurvplot(fit_PR,
           break.time.by=24,
           pval=F,
           #legend= "right",
           #legend.title="Sucrose Content (w/v)",
           main="Diets",
           legend.labs=
             c("0% Sucrose", "2% Sucrose", "4% Sucrose", 
               "8% Sucrose", "16% Sucrose", "24% Sucrose"),
           ylab="Proportion Surviving",
           xlab="Time (hrs)",
           size = 4,
           palette = six_color,
           font.y = c(25, face = "bold"),
           font.x = c(25, face = "bold"),
           font.tickslab = c(20, "plain"),
           font.legend = c(18),
           risk.table = T,
           risk.table.col = "strata",
           risk.table.height = .5,
           sur.plot.height = 1,
           #linetype= "strata"
)

plot_PR
ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  ...
)
 
#cox Proportional Mixed Effects Model 
coxme_PR = coxme(Surv(hours_til_death,censor) ~ diet + (1|block), data = infection_PR)
summary(coxme_PR)

##emmeans for comparisons between diets
emm_PR <- lsmeans(coxme_PR, ~ diet)


##cld reports letters to signify significant difference across diets
cld_PR <- cld(emm_PR,
              alpha=0.05,
              Letters=letters,
              adjust="sidak")



##Lactococcus ----------------------------------------------------------------------------------------------
fit_LL<- survfit(Surv(hours_til_death,censor==1)~ diet, data = infection_LL)
 
plot_LL<- ggsurvplot(fit_LL,
                   break.time.by=24,
                   pval=F,
                   #legend= "right",
                   #legend.title="Sucrose Content (w/v)",
                   main="Diets",
                  legend.labs=
                    c("0% Sucrose", "2% Sucrose", "4% Sucrose", 
                      "8% Sucrose", "16% Sucrose", "24% Sucrose"),
                   xlab="Time (hrs)",
                   size = 4,
                   palette = six_color,
                   font.y = c(25, face = "bold"),
                   font.x = c(25, face = "bold"),
                   font.tickslab = c(20, "plain"),
                   font.legend = c(18),
                   risk.table = T,
                   risk.table.col = "strata",
                   risk.table.height = .5,
                   sur.plot.height = 1,
                   #linetype= "strata"
 )

plot_LL

ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  ...
)

#cox Proportional Hazards model
coxme_LL = coxme(Surv(hours_til_death,censor) ~ diet + (1|block), data = infection_LL)
summary(coxme_LL)
emm_LL <- lsmeans(coxme_LL, ~ diet)
cld_LL <- cld(emm_LL,
               alpha=0.05,
               Letters=letters,
               adjust="sidak")


##Enterococcus---------------------------------------------------------------------------------------------
fit_EF<- survfit(Surv(hours_til_death,censor==1)~ diet, data = infection_EF)

plot_EF<- ggsurvplot(fit_EF,
                  break.time.by=24,
                  pval=F,
                  #legend= "right",
                  #legend.title="Sucrose Content (w/v)",
                  main="Diets",
                  legend.labs=
                    c("0% Sucrose", "2% Sucrose", "4% Sucrose", 
                      "8% Sucrose", "16% Sucrose", "24% Sucrose"),
                  xlab="Time (hrs)",
                  size = 4,
                  palette = six_color,
                  font.y = c(25, face = "bold"),
                  font.x = c(25, face = "bold"),
                  font.tickslab = c(20, "plain"),
                  font.legend = c(18),
                  risk.table = T,
                  risk.table.col = "strata",
                  risk.table.height = .5,
                  sur.plot.height = 1,
                  #linetype= "strata"
)
plot_EF
ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  ...
)

# Enterococcus cox Proportional Hazards model
coxme_EF = coxme(Surv(hours_til_death,censor) ~ diet + (1|block), data = infection_EF)
summary(coxme_EF)

emm_EF <- lsmeans(coxme_EF, ~ diet)
cld_EF <- cld(emm_EF,
               alpha=0.05,
               Letters=letters,
               adjust="sidak")

##Serratia ----------------------------------------------------------------------------------------------
fit_SM<- survfit(Surv(hours_til_death,censor==1)~ diet, data = infection_SM)

plot_SM<- ggsurvplot(fit_SM,
                     break.time.by=24,
                     pval=F,
                     #legend= "right",
                     #legend.title="Sucrose Content (w/v)",
                     main="Diets",
                     legend.labs=
                       c("0% Sucrose", "2% Sucrose", "4% Sucrose", 
                         "8% Sucrose", "16% Sucrose", "24% Sucrose"),
                     xlab="Time (hrs)",
                     size = 4,
                     palette = six_color,
                     font.y = c(25, face = "bold"),
                     font.x = c(25, face = "bold"),
                     font.tickslab = c(20, "plain"),
                     font.legend = c(18),
                     risk.table = T,
                     risk.table.col = "strata",
                     risk.table.height = .5,
                     sur.plot.height = 1,
                     #linetype= "strata"
)

plot_SM
ggsave(
  filename,
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  ...
)
  
##S. marcescens cox proportional hazards model
coxme_SM = coxme(Surv(hours_til_death,censor) ~ diet + (1|block), data =infection_SM)
summary(coxme_SM)
  
emm_SM <- lsmeans(coxme_SM, ~ diet)
cld_SM <- cld(emm_SM,
                alpha=0.05,
                Letters=letters,
                adjust="sidak")


## proportion of flies survived across all 4 pathogens on day 5

prop = read.csv("pathogens_proportion_survived.csv")

prop$diet = as.factor(prop$diet)

prop$diet <- factor(prop$diet,
       levels = c("0%", "2%", "4%", "8%", "16%", "24%"))

day_5 <- ggplot(prop, aes(x = diet, y =day_5_prop, group = pathogen))+
  geom_line(aes(color = pathogen), size = 1.25)+
  geom_point(aes(color =pathogen), size =3) +
  labs(y= "Proportion Survived", x = "Diet")+
  theme(axis.line=element_line(size=2))+ ggtitle("Proportion Survived: Day 5")+
  ylim(0,1)+
  theme_classic()

day_5

day_1<- ggplot(prop, aes(x = diet, y =day_1_proportion, group = pathogen))+
  geom_line(aes(color = pathogen), size = 1.25)+
  geom_point(aes(color =pathogen), size =3) +
  labs(y= "Proportion Survived", x = "Diet")+
  theme(axis.line=element_line(size=2))+ ggtitle("Proportion Survived: Day 1")+
  ylim(0,1)+
  theme_classic()

day_1

day_2<- ggplot(prop, aes(x = diet, y =day_2_prop, group = pathogen))+
  geom_line(aes(color = pathogen), size = 1.25)+
  geom_point(aes(color =pathogen), size =3) +
  labs(y= "Proportion Survived", x = "Diet")+
  theme(axis.line=element_line(size=2))+ ggtitle("Proportion Survived: Day 2")+
  ylim(0,1)+
  theme_classic()

day_2

