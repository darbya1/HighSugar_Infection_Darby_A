
#infection survival done on flies infected with Providencia rettgeri fed on 2% and 16% sucrose diets
df= read.csv("axenic_survivorship.csv", header = T)


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
library(emmeans)

#made diet a factor 
df$diet = as.factor(df$diet)

#use split function to split the infected flies from PBS
Allinterinfection <- split(df, f = df$treatment)

#PBS injected- sterile control
PBS <- Allinterinfection$`PBS`

#Providencia infection
Providencia<- Allinterinfection$`Providencia`


################## Figure making #######################################################

#color for plot 2%= blue 16% =orange 
two_diet = c( "#3399FF", "#FF9966" )

#sort increasing order of diets
Providencia$diet <- factor(Providencia$diet,
                            levels = c( "2%", "16%"))


### survival plots 
#Fit and plots for providencia infected axenic flies

fit<- survfit(Surv(hours_til_death,censor==1) ~  diet, data = Providencia)

ggsurvplot(fit,
           break.time.by=24,
           pval=F,
           legend="bottom",
           legend.title="Sucrose (w/v)",
           main="Diets",
           legend.labs=c("2%", "16%"),
           ylab="Proportion Surviving",
           xlab= "Time (in hours)",
           size = 2.5,
           font.y = c(25, face = "bold"),
           font.x = c(25, face = "bold"),
           font.tickslab = c(20, "plain"),
           font.legend = c(18),
           risk.table = F,
           palette = two_diet,
           risk.table.col = "strata",
           risk.table.height = .5,
           sur.plot.height = 1,
          
           
)


#Cox mixed effects modeling testing main effect of diet on mortality in axenic flies

pro.coxme = coxme(Surv(hours_til_death,censor) ~ diet + (1|block), data = Providencia)

#pairwise comparisons
pro.emm<-emmeans(pro.coxme, pairwise ~ diet)

pro.emm


