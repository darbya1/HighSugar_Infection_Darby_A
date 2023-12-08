###infection-induced AMP gene expression

# Packages
library(multcompView)
library(multcomp)
library(emmeans)
library(ggplot2)
library(scales)
library(ggrepel)
library(dplyr)
library(tidyr)

#setwd



##all serratia + providencia samples
df = read.csv("official_amp_qpcr.csv", header = T)

#
#assign levels to diet
df$diet <- factor(df$diet,
                  levels = c("2%", "16%"))

serratia <- subset(df, df$infection =="2698B" |df$infection == "uninfected")
providencia <-subset(df,df$infection =="providencia" |df$infection == "uninfected"  )

##Providencia########################################################################################################
#### Simple model for each AMP  ####

#Y = u + act_ct + Diet + infection + Diet x infection

# Dpt Model##########################################################################################################
Dpt_LS <- lm(dpt_ct ~ act_ct  + diet + infection * diet, data= providencia, na.action = na.omit)
anova(Dpt_LS)
summary(Dpt_LS)


# Extract Least Squared Means from model
marginal.1 = lsmeans(Dpt_LS, ~ diet:infection) # extracts least square means

CLD_dpt<- cld(marginal.1,
                alpha=0.05,
                Letters=letters,
                adjust="sidak") # compact letter display and correction

CLD_dpt

# Subtracts LS Means: Infected-Uninfected
two_dpt_pairs<-update(pairs(marginal.1, by= "diet"), by=NULL)
summary(two_dpt_pairs)

#extracts contrasts between infected and uninfected stored into a dataframe
two_pairs_dpt_CLD<- cld(two_dpt_pairs,
                          alpha=0.05,
                          Letters=letters,
                          adjust="sidak")
two_pairs_dpt_CLD




# AttA Model######################################################################################################################################
Att_LS <- lm(atta_ct ~ act_ct  + diet + infection * diet, data= providencia, na.action = na.omit)
anova(Att_LS)
summary(Att_LS)


# Extract Least Squared Means from model
marginal.2 = lsmeans(Att_LS, ~ diet:infection) # extracts least square means

CLD_att<- cld(marginal.2,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

CLD_att

# Subtracts LS Means: Infected-Uninfected
two_att_pairs<-update(pairs(marginal.2, by= "diet"), by=NULL)
summary(two_att_pairs)

two_pairs_att_CLD<- cld(two_att_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
two_pairs_att_CLD



###################### Drosocin Model ##################################################################
Dro_LS <- lm(dro_ct ~ act_ct  + diet + infection * diet, data= providencia, na.action = na.omit)
anova(Dro_LS)
summary(Dro_LS)

  

# Extract Least Squared Means from model
marginal.3 = emmeans(Dro_LS, ~ diet:infection) # extracts least square means

CLD_dro<- cld(marginal.3,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

CLD_dro



# Subtracts LS Means: Infected-Uninfected
two_dro_pairs<-update(pairs(marginal.3, by= "diet"), by=NULL)
summary(two_dro_pairs)

two_pairs_dro_CLD<- cld(two_dro_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
two_pairs_dro_CLD


# Defensin Model#############################################################################################
Def_LS <- lm(def_ct ~ act_ct  + diet + infection * diet, data= providencia, na.action = na.omit)
anova(Def_LS)
summary(Def_LS)



# Extract Least Squared Means from model
marginal.4 = lsmeans(Def_LS, ~ diet:infection) # extracts least square means

CLD_def<- cld(marginal.4,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

CLD_def


# Subtracts LS Means: Infected-Uninfected
two_def_pairs<-update(pairs(marginal.4, by= "diet"), by=NULL)
summary(two_def_pairs)

two_pairs_def_CLD<- cld(two_def_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
two_pairs_def_CLD




# Cecropin Model###################################################################################
Cec_LS <- lm(cec_ct ~ act_ct  + diet + infection * diet, data= providencia, na.action = na.omit)
anova(Cec_LS)
summary(Cec_LS)

# Extract Least Squared Means from model
marginal.6 = lsmeans(Cec_LS, ~ diet:infection) # extracts least square means

CLD_cec<- cld(marginal.6,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

CLD_cec


# Subtracts LS Means: Infected-Uninfected
two_cec_pairs<-update(pairs(marginal.6, by= "diet"), by=NULL)
summary(two_cec_pairs)

two_pairs_cec_CLD<- cld(two_cec_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
two_pairs_cec_CLD



# Drosmycin Model
Drs_LS <- lm(drs_ct ~ act_ct  + diet + infection * diet, data= providencia, na.action = na.omit)
anova(Drs_LS)
summary(Drs_LS)


# Extract Least Squared Means from model
marginal.5 = lsmeans(Drs_LS, ~ diet:infection) # extracts least square means

CLD_drs<- cld(marginal.5,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

CLD_drs



# Subtracts LS Means: Infected-Uninfected
two_drs_pairs<-update(pairs(marginal.5, by= "diet"), by=NULL)
summary(two_drs_pairs)

two_pairs_drs_CLD<- cld(two_drs_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
two_pairs_drs_CLD


###############Serratia ##########################################################################################################
#### Simple model for each AMP  ####

#Y = u + act_ct + Diet + infection + Diet x infection

# Dpt Model
Dpt_LS1 <- lm(dpt_ct ~ act_ct  + diet + infection * diet, data= serratia, na.action = na.omit)
anova(Dpt_LS1)
summary(Dpt_LS1)


# Extract Least Squared Means from model
marginalsm.1 = emmeans(Dpt_LS1, ~ diet:infection) # extracts least square means

smCLD_dpt<- cld(marginalsm.1,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

smCLD_dpt


# Subtracts LS1 Means: Infected-Uninfected
smtwo_dpt_pairs<-update(pairs(smmarginal.1, by= "diet"), by=NULL)
summary(smtwo_dpt_pairs)

smtwo_pairs_dpt_CLD<- cld(smtwo_dpt_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
smtwo_pairs_dpt_CLD


# AttA Model
Att_LS1 <- lm(atta_ct ~ act_ct  + diet + infection * diet, data= serratia, na.action = na.omit)
anova(Att_LS1)
summary(Att_LS1)


# Extract Least Squared Means from model
smmarginal.2 = emmeans(Att_LS1, ~ diet:infection) # extracts least square means

smCLD_att<- cld(smmarginal.2,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

smCLD_att



# Subtracts LS1 Means: Infected-Uninfected
smtwo_att_pairs<-update(pairs(smmarginal.2, by= "diet"), by=NULL)
summary(smtwo_att_pairs)

smtwo_pairs_att_CLD<- cld(smtwo_att_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
smtwo_pairs_att_CLD


# Drosocin Model
Dro_LS1 <- lm(dro_ct ~ act_ct  + diet + infection *diet, data= serratia, na.action = na.omit)
anova(Dro_LS1)
summary(Dro_LS1)


# Extract Least Squared Means from model
smmarginal.3 = emmeans(Dro_LS1, ~ diet:infection) # extracts least square means

smCLD_dro<- cld(smmarginal.3,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

smCLD_dro

# Subtracts LS1 Means: Infected-Uninfected
smtwo_dro_pairs<-update(pairs(smmarginal.3, by= "diet"), by=NULL)
summary(smtwo_dro_pairs)

smtwo_pairs_dro_CLD<- cld(smtwo_dro_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
smtwo_pairs_dro_CLD


# Defensin Model
Def_LS1 <- lm(def_ct ~ act_ct  + diet + infection * diet, data= serratia, na.action = na.omit)
anova(Def_LS1)
summary(Def_LS1)

# Extract Least Squared Means from model
smmarginal.4 = emmeans(Def_LS1, ~ diet:infection) # extracts least square means

smCLD_def<- cld(smmarginal.4,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

smCLD_def


# Subtracts LS1 Means: Infected-Uninfected
smtwo_def_pairs<-update(pairs(smmarginal.4, by= "diet"), by=NULL)
summary(two_def_pairs)

smtwo_pairs_def_CLD<- cld(smtwo_def_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
smtwo_pairs_def_CLD


############# Cecropin Model################################################################################
Cec_LS1 <- lm(cec_ct ~ act_ct  + diet + infection * diet, data= serratia, na.action = na.omit)
anova(Cec_LS1)
summary(Cec_LS1)


# Extract Least Squared Means from model
smmarginal.6 = emmeans(Cec_LS1, ~ diet:infection) # extracts least square means

smCLD_cec<- cld(smmarginal.6,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

smCLD_cec


##Subtracts LS1 Means: Infected-Uninfected
smtwo_cec_pairs<-update(pairs(smmarginal.6, by= "diet"), by=NULL)
summary(smtwo_cec_pairs)

smtwo_pairs_cec_CLD<- cld(smtwo_cec_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
smtwo_pairs_cec_CLD


# Drosmycin Model
Drs_LS1 <- lm(drs_ct ~ act_ct  + diet + infection * diet, data= serratia, na.action = na.omit)
anova(Drs_LS1)
summary(Drs_LS1)



# Extract Least Squared Means from model
smmarginal.5 = emmeans(Drs_LS1, ~ diet:infection) # extracts least square means

smCLD_drs<- cld(smmarginal.5,
              alpha=0.05,
              Letters=letters,
              adjust="sidak") # compact letter display and correction

smCLD_drs


# Subtracts LS1 Means: Infected-Uninfected
smtwo_drs_pairs<-update(pairs(smmarginal.5, by= "diet"), by=NULL)
summary(smtwo_drs_pairs)

smtwo_pairs_drs_CLD<- cld(smtwo_drs_pairs,
                        alpha=0.05,
                        Letters=letters,
                        adjust="sidak")
smtwo_pairs_drs_CLD

