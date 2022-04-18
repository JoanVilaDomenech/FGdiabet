rm(list=ls())
library(readxl)
library(Hmisc)
library(gdata)
library(tidyverse)

library(magrittr)
library(emmeans)
library(lme4) 
library(multcomp)
library(compareGroups)

dat <- data.frame(read_excel("./dat/FGR_evolutiu_DM.xlsx", sheet = 1, col_names = TRUE))

dat <- rename.vars(dat, "nº_pacie", "id")

################################################################################
# 1) L'evolució del filtrat glomerular (epi) dels 935 pacients amb diabetes mellitus
with(dat, table(DIABETES, useNA = "ifany"))
length(unique(dat$id))

################################################################################
# 2) Analitzar i comparar els diferents subgrups:
#    -Homes (1) vs Dones (0)
with(dat, table(SEX, useNA = "ifany"))
dat$SEX <- factor(with(dat, ifelse(SEX==0, "Women", ifelse(SEX==1, "Men", NA))), 
                  levels= c("Women", "Men"))
Hmisc::label(dat$SEX) <- "Sex"

#    -IMC (<18.5, 18.5-24.9, 25-29.9, >=30)
dat$IMCcat <- factor(with(dat, ifelse(IMC <18.5, "<18.5", ifelse(IMC<25, "18.5-<25",
                  ifelse(IMC<30, "25-<30", ifelse(IMC>=30, "30+", NA))))),
                  levels = c("<18.5", "18.5-<25", "25-<30", "30+"))
Hmisc::label(dat$IMCcat) <- "BMI (categories)"
with(dat, table(IMCcat, useNA = "ifany"))

#    -Bon control HbA1c (>70% de mesures ben controlades, tot i que podem agafar un altre punt de tall)
xxx <- subset(dat, !is.na(Percentatge_bon_control_HBA1c))
length(unique(xxx$id))

dat$HbA1ccat <- factor(with(dat, ifelse(Percentatge_bon_control_HBA1c <=70, "<=70",
                                 ifelse(Percentatge_bon_control_HBA1c >70, ">70", NA))),
                  levels = c("<=70", ">70"))
Hmisc::label(dat$HbA1ccat) <- "%HBA1c controlled"
with(dat, table(HbA1ccat, useNA = "ifany"))


#    -Etiologia isquèmica (1) vs Altres
dat$etiology <- factor(with(dat, ifelse(ETIOLOGIA==1, "Ischemic", ifelse(ETIOLOGIA!=1, "Non-Ischemic", NA))), 
                       levels= c("Non-Ischemic", "Ischemic"))
Hmisc::label(dat$etiology) <- "Etiology"
with(dat, table(etiology, useNA = "ifany"))

#    -Hipertensió arterial (Si/No)
dat$HTA <- factor(with(dat, ifelse(HTA==1, "Yes", ifelse(HTA==0, "No", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$HTA) <- "Hypertension"
with(dat, table(HTA, useNA = "ifany"))

#    -EF<=40%/EF 41-49%/EF>=50%
dat$EFcat <- factor(with(dat, ifelse(FE<=40, "<=40%", ifelse(FE<=49, "41-49%", 
                              ifelse(FE>=50, "50+", NA)))), 
                       levels= c("<=40%", "41-49%", "50+"))
Hmisc::label(dat$EFcat) <- "Ejection Fraction"
with(dat, table(EFcat, useNA = "ifany"))

#    -Hospitalizations (Ingressos_agrupats):0(0)/ 1-2(1)/ >=3(2)  
dat$admicat <- factor(with(dat, ifelse(n_total_ingressos==0, "0", 
                                ifelse(n_total_ingressos<=2, "1-2", 
                                ifelse(n_total_ingressos>=3, "3+", NA)))), 
                       levels= c("0", "1-2", "3+"))
Hmisc::label(dat$admicat) <- "Hospital Admissions"
with(dat, table(admicat, useNA = "ifany"))

#    -IECA/ARA II (Si/No)
dat$ACEI_ARB_s <- factor(with(dat, ifelse(ACEI_ARB_s==0, "No", ifelse(ACEI_ARB_s==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$ACEI_ARB_s) <- "ACE/ARA-II"
with(dat, table(ACEI_ARB_s, useNA = "ifany"))

#    -ARNI (Si/No)
dat$ARNI_s <- factor(with(dat, ifelse(ARNI_s==0, "No", ifelse(ARNI_s==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$ARNI_s) <- "ARNI"
with(dat, table(ARNI_s, useNA = "ifany"))



# 3)Crear quartils en funció de la pendent ("slope") de caiguda del filtrat i mortalitat global, mortalitat cardiovascular i hospitaltzacions per insuficiència cardíaca en el seguiment segons el quartil de la pendent. 
# Les causes de mort estan codificades com:
#     Cardiovascular=1,2,3,4,5,7
#     No cardiovascular=6, 
#     Desconeguda=0/missing. 
dat$status <- factor(with(dat, ifelse(EXITUS==0, "Alive", ifelse(EXITUS==1, "Death", NA))), 
                       levels= c("Alive", "Death"))
Hmisc::label(dat$status) <- "Follow-up Status"
with(dat, table(status, useNA = "ifany"))


dat$DeathStat <- factor(with(dat, ifelse(EXITUS==0, "Alive", ifelse(CAUSA_EXITUS==0, "Death unknown cause", 
                               ifelse(CAUSA_EXITUS%in%c(1,2,3,4,5,7), "Cardiovascular death",
                               ifelse(CAUSA_EXITUS==6, "Non-Cardiovascular death", NA))))), 
                       levels= c("Alive", "Non-Cardiovascular death", "Cardiovascular death", "Death unknown cause"))
Hmisc::label(dat$DeathStat) <- "Death status"
with(dat, table(DeathStat, useNA = "ifany"))

# "Creat_st"  
Hmisc::label(dat$Creat_st) <- "Creatinine"

# "Edat"
Hmisc::label(dat$Edat) <- "Age"

# "epi"  
Hmisc::label(dat$epi) <- "EPI (GFR)"

summary(subset(dat, VISITMONTH==0)$epi)
dat$EPIcat <- factor(with(dat, ifelse(epi <40, "<40", ifelse(epi<60, ">=40; <60", 
                               ifelse(epi <80, ">=60; <80",
                               ifelse(epi>=80, "80+", NA))))), 
                       levels= c("<40", ">=40; <60", ">=60; <80", "80+"))
Hmisc::label(dat$EPIcat) <- "EPI (grouped)"
with(dat, table(EPIcat, useNA = "ifany"))



# "TEMPS_EVOL"    
Hmisc::label(dat$TEMPS_EVOL) <- "Evolution Time"

# "NYHA"
dat$NYHAcat <- factor(with(dat, ifelse(NYHA%in%c(1,2), "I-II", ifelse(NYHA%in%c(3,4), "III-IV", NA))), 
                       levels= c("I-II", "III-IV"))
Hmisc::label(dat$NYHAcat) <- "NYHA"
with(dat, table(NYHAcat, useNA = "ifany"))

# "FA_FT"                        
dat$FA_FT <- factor(with(dat, ifelse(FA_FT==0, "No", ifelse(FA_FT==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$FA_FT) <- "FA/FT"
with(dat, table(FA_FT, useNA = "ifany"))

# "Beta_block_s"                  
dat$Beta_block_s <- factor(with(dat, ifelse(Beta_block_s==0, "No", ifelse(Beta_block_s==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$Beta_block_s) <- "Beta blockers"
with(dat, table(Beta_block_s, useNA = "ifany"))

# "MRA_s"
dat$MRA_s <- factor(with(dat, ifelse(MRA_s==0, "No", ifelse(MRA_s==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$MRA_s) <- "MRA"
with(dat, table(MRA_s, useNA = "ifany"))

# "Loop_diur_s"
dat$Loop_diur_s <- factor(with(dat, ifelse(Loop_diur_s==0, "No", ifelse(Loop_diur_s==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$Loop_diur_s) <- "Loop diuretics"
with(dat, table(Loop_diur_s, useNA = "ifany"))

# "Digoxin_s"                    
dat$Digoxin_s <- factor(with(dat, ifelse(Digoxin_s==0, "No", ifelse(Digoxin_s==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$Digoxin_s) <- "Digoxin"
with(dat, table(Digoxin_s, useNA = "ifany"))

# "CRT_s"      
dat$CRT_s <- factor(with(dat, ifelse(CRT_s==0, "No", ifelse(CRT_s==1, "Yes", NA))), 
                       levels= c("No", "Yes"))
Hmisc::label(dat$CRT_s) <- "CRT"
with(dat, table(CRT_s, useNA = "ifany"))

save(dat, file = "./dat/dat.rda")



