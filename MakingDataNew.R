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

dat <- data.frame(read_excel("./dat/BDD_FINAL_FINAL.xls", sheet = 1, col_names = TRUE))

dat <- rename.vars(dat, "n_pacie", "id")
sav
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
xxx <- subset(dat, !is.na(percentatge_bon_control))
length(unique(xxx$id))

dat$HbA1ccat <- factor(with(dat, ifelse(percentatge_bon_control <=70, "<=70",
                                 ifelse(percentatge_bon_control >70, ">70", NA))),
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

################################################################################
#  per comprovar que .- "numero_bon_control" i "percentatge_bon_control" de la BBDD BDD_FINAL_FINAL.xls
# hi la mateixa informació que a les variables 
# .- "percentatge_bon_control_HbA1c" i "numero_mesures" de la BBDD BDD REVISADA_25.04.22.xls
################################################################################
yyy <- dat[, c("id","VISITMONTH", "Hb_glicosilada", "HbA1C_75", "numero_mesures", "numero_bon_control", "percentatge_bon_control")]

partdat <- data.frame(read_excel("./dat/BDD REVISADA_25.04.22.xls", sheet = 1, col_names = TRUE))

partdat <- partdat[, c("n_pacie", "percentatge_bon_control_HbA1c", "numero_mesures")]
partdat <- rename.vars(partdat, "n_pacie", "id")
partdat <- rename.vars(partdat, "percentatge_bon_control_HbA1c", "percent")
partdat <- rename.vars(partdat, "numero_mesures", "n_mesures")
partdat <- partdat[order(partdat$id), ]
zzz <- unique(partdat)

(repes <- with(zzz, table(id)))[repes>1]
yyy <- merge(yyy, zzz, by="id", all.x=T)

subset(yyy, !is.na(percentatge_bon_control) & is.na(percent))
subset(yyy, !is.na(percent) & is.na(percentatge_bon_control))
subset(yyy, !is.na(percent) & !is.na(percentatge_bon_control) & (round(percent, 5)!= round(percentatge_bon_control, 5)))

subset(yyy, !is.na(numero_mesures) & is.na(n_mesures))
subset(yyy, !is.na(n_mesures) & is.na(numero_mesures))
subset(yyy, !is.na(numero_mesures) & !is.na(n_mesures) & (round(numero_mesures, 5)!= round(n_mesures, 5)))


unique(subset(yyy, is.na(percent))$id)
save(dat, file = "./dat/dat.rda")

################################################################################
################################################################################
################################################################################
rm(list=ls())
library(Hmisc)
library(gdata)

setwd("/Users/jvila/Dropbox/JLupon/FGdiabet")

load(file = "./dat/dat.rda")

# transfos to reproduce "Evolucio_FG.pdf/html"
dat$VISIT_YEARS <- dat$VISITMONTH/12

# ho salvo com "datpre" per errors varis
#########################################

dat <- dat[order(dat$id, dat$VISITMONTH), ]

##############################################
# creo idXL  
##############################################
xxx<-dat$id
xxx<-c(9999, xxx)
xxx<-format(xxx)
xxx<-gsub(" ", "0", xxx)
xxx<-xxx[-1]

yyy<-dat$VISITMONTH
yyy<-c(999, yyy)
yyy<-format(yyy)
yyy<-gsub(" ", "0", yyy)
yyy<-yyy[-1]


dat$idXL <- paste(xxx, yyy, sep="-")

# Creo la mort tal com l'han treballat a Evolucio FG.html
#########################################################
xdeath <- subset(dat, !is.na(EXITUS) & EXITUS==1)[, c("id" )]
dat$mort_rec <- with(dat, ifelse(id%in%xdeath, 1, 0))

prop.table(table(dat$mort_rec, useNA = "ifany"))

# ho salvo com a dat pre, perque ara ho refaig tot per poder fer time depending
datpre <- dat
save(datpre, file = "./dat/datpre.rda")


################################################################################
#############    Creo datok. per poder fer time depending  #####################
################################################################################
rm(list=ls())
library(Hmisc)
library(gdata)
library(chron)

setwd("/Users/jvila/Dropbox/JLupon/FGdiabet")



load(file = "./dat/datpre.rda")
dat <- datpre

##############################################
# Errors en les dates
##############################################
subset(dat, id==830)[, c("id", "VISITMONTH", "EXITUS", "DATA_EXITUS", "CAUSA_EXITUS",	"DATA_0")]
xDATA_0 <- subset(dat, id==830 & VISITMONTH==9)[, "DATA_0"]
dat$DATA_0 <- with(dat, ifelse(id==830 & VISITMONTH==0, xDATA_0, DATA_0))
dat$EXITUS <- with(dat, ifelse(id==830 & VISITMONTH!=0, NA, EXITUS))
dat$DATA_EXITUS <- with(dat, ifelse(id==830 & VISITMONTH!=0, NA, DATA_EXITUS))
dat$CAUSA_EXITUS <- with(dat, ifelse(id==830 & VISITMONTH!=0, NA, CAUSA_EXITUS))
dat$DATA_0 <- with(dat, ifelse(id==830 & VISITMONTH!=0, NA, DATA_0))


dat$todele <- with(dat, ifelse(id==1671 & VISITMONTH==6 & DATA_0=="#NULL!", 1, 0))
subset(dat, id==1671)[, c("todele", "id", "VISITMONTH", "EXITUS", "DATA_EXITUS", "CAUSA_EXITUS",	"DATA_0")]
dat <- subset(dat, todele==0)
xDATA_0 <- subset(dat, id==1671 & VISITMONTH==6)[, "DATA_0"]
dat$DATA_0 <- with(dat, ifelse(id==1671 & VISITMONTH==0, xDATA_0, DATA_0))
dat$EXITUS <- with(dat, ifelse(id==1671 & VISITMONTH!=0, NA, EXITUS))
dat$DATA_EXITUS <- with(dat, ifelse(id==1671 & VISITMONTH!=0, NA, DATA_EXITUS))
dat$CAUSA_EXITUS <- with(dat, ifelse(id==1671 & VISITMONTH!=0, NA, CAUSA_EXITUS))
dat$DATA_0 <- with(dat, ifelse(id==1671 & VISITMONTH!=0, NA, DATA_0))

subset(dat, id==1679)[, c("todele", "id", "VISITMONTH", "EXITUS", "DATA_EXITUS", "CAUSA_EXITUS",	"DATA_0")]
xDATA_0 <- subset(dat, id==1679 & VISITMONTH==66 & DATA_0!="#NULL!")[, "DATA_0"]
dat$DATA_0 <- with(dat, ifelse(id==1679 & VISITMONTH==0, xDATA_0, DATA_0))
dat$EXITUS <- with(dat, ifelse(id==1679 & VISITMONTH==0, 0, EXITUS))
dat$todele <- with(dat, ifelse(id==1679 & VISITMONTH==66 & DATA_0!="#NULL!", 1, 0))
dat <- subset(dat, todele==0)

subset(dat, id==1933)[, c("todele", "id", "VISITMONTH", "EXITUS", "DATA_EXITUS", "CAUSA_EXITUS",	"DATA_0")]
xDATA_0 <- subset(dat, id==1933 & VISITMONTH==6)[, "DATA_0"]
dat$DATA_0 <- with(dat, ifelse(id==1933 & VISITMONTH==0, xDATA_0, DATA_0))
dat$EXITUS <- with(dat, ifelse(id==1933 & VISITMONTH==0, 0, EXITUS))
dat$EXITUS <- with(dat, ifelse(id==1933 & VISITMONTH==6, NA, EXITUS))
dat$DATA_0 <- with(dat, ifelse(id==1933 & VISITMONTH!=0, NA, DATA_0))

subset(dat, id==2250)[, c("todele", "id", "VISITMONTH", "EXITUS", "DATA_EXITUS", "CAUSA_EXITUS",	"DATA_0")]
xDATA_0 <- subset(dat, id==2250 & VISITMONTH==6)[, "DATA_0"]
dat$DATA_0 <- with(dat, ifelse(id==2250 & VISITMONTH==0, xDATA_0, DATA_0))
dat$EXITUS <- with(dat, ifelse(id==2250 & VISITMONTH==0, 0, EXITUS))
dat$EXITUS <- with(dat, ifelse(id==2250 & VISITMONTH==6, NA, EXITUS))
dat$DATA_0 <- with(dat, ifelse(id==2250 & VISITMONTH!=0, NA, DATA_0))

dat <- remove.vars(dat, "todele")
##############################################
##############################################


##############################################
# elimino id repetits 
##############################################
(repes <- with(dat, table(idXL)))[repes>1]

dat <- dat[order(dat$idXL), ]             
x2<-sort(dat$idXL)
tt<-table(x2)
system.time(
ordre <- sapply(1:length(tt), function(i) 1:tt[i])
)
ordre<-c(unlist(ordre))
dat$ordre <- ordre

dat <- subset(dat, ordre==1)
(repes <- with(dat, table(idXL)))[repes>1]
##############################################
##############################################

##############################################
#########  Creo unes dades tmp on treballar dates
##############################################
tmp <- dat[, c("id", "idXL", "DATA_0", "VISITMONTH", "DATA_EXITUS", "EXITUS", "CAUSA_EXITUS")]
xvisit <- c(tmp$VISITMONTH[-1], NA)
tmp$end <- ifelse(xvisit==0, NA, xvisit)

# en el morts (calculo ultim mes)
#################################
death <- dat[!is.na(dat$DATA_EXITUS) & dat$DATA_EXITUS!= "#NULL!", c("id", "idXL", "DATA_0", "VISITMONTH", "DATA_EXITUS", "EXITUS", "CAUSA_EXITUS")]
death$day0 <- as.Date(chron(death$DATA_0, format="d-mon-Y", out.format=c(dates="day-mon-year")))
death$lastday <- as.Date(chron(death$DATA_EXITUS, format="d-mon-Y", out.format=c(dates="day-mon-year")))
death$lastmonth <- as.numeric(round((death$lastday - death$day0)/(365.25/12)))

death <- rename.vars(death, c("EXITUS", "CAUSA_EXITUS", "day0", "lastday", "lastmonth"),
            c("xEXITUS", "xCAUSA_EXITUS", "xday0", "xlastday", "xlastmonth"))
tmp <- merge(tmp, death[, c("id", "xEXITUS", "xCAUSA_EXITUS", "xday0", "xlastday", "xlastmonth")],
        by = "id", all.x =T)
tmp <- tmp[order(tmp$idXL), ]


tmp$DATA_0 <- tmp$xday0

tmp$DATA_EXITUS <- as.Date(chron(with(tmp, ifelse(is.na(end), xlastday, NA))))

tmp$EXITUS <- with(tmp, ifelse(is.na(end), xEXITUS, 0))

tmp$CAUSA_EXITUS <- with(tmp, ifelse(is.na(end), xCAUSA_EXITUS, NA))

tmp$end <- with(tmp, ifelse(is.na(end), xlastmonth, end))

tmp <- remove.vars(tmp, c("xEXITUS", "xCAUSA_EXITUS", "xday0", "xlastday", "xlastmonth"))


# els vius, elimino dades ultima visita
#######################################
sum(is.na(subset(tmp, id%nin%unique(death$id))$end))
nrow(subset(tmp, is.na(end)))
tmp <- subset(tmp, !is.na(end))

# correcting death== last(VISITMONTH)
#####################################
tmp$toend <- with(tmp, end - VISITMONTH)
tobecorrect <- subset(tmp, id%in%subset(tmp, toend ==0)$id)
tobecorrect <-  tobecorrect[order(tobecorrect$idXL), ]

tt<-table(tobecorrect$id)
ordre <- sapply(1:length(tt), function(i) 1:tt[i])
ordre<-c(unlist(ordre))
tobecorrect$ordre <- ordre

tmp$todele <- 0
xids <- unique(tobecorrect$id)
for (i in 1:length(xids)){
  xindex <- max(subset(tobecorrect, id==xids[i])$ordre)
  badid <- subset(tobecorrect, id==xids[i] & ordre==xindex)$idXL
  okid <- subset(tobecorrect, id==xids[i] & ordre==xindex-1)$idXL
  subset(tmp, idXL==badid)$DATA_EXITUS
  tmp$DATA_EXITUS <- as.Date(chron(with(tmp, ifelse(idXL==okid, 
                subset(tmp, idXL==badid)$DATA_EXITUS, DATA_EXITUS))))
  tmp$EXITUS <- with(tmp, ifelse(idXL==okid, 
                subset(tmp, idXL==badid)$EXITUS, EXITUS))
  tmp$CAUSA_EXITUS <- with(tmp, ifelse(idXL==okid, 
                subset(tmp, idXL==badid)$CAUSA_EXITUS, CAUSA_EXITUS))
  tmp$todele <- with(tmp, ifelse(idXL==badid, 1, todele))
}
tmp <- subset(tmp, todele!=1)
  
# correcting death < last(VISITMONTH)
#####################################
tobecorrect <- subset(tmp, id%in%subset(tmp, toend <=0)$id)
tobecorrect <-  tobecorrect[order(tobecorrect$idXL), ]

tt<-table(tobecorrect$id)
ordre <- sapply(1:length(tt), function(i) 1:tt[i])
ordre<-c(unlist(ordre))
tobecorrect$ordre <- ordre

tmp$todele <- 0
xids <- unique(tobecorrect$id)
for (i in 1:length(xids)){
  xindex <- max(subset(tobecorrect, id==xids[i])$ordre)
  badid <- subset(tobecorrect, id==xids[i] & ordre==xindex)$idXL
  okid <- subset(tobecorrect, id==xids[i] & ordre==xindex-1)$idXL
  subset(tmp, idXL==badid)$DATA_EXITUS
  tmp$DATA_EXITUS <- as.Date(chron(with(tmp, ifelse(idXL==okid, 
                subset(tmp, idXL==badid)$DATA_EXITUS, DATA_EXITUS))))
  tmp$EXITUS <- with(tmp, ifelse(idXL==okid, 
                subset(tmp, idXL==badid)$EXITUS, EXITUS))
  tmp$CAUSA_EXITUS <- with(tmp, ifelse(idXL==okid, 
                subset(tmp, idXL==badid)$CAUSA_EXITUS, CAUSA_EXITUS))
  tmp$todele <- with(tmp, ifelse(idXL==badid, 1, todele))
}
tmp <- subset(tmp, todele!=1)

xxx <- subset(tmp, id%in%tobecorrect$id )
xxx$xend <- as.numeric(round((xxx$DATA_EXITUS - xxx$DATA_0)/(365.25/12)))

tmp <- merge(tmp, subset(xxx, !is.na(xend))[, c("idXL", "xend")], by="idXL", all.x=T)
tmp$end <- with(tmp, ifelse(!is.na(xend), xend, end))
tmp$toend <- with(tmp, ifelse(!is.na(xend), xend- VISITMONTH, toend))

# remestons (n=2)
#####################################
subset(tmp, id%in%subset(tmp, toend<=0)$id)
tmp$end <- with(tmp, ifelse(idXL=="0617-162", 163, end))
tmp$end <- with(tmp, ifelse(idXL=="0743-030", 31, end))
tmp$toend <- with(tmp, end-VISITMONTH)
tmp <- remove.vars(tmp, c("todele", "xend", "toend"))

# base ok
##################################
dat <- remove.vars(dat, c("id","DATA_0", "VISITMONTH", "DATA_EXITUS", "EXITUS", "CAUSA_EXITUS"))
datok <- merge(tmp, dat, by="idXL", all.x = TRUE)
datok <- subset(datok, VISIT_YEARS < 15)

save(datok, file = "./dat/datok.rda")

################################################################################
################################################################################
################################################################################
rm(list=ls())
library(readxl)
library(Hmisc)
library(gdata)

setwd("/Users/jvila/Dropbox/JLupon/FGdiabet/")

load("./dat/datok.rda")

# -% de Bon control de HbA1c. La endocrí em va comentar que ella utilitzaria 
# el punt del tall de 50% per crear els dos grups. 
# És el punt de tall que vau utilitzar a l'últim article de diabètics amb en Josep. Tenim 171 missings, però ella ha insistit molt en l'interès d'aquest subanàlisi. 
datok$HbA1ccat50 <- factor(with(datok, ifelse(percentatge_bon_control <=50, "<=50",
                                 ifelse(percentatge_bon_control >50, ">50", NA))),
                  levels = c("<=50", ">50"))
Hmisc::label(datok$HbA1ccat50) <- "%HBA1c controlled (<=50)"
with(datok, table(HbA1ccat50, useNA = "ifany"))

save(datok, file = "./dat/datok.rda")


################################################################################
################################################################################
################################################################################
# Arranging data to analyse recurrent events

# NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO!!!!!!! VISITMONTH no es la visita de reingres
rm(list=ls())
setwd("/Users/jvila/Dropbox/JLupon/FGdiabet/")

require(Hmisc)
require(gdata)

load("./dat/datok.rda")

datok$AdmiRecur <- with(datok, ifelse(EXITUS==0, 1, 0))

# stating sequence
datok <- datok[order(datok$id, datok$end), ]
kk <- table(sort(datok$id))
orden <- sapply(1:length(kk), function(i) 1:kk[i])
datok$orden <- unlist(orden)

# adding times an specific id is repeated 
tt <- table(datok$id)
dat2 <- data.frame(id=names(tt),freq=as.integer(tt))
datok <- merge(datok, dat2,by="id",all.x=TRUE)

# indicator for the last
datok$TheLast <- with(datok, ifelse(orden==freq, 1, 0))

# AdmiRecur == 0 if last entry and not dead
datok$AdmiRecur <- with(datok, ifelse(TheLast==1 & EXITUS==0, 0, AdmiRecur))
datok <- remove.vars(datok, c("freq", "TheLast"))

################################################################################
datok$IMCcatB <- factor(as.character(datok$IMCcat), levels = c("18.5-<25", "<18.5", "25-<30", "30+"))
Hmisc::label(datok$IMCcatB) <- Hmisc::label(datok$IMCcat)

save(datok, file = "./dat/datok.rda")


################################################################################
####### Noves transfos sobre la base original ##################################
################################################################################
rm(list=ls())
library(Hmisc)
library(gdata)

setwd("/Users/jvila/Dropbox/JLupon/FGdiabet")
load("./dat/datpre.rda")

datpre$ingressos_agrupats <- factor(with(datpre, ifelse(ingressos_agrupats==0, "0",
             ifelse(ingressos_agrupats==1, "1-2",
             ifelse(ingressos_agrupats==2,  ">=3", NA)))), 
       levels= c("0", "1-2", ">=3"))

Hmisc::label(datpre$ingressos_agrupats) <- "Admissions (grouped)"

datpre$evolMoths <- factor(with(datpre, ifelse(TEMPS_EVOL <12, "<12",
             ifelse(TEMPS_EVOL>=12, "12+", NA))),
       levels = c("<12", "12+"))
Hmisc::label(datpre$evolMoths) <- "Months evolution HF"

datpre$HbA1ccat50 <- factor(with(datpre, ifelse(percentatge_bon_control <=50, "<=50",
                                 ifelse(percentatge_bon_control >50, ">50", NA))),
                  levels = c("<=50", ">50"))
Hmisc::label(datpre$HbA1ccat50) <- "%HBA1c controlled"

xxx <- subset(datpre, VISITMONTH==0)[, c("id", "epi")]
xxx$FGbasal <- factor(with(xxx, ifelse(epi <15, "<15", 
          ifelse(epi <30, "15-<30",
          ifelse(epi <45, "30-<45",
          ifelse(epi <60, "45-<60",
          ifelse(epi <90, "60-90", 
          ifelse(epi >=90, "90+", NA))))))),
       levels = c("<15", "15-<30", "30-<45", "45-<60", "60-90", "90+"))
datpre <- merge(datpre, xxx[, c("id", "FGbasal")], by="id", all.x=T)
Hmisc::label(datpre$FGbasal) <- "Baseline Epi"

datpre$SexAge <- factor(with(datpre, ifelse(SEX=="Women" & Edat <50, "Women/<50",
             ifelse(SEX=="Women" & Edat>=50, "Women/50+",
             ifelse(SEX=="Men" & Edat <50, "Men/<50",
             ifelse(SEX=="Men" & Edat>=50, "Men/50+", NA))))),
       levels = c("Women/<50", "Women/50+", "Men/<50", "Men/50+"))

Hmisc::label(datpre$SexAge) <- "Age groups / Sex"

with(datpre, table(FGbasal, useNA = "ifany"))

save(datpre, file = "./dat/datpre.rda")
