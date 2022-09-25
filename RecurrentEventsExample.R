rm(list=ls())
library("reReg")
data(simDat) 
head(simDat)

subset(simDat, event==0 & status==0)
subset(simDat, id==199)
simDat <- rbind(c(201, 0, 2.6, 0, 1, 1, -0.078),simDat)


## Nonparametric estimate 
plot(reReg(Recur(t.start %to% t.stop, id, event, status) ~ 1, data = simDat, B = 50)) 
fm <- Recur(t.start %to% t.stop, id, event, status) ~ x1 + x2 

## Fit the Cox rate model 
summary(reReg(fm, data = simDat, model = "cox", B = 50))


# A simulated data frame with the following variables 
# id subjects identification 
# t.start start of the interval 
# t.stop endpoint of the interval; when time origin is 0 this variable also marks 
#     the recurrence or terminal/censoring time 
# status terminal event indicator; 1 if a terminal event is recorded 
# event recurrent event indicator; 1 if a recurrent event is recorded 
# x1 baseline covariate generated from a standard uniform distribution 
# x2 baseline covariate generated from a standard uniform distribution (independent from x1)


setwd("/Users/jvila/Dropbox/JLupon/FGdiabet/")
load("./dat/datok.rda")
library(Hmisc)
library(gdata)
require("Epi")

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
################################################################################
################################################################################
datok$IMCcatB <- factor(as.character(datok$IMCcat), levels = c("18.5-<25", "<18.5", "25-<30", "30+"))
Hmisc::label(datok$IMCcatB) <- Hmisc::label(datok$IMCcat)

xvaris <- c("epi", "SEX", "Edat", "IMCcatB", "etiology", "HTA", "EFcat", "ACEI_ARB_s", "ARNI_s", "MRA_s")
#  HbA1ccat50
datNoMiss <- subset(datok, complete.cases(datok[, xvaris]))

formul <- paste("Recur(VISITMONTH %to% end, id, AdmiRecur, EXITUS) ~ ", paste(xvaris, collapse= " + "), sep="")
fm <- as.formula(noquote(formul))
resu1 <- summary(reReg(fm, data = datNoMiss, model = "cox", B = 50))

fm <- as.formula(noquote(sub("ACEI_ARB_s", "", formul)))
resu2 <- summary(reReg(fm, data = datNoMiss, model = "cox", B = 50))

resuPre <- resu1
resuPost <- resu2


compareModels <- function(resuPre, resuPost){
    xxx1 <- data.frame(id = seq(1, length(resuPre$coefficients[, 1])), 
                       coeffPre=resuPre$coefficients[, 1], 
                       pValuePre=round(resuPre$coefficients[, 4], 4))
    xxx1$xvari <- rownames(xxx1)
    
    xxx2 <- data.frame(coeffPost=resuPost$coefficients[, 1], 
                       pValuePost=round(resuPost$coefficients[, 4], 4))
    xxx2$xvari <- rownames(xxx2)
    
    xxx <- merge(xxx1, xxx2, by= "xvari", all.x=TRUE)
    xxx <- xxx[order(xxx$id), ]
    xxx <- remove.vars(xxx, "id")
    xxx$percent <- with(xxx, round((coeffPre-coeffPost)/coeffPre*100, 2))
    return(xxx)
}

compareModels(resu1, resu2)




stat.table(list(Sex=SEX, ACE_ARB=ACEI_ARB_s), list(N=count(),
      '%'=percent(SEX)),data=subset(datNoMiss, orden==1), margins=T)

stat.table(index= list(ACEI_ARB_s), 
           contents= list(mean(Edat), sd(Edat), count()),
           data=subset(datNoMiss, orden==1),
           margins=c(TRUE))

stat.table(list(MRA=MRA_s, ACE_ARB=ACEI_ARB_s), list(N=count(),
      '%'=percent(MRA_s)),data=subset(datNoMiss, orden==1), margins=T)


xxx <- datNoMiss[, c("id", "orden", "SEX", "ACEI_ARB_s")]



fm <- Recur(VISITMONTH %to% end, id, AdmiRecur, EXITUS) ~  epi 
(x1 <- summary(reReg(fm, data = datNoMiss, model = "cox", B = 50)))

fm <- Recur(VISITMONTH %to% end, id, AdmiRecur, EXITUS) ~  epi + SEX
(x2 <- summary(reReg(fm, data = datNoMiss, model = "cox", B = 50)))

compareModels(x2, x1)

fm <- Recur(VISITMONTH %to% end, id, AdmiRecur, EXITUS) ~  epi + SEX + Edat
(x3 <- summary(reReg(fm, data = datNoMiss, model = "cox", B = 50)))

compareModels(x3, x2)

stat.table(index= list(SEX), 
           contents= list(mean(Edat), sd(Edat), count()),
           data=subset(datNoMiss, orden==1),
           margins=c(TRUE))

xvaris <- c("epi", "SEX", "Edat", "IMCcatB", "etiology", "HTA", "EFcat", "ACEI_ARB_s", "ARNI_s", "MRA_s")


# no change example: ARNI_s
fm <- Recur(VISITMONTH %to% end, id, AdmiRecur, EXITUS) ~  ARNI_s
(x4 <- summary(reReg(fm, data = datNoMiss, model = "cox", B = 50)))

stat.table(list(ARNI=ARNI_s, ACE_ARB=ACEI_ARB_s), list(N=count(),
      '%'=percent(ARNI_s)),data=subset(datNoMiss, orden==1), margins=T)
chisq.test(with(subset(datNoMiss, orden==1), table(ARNI_s, ACEI_ARB_s)))

# strong change example: SEX
fm <- Recur(VISITMONTH %to% end, id, AdmiRecur, EXITUS) ~  SEX
(x5 <- summary(reReg(fm, data = datNoMiss, model = "cox", B = 50)))

stat.table(list(Sex=SEX, ACE_ARB=ACEI_ARB_s), list(N=count(),
      '%'=percent(SEX)),data=subset(datNoMiss, orden==1), margins=T)
chisq.test(with(subset(datNoMiss, orden==1), table(SEX, ACEI_ARB_s)))


