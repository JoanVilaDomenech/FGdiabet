rm(list=ls())
library(readxl)
library(Hmisc)
library(gdata)
library(tidyverse)

library(magrittr)
library(emmeans)
library(lme4) 
library(multcomp)

setwd("/Users/jvila/Dropbox/JLupon/FGdiabet/")
dades <- data.frame(read_excel("./dat/BBD_Evolucio_eGFR_FINAL.xlsx", sheet = 1, col_names = TRUE))


dades %>% mutate_at(c("SEX","ETIOLOGIA","DIABETES","HTA","NYHA","FA_FT"), list(~factor(.))) -> dades
summary(dades)

dades$VISIT_YEARS <- dades$VISITMONTH/12
dades <- subset(dades, VISIT_YEARS < 15)
dades$VISIT_YEARS %>% summary

dades$VISIT_YEARS_Cat <- cut(dades$VISIT_YEARS, breaks = 0:19)
dades$VISIT_YEARS_Cat %>% table -> table_years_follow_up
table_years_follow_up

table_years_follow_up %>% plot

Mixed_models_FG <- function(x, y){
  
  #Model time continuous
  lmer(y ~ x + x:VISIT_YEARS + (1|ID), data=dades) -> model_temps_num
  model_temps_num %>% cftest -> cftest_temps_num 
  model_temps_num %>% anova -> anova_temps_num
  
  #Model time categorical
  lmer(y ~ x + x:VISIT_YEARS_Cat + (1|ID), data=dades) -> model_temps_cat
  model_temps_cat %>% cftest -> cftest_temps_cat
  model_temps_cat %>% anova -> anova_temps_cat

  emmeans(model_temps_cat, ~ VISIT_YEARS_Cat*x) -> emmeans_model_temps_cat
  
  plot(emmeans_model_temps_cat) + coord_flip() -> plot_emmeans
  
  return(list(model_tnum = model_temps_num, anova_tnum = anova_temps_num, cftest_tnum = cftest_temps_num, 
              model_tcat = model_temps_cat, anova_tcat = anova_temps_cat, cftest_tcat = cftest_temps_cat, 
              emmeans_model_tcat = emmeans_model_temps_cat, plot_marginal_means = plot_emmeans))
}

logitudinal_plot <- function(dades_plot){

    dades_plot %>% as.data.frame -> dades_plot

    figura <- ggplot(data = dades_plot, aes(x =  VISIT_YEARS_Cat, y = emmean, group=x, fill=x)) + 
      # geom_point(size=5, col=I("black")) + 
      geom_line(aes(x = VISIT_YEARS_Cat, y = emmean, col=I("black")), lwd=1.3) + 
      geom_ribbon(aes(ymin = asymp.LCL,ymax = asymp.UCL), lwd=1.5, width=0.5, alpha = 0.5) + 
      theme_grey(base_size = 20) + xlab("Years") + ylab("...") + # facet_grid(x~.) + 
      theme(axis.text.x=element_text(angle=90, hjust=1)) + theme(legend.title = element_blank()) 
  
  return(figura)
}

tapply(dades$epi, dades$VISIT_YEARS_Cat, mean)

model_1 <- lm(epi ~ VISIT_YEARS, data=dades)
model_1 %>% summary

model_1 %>% plot


model_2 <- lm(epi ~ VISIT_YEARS_Cat, data=dades)
emmeans(model_2, ~ VISIT_YEARS_Cat) 


model_3 <- lmer(epi ~ VISIT_YEARS + (1|ID), data=dades)
model_3 %>% cftest

model_3 %>% plot

model_4 <- lmer(epi ~ VISIT_YEARS_Cat + (1|ID), data=dades)
model_4 %>% cftest


emmeans_model_temps_cat <- emmeans(model_4, ~ VISIT_YEARS_Cat)

emmeans_model_temps_cat

dades$SEX_rec <-  ifelse(dades$SEX == 1, 'Male', 'Female')
model_sex <- Mixed_models_FG(dades$SEX_rec, dades$epi)
model_sex$anova_tnum

model_sex$cftest_tnum

model_sex$anova_tcat


model_sex$cftest_tcat

model_sex$emmeans_model_tcat

model_sex$plot_marginal_means

fiber.emt <- emtrends(model_sex$model_tnum, "x", var = "VISIT_YEARS")

fiber.emt

pairs(fiber.emt)

model_sex$emmeans_model_tcat %>% logitudinal_plot + ylab('eFGR')
