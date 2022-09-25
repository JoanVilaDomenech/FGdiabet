model_3 <- lmer(epi ~ VISIT_YEARS + (1|id), data=dades)
model_3 %>% cftest

model_3 %>% plot

model_4 <- lmer(epi ~ VISIT_YEARS_Cat + (1|id), data=dades)
model_4 %>% cftest

emmeans_model_temps_cat <- emmeans(model_4, ~ VISIT_YEARS_Cat)

emmeans_model_temps_cat

plot(emmeans_model_temps_cat) + coord_flip()

# SEX
###########################
model_sex <- Mixed_models_FG(dades$SEX, dades$epi)
model_sex$anova_tnum
model_sex$anova_tcat
model_sex$cftest_tcat
model_sex$emmeans_model_tcat
model_sex$plot_marginal_means
fiber.emt <- emtrends(model_sex$model_tnum, "x", var = "VISIT_YEARS")
fiber.emt
pairs(fiber.emt)
model_sex$emmeans_model_tcat %>% logitudinal_plot + ylab('eFGR')

# DIABETES (makes no sense)
###########################

# FE cat
###########################
dades$FE_cat <- cut(dades$FE, breaks = c(4,40,49,86))
dades$FE_cat_rec <- as.factor(dades$FE_cat)
levels(dades$FE_cat_rec) <- c('EF<=40%', 'EF 41-49%', 'EF>=50')
table(dades$FE_cat_rec, dades$FE_cat)
model_FE <- Mixed_models_FG(dades$FE_cat_rec, dades$epi)
model_FE$anova_tnum
model_FE$cftest_tnum
model_FE$anova_tcat
model_FE$cftest_tcat
model_FE$emmeans_model_tcat
model_FE$plot_marginal_means
fiber.emt <- emtrends(model_FE$model_tnum, "x", var = "VISIT_YEARS")
fiber.emt
pairs(fiber.emt)
model_FE$emmeans_model_tcat %>% logitudinal_plot + ylab('eFGR')

# ETIOLOGIA
###########################
dades$ETIOLOGIA_rec <- ifelse(dades$ETIOLOGIA == 1, 'Ischemic', 'Non ischemic')
model_ETIOLOGIA <- Mixed_models_FG(dades$ETIOLOGIA_rec, dades$epi)
model_ETIOLOGIA$anova_tnum
model_ETIOLOGIA$cftest_tnum
model_ETIOLOGIA$anova_tcat
model_ETIOLOGIA$cftest_tcat
model_ETIOLOGIA$emmeans_model_tcat
model_ETIOLOGIA$plot_marginal_means
fiber.emt <- emtrends(model_ETIOLOGIA$model_tnum, "x", var = "VISIT_YEARS")
fiber.emt
pairs(fiber.emt)
model_ETIOLOGIA$emmeans_model_tcat %>% logitudinal_plot + ylab('eFGR')

# MORT
###########################
tmp <- dades[, c("id", "EXITUS")]
tmp <- subset(tmp, !is.na(EXITUS))
tmp <- unique(tmp)
(repes <- with(tmp, table(id)))[repes>1]
unique(dades$id)[unique(dades$id)%nin%tmp$id]
table(tmp$EXITUS)
tmp <- rename.vars(tmp, "EXITUS", "Mort")
dades <- merge(dades, tmp[, c("id", "Mort")], by = "id", all.x = TRUE)

dades$Mort <- ifelse(dades$Mort == 0, 'Alive', 'Dead')
model_MORT <- Mixed_models_FG(dades$Mort, dades$epi)
model_MORT$anova_tnum
model_MORT$cftest_tnum
model_MORT$anova_tcat
model_MORT$cftest_tcat
model_MORT$emmeans_model_tcat
model_MORT$plot_marginal_means
fiber.emt <- emtrends(model_MORT$model_tnum, "x", var = "VISIT_YEARS")
fiber.emt
pairs(fiber.emt)
model_MORT$emmeans_model_tcat %>% logitudinal_plot + ylab('eFGR')
