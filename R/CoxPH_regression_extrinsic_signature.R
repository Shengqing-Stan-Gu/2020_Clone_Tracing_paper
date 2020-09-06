getwd()
library(survival)


## Analyze Mariathasan dataset
mariathasan_surv = read.csv("Exploratory_analysis/Mariathasan/2019-07-17_Mariathasan_survival.csv")
mariathasan_model = coxph(Surv(OS, Event) ~ cor_invivo.cor_order., data = mariathasan_surv)
summary(mariathasan_model)
write.csv(summary(mariathasan_model)$coefficients, "Exploratory_analysis/Mariathasan/2020-07-09_CoxPH_Mariathasan.csv")

## Analyze Van Allen dataset
vanallen_surv = read.csv("Exploratory_analysis/VanAllen/2019-07-17_VanAllen_survival.csv")
vanallen_model = coxph(Surv(OS, OS.Event) ~ cor_invivo.cor_order., data = vanallen_surv)
summary(vanallen_model)
write.csv(summary(vanallen_model)$coefficients, "Exploratory_analysis/VanAllen/2020-07-09_CoxPH_VanAllen.csv")

## Analyze Riaz dataset
riaz_on_surv = read.csv("Exploratory_analysis/Riaz/2019-08-20_Riaz_survival_On.csv")
riaz_on_model = coxph(Surv(OS, OS.Event) ~ cor_invivo.cor_order., data = riaz_on_surv)
summary(riaz_on_model)
write.csv(summary(riaz_on_model)$coefficients, "Exploratory_analysis/Riaz/2020-07-10_CoxPH_Riaz_on.csv")

riaz_pre_surv = read.csv("Exploratory_analysis/Riaz/2019-07-18_Riaz_survival.csv")
riaz_pre_model = coxph(Surv(OS, OS.Event) ~ cor_invivo.cor_order., data = riaz_pre_surv)
summary(riaz_pre_model)
write.csv(summary(riaz_pre_model)$coefficients, "Exploratory_analysis/Riaz/2020-07-10_CoxPH_Riaz_pre.csv")

## Analyze Snyder dataset
snyder_surv = read.csv("Exploratory_analysis/Snyder/2019-07-18_Snyder_survival.csv")
snyder_model = coxph(Surv(OS, OS.Event) ~ cor_invivo.cor_order., data = snyder_surv)
summary(snyder_model)
write.csv(summary(snyder_model)$coefficients, "Exploratory_analysis/Snyder/2020-07-10_CoxPH_Snyder.csv")

## Analyze Gide dataset
gide_surv = read.csv("Exploratory_analysis/Gide/2019-08-20_Gide_survival.csv")
gide_model = coxph(Surv(PFS, PFS.Event) ~ cor_invivo.cor_order., data = gide_surv)
summary(gide_model)
write.csv(summary(gide_model)$coefficients, "Exploratory_analysis/Gide/2020-07-10_CoxPH_Gide.csv")

## Analyze Hugo dataset
hugo_surv = read.csv("Exploratory_analysis/Hugo/2019-07-18_Hugo_survival.csv")
hugo_model = coxph(Surv(OS, Event) ~ cor_invivo.cor_order., data = hugo_surv)
summary(hugo_model)
write.csv(summary(hugo_model)$coefficients, "Exploratory_analysis/Hugo/2020-07-10_CoxPH_Hugo.csv")










