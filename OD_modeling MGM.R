library(tidyverse)

OD_info <- read_csv("Downloads/Overdose_Information_Network_Data_CY_January_2018_-_Current_Monthly_County_State_Police.csv")

OD_unique <- OD_info %>%
  dplyr::select(
    `Incident ID`:`Victim County`, `Incident County FIPS Code`:`Victim County Latitude and Longitude`, `Naloxone Administered`, `Response Time Desc`,
    Survive
  ) %>%
  distinct()

OD_drug <- OD_info %>%
  dplyr::select(`Accidental Exposure`:`Susp OD Drug Desc`, `Incident ID`, `Victim ID`) %>%
  distinct() %>%
  mutate(tmp = 1) %>%
  pivot_wider(id_cols = c(`Incident ID`, `Victim ID`), values_from = tmp, values_fill = list(tmp = 0), names_from = `Susp OD Drug Desc`)

OD_all <- inner_join(OD_unique, OD_drug)

names(OD_all) <- make.names(names(OD_all)) # for convenience

OD_all$Died <- ifelse(OD_all$Survive == "N", 1, 0)
OD_all$Lived <- ifelse(OD_all$Survive == "Y", 1, 0)

OD_all$Age2 <- ifelse(OD_all$Age.Range %in% c("0 - 9", "10 - 14"), "15 - 19", OD_all$Age.Range)
OD_all$Age2 <- ifelse(OD_all$Age2 %in% c("70 - 79", "80 - *"), "60 - 69", OD_all$Age2)
OD_all$Race <- relevel(as.factor(OD_all$Race), "White")
OD_all$Race2 <- OD_all$Race
levels(OD_all$Race2) <- c("White", "Other or Unknown", "Other or Unknown", "Black", "Other or Unknown")
OD_all$Ethnicity.Desc <- relevel(as.factor(OD_all$Ethnicity.Desc), "Not Hispanic")
OD_all$HEROIN <- as.factor(OD_all$HEROIN)
OD_all$FENTANYL <- as.factor(OD_all$FENTANYL)
OD_all$UNKNOWN <- as.factor(OD_all$UNKNOWN)

OD_all$Female <- OD_all$Gender.Desc == "Female"

library(rms)
frm_died <- Died ~ Naloxone.Administered * (HEROIN * FENTANYL + UNKNOWN + Age2) + Race2  + Female
frm_lived <- Lived ~ Naloxone.Administered * (HEROIN * FENTANYL + UNKNOWN + Age2) + Race2  + Female

lrm_died <- lrm(frm_died, OD_all, x = TRUE, y = TRUE)
lrm_lived <- lrm(frm_lived, OD_all, x = TRUE, y = TRUE)

val_died <- rms::validate(lrm_died, B = 200)
val_lived <- rms::validate(lrm_lived, B = 200)

cal_died <- calibrate(lrm_died, B = 200)
cal_lived <- calibrate(lrm_lived, B = 200)

glm_died_1 <- glm(frm_died, data = OD_all, family = binomial)
glm_died_2 <- glmer(update(frm_died, . ~ . + (1 | Incident.County.FIPS.Code)), data = OD_all, family = binomial)
glm_died_3 <- glmer(update(frm_died, . ~ . + (Naloxone.Administered | Incident.County.FIPS.Code)), data = OD_all, family = binomial)


anova(glm_died_3, glm_died_2, glm_died_1)

drop1(glm_died_3, test = "Chisq")


OD_all %>%
  group_by(HEROIN, FENTANYL, Naloxone.Administered) %>%
  summarise(n = n(), Died = sum(Died), rate = mean(Died)) %>%
  arrange(Naloxone.Administered, HEROIN, FENTANYL)

OD_all %>%
  group_by(UNKNOWN, Naloxone.Administered) %>%
  summarise(n = n(), Died = sum(Died), rate = mean(Died)) %>%
  arrange(Naloxone.Administered, UNKNOWN)



library(ggeffects)
library(sjPlot)
library(scales)





OD_all %>%
  group_by(Naloxone.Administered, Age2) %>%
  summarise(n = n(), Died = sum(Died), Lived = sum(Lived)) %>%
  arrange(Naloxone.Administered, Age2)

library(tableone)

table_vars <- c("Age.Range", "Race", "Ethnicity.Desc", "Gender.Desc", "HEROIN", "FENTANYL", "UNKNOWN", "Survive")
factor_vars <- table_vars
table1 <- CreateTableOne(table_vars, "Naloxone.Administered", OD_all, factor_vars)

table_vars <- c("Age.Range", "Race", "Ethnicity.Desc", "Gender.Desc", "Naloxone.Administered", "UNKNOWN", "Survive")
factor_vars <- table_vars
table2 <- CreateTableOne(table_vars, c("HEROIN", "FENTANYL"), OD_all, factor_vars)

table_vars <- c("Age.Range","Race", "Ethnicity.Desc", "Gender.Desc", "HEROIN", "FENTANYL", "Naloxone.Administered", "Survive")
factor_vars <- table_vars
table3 <- CreateTableOne(table_vars, "UNKNOWN", OD_all, factor_vars)

#propensity score matching was the intent, but also useful to identify patterns in useage of Naxalone

OD_all$treat <- as.numeric(OD_all$Naloxone.Administered=="Y")

frm_ps <- treat ~ HEROIN * FENTANYL + UNKNOWN +  Ethnicity.Desc +
  Age2 + Race2 + Female + (1|Incident.County.FIPS.Code)

glm_ps <-glmer(frm_ps,family=binomial,data=OD_all)
drop1(glm_ps,test="Chisq")

library(survey)

ps=fitted(glm_ps)

OD_all$psw <- ifelse(OD_all$Naloxone.Administered=="Y",1/ps,1/(1-ps))

design_ps <-svydesign(ids=~Incident.County.FIPS.Code,weights=~psw,data=OD_all)

ps_model <-svyglm(Died~Naloxone.Administered*(HEROIN * FENTANYL+UNKNOWN+Age2),family=quasibinomial,design=design_ps)

table_vars <- c("Age.Range", "Race", "Ethnicity.Desc", "Gender.Desc", "HEROIN", "FENTANYL", "UNKNOWN", "Survive")
factor_vars <- table_vars
svyCreateTableOne(table_vars,"Naloxone.Administered" ,design_ps,factor_vars)


### county and age effects


age_effect_raw <- OD_all %>%
  mutate(Age2 = case_when(Age2 == "15 - 19" ~ "19 or less", Age2 == "60 - 69" ~ "60 or more", TRUE ~ Age2)) %>%
  group_by(Age2, Naloxone.Administered) %>%
  summarise(N = n(), r_died = mean(Died), r_lived = mean(Lived)) %>%
  mutate(se_died = sqrt(r_died * (1 - r_died) / N), se_lived = sqrt(r_lived * (1 - r_lived) / N))

age_effect_raw %>%
  mutate(lower = r_died - 1.96 * se_died, upper = r_died + 1.96 * se_died) %>%
  ggplot(aes(x = Age2, y = r_died, color = Naloxone.Administered)) + geom_point(size = 3) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 1) + theme_classic() + scale_y_continuous(label = percent) + xlab("Age Bracket") + ylab("% Died") + theme(axis.text = element_text(size = rel(1)))

age_effect_raw %>%
  mutate(lower = r_lived - 1.96 * se_lived, upper = r_lived + 1.96 * se_lived) %>%
  ggplot(aes(x = Age2, y = r_lived, color = Naloxone.Administered)) + geom_point(size = 3) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 1) + theme_classic() + scale_y_continuous(label = percent) + xlab("Age Bracket") + ylab("% Lived") + theme(axis.text = element_text(size = rel(1)))


county_effect_raw <- OD_all %>%
  group_by(Naloxone.Administered, Incident.County.FIPS.Code) %>%
  summarise(N = n(), r_died = mean(Died), r_lived = mean(Lived)) %>%
  mutate(se_died = sqrt(r_died * (1 - r_died) / N), se_lived = sqrt(r_lived * (1 - r_lived) / N))

county_effect_raw %>%
  ggplot(aes(x = r_died)) + geom_histogram() + facet_wrap(Naloxone.Administered ~ ., ncol = 1, labeller = label_both) + theme_classic() + xlab("% died") + ylab("# of counties") + scale_x_continuous(label = percent) + theme(axis.text = element_text(size = rel(1.0)), strip.text = element_text(size = rel(1.0)))

county_effect_raw %>%
  pivot_wider(id_cols = Incident.County.FIPS.Code, names_from = Naloxone.Administered, values_from = r_died) %>%
  ggplot(aes(x = N, y = Y)) + geom_point() + stat_smooth() + theme_classic() + xlab("% died, no Narcan") + ylab("% died, with Narcan") + scale_x_continuous(label = percent) + scale_y_continuous(label = percent) + theme(axis.text = element_text(size = rel(1.0))) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

county_narcan <- OD_all %>%
  group_by(Incident.County.FIPS.Code) %>%
  summarise(N = n(), r_treat = mean(treat)) %>%
  mutate(se_treat = sqrt(r_treat * (1 - r_treat) / N))
              

county_narcan %>%
  ggplot(aes(x = r_treat)) + geom_histogram()  + theme_classic() + xlab("% cases treated") + ylab("# of counties") + scale_x_continuous(label = percent) + theme(axis.text = element_text(size = rel(1.0)), strip.text = element_text(size = rel(1.0)))


county_narcan %>%
  ggplot(aes(x = N,y=r_treat)) + geom_point() +geom_smooth()  + theme_classic() + xlab("# of cases per county") + ylab("% treated") + scale_y_continuous(label = percent) + theme(axis.text = element_text(size = rel(1.0)), strip.text = element_text(size = rel(1.0)))+scale_x_log10()
    
  )
