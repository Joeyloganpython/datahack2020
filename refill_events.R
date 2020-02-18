library(tidyverse)

refill <- read_csv("pp_refill_events.csv")

#Format dates
refill$refill_date <- as.Date(refill$refill_date, format = "%m/%d/%y")

#Number visits per participant
id_unknown <- filter(refill, h_id == 999999)
id_known <- filter(refill, h_id != 999999)
refill <- id_known %>% group_by(h_id) %>% mutate(visit_number = paste("Visit", min_rank(refill_date))) %>% full_join(id_unknown)

#Identify repeat visits in 2018
refill$repeated_visit <- as.factor(ifelse(refill$visit_number != "Visit 1", "yes", "no"))

#Survival as factor
refill$survival <- as.factor(ifelse(refill$outcome_death == 1, "no", ifelse(refill$outcome_ok == 1, "yes", "unknown")))

#Clean Narcan type
refill$used_narcan_type <- as.factor(refill$used_narcan_type)
refill <- refill %>% mutate(used_narcan_type.collapse = fct_recode(used_narcan_type, "intranasal" = "nose", "intranasal" = "intransal", "intranasal" = "instransal", "intranasal" = "intrtanasal", "intranasal" = "2",
                                                     "injectable" = "muscle", "injectable" = "muslce", "injectable" = "1",
                                                     "Both" = "1 and 2", "Both" = "injectable, intranasal", "Both" = "nose, muscle", "Both" = "1,2", "Both" = "1 & 2", "Both" = "1& 2", "Both" = "Injectable and Intranasal", "Both" = "1&2", "Both" = "1,2",
                                                     NULL = "<NA>", NULL = "F and Clearfield", NULL = "heroin, cocaine", NULL = "999", NULL = "4"))
#Clean Med Program
refill$med_program_used <- as.factor(refill$med_program_used)
refill <- refill %>% mutate(med_program_used.collapse = fct_recode(med_program_used, "sep-onsite" = "sep onsite", "sep-onsite" = "sep-oniste",
                                                          "sep-mobile" = "sep mobile", "drop-in" = "drop in/testing",
                                                          "other" = "cm", "other" = "wn", "other" = "external", "other" = "testing", "other" = "shp", "other" = "outreach", "other" = "wcc", "other" = "wound care", "other" = "outreach/homeless",
                                                          "other" = "special projects", "other" = "women's night", "other" = "sunray", "other" = "admin", "other" = "administration",
                                                          NULL = "2", NULL = "7", NULL = "6", NULL = "1", NULL = "999"))

#Clean Train Program
refill$train_program <- as.factor(refill$train_program)
refill <- refill %>% mutate(train_program.collapse = fct_recode(train_program, "drop-in" = "drop in/testing",
                                                                   "other" = "cm", "other" = "wn", "other" = "external", "other" = "testing", "other" = "shp", "other" = "outreach", "other" = "wcc", "other" = "wound care", "other" = "outreach/homeless",
                                                                   "other" = "special projects", "other" = "women's night", "other" = "sunray", "other" = "admin", "other" = "streetside health project", "other" = "staff", "other" = "d-o", "other" = "case management",
                                                                   NULL = "2", NULL = "7", NULL = "6", NULL = "1"))

#Remove implausible values
is.na(refill$n_pp_refill) <- which(refill$n_pp_refill > 1000)
is.na(refill$n_o_refill) <- which(refill$n_o_refill > 1000)
is.na(refill$n_ods) <- which(refill$n_ods > 1000)
is.na(refill$n_admins) <- which(refill$n_admins > 1000)
is.na(refill$n_revivals) <- which(refill$n_revivals > 1000)
is.na(refill$ml_naloxone) <- which(refill$ml_naloxone > 100)
is.na(refill$outcome_ed) <- which(refill$outcome_ed == 3)
is.na(refill$od_drug_heroin) <- which(refill$od_drug_heroin == 2)
is.na(refill$od_drug_cocaine) <- which(refill$od_drug_cocaine == 2 | refill$od_drug_cocaine == 99)
is.na(refill$od_drug_fentanyl) <- which(refill$od_drug_fentanyl == 99)

#Clean Other Factors
cols <- c("gender", "race", "cpr_used", "outcome_ems", "outcome_police", "outcome_ed", "od_drug_heroin", "od_drug_cocaine", "od_drug_fentanyl")
refill[cols] <- lapply(refill[cols], factor)
refill <- refill %>% mutate(race.collapse = fct_recode(race, "hispanic" = "latino", "multiracial" = "mixed"))
