################################################################
# Codes to get demo table in JITAI group
# Jitao Wang
# 4/20/2022
################################################################

library(tidyverse)
library(ggplot2)
library(geepack)
library(sjPlot)
library(haven)

if (.Platform$OS.type == "unix") {
  prefix = "~/Dropbox (University of Michigan)/CODA Accelerating Synergy/Standard Analytic Files/Accelerating Syngergy Caregiver SAF 202111/Standard Analytic Files/SAS/"
} else {
  prefix = "C:/Users/jtwan/Dropbox (University of Michigan)/CODA Accelerating Synergy/Standard Analytic Files/Accelerating Syngergy Caregiver SAF 202111/Standard Analytic Files/SAS/"
}

JITAI_participants <- read_sas(paste(prefix,"fitbit_baseline.sas7bdat",sep="")) %>%
  select(participant_ID,subjectid,elig_cohortdefine,group) %>%
  filter(group==2) %>% mutate(subgroup = recode(elig_cohortdefine,`1`="HD",`2`="SCI",`3`="HCT")) %>%
  select(participant_ID,subgroup)

dat_demo <- read_sas(paste(prefix,"participant.sas7bdat",sep="")) %>%
  select(participant_ID,access_code,subjectid,dem_age,dem_sex,dem_race,dem_ethnicity,dem_relationship,dem_cgrelation,dem_age_ppcf,dem_workstatus,dem_timespent,dem_hd2,dem_sci2,dem_hct1a) %>%
  inner_join(JITAI_participants) %>%
  mutate(dem_careduration = case_when(subgroup == "HD" ~ dem_hd2,
                                      subgroup == "SCI" ~ dem_sci2,
                                      subgroup == "HCT" ~ dem_hct1a)) %>% 
  mutate(dem_sex = recode(dem_sex,`1`="male",`2`="female")) %>%
  mutate(dem_race = recode(dem_race,`5`="Caucasian",`4`="African American",`2`="Asian",`6`="More than 1")) %>%
  mutate(dem_ethnicity = recode(dem_ethnicity,`1`="Non-Hispanic",`2`="Hispanic",.missing = "NA")) %>%
  mutate(dem_relationship = recode(dem_relationship,`0`="Married/Cohabitating",`1`="Married/Cohabitating",`2`="Partenered,living apart",`4`="Single & divorced",`5`="Single & never married",`7`="Missing",.missing ="NA")) %>%
  mutate(dem_workstatus = recode(dem_workstatus,`0`="Full-time",`1`="Part-time",`2`="Homemaker",`3`="Student",`4`="Retired",`5`="Retired early, disability",
                                 `6`="Unemployed < 1yr, LOOKING for work",`8`="Unemployed > 1yr, LOOKING for work", `9`="Unemployed > 1yr, NOT LOOKING for work",`11`="Other")) %>%
  mutate(dem_cgrelation = recode(dem_cgrelation,`1`="Spouse/Partner",`2`="Spouse/Partner",`3`="Spouse/Partner",`4`="Spouse/Partner",`5`="Child",`6`="Child",
                                 `7`="Parent",`8`="Parent",`9`="Sibling",`10`="Sibling",`11`="Friend",`12`="In-law",`13`="In-law",`14`="In-law",`15`="In-law",
                                 `16`="In-law",`17`="In-law",`18`="In-law",`19`="In-law", `20`="In-law",`21`="Other")) %>%
  mutate(dem_timespent = recode(dem_timespent,`1`="1-2 hrs/day or less",`2`="3-4 hrs/day (half of a working day)",`3`="5-8 hrs/day (full working day)",`4`="9-12 hrs/day",
                                `5`="> 12 hrs/day or round-the-clock care"))
  

# for all JITAI groups
tmp = dat_demo
describe(tmp,exclude.missing = FALSE)
print(tmp %>% pull(dem_careduration) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age_ppcf) %>% sd(na.rm=TRUE),digits=5)

# for HD groups
tmp = dat_demo %>% filter(subgroup == "HD")
describe(tmp,exclude.missing = FALSE)
print(tmp %>% pull(dem_careduration) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age_ppcf) %>% sd(na.rm=TRUE),digits=5)

# for SCI groups
tmp = dat_demo %>% filter(subgroup == "SCI")
describe(tmp,exclude.missing = FALSE)
print(tmp %>% pull(dem_careduration) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age_ppcf) %>% sd(na.rm=TRUE),digits=5)

# for HCT groups
tmp = dat_demo %>% filter(subgroup == "HCT")
describe(tmp,exclude.missing = FALSE)
print(tmp %>% pull(dem_careduration) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age) %>% sd(na.rm=TRUE),digits=5)
print(tmp %>% pull(dem_age_ppcf) %>% sd(na.rm=TRUE),digits=5)