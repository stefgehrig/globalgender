library(tidyverse)
library(haven)
library(countrycode)
library(labelled)

##########################
#### import WEVS data ####
##########################
df_wvs <- readRDS("data/raw/WVS_Time_Series_1981-2022_rds_v5_0.rds")
df_evs <- read_dta("data/raw/ZA7503_v3-0-0.dta")

# set/harmonize some basic variables before merging
# ... in WVS
df_wvs$srv <- "WVS"
df_wvs$wave <- as.character(df_wvs$S002VS)
df_wvs$id   <- as.character(remove_labels(df_wvs$S007))

# ... in EVS
df_evs$srv <- "EVS"
df_evs$wave <- as.character(df_evs$S002EVS)
df_evs$X049 <- df_evs$x049a # townsize

###########################
#### merge WVS and EVS ####
###########################
df_merged <- map_dfr(list(df_wvs, df_evs), function(x){

    x %>% 
      transmute(id = as.character(remove_labels(S007)),
                country = to_character(S003),
                srv,
                year = remove_val_labels(S020),
                wave = remove_val_labels(wave),
                age = ifelse(X003<0, NA_integer_, X003),
                sex = to_character(X001),
                townsize = case_when(
                  x$srv == "EVS" ~ remove_val_labels(X049),
                  x$srv == "WVS" & remove_val_labels(X049) %in% c(1, 2) ~ 1, # under 2000, 2-5000
                  x$srv == "WVS" & remove_val_labels(X049) %in% c(3, 4) ~ 2, # 5-10000, 10-20000
                  x$srv == "WVS" & remove_val_labels(X049) %in% c(5, 6) ~ 3, # 20-50000, 50-100000
                  x$srv == "WVS" & remove_val_labels(X049) %in% c(7) ~ 4,    # 100-500000
                  x$srv == "WVS" & remove_val_labels(X049) %in% c(8) ~ 5,    # 500000 and more
                  TRUE ~ remove_val_labels(X049)
                ),
                income_lvl = ifelse(x$srv == "WVS", 
                                    to_character(X047R_WVS), 
                                    to_character(X047R_EVS)),
                education = case_when(
                  x$srv == "EVS" ~ remove_val_labels(X025), # non-isced
                  x$srv == "WVS" & x$wave != "7" ~ remove_val_labels(X025), # non-isced
                  x$srv == "WVS" & x$wave == "7" ~ remove_val_labels(X025A_01) # isced
                ),
                gender_lead = remove_val_labels(D059),
                gender_univ = remove_val_labels(D060),
                gender_jobs = remove_val_labels(C001),
                relig_child = remove_val_labels(A040),
                relig_persn = remove_val_labels(F034),
                relig_attnd = remove_val_labels(F028),
                relig_imprt = remove_val_labels(A006),
                children_nr = pmin(remove_val_labels(X011), 5), # in some survey (waves), the scale was limited at "5 or more", in others it wasn't
                )
}) 

# set India's wave 7 of WVS to 2022, the year which officially closes wave 7
df_merged$year[df_merged$year == 2023] <- 2022

################
#### recode ####
################
df <- df_merged %>% 
  mutate(sex = ifelse(sex %in% c("Male", "Female"), sex, NA), # recode labels for missingness
         age = ifelse(age %in% c(998,999) | age < 0, NA, age), # recode numeric values for missingness
         income_lvl = ifelse(!income_lvl %in% c("Low", "Medium", "High"), NA, income_lvl), # recode for missingness
         children_nr = ifelse(children_nr < 0, NA, children_nr), # recode numeric value for missingness
         townsize = case_when( # recode msisingness and give harmonized names
           townsize < 0 ~ NA_character_,
           townsize == 1 ~ "<5,000",
           townsize == 2 ~ "5,000-20,000",
           townsize == 3 ~ "20,000-100,000",
           townsize == 4 ~ "100,000-500,000",
           townsize == 5 ~ ">500,000"
         ),
         education = case_when(  # recode msisingness and bring all used variables on a same scale
           education < 0  ~ NA_character_,
           # isced
           srv == "WVS" & wave == 7 & education %in% c(0,1)     ~ "Primary or less",
           srv == "WVS" & wave == 7 & education %in% c(2)       ~ "Lower secondary",
           srv == "WVS" & wave == 7 & education %in% c(3,4)     ~ "Upper secondary",
           srv == "WVS" & wave == 7 & education %in% c(5,6,7,8) ~ "Tertiary",
           # non-isced
           education %in% c(1,2)   ~ "Primary or less",
           education %in% c(3)     ~ "Lower secondary",
           education %in% c(4,5,6) ~ "Upper secondary",
           education %in% c(7,8)   ~ "Tertiary"
         ),
         # scale all gender role items between 0 (least gender equal) and 1 (most gender equal)
         # equally spaced and recode the values for missingness
         gender_lead = case_when(
           gender_lead < 0 ~ NA,
           gender_lead == 4 ~ 1, # Strongly disagree
           gender_lead == 3 ~ 2/3, # Disagree
           gender_lead == 2 ~ 1/3, # Agree
           gender_lead == 1 ~ 0 # Agree strongly
         ),
         gender_univ = case_when(
           gender_univ < 0 ~ NA,
           gender_univ == 4 ~ 1, # Strongly disagree
           gender_univ == 3 ~ 2/3, # Disagree
           gender_univ == 2 ~ 1/3, # Agree
           gender_univ == 1 ~ 0 # Agree strongly
         ),
         gender_jobs = case_when(
           gender_jobs < 0 ~ NA,
           gender_jobs %in% c(8,9) ~ NA,
           gender_jobs == 1 ~ 0, # Agree
           gender_jobs == 3 ~ 1/2, # Neither
           gender_jobs == 2 ~ 1, # Disagree
         ),
         
         # scale religiosity items
         relig_child = case_when(
           relig_child < 0 ~ NA, 
           relig_child == 1 ~ 1, # 1 Important;
           TRUE ~ 0
         ),
         relig_persn = case_when(
           relig_persn < 0 ~ NA,
           relig_persn == 1 ~ 1, # 1 Religious person
           relig_persn == 2 ~ 1/2, # 2 Not religious person
           relig_persn == 3 ~ 0 # 3 Convinced atheist
         ),
         relig_attnd = case_when(
           relig_attnd < 0 ~ NA,
           TRUE ~ (7-ifelse(relig_attnd %in% 5:8, relig_attnd - 1, relig_attnd))/6 
           # 1 More than once a week; 7 Never, practically never 
         ),
         relig_imprt = case_when(
           relig_imprt < 0 ~ NA,
           TRUE ~ (4-relig_imprt)/3 , # 1 Very important; 4 Not at all important
         )) %>% 
  # apply 0-100 scale
  mutate(across(contains("gender") | contains("relig"), ~.x*100)) %>% 
  # compute indices
  mutate(
    grb = (gender_univ + gender_lead + gender_jobs) / 3,
    religiosity = (relig_child + relig_persn + relig_attnd + relig_imprt) / 4
  ) %>% 
  # apply harmonized country names
  mutate(
   country =  countrycode(country, "country.name", "country.name",
                          custom_match = c("Northern Ireland" = "Northern Ireland", 
                                           " Northern Ireland" = "Northern Ireland"))
  )

###############################
#### merge other data sets ####
###############################
# human development index
df_hdi <- read_csv("data/raw/HDR23-24_Composite_indices_complete_time_series.csv") %>% 
  filter(!str_length(iso3) > 3) %>% 
  mutate(country = countrycode(iso3, "iso3c", "country.name")) %>% 
  # use hdi in 2009, because it's the midpoint of the observation period
  select(country, hdi_2009) %>% 
  mutate(
    # for cutoff points, see https://hdr.undp.org/sites/default/files/2023-24_HDR/hdr2023-24_technical_notes.pdf
    hdicode = case_when(
      is.na(hdi_2009) ~ NA_character_,
      hdi_2009 < 0.55 ~ "Low",
      hdi_2009 < 0.7 ~ "Medium",
      hdi_2009 < 0.8 ~ "High",
      TRUE ~ "Very High"
    )
  ) %>% 
  select(country, hdicode) %>% 
  # assign missing levels based on contextual knowledge
  bind_rows(
    tibble(
      country = c("Puerto Rico", "Taiwan"),
      hdicode = rep("Very High", 2)
    )
  )

# cultural zones
df_cul <- read_csv2("data/raw/Cultural_Map.csv") %>% 
  transmute(
    country = countrycode(Country, "country.name", "country.name",
                          custom_match = c("Northern Ireland" = "Northern Ireland")),
    culzone = Zone
  )


# join all
df <- df %>% 
  left_join(df_hdi, by = join_by(country)) %>% 
  left_join(df_cul, by = join_by(country))

saveRDS(df, "data/wevsdata.rds")

####################################
#### prepare UN population data ####
####################################
# import adult population data
dfpp <- read_csv("data/raw/unpopulation_dataportal_20250118100557.csv")

# prepare data for aggregation
df_pop <- dfpp %>% 
  select(Location, year = Time, Sex, age = Age, Value) %>% 
  mutate(Location = ifelse(
    Location == "Micronesia", "Federated States of Micronesia", Location
  )) %>% 
  filter(
    !Location %in% c(
      "Africa", 
      "Americas", 
      "Asia", 
      "Caribbean", 
      "Least developed countries", 
      "Melanesia", 
      "Middle-income countries", 
      "Middle Africa", 
      "Oceania", 
      "Oceania (excluding Australia and New Zealand)", 
      "Polynesia", 
      "World"
    )
  ) %>% 
  pivot_wider(names_from = "Sex", values_from = "Value") %>% 
  mutate(age = ifelse(age == "100+", "100", age)) %>% 
  mutate(country = countrycode(Location, "country.name", "country.name")) %>% 
  select(-Location) %>% 
  mutate(age = as.numeric(age)) %>% 
  mutate(age_grp = case_when(
    is.na(age) ~ NA_character_,
    age <= 34  ~ "18 - 34",
    age <= 59  ~ "35 - 59",
    age >= 60  ~ "60 +",
  ))

# create adult population share by age and sex
df_pop_shares <- df_pop  %>% 
  group_by(country, year, age_grp) %>%
  summarise(across(.cols = c(Male, Female), sum), .groups = "drop_last") %>%
  pivot_longer(cols = c(Male, Female), names_to = "sex", values_to = "population") %>% 
  mutate(pr_stratum = population / sum(population)) %>% ungroup # adult population share
saveRDS(df_pop_shares, "data/UN_adultpop_bydemogr.rds")

# get global adult population
df_pop_world <- dfpp  %>% 
  filter(Location == "World") %>% 
  select(Location, year = Time, Sex, age = Age, Value) %>% 
  group_by(Location, year) %>% 
  summarise(population = sum(Value), .groups = "drop")
saveRDS(df_pop_world, "data/UN_adultpop_world.rds")
