# 
# Description: Prepare data for sample size calculation
# 
# This script prepares the dataset needed for Alain Amstutz to calculate the 
# sample size for the NUDGE project. 
# SHCS database used: 01/2026
# 
# Last updated: 10.2.2026
#

# Libraries and functions -------------------------------------------------

library(haven)
library(tidyverse)


# Reads data from the right folder
shcs_dta <- function(file) {
  path <- "data/2601stata/"
  read_dta(here::here(path, file))
}


# Read data ---------------------------------------------------------------

fup <- shcs_dta("fup.dta")
stop <- shcs_dta("stop.dta")
pat <- shcs_dta("pat.dta")
lab2 <- vroom::vroom(file = "data/2601stata/lab2.csv")
clinical <- shcs_dta("clinical.dta")
all_drugs <- read_rds("data/2601stata/00-all_drugs.rds")
art_drugs <- read_rds("data/2601stata/00-art_drugs.rds")



# Prepare data ------------------------------------------------------------


# Determine whether people were on a statin
on_statin <- all_drugs |> 
  # Filter ATC Code for statin
  filter(str_detect(atc, "C10AA")) |> 
  # Smaller dataset
  select(id, start_date, stop_date, brand_name) |> 
  # Only keep last statin entry per ID
  arrange(id, desc(start_date)) |> 
  group_by(id) |> 
  slice(1) |> 
  ungroup() |> 
  # Only select those that are currently on statin (i.e. no stop date)
  filter(is.na(stop_date)) 

# Determine whether people are on ART
on_art <- art_drugs |> 
  # Only keep latest entry 
  arrange(id, desc(start_date)) |> 
  group_by(id) |>
  slice(1) |> 
  ungroup() |> 
  # Only select those that are currently on ART (i.e. no stop date)
  filter(is.na(stop_date)) |> 
  select(id, art_start = start_date, art = brand_name)



# Information on DM: clinical diagnosis
dia_clinical <- clinical |> 
  filter(clin_id == "DIA") |> 
  select(id, diabetes = clin_id, diab_clin_date = clin_date) |> 
  # only keep first clinical diagnosis of DM
  arrange(id, diab_clin_date) |>
  group_by(id) |>
  slice(1) |>
  ungroup()

# Lab HBA1c above 6.5
hba1c <- lab2 |> 
  filter(item == "HBA1C") |> 
  # Only keep highest HbA1c
  arrange(id, desc(value)) |> 
  group_by(id) |> 
  slice(1) |> 
  ungroup() |> 
  # Some centers report HbA1c in mmol/mol, convert using formula
  mutate(value = if_else(value >= 24, (value * 0.0915) + 2.15, value)) |> 
  # Keep those with HbA1c >= 6.5%
  filter(value >= 6.5) |> 
  select(id, hba1c = value)


# CVD events: 
cvd_event <- clinical |> 
  # Select CVD events
  filter(clin_id %in% c("AMI", "ANG", "CEI", "END", "PRO")) |> 
  # If there were more than one CVD events per individual, take the first
  arrange(id, clin_date) |> 
  group_by(id) |>
  slice(1) |> 
  ungroup() |> 
  select(id, first_cvd_date = clin_date)



# Select participants and physicians --------------------------------------

# Step 1: Get all visits from 2025 and add relevant data
visits_25 <- fup |> 
  select(id, fupdate, physician, center, source) |> 
  # Only consider visits from 2025
  filter(fupdate >= ymd("2025-01-01"),
         fupdate <= ymd("2025-12-31")) |> 
  # Remove individuals who were stopped in the SHCS
  anti_join(stop, by = "id") |> 
  # Clean center names
  mutate(center = case_when(
    center == 10 ~ "Zürich",
    center == 20 ~ "Basel",
    center == 30 ~ "Bern",
    center == 40 ~ "Geneva",
    center == 50 ~ "Lausanne",
    center == 60 ~ "Lugano",
    center == 70 ~ "Sankt Gallen"
  )) |> 
  # Get info on sex and age
  left_join(pat |> select(id, sex, born)) |> 
  mutate(age = year(fupdate) - born) |> 
  # Get info on statin
  left_join(on_statin, by = "id") |> 
  rename(statin_start = start_date, 
         statin_stop = stop_date,
         statin_name = brand_name) |>
  # Determine if the patient was on statin at the time of the visit
  mutate(on_statin = case_when(
    statin_start <= fupdate & statin_stop > fupdate ~ TRUE,
    statin_start <= fupdate & is.na(statin_stop) ~ TRUE,
    TRUE ~ FALSE)) |> 
  select(-starts_with("statin")) |>
  # Get info on ART (at last visit, not very elegant but unlikely to change much
  left_join(on_art, by = "id") |> 
  # Get info on Diabetes
  left_join(dia_clinical, by = "id") |> 
  left_join(hba1c, by = "id") |> 
  # Determine if the patient had diabetes at the time of the visit
  mutate(diab = diab_clin_date <= fupdate | hba1c >= 6.5) |> 
  mutate(diab = replace_na(diab, FALSE)) |> 
  select(-diabetes, -diab_clin_date, -hba1c) |> 
  # Get date of first CVD event
  left_join(cvd_event, by = "id") |> 
  mutate(cvd = fupdate >= first_cvd_date, 
         cvd = replace_na(cvd, FALSE)) |>
  select(-first_cvd_date)


# Get eligible individuals (people on statin still included for later)
elig_25 <- visits_25 |> # 15'702
  # Remove people not on ART
  filter(!is.na(art)) |> # 15'427
  # Only consider people between 40 and 75 years of age
  filter(age >= 40, age <= 75) |> # 13'023
  # Remove people with a history of CVD
  filter(!cvd) |> # 11'950
  # Remove people with diabetes
  filter(!diab) # 10'724

  
# Get overall statin prescription rate by physician and center
tab <- elig_25 |> 
  # Get to distinct participant-id pairs
  # NOTE: individuals who had a visit without and one with statin will be 
  # counted as two separate individuals # Reduces visits by 131 (~65 patients)
  distinct(id, physician, center, source, on_statin) |> 
  group_by(physician, center, source) |> 
  # Get the number of peoplel overal, and the number of people not on statin
  summarise(n_overall = n(),
            n_without_statin = sum(!on_statin), 
            statin_presc_rate = 1 - n_without_statin/n_overall) |> 
  arrange(center, source, desc(n_without_statin))




# Write out data
write_rds(elig_25, "processed/01-eligible_visits25.rds")
write_rds(tab, "processed/01-visits_per_physician25.rds")    






