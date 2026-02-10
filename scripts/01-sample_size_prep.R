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
lab2 <- shcs_dta("lab2.dta")
clinical <- shcs_dta("clinical.dta")
all_drugs <- read_rds("data/2601stata/00-all_drugs.rds")
art_drugs <- read_rds("data/2601stata/00-art_drugs.rds")



# Prepared data -----------------------------------------------------------
# Select latest follow-up of all partipicants and keep physician, 
# center, and source of follow-up. Remove those that are stopped in the DB.

active <- fup |> 
  select(id, fupdate, physician, source, center) |> 
  arrange(id, desc(fupdate)) |> 
  group_by(id) |> 
  # Subset to latest follow-up per ID
  slice(1) |>  
  ungroup() |> # N = 22'390
  # Remove participants that are stopped in the DB
  anti_join(stop, by = "id") |> 
  mutate(center = case_when(
    center == 10 ~ "Zürich",
    center == 20 ~ "Basel",
    center == 30 ~ "Bern",
    center == 40 ~ "Geneva",
    center == 50 ~ "Lausanne",
    center == 60 ~ "Lugano",
    center == 70 ~ "Sankt Gallen"
  ))


# Add age and sex to list
active <- active |> 
  left_join(pat |> select(id, sex, born)) |> 
  mutate(age = 2026 - born)

# ATC code for statin: C10AA
on_statin <- all_drugs |> 
  filter(str_detect(atc, "C10AA")) |> 
  select(id, start_date, stop_date, brand_name) |> 
  arrange(id, desc(start_date)) |> 
  group_by(id) |> 
  slice(1) |> # Only take last statin entry
  ungroup() |> 
  filter(is.na(stop_date)) # Only select those that are currently on statin (i.e. no stop date)


# Add statins to list
active <- active |> 
  left_join(on_statin, by = "id") |> 
  rename(statin_start = start_date, 
         statin_name = brand_name) |> 
  select(-stop_date)

# On ART
on_art <- art_drugs |> 
  arrange(id, desc(start_date)) |> 
  group_by(id) |>
  slice(1) |> 
  ungroup() |> 
  filter(is.na(stop_date)) |> 
  select(id, art_start = start_date, art = brand_name)

# Add current ART to list
active <- active |> 
  left_join(on_art, by = "id")

# Clinical diagnosis of diabetes in DB
dia_clinical <- clinical |> 
  filter(clin_id == "DIA") |> 
  select(id, diabetes = clin_id)

# HBA1c above 6.5
hba1c <- lab2 |> 
  filter(item == "HBA1C") |> 
  arrange(id, desc(value)) |> 
  group_by(id) |> 
  slice(1) |> 
  arrange(desc(value)) |> 
  # Some centers report HbA1c in mmol/mol, convert using formula
  mutate(value = if_else(value >= 24, (value * 0.0915) + 2.15, value)) |> 
  filter(value >= 6.5) |> 
  select(id, hba1c = value)

# Add info on diabetes
active <- active |> 
  left_join(dia_clinical, by = "id") |> 
  left_join(hba1c, by = "id") |> 
  mutate(diab = diabetes == "DIA" | hba1c >= 6.5) |> 
  mutate(diab = replace_na(diab, FALSE)) |> 
  select(-diabetes, -hba1c)

# CVD events: 
cvd_event <- clinical |> 
  filter(clin_id %in% c("AMI", "ANG", "CEI", "END", "PRO")) |> 
  arrange(id, desc(clin_date)) |> 
  group_by(id) |>
  slice(1) |> 
  ungroup() |> 
  select(id, last_cvd_date = clin_date)

# Add info on CVD events
active <- active |> 
  left_join(cvd_event, by = "id") |> 
  mutate(cvd = if_else(is.na(last_cvd_date), FALSE, TRUE))


# Select study population -------------------------------------------------

# I will create two study populations: 
# (1) strict version where we exclude participants with DM and CVD.
# (2) liberal version where we keep them as they are also eligible for statins.


# Select study population (strict = exclude diabetes and CVD)
strict <- active |> # N = 9369
  filter(!is.na(art_start)) |> # N = 9194
  filter(age >= 40, age <= 75) |> # N = 7803
  # filter(is.na(statin_start)) |> # N = 5612
  filter(!cvd) |> # N = 5320
  filter(!diab) # N = 4946

# Create table with number of participants and statin prescription rate by physician
tab_strict <- strict |> 
  group_by(center, source, physician) |> 
  summarise(n = n(), 
            n_statin = sum(!is.na(statin_name)), 
            n_elig_study = sum(is.na(statin_name))) |> 
  mutate(statin_presc_rate = n_statin / n)


# Select study population (liberal = keep diabetes and CVD)
liberal <- active |> # N = 9369
  filter(!is.na(art_start)) |> # N = 9194
  filter(age >= 40, age <= 75) # N = 7803

# Create table with number of participants and statin prescription rate by physician
tab_liberal <- liberal |> 
  group_by(center, source, physician) |> 
  summarise(n = n(), 
            n_statin = sum(!is.na(statin_name)), 
            n_elig_study = sum(is.na(statin_name))) |> 
  mutate(statin_presc_rate = n_statin / n)


# Write out data
write_rds(active, "processed/01-active_participants.rds")
write_rds(tab_strict, "processed/01-strict_selection.rds")    
write_rds(tab_liberal, "processed/01-liberal_selection.rds")    





