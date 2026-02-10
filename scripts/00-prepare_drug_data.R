library(bernr)
library(tidyverse)

# Gather a drug and ART dataset

# read current SHCS data download
shcs_read <- function(file) {
  haven::read_dta(here::here("data", "2601stata", file))
}

# Load data
med_treatment <- shcs_read("med_treatment.dta")
med_product <- shcs_read("med_product.dta")
med_drug_code <- shcs_read("med_drug_code.dta")
med_category <- shcs_read("med_category.dta")
med_substance_in_product <- shcs_read("med_substance_in_product.dta")
med_substance <- shcs_read("med_substance.dta")


# Put them together
step1 <- med_treatment %>% 
  select(id, product_id, start_date, start_acc = start_date_accuracy, 
         stop_date, stopp_add = stop_date_accuracy)

step2 <- step1 %>% 
  left_join(
    med_substance_in_product %>% 
      select(product_id = containing_product_id, 
             substance_id = contained_substance_id), 
    by = "product_id") %>% 
  left_join(
    med_substance %>% select(substance_id, drug_code), 
    by = "substance_id"
  )


all_drugs <- step2 %>% 
  left_join(med_product %>% select(product_id, brand_name, atc)) %>% 
  left_join(med_drug_code %>% select(-description)) %>% 
  left_join(med_category %>% rename(category_id = code))

art_drugs <- all_drugs %>% 
  filter(indication == "A") 


# Write out
write_rds(all_drugs, "data/2601stata/00-all_drugs.rds")
write_rds(art_drugs, "data/2601stata/00-art_drugs.rds")
