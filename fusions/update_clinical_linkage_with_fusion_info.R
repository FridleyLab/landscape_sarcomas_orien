library(tidyverse)
library(openxlsx)

#andrews fusion list histology reassignment
reassigned = read.xlsx("Analysis_folders/Fusions/MSK_fusion_reclassification/Strict_Fusion_intersect_AB_06_14_23.xlsx", 
                       sheet = 7, check.names = F)
#previously explored and ironed out clinical file with Dale
nice_name_clinical = read.xlsx("Analysis_folders/M2Gen_Final_Clinical/ClinicalLinkagewithFiles_20221103_niceNames.xlsx")

#now match those samples that are in the fusion list being reassigned
reassigned2 = reassigned %>%
  select(orien_avatar_key, tsb, wes, sarcoma_collapsed, changed_diagnosis) %>%
  distinct() %>% 
  mutate(changed_diagnosis = ifelse(sarcoma_collapsed == "gastrointestinal_stromal_tumor", 
                                    "gastrointestinal_stromal_tumor",
                                    changed_diagnosis)) %>%
  mutate(changed_diagnosis_clean = changed_diagnosis %>% gsub(" \\(.*", "", .) %>% gsub(" ", "", .)) %>%
  filter(changed_diagnosis_clean != "n/a") %>%
  rename("rna_seq" = tsb)

#match the reassigned samples wtth the origian
nice_name_clinical2 = nice_name_clinical %>%
  full_join(reassigned2) %>%
  select(-X1) %>%
  mutate(changed_diagnosis_clean = ifelse(is.na(changed_diagnosis_clean), sarcoma_collapsed, changed_diagnosis_clean)) %>%
  mutate(nice_name_reassigned = ifelse(is.na(changed_diagnosis), nice_name, NA))
#collect the histology subtypes that are in the data and list them between the collapsed name and the nice name
nice_name_key = nice_name_clinical2 %>%
  filter(!is.na(nice_name_reassigned)) %>%
  select(changed_diagnosis_clean, nice_name_reassigned) %>%
  distinct() %>%
  pull(changed_diagnosis_clean, nice_name_reassigned)
#Add the nice name of the reassigned samples if other samples in the cohort have the same 'collapsed' name
nice_name_clinical3 = nice_name_clinical2 %>%
  mutate(nice_name_reassigned = ifelse(is.na(nice_name_reassigned), 
                                       names(nice_name_key)[match(changed_diagnosis_clean, nice_name_key)],
                                       nice_name_reassigned))
#manually add in the 2 histology subtypes that are new after reassigning samples due to fusion information
nice_name_clinical3$nice_name_reassigned[nice_name_clinical3$changed_diagnosis_clean == "ossifying_fibromyxoid_tumor"] = "Ossifying Fibromyxoid Tumor"
nice_name_clinical3$nice_name_reassigned[nice_name_clinical3$changed_diagnosis_clean == "round_cell_sarcoma"] = "Round Cell Sarcoma"

#save the file for use with all downstream analyses
wb = createWorkbook()
addWorksheet(wb, "fusion_reassigned")
writeData(wb, 1, nice_name_clinical3)
saveWorkbook(wb, "Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20230731_niceNames.xlsx", overwrite = F)



