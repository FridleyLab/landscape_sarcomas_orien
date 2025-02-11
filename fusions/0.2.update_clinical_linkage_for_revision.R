library(dplyr)
library(openxlsx)

clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/archive/ClinicalLinkagewithFiles_20230731_niceNames.xlsx") %>%
  filter(is.na(sarcoma))

clinical_key = clinical %>% 
  select(changed_diagnosis_clean, nice_name_reassigned) %>% 
  distinct()
clinical_key2 = clinical_key %>%
  mutate(cdc_revision = case_when(changed_diagnosis_clean == "liposarcoma" ~ "liposarcoma_NOS",
                                  changed_diagnosis_clean == "fibromyxosarcoma" ~ "myxofibrosarcoma",
                                  changed_diagnosis_clean == "fibrosarcoma" ~ "fibrosarcoma_NOS",
                                  changed_diagnosis_clean == "sarcoma" ~ "sarcoma_NOS",
                                  changed_diagnosis_clean == "myxosarcoma" ~ "sarcoma_NOS",
                                  changed_diagnosis_clean == "stromal_sarcoma" ~ "sarcoma_NOS",
                                  changed_diagnosis_clean == "stromal_tumor" ~ "sarcoma_NOS",
                                  changed_diagnosis_clean == "small_cell_sarcoma" ~ "sarcoma_NOS",
                                  changed_diagnosis_clean == "RMS" ~ "rhabdomyosarcoma",
                                  changed_diagnosis_clean == "myofibroblastic_tumor" ~ "inflammatory_myofibroblastic_tumor",
                                  TRUE ~ changed_diagnosis_clean)) %>%
  mutate(cdc_nice_name = case_when(cdc_revision == "sarcoma_NOS" ~ "Sarcoma, NOS",
                                   cdc_revision == "liposarcoma_well_differentiated" ~ "Well Differentiated Liposarcoma",
                                   cdc_revision == "myxofibrosarcoma" ~ "Myxofibrosarcoma",
                                   cdc_revision == "liposarcoma_NOS" ~ "Liposarcoma, NOS",
                                   cdc_revision == "fibrosarcoma_NOS" ~ "Fibrosarcoma, NOS",
                                   cdc_revision == "rhabdomyosarcoma" ~ "Rhabdomyosarcoma",
                                   cdc_revision == "inflammatory_myofibroblastic_tumor" ~ "Inflammatory Myofibroblastic Tumor",
                                   changed_diagnosis_clean == cdc_revision ~ nice_name_reassigned,
                                   TRUE ~ NA))
#save out for andrew to assess
#write.csv(clinical_key2, "Analysis_folders/1.0.Histology_Reassignment/andrew_histology_revision.csv")
#andrew check off on histology reclasses July 17, 2024
#ensure the noted histology subtypes are also added to the remove list
#atypical_fibrous_histiocytoma
#atypical_lipoma
#giant_cell_tumor_of_bone
#myxoma
#tenosynovial_giant_cell_tumor
#ossifying_fibromyxoid_tumor


clinical_key3 = clinical_key2 %>%
  mutate(reviewer_remove = ifelse(changed_diagnosis_clean %in% c("atypical_fibrous_histiocytoma",
                                                                 "atypical_lipoma",
                                                                 "giant_cell_tumor_of_bone",
                                                                 "myxoma",
                                                                 "tenosynovial_giant_cell_tumor",
                                                                 "ossifying_fibromyxoid_tumor"),
                                  TRUE,
                                  FALSE))
#merge back to the other original clinical information for downstream analyses
clinical_final = clinical %>%
  full_join(clinical_key3)
#save the file for use with all downstream analyses
wb = createWorkbook()
addWorksheet(wb, "reviewer_requested")
writeData(wb, 1, clinical_final)
saveWorkbook(wb, "Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20240716_niceNames.xlsx", overwrite = F)

