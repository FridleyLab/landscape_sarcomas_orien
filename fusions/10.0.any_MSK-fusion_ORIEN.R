rm(list=ls())
library(tidyverse)
library(openxlsx)

#clinical data
clinical = read.xlsx("Analysis_folders/M2Gen_Final_Clinical/ClinicalLinkagewithFiles_20221103_niceNames.xlsx")
clin2 = clinical %>% 
  filter(!is.na(somatic_file), #get only samples that we have mutations for
         tumor_germline == "Tumor", #remove the germline samples
         is.na(sarcoma)) %>% #make sure only to keep those that don't have a flag for not being sarcoma
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>% #extract the sample ID from the file name
  group_by(sarcoma_collapsed) %>% #group by sarcoma histology that Andrew Collapsed
  mutate(new_collapsed = ifelse(n() < 5, "other", sarcoma_collapsed), #if there are less than 5 samples, make new collapsed "other"
         nice_name_collapsed = ifelse(n() < 5, "Other", nice_name)) %>% #if there are less than 5 samples, group them to Other
  ungroup()
#msk fusion list
msk_fusion_list = read.xlsx("Analysis_folders/Fusions/data/MSK_fusion_reclassification.xlsx")

## Sample Fusion Data
fusion_data = readRDS("Analysis_folders/Fusions/data/All_fusion_table.rds")

#unique msk fusions
msk_fusions = lapply(msk_fusion_list$CuratedFusions, function(fs){
  fs %>% 
    gsub("\\;.*", "", .) %>% #remove the fusions in the "(A)" position since not using missing as evidence
    str_split(., ",") %>% #make vector for search
    unlist() %>%
    gsub("-", "--", .)
}) %>% unlist() %>% unique()

df = fusion_data[fusion_data$FusionName %in% msk_fusions,]
