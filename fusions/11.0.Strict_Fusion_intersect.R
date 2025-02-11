rm(list=ls())
library(tidyverse)
library(openxlsx)

#clinical data
clinical = read.xlsx("Analysis_folders/M2Gen_Final_Clinical/ClinicalLinkagewithFiles_20221103_niceNames.xlsx")
clin2 = clinical %>% 
  filter(!is.na(somatic_file), #get only samples that we have mutations for - doesn't matter here as it's just for histology namess
         tumor_germline == "Tumor", #remove the germline samples
         is.na(sarcoma)) %>% #make sure only to keep those that don't have a flag for not being sarcoma
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>% #extract the sample ID from the file name
  group_by(sarcoma_collapsed) %>% #group by sarcoma histology that Andrew Collapsed
  mutate(new_collapsed = ifelse(n() < 5, "other", sarcoma_collapsed), #if there are less than 5 samples, make new collapsed "other"
         nice_name_collapsed = ifelse(n() < 5, "Other", nice_name)) %>% #if there are less than 5 samples, group them to Other
  ungroup()

msk_fusion_list = read.xlsx("Analysis_folders/Fusions/data/MSK_fusion_reclassification.xlsx")
msk_fusion_list = msk_fusion_list %>%
  mutate(From_filled = zoo::na.locf(zoo::zoo(msk_fusion_list$From)) %>% as.character(),
         .after = From) %>%
  #make the names matching between MSK and ORIEN
  mutate(
    nice_name_collapsed_from = case_when(
      From_filled == "Well/de-differentiated liposarcoma" ~ "Dedifferentiated Liposarcoma",
      From_filled == "Dermatofibrosarcoma protuberans" ~ "Dermatofibrosarcoma",
      From_filled == "Desmoplastic small round cell tumor" ~ "Desmoplastic Small Round Cell Tumor",
      From_filled == "Epithelioid sarcoma" ~ "Epithelioid Sarcoma",
      From_filled == "Ewing" ~ "Ewing Sarcoma",
      From_filled == "Liposarcoma, NOS" ~ "Liposarcoma",
      From_filled == "Malignant peripheral nerve sheath tumor (mpnst)" ~ "Malignant Peripheral Nerve Sheath Tumor",
      From_filled == "Myxoid liposarcoma" ~ "Myxoid Liposarcoma",
      From_filled == "Pleomorphic liposarcoma" ~ "Pleomorphic Liposarcoma",
      From_filled == "Rhabdomyosarcoma, NOS" ~ "RMS",
      From_filled == "Sarcoma, NOS" ~ "Sarcoma",
      From_filled == "Solitary fibrous tumor" ~ "Solitary Fibrous Tumor Malignant",
      From_filled == "Synovial sarcoma" ~ "Synovial Sarcoma",
      From_filled == "Undifferentiated pleomorphic sarcoma" ~ "Undifferentiated Pleomorphic Sarcoma",
      TRUE ~ From_filled
    ),
    nice_name_collapsed_to = case_when(
      To == "Well/de-differentiated liposarcoma" ~ "Dedifferentiated Liposarcoma",
      To == "Dermatofibrosarcoma protuberans" ~ "Dermatofibrosarcoma",
      To == "Ewing" ~ "Ewing Sarcoma",
      To == "Myxoid liposarcoma" ~ "Myxoid Liposarcoma",
      To == "Pleomorphic liposarcoma" ~ "Pleomorphic Liposarcoma",
      To == "Rhabdomyosarcoma, NOS" ~ "RMS",
      To == "Sarcoma, NOS" ~ "Sarcoma",
      To == "Synovial sarcoma" ~ "Synovial Sarcoma",
      TRUE ~ To
    )
  ) 
#add rows for liposarcoma well differentiated
#MSK groups these together, but this is "from" and are reclassified to something else.
#the issue might arise when something needs to be reclassified TO either dedifferentiated or well differentiated liposarcoma
msk_fusion_list = bind_rows(
  msk_fusion_list,
  msk_fusion_list[msk_fusion_list$nice_name_collapsed_from == "Dedifferentiated Liposarcoma",] %>%
    mutate(nice_name_collapsed_from = "Liposarcoma (Well Differentiated)")
)

## Sample Fusion Data
fusion_data = readRDS("Analysis_folders/Fusions/data/All_fusion_table.rds")

## Clean fusion list to only those that are in the ORIEN data
msk_fusion_list = msk_fusion_list %>%
  filter(nice_name_collapsed_from %in% unique(clin2$nice_name_collapsed),
         nice_name_collapsed_to %in% unique(clin2$nice_name_collapsed)) #remove this line for broader list

## comb through data
samples_to_be_reclassified = lapply(unique(msk_fusion_list$nice_name_collapsed_from), function(from){
  wkn_dat = msk_fusion_list %>% filter(nice_name_collapsed_from == from)
  out = lapply(unique(wkn_dat$nice_name_collapsed_to), function(to){
    fusions = wkn_dat$CuratedFusions[wkn_dat$nice_name_collapsed_to == to] %>% 
      gsub("\\;.*", "", .) %>% #remove the fusions in the "(A)" position since not using missing as evidence
      str_split(., ",") %>% #make vector for search
      unlist() %>%
      gsub("-", "--", .)
    
    sample_dat = fusion_data %>%
      filter(nice_name_collapsed == from)
    #if is 'from', and has fusion, belongs with 'to'
    data.frame(from, to, Samples_with_fusion = TRUE %in% (sample_dat$FusionName %in% fusions))
  }) %>%
    do.call(bind_rows, .)
  return(out)
}) %>% 
  do.call(bind_rows, .)

#which samples to reclassify
from_to_combinations = samples_to_be_reclassified %>%
  filter(Samples_with_fusion == TRUE) %>%
  mutate(to = ifelse(to == "Dedifferentiated Liposarcoma", "Well/De-differentiated Liposarcoma", to))
#the 'to' of dedifferentiated liposarcoma is well/de-differentiated liposarcoma
reclassify_fusions = lapply(1:nrow(from_to_combinations), function(r){
  fusions = msk_fusion_list %>% filter(nice_name_collapsed_from == from_to_combinations[r, "from"],
                                       nice_name_collapsed_to == from_to_combinations[r, "to"]) %>%
    pull(CuratedFusions) %>% 
    gsub("\\;.*", "", .) %>% #remove the fusions in the "(A)" position since not using missing as evidence
    str_split(., ",") %>% #make vector for search
    unlist() %>%
    gsub("-", "--", .)
  
  sample_dat = fusion_data %>%
    filter(nice_name_collapsed == from_to_combinations[r, "from"]) %>%
    filter(FusionName %in% fusions) %>%
    mutate(Proposed_Reclassification = from_to_combinations[r, "to"])
  
  return(sample_dat)
}) %>% 
  do.call(bind_rows, .)

wb = createWorkbook()
addWorksheet(wb, sheetName = "ORIEN MSK from intersect")
writeData(wb, 1, reclassify_fusions)
saveWorkbook(wb, "Analysis_folders/Fusions/MSK_fusion_reclassification/Strict_Fusion_intersect.xlsx")