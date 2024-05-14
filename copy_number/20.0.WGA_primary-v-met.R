rm(list=ls())
library(tidyverse)
library(openxlsx)

#data
cn_data = readRDS("Analysis_folders/CopyNumberVariation/segmentationFiles/data/copy_number_level-by-length-noy.rds")
clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20230731_niceNames.xlsx") %>%
  filter(is.na(sarcoma), !is.na(somatic_file))
clinical$segment_file = str_split(clinical$somatic_file, pattern = "_") %>% 
  do.call(rbind, .) %>% data.frame() %>%
  rowwise() %>% mutate(F = paste0(X1, "_", X4, "_segments.txt")) %>% pull(F)
clin2 = clinical %>%
  filter(!is.na(somatic_file), #make sure that we have the somatic mutaiton file for the sample
         tumor_germline == "Tumor", #remove germline from clinical file
         is.na(sarcoma)) %>% #remove 'non-sarcoma' samples
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>% #remove everything after "." to get sample ID
  group_by(changed_diagnosis_clean) %>% #group by the sarcoma histology Andrew collapsed
  mutate(new_collapsed = ifelse(n() < 5, "other", changed_diagnosis_clean), #if there are less than 5 samples, make new collapsed "other" using the new histology subtype categoties
         nice_name_collapsed = ifelse(n() < 5, "Other", nice_name_reassigned)) %>% #with a nice name called "Other" for plotting
  ungroup()

#summarising WGA v not
cn_data2 = cn_data %>%
  filter(Tumor_Sample_Barcode %in% clin2$Tumor_Sample_Barcode) %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(All_Segments_Length = sum(Length),
         cn_proportion = Length / All_Segments_Length) %>%
  filter(cn == "AMP") %>%
  mutate(WGA = ifelse(cn_proportion > 0.50, "Yes", "No")) %>%
  left_join(clin2)
table(cn_data2$WGA, cn_data2$primary_met)

#chi squared test
tab = cn_data2 %>%
  janitor::tabyl(WGA, primary_met)
res = tab %>%
  janitor::chisq.test()
tab %>%
  gather("Site", "Samples", -WGA) %>%
  group_by(Site) %>%
  mutate(Freq = Samples / sum(Samples) * 100)

#contingency tables for top histology subtypes
#all samples
cn_data2 %>%
  janitor::tabyl(WGA, primary_met)
#find most sample histology subtypes
top_hists = cn_data2 %>%
  group_by(nice_name_collapsed) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  slice(1:5) %>% pull(nice_name_collapsed)

top_hists_res = lapply(setNames(top_hists, top_hists), function(x){
  tab = cn_data2 %>%
    filter(nice_name_collapsed == x) %>%
    janitor::tabyl(WGA, primary_met)
  res = tab %>%
    janitor::chisq.test()
  return(list(Table = tab,
              Results = res))
})
#results
top_hists_res[1]
top_hists_res[2]
top_hists_res[3]
top_hists_res[4]
top_hists_res[5]

#gist frequencies
top_hists_res$`Gastrointestinal Stromal Tumor`$Table %>%
  gather("Site", "Samples", -WGA) %>%
  group_by(Site) %>%
  mutate(Freq = Samples / sum(Samples) * 100)


