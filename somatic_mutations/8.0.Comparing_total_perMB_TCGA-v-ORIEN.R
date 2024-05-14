# libraries ---------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(openxlsx)


# data --------------------------------------------------------------------

clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20230731_niceNames.xlsx")
clin2 = clinical %>%
  filter(!is.na(somatic_file),
         tumor_germline == "Tumor",
         is.na(sarcoma)) %>%
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>%
  group_by(changed_diagnosis_clean) %>%
  mutate(changed_diagnosis_collapsed = ifelse(n() < 5, "other", changed_diagnosis_clean),
         nice_name_reassigned_collapsed = ifelse(n() < 5, "other", nice_name_reassigned)) %>%
  ungroup()
all_studies_total_perMB = lapply(list.files("Analysis_folders/FilteringVCFfiles/MutationPerMB10", full.names = T, pattern = "csv$"), 
                                 read.csv, check.names = F) %>%
  do.call(bind_rows, .) %>% distinct() #since all will have the TCGA, remove duplicates
#get the updated histology subtype assignments
sum(all_studies_total_perMB$Tumor_Sample_Barcode %in% clin2$Tumor_Sample_Barcode)
#result is 1168 because there are 2 samples with none

#merging the clinical (all samples) and the tumor mutation burden
sarc_dat = clin2 %>%
  #rename("cohort" = 2) %>%
  full_join(all_studies_total_perMB)  %>% 
  filter(!is.na(wes) | cohort == "SARC") %>%
  select(wes, Tumor_Sample_Barcode, total, cohort, total_perMB, Capture_Region)
#number of avatar members
sum(sarc_dat$Tumor_Sample_Barcode %in% clin2$Tumor_Sample_Barcode)
#1170 - all samples

#calculate burden
sarc_summary = sarc_dat %>%
  mutate(total_perMB = ifelse(is.na(total_perMB), 0, total_perMB),
         burden = case_when(total_perMB < 5 ~ "Low",
                            total_perMB >= 5 & total_perMB < 10 ~ "Intermediate",
                            total_perMB >= 10 ~ "High"),
         burden = factor(burden, levels = c('Low', "Intermediate", "High")),
         Study = ifelse(is.na(wes), "TCGA", "ORIEN")) %>%
  group_by(Study, burden) %>%
  summarise(samples = n()) %>%
  spread(Study, samples) %>%
  column_to_rownames('burden')
#print the table
sarc_summary
#low inter %
(1134 + 14) / (1134 + 14 + 22)
(226 + 5) / (226 + 5 + 8)
#do fishers exact
fisher.test(sarc_summary)

#TMB range in ORIEN
sarc_dat %>%
  filter(cohort != "SARC" | is.na(cohort)) %>% 
  summarise(range(total_perMB, na.rm = TRUE))
#minimum is 0 from the two samples, max is 56.6

#percent in each group
sapply(colnames(sarc_summary), function(cn){
  v = sarc_summary[,cn]
  names(v) = rownames(sarc_summary)
  round(v/sum(v)*100, 2)
})

#median by histology
sarc_dat %>%
  filter(cohort != "SARC") %>% #1168 - 2 with no somatic mutations
  select(Tumor_Sample_Barcode, total_perMB) %>%
  right_join(clin2 %>%
              select(Tumor_Sample_Barcode, nice_name_reassigned_collapsed)) %>%
  group_by(nice_name_reassigned_collapsed) %>%
  mutate(total_perMB = ifelse(is.na(total_perMB), 0, total_perMB)) %>%
  summarise(median = median(total_perMB),
            min = min(total_perMB),
            max = max(total_perMB)) %>%
  arrange(desc(median))


