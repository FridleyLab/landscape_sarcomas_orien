rm(list=ls())

#libraries
library(tidyverse)
library(openxlsx)
library(ggpubr)

#clinical data that we are working with
clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20240716_niceNames.xlsx") %>%
  mutate(preservation_method = case_when(preservation_method %in% c("Snap Frozen", "Frozen") ~ "Fresh Frozen",
                                         TRUE ~ preservation_method))
survival_table = read.xlsx("data/RNAseq/oo/SAR.RIG.dx.table.xlsx") %>%
  mutate(Survival = RightAnchor - LeftAnchor.p)
diagnosis_table = read.csv("data/clinical/19PRJ024RIG-SAR_NormalizedFiles/19PRJ024RIG-SAR_20220718_Diagnosis_V4.csv")

#subset out the rna and wes samples that we used for analysis
rna = clinical %>%
  filter(!is.na(rna_file), is.na(sarcoma), !is.na(preservation_method), !reviewer_remove) %>%
  group_by(cdc_nice_name) %>%
  mutate(nice_name_reassigned_collapsed = ifelse(n() < 5, "Other", cdc_nice_name)) %>%
  group_by(nice_name) %>%
  mutate(nice_name_collapsed = ifelse(n() < 5, "Other", nice_name)) %>% ungroup()

wes = clinical %>%
  filter(is.na(sarcoma), !is.na(somatic_file), !reviewer_remove) %>% #remove 'non-sarcoma' samples
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>% #remove everything after "." to get sample ID
  group_by(cdc_nice_name) %>% #group by the sarcoma histology Andrew collapsed
  mutate(nice_name_reassigned_collapsed = ifelse(n() < 5, "Other", cdc_nice_name)) %>% #with a nice name called "Other" for plotting
  ungroup()

#combine data back together
clin2 = full_join(
  rna,
  wes
)

total_tumors = clin2 = full_join(
  rna %>% select(-contains("collapsed")),
  wes %>% select(-contains("collapsed"))
)
#unique avatars
length(unique(total_tumors$orien_avatar_key))

# # Begin Summarizing -------------------------------------------------------
# #first total cohort level
# #when moving to histology then will have to go to individual data types because 'other' had different samples depending on wes and rna
# #unique avatar IDs
# length(unique(clin2$orien_avatar_key))
# #avatars with primary/met/both
# clin2 %>%
#   select(orien_avatar_key, primary_met) %>%
#   distinct() %>%
#   mutate(has = "Yes") %>%
#   spread("primary_met", "has", fill = "No") %>%
#   mutate(primary_met = case_when(Metastatic == "No" & Primary == "Yes" ~ "Primary",
#                                  Metastatic == "Yes" & Primary == "No" ~ "Metastatic",
#                                  Metastatic == "Yes" & Primary == "Yes" ~ "Both")) %>%
#   group_by(primary_met) %>%
#   summarise(n = n())
# #avatars with wes/rna/both on a sample
# #avatar might have multiple samples
# tmp = clin2 %>%
#   select(orien_avatar_key, somatic_file, rna_file, preservation_method) %>%
#   mutate(wes = ifelse(is.na(somatic_file), "No", "Yes"),
#          rna = ifelse(is.na(rna_file), "No", "Yes")) %>% 
#   select(orien_avatar_key, wes, rna, preservation_method) %>% distinct() %>%
#   mutate(both = ifelse(wes == "Yes" & rna == "Yes", "Yes", "No"))
# #wes ids
# wes_ids = tmp %>% filter(wes == "Yes") %>% distinct() %>% pull(orien_avatar_key) %>% unique()
# #rna ids
# rna_ids = tmp %>% filter(rna == "Yes", !is.na(preservation_method)) %>% distinct() %>% pull(orien_avatar_key) %>% unique()#remove the unknown preservation methos
# #number of patients with at least 1 samples with wes and 1 sample with rna, don't have to be same sample because patient level
# sum(wes_ids %in% rna_ids)
# #number of patients no wes in any sample and at least 1 sample with rna
# length(rna_ids[!(rna_ids %in% wes_ids)])
# #number of patients with no rna in any sample and at least 1 sample with wes
# length(wes_ids[!(wes_ids %in% rna_ids)])
# 
# #wes samples per avatar range
# clin2 %>%
#   filter(!is.na(somatic_file)) %>%
#   select(orien_avatar_key, somatic_file) %>%
#   distinct() %>%
#   group_by(orien_avatar_key) %>%
#   summarise(n = n()) %>%
#   summarise(min = min(n),
#             max = max(n))
# #rna samples per avatar range
# clin2 %>%
#   filter(!is.na(rna_file), !is.na(preservation_method)) %>% #have to remove the rna sample again with no preservation method kept for WES
#   select(orien_avatar_key, rna_file) %>%
#   distinct() %>%
#   group_by(orien_avatar_key) %>%
#   summarise(n = n()) %>%
#   summarise(min = min(n),
#             max = max(n))
# 
# #wes preservation methods
# clin2 %>%
#   filter(!is.na(somatic_file)) %>%
#   select(orien_avatar_key, somatic_file, preservation_method) %>%
#   distinct() %>%
#   group_by(preservation_method) %>%
#   summarise(n = n())
# #rna preservation methods
# clin2 %>%
#   filter(!is.na(rna_file), !is.na(preservation_method)) %>% #have to remove that NA preservation method, row kept because WES
#   select(orien_avatar_key, rna_file, preservation_method) %>%
#   distinct() %>%
#   group_by(preservation_method) %>%
#   summarise(n = n())
# 
# sum(wes$wes %in% survival_table$WES)
# 
# 

# WES histology level -----------------------------------------------------

#primary_met per hist
wes %>%
  select(orien_avatar_key, nice_name_reassigned_collapsed, primary_met) %>%
  group_by(nice_name_reassigned_collapsed) %>%
  mutate(Avatars = length(unique(orien_avatar_key))) %>%
  group_by(Avatars, primary_met, .add = TRUE) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = 'primary_met', values_from = 'n', values_fill = 0) %>% View()

#samples per preservation method per hist
wes %>%
  select(orien_avatar_key, nice_name_reassigned_collapsed, preservation_method) %>%
  group_by(nice_name_reassigned_collapsed, preservation_method) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = 'preservation_method', values_from = 'n', values_fill = 0) %>% View()

tmp = survival_table %>%
  filter(WES %in% wes$wes)

tmp2 = diagnosis_table %>%
  filter(!grepl("Age|Unk", AgeAtDiagnosis))

tmp3 = tmp %>%
  left_join(tmp2 %>%
              mutate(AgeAtDiagnosis = as.numeric(AgeAtDiagnosis)),
            c("AvatarKey", 
              "PrimaryDiagnosisSiteCode.p" = "PrimaryDiagnosisSiteCode", 
              "AgeAtDiagnosis.p" = "AgeAtDiagnosis")) %>%
  select(-contains(".e")) %>%
  filter(dx2sc.p < 5) %>%
  mutate(stage_clean=case_when(str_detect(PathGroupStage, "Unknown|TNM|Occult") ~ 'unknown', # Clean and collapse stage information
                               str_detect(PathGroupStage, "^I{1}[ABCS]") ~ 'I', 
                               str_detect(PathGroupStage, "^I{2}[ABCS]") ~ 'II',
                               str_detect(PathGroupStage, "^I{3}[ABCS]") ~ 'III',
                               str_detect(PathGroupStage, "^IV[ABCS]") ~ 'IV',
                               str_detect(PathGroupStage, "0") ~ '0',
                               TRUE ~ PathGroupStage)) %>%
  mutate(stage_clean_recode=recode(stage_clean, '0'=0, 'I'=1, 'II'=2, 'III'=3, 'IV'=4, 'unknown'=-1)) %>% # Recode stages 
  right_join(wes %>% 
               select(wes, nice_name_reassigned_collapsed), 
             ., c('wes' = "WES"))

multisensored = tmp3 %>%
  group_by(AvatarKey) %>%
  summarise(Censorship = length(unique(Censorship)))
unique(multisensored$Censorship)
#whoohoo all have one censorship

tmp3 %>%
  group_by(AvatarKey) %>%
  arrange(desc(stage_clean_recode)) %>%
  slice(1) %>%
  group_by(nice_name_reassigned_collapsed) %>%
  summarise(average = mean(Survival),
            std = sd(Survival)) %>%
  mutate(values = paste0(round(average, 2), " (", round(std, 2), ")")) %>% View()


# RNA histology level -----------------------------------------------------

#primary_met per hist
rna %>%
  select(orien_avatar_key, nice_name_reassigned_collapsed, primary_met) %>%
  group_by(nice_name_reassigned_collapsed) %>%
  mutate(Avatars = length(unique(orien_avatar_key))) %>%
  group_by(Avatars, primary_met, .add = TRUE) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = 'primary_met', values_from = 'n', values_fill = 0) %>% View()

#samples per preservation method per hist
rna %>%
  select(orien_avatar_key, nice_name_reassigned_collapsed, preservation_method) %>%
  group_by(nice_name_reassigned_collapsed, preservation_method) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = 'preservation_method', values_from = 'n', values_fill = 0) %>% View()

tmp = survival_table %>%
  filter(RNASeq %in% rna$rna_seq)

tmp2 = diagnosis_table %>%
  filter(!grepl("Age|Unk", AgeAtDiagnosis))

tmp3 = tmp %>%
  left_join(tmp2 %>%
              mutate(AgeAtDiagnosis = as.numeric(AgeAtDiagnosis)),
            c("AvatarKey", 
              "PrimaryDiagnosisSiteCode.p" = "PrimaryDiagnosisSiteCode", 
              "AgeAtDiagnosis.p" = "AgeAtDiagnosis")) %>%
  select(-contains(".e")) %>%
  filter(dx2sc.p < 5) %>%
  mutate(stage_clean=case_when(str_detect(PathGroupStage, "Unknown|TNM|Occult") ~ 'unknown', # Clean and collapse stage information
                               str_detect(PathGroupStage, "^I{1}[ABCS]") ~ 'I', 
                               str_detect(PathGroupStage, "^I{2}[ABCS]") ~ 'II',
                               str_detect(PathGroupStage, "^I{3}[ABCS]") ~ 'III',
                               str_detect(PathGroupStage, "^IV[ABCS]") ~ 'IV',
                               str_detect(PathGroupStage, "0") ~ '0',
                               TRUE ~ PathGroupStage)) %>%
  mutate(stage_clean_recode=recode(stage_clean, '0'=0, 'I'=1, 'II'=2, 'III'=3, 'IV'=4, 'unknown'=-1)) %>% # Recode stages 
  right_join(rna %>% 
               select(rna_seq, nice_name_reassigned_collapsed), 
             ., c('rna_seq' = "RNASeq"))

multisensored = tmp3 %>%
  group_by(AvatarKey) %>%
  summarise(Censorship = length(unique(Censorship)))
unique(multisensored$Censorship)
#whoohoo all have one censorship

tmp3 %>%
  group_by(AvatarKey) %>%
  arrange(desc(stage_clean_recode)) %>%
  slice(1) %>%
  group_by(nice_name_reassigned_collapsed) %>%
  summarise(average = mean(Survival),
            std = sd(Survival)) %>%
  mutate(values = paste0(round(average, 2), " (", round(std, 2), ")")) %>% View()



# samples per avatar histograms -------------------------------------------

bind_rows(wes %>%
            group_by(orien_avatar_key) %>%
            summarise(samples = n()) %>%
            mutate(Data = "WES"),
          rna %>%
            group_by(orien_avatar_key) %>%
            summarise(samples = n()) %>%
            mutate(Data = "RNA")
) %>%
  ggplot() +
  geom_bar(aes(x = samples, fill = Data), color = "black", position = "dodge") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(x = "Samples", y = "Avatars", title = "Samples per Avatar Patient") +
  scale_fill_manual(values = c("white", "gray"))



# samples per histology (not collapsed) -----------------------------------

clin3 = full_join(
  rna %>% select(-nice_name_reassigned_collapsed, -nice_name_collapsed),
  wes %>% select(-nice_name_reassigned_collapsed)
)

clin3 %>% 
  select(orien_avatar_key, preservation_method) %>%
  distinct() %>%
  group_by(preservation_method) %>%
  summarise(n())

clin3 %>% 
  select(orien_avatar_key, cdc_nice_name) %>%
  distinct() %>%
  group_by(cdc_nice_name) %>%
  summarise(n()) %>%
  View()

wes_other = wes %>%
  filter(nice_name_reassigned_collapsed == "Other") %>%
  group_by(cdc_nice_name) %>%
  summarise(n = n())
wes_other %>%
  unite("cleaner", cdc_nice_name:n, sep=" (n = ") %>%
  pull(cleaner) %>%
  paste0(., collapse = "), ")
sum(wes_other$n)

rna_other = rna %>%
  filter(nice_name_reassigned_collapsed == "Other") %>%
  group_by(cdc_nice_name) %>%
  summarise(n = n())
rna_other %>%
  unite("cleaner", cdc_nice_name:n, sep=" (n = ") %>%
  pull(cleaner) %>%
  paste0(., collapse = "), ")
sum(rna_other$n)





# #patients by histology subtype
# tab1 = clin2 %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed) %>%
#   group_by(nice_name_reassigned_collapsed) %>%
#   distinct() %>%
#   summarise(`Unique Avatars by histology subtype` = n())
# #patients with both wes and rna on a sample by histology subtype
# tab2 = clin2 %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed, wes, rna_seq, somatic_file, rna_file) %>%
#   filter(!is.na(somatic_file), !is.na(rna_file)) %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed) %>%
#   group_by(nice_name_reassigned_collapsed) %>%
#   distinct() %>%
#   summarise(`Avatars with samples that have both WES and RNA` = n())
# #patients with wes but not rna on a sample by histology subtype
# tab3 = clin2 %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed, wes, rna_seq, somatic_file, rna_file) %>%
#   filter(!is.na(somatic_file), is.na(rna_file)) %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed) %>%
#   group_by(nice_name_reassigned_collapsed) %>%
#   distinct() %>%
#   summarise(`Avatars with samples that have WES but not RNA` = n())
# #patients with wes but not rna on a sample by histology subtype
# tab4 = clin2 %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed, wes, rna_seq, somatic_file, rna_file) %>%
#   filter(is.na(somatic_file), !is.na(rna_file)) %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed) %>%
#   group_by(nice_name_reassigned_collapsed) %>%
#   distinct() %>%
#   summarise(`Avatars with samples that have RNA but not WES` = n())
# 
# #samples of wes by histology subtype
# #gets fishy with assignments because other may be different between wes and RNA
# tab5 = clin2 %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed, wes, somatic_file) %>%
#   filter(!is.na(somatic_file),) %>%
#   select(orien_avatar_key, nice_name_reassigned_collapsed) %>%
#   group_by(nice_name_reassigned_collapsed) %>%
#   distinct() %>%
#   summarise(`Avatars with samples that have RNA but not WES` = n())
# 
# 
# out = Reduce(full_join, lapply(grep("tab[0-9]", ls(), value = TRUE), get)); View(out)
# 
# 


# Apr 3, 2024 - Andrews update to tumor level instead of avatar le --------
#have to set the rna file without preservation method to NA
total_tumors2 = total_tumors %>%
  mutate(rna_file = ifelse(is.na(preservation_method), NA, rna_file),
         preservation_method = ifelse(is.na(preservation_method), "Formalin-Fixed", preservation_method)) %>%
  group_by(cdc_nice_name) %>%
  mutate(nice_name_reassigned_collapsed = ifelse(n() < 5, 'Other', cdc_nice_name))

total_tumors2 %>% janitor::tabyl(primary_met) %>% mutate(percent = percent * 100) %>% View()
total_tumors2 %>% janitor::tabyl(preservation_method) %>% mutate(percent = percent * 100) %>% View()
total_tumors2 %>%
  mutate(molecular = case_when(!is.na(rna_file) & !is.na(somatic_file) ~ "Both",
                               !is.na(rna_file) & is.na(somatic_file) ~ "RNA",
                               is.na(rna_file) & !is.na(somatic_file) ~ "WES")) %>%
  janitor::tabyl(molecular) %>% mutate(percent = percent * 100) %>% View()
total_tumors2 %>% janitor::tabyl(nice_name_reassigned_collapsed) %>% mutate(percent = percent * 100) %>% View()

total_tumors2 %>% 
  janitor::tabyl(cdc_nice_name) %>%
  filter(n<5) %>%
  unite("cleaner", cdc_nice_name:n, sep=" (n = ") %>%
  pull(cleaner) %>%
  paste0(., collapse = "), ")
