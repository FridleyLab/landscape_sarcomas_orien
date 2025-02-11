rm(list=ls())
gc(full=TRUE)
# libraries ---------------------------------------------------------------
library(tidyverse)
library(openxlsx)

# data --------------------------------------------------------------------

clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20240716_niceNames.xlsx")
clin2 = clinical %>%
  filter(!is.na(somatic_file),
         tumor_germline == "Tumor",
         is.na(sarcoma),
         !reviewer_remove) %>%
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>%
  group_by(cdc_revision) %>%
  mutate(changed_diagnosis_collapsed = ifelse(n() < 5, "other", cdc_revision),
         nice_name_reassigned_collapsed = ifelse(n() < 5, "other", cdc_nice_name)) %>%
  ungroup()
all_studies_total_perMB = lapply(list.files("Analysis_folders/FilteringVCFfiles/MutationPerMB10", full.names = T, pattern = "csv$"), 
                                 read.csv, check.names = F) %>%
  do.call(bind_rows, .) %>% distinct() #since all will have the TCGA, remove duplicates
#get the updated histology subtype assignments
sum(all_studies_total_perMB$Tumor_Sample_Barcode %in% clin2$Tumor_Sample_Barcode)
#result is 1160 because there are 2 samples with none

#merging the clinical (all samples) and the tumor mutation burden
all_studies_total_perMB = clin2 %>%
  rename("cohort" = 2) %>%
  left_join(all_studies_total_perMB %>%
              filter(!grepl("TCGA", Tumor_Sample_Barcode)) %>%
              select(-cohort)) %>%
  relocate(cohort, .after = total)

#identify high/middle/low thresholds 
vec = log2(all_studies_total_perMB$total_perMB)
quants = quantile(vec, c(0.33, 0.66), na.rm = TRUE)

#add back to the sample data frmae
analysis_data = all_studies_total_perMB %>%
  mutate(total_perMB = ifelse(is.na(total_perMB), 0.000001, total_perMB)) %>% #for 2 samples with none, add tiny tiny value - will be LOW TMB load
  mutate(TMB_load = case_when(log2(total_perMB) < quants['33%'] ~ "LOW",
                              log2(total_perMB) >= quants['33%'] & log2(total_perMB) < quants['66%'] ~ "MEDIUM",
                              log2(total_perMB) >= quants['66%'] ~ "HIGH"),
         TMB_load = factor(TMB_load, levels = c("LOW", "MEDIUM", "HIGH")),
         primary_met = factor(primary_met, level = c("Primary", "Metastatic"))) %>%
  select(orien_avatar_key, primary_met, total_perMB, TMB_load, nice_name_reassigned_collapsed)


# plotting ----------------------------------------------------------------

analysis_data %>%
  ggplot() +
  geom_violin(aes(x = TMB_load, y = total_perMB, fill = primary_met)) +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Tumor Mutation Burden by Site",
       y = "Tumor Mutation Burden\n(per MB)",
       x = "Mutation Burden Load") +
  guides(fill = guide_legend(title = "Primary/\nMetastatic"))

#looking at patients with primary and met from same disease type
analysis_data %>%
  group_by(orien_avatar_key) %>% 
  filter(n() == 2, 
         length(unique(nice_name_reassigned_collapsed)) == 1, 
         length(unique(primary_met)) == 2) %>%
  select(-TMB_load) %>%
  spread(primary_met, total_perMB) %>%
  ggplot() +
  geom_point(aes(x = Primary, y = Metastatic)) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans = "log10") +
  coord_equal() +
  lims(x = c(0, NA),
       y = c(0, NA))


# comparing ---------------------------------------------------------------

#how many patients have multiple samples
analysis_data %>%
  group_by(orien_avatar_key) %>%
  filter(n() > 1) %>%
  pull(orien_avatar_key) %>%
  unique() %>%
  length()
#51 patients have multiple samples

tmp = analysis_data %>%
  group_by(primary_met, TMB_load) %>%
  summarise(samples = n()) %>%
  spread(TMB_load, samples) %>%
  column_to_rownames("primary_met")
mosaicplot(tmp,
           main = "TMB Load")  
chisq.test(tmp)
chisq.test(t(tmp))$expected
#looks like the composition of TMB load is not independent of tumor site

wilcox.test(total_perMB ~ primary_met, data = analysis_data)

#identify high vs low based on TCGA classification per histology per tumor site
interest_hists = all_studies_total_perMB %>%
  mutate(total_perMB = ifelse(is.na(total_perMB), 0.000001, total_perMB)) %>%
  select(nice_name_reassigned_collapsed, primary_met, total_perMB) %>%
  group_by(nice_name_reassigned_collapsed, primary_met) %>%
  summarise(n = n()) %>%
  filter(n >= 20) %>%
  filter(length(unique(primary_met)) == 2) %>%
  pull(nice_name_reassigned_collapsed) %>% 
  unique()

all_studies_total_perMB %>%
  mutate(total_perMB = ifelse(is.na(total_perMB), 0.000001, total_perMB)) %>%
  select(nice_name_reassigned_collapsed, primary_met, total_perMB) %>%
  filter(nice_name_reassigned_collapsed %in% interest_hists) %>%
  arrange(nice_name_reassigned_collapsed, primary_met) %>%
  group_by(nice_name_reassigned_collapsed) %>%
  group_map(~ {
    tmp = wilcox.test(total_perMB ~ primary_met, data = .x)
    data.frame('nice_name_reassigned_collapsed' = unique(.y$nice_name_reassigned_collapsed),
               'statistic' = tmp$statistic,
               'p-value' = tmp$p.value)
  }) %>%
  do.call(bind_rows, .)

all_studies_total_perMB %>%
  mutate(total_perMB = ifelse(is.na(total_perMB), 0.000001, total_perMB)) %>%
  select(nice_name_reassigned_collapsed, primary_met, total_perMB) %>%
  filter(nice_name_reassigned_collapsed %in% interest_hists) %>%
  arrange(nice_name_reassigned_collapsed, primary_met) %>%
  ggplot() + 
  geom_violin(aes(x = nice_name_reassigned_collapsed, color = primary_met, y = total_perMB))










