rm(list=ls())
library(tidyverse)
library(openxlsx)
library(maftools)

# data --------------------------------------------------------------------
# External Files ----------------------------------------------------------
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
filtered_annotated = list.files("Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/",
                                pattern = "vcf.gz$", full.names = T)
ordered_filtered_annotated = lapply(clin2$wes, function(v) grep(v, filtered_annotated, value = T)) %>% unlist()
# maf = annovarToMaf(ordered_filtered_annotated)

cosmic_genes = read.csv("Analysis_folders/Significantly_Mutated_Genes/Census_all_COSMIC.csv")
robust_regression = readRDS("Analysis_folders/Significantly_Mutated_Genes/Robust_Regression_object_reassigned.rds")
#standardized residuals are standard deviations so can use to filter genes
tumor_mutation_burden = read.csv("Analysis_folders/TumorMutationBurden/archive/Somatic_MutationsPerMB_10rFiltered.csv",
                                 row.names = NULL)

# Finding Genes to Plot ---------------------------------------------------

sig_mutated = lapply(names(robust_regression[1:30]), function(nam){
  robust_regression[[nam]][["Data"]] %>% data.frame(check.names = F) %>%
    mutate_all(unname) %>%
    filter(`Standardized Residuals` > 4) %>%
    arrange(desc(`Standardized Residuals`)) %>%
    mutate(!!paste0(nam, " Weight") := log2(1:n() * `Standardized Residuals`))
})

from_all = robust_regression$All_Sarcomas$Data %>%
  data.frame(check.names = F) %>%
  arrange(desc(`Standardized Residuals`)) %>% 
  slice(1:20) %>% pull(`gene symbol`)

names = lapply(sig_mutated, function(x){
  x$`gene symbol`[1:10]
})
filtd = sort(table(unlist(names)), decreasing = T)
cosmic_sarcoma_genes = intersect(names(filtd), cosmic_genes$Gene.Symbol)
filtd = names(filtd[filtd>1])

# MAF ---------------------------------------------------------------------

maf = readRDS("Analysis_folders/Oncoplot/somatic_AF4_FR11_10/1170_somatic_maf.rds")
maf = read.maf(maf = maf, clinicalData = clin2)

#subset maf to only COSMIC genes and those that are kept by maftools filtering - removes location types as well as unknown and silent
maf_dat = maf@data %>%
  filter(Hugo_Symbol %in% cosmic_genes$Gene.Symbol) %>%
  full_join(clin2)


# summarise ---------------------------------------------------------------

tmp = maf_dat %>%
  group_by(primary_met, Tumor_Sample_Barcode) %>%
  summarise(cosmic_mutations = sum(!is.na(Variant_Classification)),
            log_cosmic_mutations = log2(cosmic_mutations+1))
tmp %>%
  ggplot() + 
  geom_density(aes(x = log_cosmic_mutations, color = primary_met))
tmp %>%
  ggplot() +
  geom_boxplot(aes(x = primary_met, y = log_cosmic_mutations))

#compare between primary met
wilcox.test(cosmic_mutations ~ primary_met, data = tmp)
tmp %>%
  group_by(primary_met) %>%
  summarise(mean = mean(cosmic_mutations),
            median = median(cosmic_mutations),
            stdev = sd(cosmic_mutations))

tp53_muts = maf_dat %>% filter(Hugo_Symbol == "TP53") %>% pull(Tumor_Sample_Barcode) %>% unique()
atrx_muts = maf_dat %>% filter(Hugo_Symbol == "ATRX") %>% pull(Tumor_Sample_Barcode) %>% unique()
kit_muts = maf_dat %>% filter(Hugo_Symbol == "KIT") %>% pull(Tumor_Sample_Barcode) %>% unique()
#tp53
tmp %>% #whole cohort
  mutate(TP53_mutation = ifelse(Tumor_Sample_Barcode %in% tp53_muts, "yes", "no")) %>%
  ungroup() %>%
  summarise(sum(TP53_mutation == 'yes')/n() * 100)
tmp %>% #primary and met
  mutate(TP53_mutation = ifelse(Tumor_Sample_Barcode %in% tp53_muts, "yes", "no")) %>%
  group_by(primary_met) %>%
  summarise(sum(TP53_mutation == 'yes')/n() * 100)
tmp %>%
  mutate(TP53_mutation = ifelse(Tumor_Sample_Barcode %in% tp53_muts, "yes", "no")) %>%
  group_by(primary_met, TP53_mutation) %>%
  summarise(samps = n()) %>%
  spread(TP53_mutation, samps) %>%
  column_to_rownames("primary_met") %>%
  chisq.test()
#atrx
tmp %>%
  mutate(ATRX_mutation = ifelse(Tumor_Sample_Barcode %in% atrx_muts, "yes", "no")) %>%
  group_by(primary_met) %>%
  summarise(sum(ATRX_mutation == 'yes')/n() * 100)
tmp %>%
  mutate(ATRX_mutation = ifelse(Tumor_Sample_Barcode %in% atrx_muts, "yes", "no")) %>%
  group_by(primary_met, ATRX_mutation) %>%
  summarise(samps = n()) %>%
  spread(ATRX_mutation, samps) %>%
  column_to_rownames("primary_met") %>%
  chisq.test()
#kit
tmp %>%
  mutate(KIT_mutation = ifelse(Tumor_Sample_Barcode %in% kit_muts, "yes", "no")) %>%
  group_by(primary_met) %>%
  summarise(sum(KIT_mutation == 'yes')/n() * 100)
tmp %>%
  mutate(KIT_mutation = ifelse(Tumor_Sample_Barcode %in% kit_muts, "yes", "no")) %>%
  group_by(primary_met, KIT_mutation) %>%
  summarise(samps = n()) %>%
  spread(KIT_mutation, samps) %>%
  column_to_rownames("primary_met") %>%
  chisq.test()





