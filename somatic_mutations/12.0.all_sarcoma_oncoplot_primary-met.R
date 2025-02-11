rm(list=ls())
library(maftools)
library(tidyverse)
library(openxlsx)
library(pals)
library(ComplexHeatmap)

`%notin%` = negate(`%in%`)

# External Files ----------------------------------------------------------
clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20240716_niceNames.xlsx")
clin2 = clinical %>%
  filter(!is.na(somatic_file),
         tumor_germline == "Tumor",
         is.na(sarcoma), 
         !reviewer_remove) %>%
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>%
  group_by(cdc_revision) %>%
  mutate(changed_diagnosis_collapsed = ifelse(n() < 5, "other", cdc_revision),
         nice_name_reassigned_collapsed = ifelse(n() < 5, "Other", cdc_nice_name)) %>%
  ungroup()
filtered_annotated = list.files("Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/",
                                pattern = "vcf.gz$", full.names = T)
ordered_filtered_annotated = lapply(clin2$wes, function(v) grep(v, filtered_annotated, value = T)) %>% unlist()
# maf = annovarToMaf(ordered_filtered_annotated)

cosmic_genes = read.csv("Analysis_folders/Significantly_Mutated_Genes/Census_all_COSMIC.csv")
robust_regression = readRDS("Analysis_folders/Significantly_Mutated_Genes/Robust_Regression_object_reassigned.rds")
#standardized residuals are standard deviations so can use to filter genes
tumor_mutation_burden = read.csv("Analysis_folders/TumorMutationBurden/Somatic_MutationsPerMB_10rFiltered.csv",
                                 row.names = NULL)

# Finding Genes to Plot ---------------------------------------------------

sig_mutated = lapply(names(robust_regression[1:27]), function(nam){
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

#TSB to remove
remove_tsb = unique(maf@data$Tumor_Sample_Barcode)[unique(maf@data$Tumor_Sample_Barcode) %notin% clin2$Tumor_Sample_Barcode]
maf = filterMaf(maf, tsb = as.character(remove_tsb))

#make new columns in MAF for oncoplot
maf@clinical.data = maf@clinical.data %>%
  data.frame() %>% 
  group_by(nice_name_reassigned_collapsed) %>%
  mutate(num_samples = n()) %>%
  data.table::data.table()
maf@clinical.data$nice_name_oncoplot = paste0(maf@clinical.data$nice_name_reassigned_collapsed, " (", maf@clinical.data$num_samples, ")")
maf@clinical.data$nice_name_oncoplot = factor(maf@clinical.data$nice_name_oncoplot,
                                              levels = unique(maf@clinical.data$nice_name_oncoplot[order(maf@clinical.data$num_samples, decreasing = TRUE)]))

maf_primary = subsetMaf(maf = maf,
                        tsb = clin2$Tumor_Sample_Barcode[clin2$primary_met == "Primary"])
maf_metastatic = subsetMaf(maf = maf,
                           tsb = clin2$Tumor_Sample_Barcode[clin2$primary_met == "Metastatic"])

#pre colors and things
#make colors
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#reproducible colors
set.seed(333)
#select colors
hist_sub_colors = sample(color, length(unique(maf@clinical.data$nice_name_oncoplot)))
pri_met_colors = sample(color, length(unique(maf@clinical.data$primary_met)))
#give colors names
names(hist_sub_colors) = unique(maf@clinical.data$nice_name_oncoplot)
names(pri_met_colors) = unique(maf@clinical.data$primary_met)
#create index of colors for variables/levels
onco_cols = list(nice_name_oncoplot = hist_sub_colors,
                 primary_met = pri_met_colors)

annotation_order = sort(table(maf@clinical.data$nice_name_oncoplot), decreasing = T) %>% names()
gene_order = maf@data %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  distinct() %>%
  full_join(clin2 %>% select(Tumor_Sample_Barcode, nice_name_reassigned_collapsed)) %>% 
  group_by(Hugo_Symbol) %>%
  summarise(samples = n()) %>%
  arrange(desc(samples)) %>%
  filter(Hugo_Symbol %in% cosmic_sarcoma_genes) %>% pull(Hugo_Symbol)
  

pdf("Analysis_folders/Oncoplot/primary_met/All_sarcomas_top_mutated_cosmix_sarcoma_intersect_primary.pdf")
oncoplot(maf = maf_primary, 
         clinicalFeatures = "nice_name_oncoplot", 
         sortByAnnotation = TRUE, removeNonMutated = FALSE, 
         genes = gene_order,
         annotationColor = onco_cols,
         annotationOrder = annotation_order,
         GeneOrderSort = FALSE, keepGeneOrder = TRUE,
         topBarData = tumor_mutation_burden %>% 
           mutate(`log2(TMB + 1)` = log2(total_perMB + 1)) %>% 
           select(Tumor_Sample_Barcode, `log2(TMB + 1)`))
dev.off()
pdf("Analysis_folders/Oncoplot/primary_met/All_sarcomas_top_mutated_cosmix_sarcoma_intersect_metastatic.pdf")
oncoplot(maf = maf_metastatic, 
         clinicalFeatures = "nice_name_oncoplot", 
         sortByAnnotation = TRUE, removeNonMutated = FALSE, 
         genes = gene_order,
         annotationColor = onco_cols,
         annotationOrder = annotation_order,
         GeneOrderSort = FALSE, keepGeneOrder = TRUE,
         topBarData = tumor_mutation_burden %>% 
           mutate(`log2(TMB + 1)` = log2(total_perMB + 1)) %>% 
           select(Tumor_Sample_Barcode, `log2(TMB + 1)`))
dev.off()

# sort(table(maf@clinical.data$nice_name_oncoplot), decreasing = T)
# 
# #SARC genes
# 
# TCGA_sarc_genes = c("TP53", "ATRX", "RB1", "KRAS", "MAP4K1", "NF1", "NF2", "ATR", "ATM", "BRCA2", "DICER1", "FANCC", "PRKDC", "ARID1A", "ARID1B", "KMT2B", "KMT2C", "KMT2D", "PBRM1", "ARID2", "CREBBP", "KMT2A", "PIK3CA", "PTEN", "DNMT3A", "FGFR1", "KIT", "RET", "APC", "IDH1")
# TCGA_sarc_genes[!(TCGA_sarc_genes %in% cosmic_genes$Gene.Symbol)]
# 
# pdf("Analysis_folders/Oncoplot/primary_met/All_sarcomas_top_mutated_TCGA_genes_primary.pdf")
# oncoplot(maf = maf_primary, clinicalFeatures = "nice_name_collapsed", sortByAnnotation = TRUE, removeNonMutated = FALSE, genes = TCGA_sarc_genes)
# dev.off()
# pdf("Analysis_folders/Oncoplot/primary_met/All_sarcomas_top_mutated_TCGA_genes_metastatic.pdf")
# oncoplot(maf = maf_metastatic, clinicalFeatures = "nice_name_collapsed", sortByAnnotation = TRUE, removeNonMutated = FALSE, genes = TCGA_sarc_genes)
# dev.off()
# 
# # Lollipop Plots ----------------------------------------------------------
# #come back to. Don't understand why some mutations are being dropped
# GIST_maf = subsetMaf(maf_primary, tsb = clin2$Tumor_Sample_Barcode[clin2$new_collapsed == "gastrointestinal_stromal_tumor"])
# lollipopPlot(maf = GIST_maf,
#              gene = "KIT", 
#              AACol = "AAChange.refGene", 
#              showMutationRate = T)
# LEIO_maf = subsetMaf(maf_primary, tsb = clin2$Tumor_Sample_Barcode[clin2$new_collapsed == "leiomyosarcoma"])
# lollipopPlot(maf = LEIO_maf,
#              gene = "TP53", 
#              AACol = "AAChange.refGene", 
#              showMutationRate = T)
# 
# 
# # heatmap by histology ----------------------------------------------------
# clin2_primary = clin2 %>%
#   filter(primary_met == "Primary")
# 
# mat_primary = sapply(sort(unique(maf_primary@clinical.data$nice_name_collapsed)), function(histology){
#   sub_maf = subsetMaf(maf_primary, tsb = clin2_primary$Tumor_Sample_Barcode[clin2_primary$nice_name_collapsed == histology])
#   total = length(unique(clin2_primary$Tumor_Sample_Barcode[clin2_primary$nice_name_collapsed == histology]))
#   fracs = sapply(cosmic_sarcoma_genes, function(gene){
#     sub_maf@data %>%
#       group_by(Tumor_Sample_Barcode) %>%
#       summarise(Mut = ifelse(TRUE %in% (gene %in% Hugo_Symbol), 1, 0)) %>%
#       pull(Mut) %>% sum()
#   })
#   fracs / total * 100
# })
# library(circlize)
# col_fun = colorRamp2(c(0, 100), c("white", "darkorange"))
# Heatmap(mat_primary, cluster_rows = F, col = col_fun)
# 
# mut_rate_primary = maf_primary@data %>%
#   group_by(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
#   filter(Hugo_Symbol %in% cosmic_sarcoma_genes) %>%
#   distinct() %>%
#   summarise(Mut_Counts = length(unique(Variant_Classification))) %>%
#   spread(Variant_Classification, Mut_Counts) %>%
#   ungroup() %>%
#   mutate(Multi_Hit = ifelse(rowSums(.[3:10], na.rm=T) > 1, 1, NA)) %>%
#   rowwise() %>%
#   mutate(across(Frame_Shift_Del:Translation_Start_Site, ~ ifelse(is.na(Multi_Hit), .x, NA))) %>%
#   gather(Variant_Classification, Mutated, -Tumor_Sample_Barcode, -Hugo_Symbol) %>%
#   filter(!is.na(Mutated)) %>%
#   group_by(Hugo_Symbol, Variant_Classification) %>% 
#   summarise(Mut = n()) %>% 
#   mutate(Mut = Mut/922 * 100)
# ord = names(sort(mut_rate_primary %>% group_by(Hugo_Symbol) %>%
#                    summarise(Mut = sum(Mut)) %>%
#                    pull(Mut, Hugo_Symbol), decreasing = T))
# mut_rate2_primary = mut_rate_primary %>%
#   mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(ord)))
# 
# primary_barplot = mut_rate2_primary %>% 
#   ggplot() +
#   geom_bar(aes(fill = Variant_Classification, y = Mut, x = Hugo_Symbol), position = "stack", stat = "identity") +
#   coord_flip() +
#   theme_bw() + 
#   labs(y = "Percent of Samples Mutated", x = NULL)
# 
# col_fun = colorRamp2(c(0, 100), c("white", "darkorange"))
# primary_heatmap = Heatmap(mat_primary[match(ord, row.names(mat_primary)),
#             match(names(sort(mat_primary["TP53",], decreasing = T)), colnames(mat_primary))],
#             name = "Primary\nFrequency",
#         cluster_rows = F, cluster_columns = F, col = col_fun)
# #save plots
# pdf("Analysis_folders/Oncoplot/primary_met/Primary_Heatmap_cosmic_intersect.pdf")
# primary_heatmap
# dev.off()
# pdf("Analysis_folders/Oncoplot/primary_met/Primary_Barplot_cosmic_intersect.pdf")
# primary_barplot
# dev.off()
# 
# # heatmap by histology - Metastatic samples ----------------------------------------------------
# clin2_metastatic = clin2 %>%
#   filter(primary_met == "Metastatic")
# 
# mat_metastatic = sapply(sort(unique(maf_metastatic@clinical.data$nice_name_collapsed)), function(histology){
#   sub_maf = subsetMaf(maf_metastatic, tsb = clin2_metastatic$Tumor_Sample_Barcode[clin2_metastatic$nice_name_collapsed == histology])
#   total = length(unique(clin2_metastatic$Tumor_Sample_Barcode[clin2_metastatic$nice_name_collapsed == histology]))
#   fracs = sapply(cosmic_sarcoma_genes, function(gene){
#     sub_maf@data %>%
#       group_by(Tumor_Sample_Barcode) %>%
#       summarise(Mut = ifelse(TRUE %in% (gene %in% Hugo_Symbol), 1, 0)) %>%
#       pull(Mut) %>% sum()
#   })
#   fracs / total * 100
# })
# colnames(mat_metastatic) = clin2_metastatic$nice_name_collapsed[match(colnames(mat_metastatic), clin2_metastatic$nice_name_collapsed)]
# library(circlize)
# col_fun = colorRamp2(c(0, 100), c("white", "darkorange"))
# Heatmap(mat_metastatic, cluster_rows = F, col = col_fun)
# 
# mut_rate_metastatic = maf_metastatic@data %>%
#   group_by(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
#   filter(Hugo_Symbol %in% cosmic_sarcoma_genes) %>%
#   distinct() %>%
#   summarise(Mut_Counts = length(unique(Variant_Classification))) %>%
#   spread(Variant_Classification, Mut_Counts) %>%
#   ungroup() %>%
#   mutate(Multi_Hit = ifelse(rowSums(.[3:10], na.rm=T) > 1, 1, NA)) %>%
#   rowwise() %>%
#   mutate(across(Frame_Shift_Del:Translation_Start_Site, ~ ifelse(is.na(Multi_Hit), .x, NA))) %>%
#   gather(Variant_Classification, Mutated, -Tumor_Sample_Barcode, -Hugo_Symbol) %>%
#   filter(!is.na(Mutated)) %>%
#   group_by(Hugo_Symbol, Variant_Classification) %>% 
#   summarise(Mut = n()) %>% 
#   mutate(Mut = Mut/248 * 100)
# ord = names(sort(mut_rate_metastatic %>% group_by(Hugo_Symbol) %>%
#                    summarise(Mut = sum(Mut)) %>%
#                    pull(Mut, Hugo_Symbol), decreasing = T))
# mut_rate2_metastatic = mut_rate_metastatic %>%
#   mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(ord)))
# 
# metastatic_barplot = mut_rate2_metastatic %>% 
#   ggplot() +
#   geom_bar(aes(fill = Variant_Classification, y = Mut, x = Hugo_Symbol), position = "stack", stat = "identity") +
#   coord_flip() +
#   theme_bw() + 
#   labs(y = "Percent of Samples Mutated", x = NULL)
# 
# col_fun = colorRamp2(c(0, 100), c("white", "darkorange"))
# metastatic_heatmap = Heatmap(mat_metastatic[match(ord, row.names(mat_metastatic)),
#                                       match(names(sort(mat_metastatic["TP53",], decreasing = T)), colnames(mat_metastatic))],
#                              name = "Metastatic\nFrequency",
#                           cluster_rows = F, cluster_columns = F, col = col_fun)
# #save plots
# pdf("Analysis_folders/Oncoplot/primary_met/Metastatic_Heatmap_cosmic_intersect.pdf")
# metastatic_heatmap
# dev.off()
# pdf("Analysis_folders/Oncoplot/primary_met/Metastatic_Barplot_cosmic_intersect.pdf")
# metastatic_barplot
# dev.off()
# 
# 
# 



