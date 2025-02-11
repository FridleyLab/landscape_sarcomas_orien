rm(list=ls())
#calculate signature?

#libraries
library(mutSignatures)
library(data.table)
library(tidyverse)
library(openxlsx)

#data
IDT_matrix = fread("Analysis_folders/Mutation_Signature/results/merge.IDT-samples_TCGA_annotated")
NIM_matrix = fread("Analysis_folders/Mutation_Signature/results/merge.NIM-samples_TCGA_annotated")
#cosmic signature file from https://cancer.sanger.ac.uk/signatures/downloads/
cosmic_sigs = fread("Analysis_folders/Mutation_Signature/data/COSMIC_v3.4_SBS_GRCh38.txt") %>%
  column_to_rownames("Type")

clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20240716_niceNames.xlsx")
clin2 = clinical %>%
  filter(!is.na(somatic_file),
         tumor_germline == "Tumor",
         is.na(sarcoma), 
         !reviewer_remove) %>%
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>%
  group_by(cdc_revision) %>%
  mutate(changed_diagnosis_collapsed = ifelse(n() < 5, "other", cdc_revision),
         nice_name_reassigned_collapsed = ifelse(n() < 5, "Other", cdc_nice_name),
         Tumor_ID = gsub("^T", "", Tumor_Sample_Barcode) %>%
           gsub("\\_t.*", "\\_t", .) %>%
           gsub("FT\\-", "FT\\.", .)) %>%
  ungroup()

#bind IDT and NIM together
#they have the first column as sample names which needs to be moved to rownames
count_matrix = bind_rows(
  IDT_matrix %>% 
    column_to_rownames("V1"),
  NIM_matrix %>%
    column_to_rownames("V1")
) %>%
  t() %>%
  data.frame(check.names = FALSE)

count_matrix = as.mutation.counts(count_matrix)

#get the signature matrix of COSMIC
# cosmx <- mutSigData$blcaSIGS %>% dplyr::select(starts_with("COSMIC"))
# cosmx <- as.mutation.signatures(cosmx)
cosmx = as.mutation.signatures(cosmic_sigs)

#calculate the signatures of samples
sarc_res = resolveMutSignatures(mutCountData = count_matrix,
                                signFreqData = cosmx)
sarc_res2 = sarc_res$results$count.result
msigPlot(sarc_res2)

output = as.data.frame(sarc_res2) %>%
  select(!!clin2$Tumor_ID)
#save output
saveRDS(output,"Analysis_folders/Mutation_Signature/results/1162_NMF-docomposition_86-cosmic-sigs.rds")

plot_dat = output %>%
  rownames_to_column("signature") %>%
  filter(signature %in% paste0("SBS", c(1, 2, 5, 6, 13))) %>%
  pivot_longer(cols = contains("_st_t"))
samp_order = plot_dat %>%
  group_by(name) %>%
  summarise(val = sum(value)) %>%
  arrange(desc(val)) %>%
  pull(name)
plot_dat %>%
  mutate(name = factor(name, levels = samp_order)) %>%
  ggplot() + 
  geom_bar(aes(x = name, y = value, fill = signature),
            position = "stack",
            stat = "identity") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#select same sigantures as TCGA and sum 2 and 13
TCGA_signatures = output %>%
  rownames_to_column("signature") %>%
  filter(signature %in% paste0("SBS", c(1, 2, 5, 13))) %>%
  pivot_longer(cols = contains("_st_t"))
TCGA_signatures = TCGA_signatures %>%
  filter(signature %in% c("SBS2", "SBS13")) %>%
  group_by(name) %>%
  summarise(value = sum(value, na.rm= TRUE)) %>%
  mutate(signature = "APOBEC (SBS2, 13)") %>%
  bind_rows(TCGA_signatures %>%
              filter(signature %in% c("SBS1", "SBS5")))
samp_order = TCGA_signatures %>%
  group_by(name) %>%
  summarise(val = sum(value)) %>%
  arrange(desc(val)) %>%
  pull(name)
tcga_waterfall = TCGA_signatures %>%
  mutate(name = factor(name, levels = samp_order)) %>%
  ggplot() + 
  geom_bar(aes(x = name, y = value, fill = signature),
           position = "stack",
           stat = "identity") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Projected Activity onto COSMIC SBS Signatures 1, 5, and 2+13",
       x = "Samples", y = "Signature Activity")
TCGA_sigs_hist = TCGA_signatures %>%
  full_join(clin2 %>%
              select(Tumor_ID, nice_name_reassigned_collapsed),
            by = c("name" = "Tumor_ID"))
tcga_hist_waterfalls = lapply(sort(unique(TCGA_sigs_hist$nice_name_reassigned_collapsed)), function(histo){
  pl_dat = TCGA_sigs_hist %>%
    filter(nice_name_reassigned_collapsed == histo)
  samp_order = pl_dat %>%
    group_by(name) %>%
    summarise(score = sum(value)) %>%
    arrange(desc(score)) %>% pull(name)
  pl_dat %>%
    mutate(name = factor(name, levels = samp_order)) %>%
    ggplot() + 
    geom_bar(aes(x = name, y = value, fill = signature),
             position = "stack",
             stat = "identity") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = paste0("Projected Activity onto COSMIC SBS Signatures 1, 5, and 2+13\n", histo),
         x = "Samples", y = "Signature Activity")
})
#save full study waterfall
pdf("Analysis_folders/Mutation_Signature/results/1162_SBS1-5-2+13_activity.pdf", height = 8, width = 15)
tcga_waterfall
dev.off()
#per histology
pdf("Analysis_folders/Mutation_Signature/results/histology_level_SBS1-5-2+13_activity.pdf", height = 5, width = 8)
tcga_hist_waterfalls
dev.off()

#histology summary kind
#calculate the signatures of samples
count_matrix = as.data.frame(count_matrix)
count_matrix_summary = sapply(unique(clin2$nice_name_reassigned_collapsed), function(histo){
  samps = clin2$Tumor_ID[clin2$nice_name_reassigned_collapsed == histo]
  count_matrix[,samps] %>% rowSums()
}) %>%
  data.frame(check.names = FALSE)
count_matrix_summary = as.mutation.counts(count_matrix_summary)
hist_res = resolveMutSignatures(mutCountData = count_matrix_summary,
                                signFreqData = cosmx)
hist_res2 = hist_res$results$count.result
msigPlot(hist_res2)

output2 = as.data.frame(hist_res2)

pl_dat_hist = output2 %>%
  rownames_to_column("signature") %>%
  filter(signature %in% paste0("SBS", c(1, 2, 5, 13))) %>%
  pivot_longer(cols = -signature)
pl_dat_hist = pl_dat_hist %>%
  filter(signature %in% c("SBS2", "SBS13")) %>%
  group_by(name) %>%
  summarise(value = sum(value, na.rm= TRUE)) %>%
  mutate(signature = "APOBEC (SBS2, 13)") %>%
  bind_rows(pl_dat_hist %>%
              filter(signature %in% c("SBS1", "SBS5")))
samp_order = pl_dat_hist %>%
  group_by(name) %>%
  summarise(val = sum(value)) %>%
  arrange(desc(val)) %>%
  pull(name)
tcga_hist_waterfall = pl_dat_hist %>%
  mutate(name = factor(name, levels = samp_order)) %>%
  ggplot() + 
  geom_bar(aes(x = name, y = value, fill = signature),
           position = "stack",
           stat = "identity") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) +
  labs(title = "Projected Activity onto COSMIC SBS Signatures 1, 5, and 2+13",
       x = "Histology Subtype", y = "Signature Activity")
pdf("Analysis_folders/Mutation_Signature/results/NMF_signature_by_histology-subtype.pdf", height = 5, width = 8)
tcga_hist_waterfall
dev.off()

#jamie's signature estimation packages
# remotes::install_url("https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz")
library(SignatureEstimation)

#make count matrix again
count_matrix = bind_rows(
  IDT_matrix %>% 
    column_to_rownames("V1"),
  NIM_matrix %>%
    column_to_rownames("V1")
) %>%
  t() %>%
  data.frame(check.names = FALSE)

QP = findSigExposures(count_matrix, signaturesCOSMIC, 
                      decomposition.method = decomposeQP)

#make sure that the rows are in the same order
rownames(count_matrix) ==row.names(cosmic_sigs)

#run decomposition
QP = apply(count_matrix, 2, function(x){
  decomposeQP(x, as.matrix(cosmic_sigs))
}) %>%
  data.frame(check.names = FALSE)
#sa is much more slow, but i think more exact
#setting up parallel
samples = colnames(count_matrix) %>% setNames(.,.)
SA = parallel::mclapply(samples, function(s){
  decomposeSA(count_matrix[,s], as.matrix(cosmic_sigs))
}, mc.cores = 23)
SA = SA %>%
  do.call(bind_cols, .)
rownames(QP) = colnames(cosmic_sigs)
rownames(SA) = colnames(cosmic_sigs)

saveRDS(QP, "Analysis_folders/Mutation_Signature/results/1170_QP-docomposition_86-cosmic-sigs.rds")
saveRDS(SA, "Analysis_folders/Mutation_Signature/results/1170_SA-docomposition_86-cosmic-sigs.rds")
