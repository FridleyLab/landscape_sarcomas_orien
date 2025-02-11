rm(list=ls())

#libraries
library(dplyr)
library(tidyr)

#list the files that are filtered and Annovar annotated
kits = list.files("Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/", pattern = "gz$")
kits2 = data.frame(VCF_files = kits) %>%
  separate(col = VCF_files, 
           into = c("Somatic", "Filter", "Kit", "Genome_Version", "Annotated", "Extension", "Compressed"), 
           remove = FALSE, 
           sep = "\\.") %>%
  mutate(Tumor_ID = gsub("^T", "", Somatic) %>%
           gsub("\\_t.*", "\\_t", .) %>%
           gsub("FT\\-", "FT\\.", .))
IDT_samples = kits2 %>%
  filter(Kit == "IDT")
NIM_samples = kits2 %>%
  filter(Kit == "NIM")

write.table(IDT_samples$Tumor_ID, file = "Analysis_folders/Mutation_Signature/data/target_regions/WES_IDT_Sample-IDs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(NIM_samples$Tumor_ID, file = "Analysis_folders/Mutation_Signature/data/target_regions/WES_NIM_Sample-IDs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
