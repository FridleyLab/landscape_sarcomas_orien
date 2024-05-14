rm(list=ls())
#libraries
library(tidyverse)
#data
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
#cosmic_genes = read.csv("Analysis_folders/Significantly_Mutated_Genes/Census_all_COSMIC.csv") #retrieve the COSMIC Tier 1 gene list
#created cnTable with the creating_cnTable.R script on the cluster
#cnTable = readRDS("Analysis_folders/CopyNumberVariation/cnTable_maftools_wNormal.rds") 

#create big data of arm-level copy number changes
out = mclapply(1:nrow(clin2), function(seg){
  fread(paste0("data/WES/somatic_CNV/",clin2$segment_file[seg]), data.table = F) %>%
    filter(chromosome != "chrY") %>%
    select(chromosome, start.pos, end.pos, CNt) %>%
    full_join(centromere_simplified, by = c("chromosome" = "chrom")) %>%
    mutate(SegmentLength = end.pos - start.pos) %>%
    filter(!is.na(CNt)) %>%
    mutate(cn = case_when(`CNt` < 1 ~ "HOMDEL", 
                          `CNt` >= 1 & `CNt` < 2 ~ "HETLOSS",
                          `CNt` == 2 ~ "NORMAL",
                          `CNt` > 2 & `CNt` <= 3 ~ "GAIN",
                          `CNt` > 3 ~ "AMP"),
           cn = factor(cn, levels = c("HOMDEL", "HETLOSS", "NORMAL", "GAIN", "AMP"))) %>%
    group_by(cn) %>%
    summarise(Length = sum(SegmentLength)) %>%
    mutate(Tumor_Sample_Barcode = clin2$Tumor_Sample_Barcode[seg], .before = 1)
  
}, mc.cores = 12) %>%
  do.call(bind_rows, .)
#save output
#saveRDS(out, "Analysis_folders/CopyNumberVariation/segmentationFiles/data/copy_number_level-by-length-noy.rds")

#explore
out %>%
  full_join(clin2) %>%
  ggplot() +
  geom_violin(aes(x = cn, y = log2(Length), color = primary_met))
