library(maftools)
library(tidyverse)
library(parallel)
library(openxlsx)
clinical = read.xlsx("data/clinical/ClinicalLinkagewithFiles_20221103.xlsx")
wes = clinical %>%
  filter(wes != "")

top = clinical %>%
  filter(tumor_germline != "Germline",
         wes != "",
         !is.na(somatic_file)) %>% 
  group_by(sarcoma_collapsed) %>% 
  summarise(n=n()) %>%
  arrange(desc(n))
top

vcfs = list.files("Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/", pattern = ".vcf", full.names = T)
tmp = lapply(top$sarcoma_collapsed[1:10], function(type){
  wes2 = wes %>%
    filter(tumor_germline != "Germline",
           wes != "",
           sarcoma_collapsed == type)
  wes2$VCF = NA
  for(i in 1:nrow(wes2)){
    wes2$VCF[i] = ifelse(length(grep(wes2$wes[i], vcfs, value = T)) == 0,
                         NA,
                         grep(wes2$wes[i], vcfs, value = T))
    }
  wes2 = wes2 %>% filter(!is.na(VCF))
  
  maf = annovarToMaf(wes2$VCF, refBuild = "hg38", MAFobj = T)
  write.mafSummary(maf = maf, basename = paste0('Analysis_folders/FilteringVCFfiles/MAFs10/',type))
})