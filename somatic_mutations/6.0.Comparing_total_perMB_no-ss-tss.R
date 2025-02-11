rm(list=ls())
gc(full=TRUE)
#libraries
library(tidyverse)
library(data.table)
library(openxlsx)
library(maftools)
library(parallel)
#get the maf objects previously created for each histology subtype - not collapsed
mafs = list.files('Analysis_folders/FilteringVCFfiles/All_Sarcomas_MAFs//', pattern = ".maf")
#genome coverage for the different kits
IDT = 39.671448
NIM = 63.377915
#list the files that are filtered and Annovar annotated
kits = list.files("Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated//")
#updated_clinical
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
#loop over the mafs and calculate the real mutations per MB using appropriate kit capture size
res = lapply(mafs, function(maf){
  laml = read.maf(paste0("Analysis_folders/FilteringVCFfiles/All_Sarcomas_MAFs/", maf))
  
  ##code to make the plots
  ##
  
  cohort = maf %>%
    gsub("pass_tumor_", "", .) %>%
    gsub("_maftools.maf", "", .)
  laml.mutload = tcgaCompare(maf = laml, cohortName = cohort, logscale = TRUE, capture_size = 39.671448)
  dat = laml@data %>%
    filter(!(Variant_Classification %in% c("Splice_Site", 'Translation_Start_Site'))) %>%
    select(Tumor_Sample_Barcode, Variant_Classification) %>%
    rowwise() %>%
    mutate(Capture_Region = ifelse(grep(Tumor_Sample_Barcode, kits, value = T) %>%
                                     grepl("IDT", .), IDT, NIM)) %>%
    ungroup() %>%
    group_by(Tumor_Sample_Barcode) %>%
    mutate(total_perMB = n() / Capture_Region) %>%
    group_by(Tumor_Sample_Barcode, Variant_Classification) %>%
    mutate(Variant_Classification_perMB = n() / Capture_Region)
  return(dat)
}) %>%
  do.call(bind_rows, .)
#merge with reassigned, collapsed, histology subtypes
res2 = res %>%
  right_join(clin2 %>%
              select(changed_diagnosis_clean:nice_name_reassigned_collapsed))
#save the data
write.csv(res2, "Analysis_folders/FilteringVCFfiles/TMB_no-ss-tss/TMB_expanded_no-ss-tss.csv")
