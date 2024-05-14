#Alex Soupir
#making cnTable for use with mafTools copy number that isn't GISTIC

library(tidyverse) #piping and data processing
library(pbmcapply) #progress bar multicore apply functions
library(data.table) #fast reading of data frames
library(openxlsx) #reading excel files into R

# all.lesions = system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
# amp.genes = system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
# del.genes = system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
# scores.gis = system.file("extdata", "scores.gistic", package = "maftools")
# 
# previous MAF with all mutations for the 1170 samples
# maf = readRDS("Analysis_folders/Oncoplot/somatic_AF4_FR11_10/1170_somatic_maf.rds")

#clinical file created with Dale
clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20230731_niceNames.xlsx")
clin2 = clinical %>% 
  filter(!is.na(somatic_file), #get only samples that we have mutations for
         tumor_germline == "Tumor", #remove the germline samples
         is.na(sarcoma)) %>% #make sure only to keep those that don't have a flag for not being sarcoma
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>% #extract the sample ID from the file name
  group_by(sarcoma_collapsed) %>% #group by sarcoma histology that Andrew Collapsed
  mutate(new_collapsed = ifelse(n() < 5, "other", changed_diagnosis_clean)) %>% #if there are less than 5 samples, group them to Other
  ungroup()


# 'cnTable'
# genes   tsb     cn
# TTN     TCGA.AB.2814    Amp
# NRAS    TCGA.AB.2872    Del
# TP53    TCGA.AB.2957    Amp
# NPM1    TCGA.AB.2904    Del
# SMC3    TCGA.AB.2824    Amp
# SMC1A   TCGA.AB.2898    Amp
# TTN     TCGA.AB.2917    Amp
# TP53    TCGA.AB.2967    Del

#per oliver
# 0 = HOMDEL
# 1 = HETLOSS
# 2 = normal
# 3 & 4 = GAIN
# 5 or more = AMP

#This is just a text file that was created that contains the cntools files in the somatic copy number folder
#a few checks were done to make sure that the 1170 samples we settled on for mutations are also present for copy number
#that way we are working with the same samples
#the reason that this was done was because getting the files with list.files in R takes a REALLY long time here. it was made in bash with ls, grep, and '> file.txt' for the cntools files
cnTools_files = read.table("data/WES/somatic_CNV/cntools_files.txt") %>% #imports as a data frame
  unlist() #create character vector

#need 3 columns for the cnTable of mafTools
# genes = gene symbol that matches Hugo Symbol
# tsb = sample ID that matches that in the maf
# cn = category of copy number change listed per oliver above
cntools_tables = pbmclapply(clin2$wes, function(samp){ #for the WES sample IDs in the clinical file
  cnt_id = grep(samp, cnTools_files) #grab the copy number for that sample
  tmp_file = fread(paste0("data/WES/somatic_CNV/", cnTools_files[cnt_id]), data.table = F) %>% #fread in sample (data.table) 
    filter(!grepl("rRNA", genename)) %>% #remove those genes that are normal with 2 copies and rRNA genes
    mutate(ma = abs(2 - .[[6]]), length = end - start) %>% group_by(genename) %>% select(genename:length) %>% distinct() %>% #select largest variation from 2 cn
    #this was done because there are some genes (HUGO) that have multiple entries, less so after the removal of rRNA genes but still just in case
    filter(length == max(length)) %>% # selecting the longest entry if there are multiple rows for a gene that only differ because of length
    rename("genes" = genename, "cn" = 2) %>% #make names compatible with 
    mutate(tsb = gsub(".cntools.tsv", "", cnTools_files[cnt_id]) %>% #remove file extension for tumor sample barcode column
             str_split("_") %>% unlist() %>% paste0(., "_st") %>% #split the ID and add in _st
             paste0(., collapse = "_t_") %>% paste0(., "_g"), #add the _t_ and _g
           cn = case_when(cn == 0 ~ "HOMDEL",#categorize the copy number value into levels mentioned by Oliver
                          cn == 1 ~ "HETLOSS",
                          cn == 2 ~ "NORMAL",
                          cn %in% c(3, 4) ~ "GAIN",
                          T ~ "AMP")) %>%
    select(genes, tsb, cn, length) %>% #get only columns of interest
    distinct() #double check they are unique
  return(tmp_file) #return sample cleaned CN data
}, mc.cores = 32) #parallel for 32 cores
#bind all tables together
dat = do.call(bind_rows, cntools_tables)
#save file
saveRDS(dat, "Analysis_folders/CopyNumberVariation/cnTable_maftools_wNormal.rds")
