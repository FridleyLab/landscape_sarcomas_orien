#Alex Soupir
#Filtering AVATAR Sarcoma Somatic VCF Files
#Recommendations for filters provided by:
#   Dale Hedges
#   Oliver Hampton
library(tidyverse)
library(vcfR)
library(openxlsx)
library(parallel)

clinical = read.xlsx("data/clinical/ClinicalLinkagewithFiles_20221103.xlsx")
#remove those that andrew said aren't sarcomas
sarcoma = clinical %>%
  filter(is.na(sarcoma))
sum(!is.na(sarcoma$somatic_file))

#pon files
old = list.files("data/WES/annotated_somatic_vcfs/", "PoN.v2.vcf.gz$")
new = list.files("data/WES/WES_annotated_somatic_vcfs/", "PoN.v2.vcf.gz$")
pon_files = c(old, new)

grepl(paste0(sarcoma$somatic_file[!is.na(sarcoma$somatic_file)] %>% 
        gsub(".tnhaplotyper2.vcf.gz", "", .), collapse = "|"), 
      pon_files %>% gsub(".ft.M2GEN.PoN.v2.vcf.gz", "", .)) %>% 
  sum()
#whoohoo all of the 1170 are there which makes sense
#now time to count multiple alt alleles
sarcoma_files = sarcoma$somatic_file[!is.na(sarcoma$somatic_file)]
#looping over all files and filtering
empty_variable = mclapply(sarcoma_files, function(f1){
  #grab sample id
  f2 = gsub("tnha.*", "", f1)
  #grab PoN sample
  f = grep(f2, list.files("data/WES/WES_annotated_somatic_vcfs/", ".vcf.gz$", full.names = T), value =T)
  if(0 == length(f)){
    f = grep(f2, list.files("data/WES/annotated_somatic_vcfs/", ".vcf.gz$", full.names = T), value =T)
  }
  #reading vcf file
  vcf_file = read.vcfR(f)
  #separate out the genotype and common columns
  gt = vcf_file@gt %>% data.frame()
  fix = vcf_file@fix %>% data.frame()
  #bind together to make single dataframe while removing the germline
  df = cbind(fix, gt) %>% 
    filter(FILTER == "PASS") %>%
    select(-grep("st_g", colnames(.)))
  
  #seq along each row in the data frame that is pass for tumor information
  out = lapply(1:nrow(df), function(df_r){
    #extract row 
    tmp = as.list(df[df_r,])
    #split info into attributes
    info_el = str_split(tmp$INFO, ";") %>%
      unlist()
    t = lapply(info_el, function(l){
      unlist(str_split(l, pattern = "="))[-1] %>%
        str_split(., ",") %>% unlist()
    })
    names(t) = lapply(info_el, function(l){
      unlist(str_split(l, pattern = "="))[1]
    }) %>% unlist()
    #alt alleles
    alt_allele = unlist(str_split(tmp$ALT, ","))
    #format anmes
    att_nam = unlist(str_split(tmp$FORMAT, ":"))
    #sample information
    s_info = unlist(str_split(tmp[[length(tmp)]], ":")) %>%
      lapply(., function(x){
        str_split(x, "/|,") %>% unlist()
      })
    names(s_info) = att_nam
    #filters for M2Gen and moffitt
    AF_04 = c(TRUE, s_info$AF >= 0.04)
    F1R2_1 = s_info$F1R2 >= 1
    F2R1_1 = s_info$F2R1 >= 1
    min20 = as.numeric(s_info$F1R2) + as.numeric(s_info$F2R1) >= 10
    #get vector of TRUE/FALSE to use for vector filtering
    keep = AF_04 & F2R1_1 & F1R2_1 & min20
    #if the only thing that passed was the filtering criteria was the reference, return null to remove
    if(sum(keep) < 2){
      return(NULL)
    }
    #if reference doesn't meet filtering criteria, remove mutation row
    if(!keep[1]){
      return(NULL)
    }
    #for the last column that should be the tumor mutation information
    tmp[length(tmp)] = lapply(names(s_info), function(n){
      #if we are looking at just the genotype information
      #collapse to the first however many passed
      #if there were 3 but the first and last passed, keep 0/1 instead of 0/2
      if(n == "GT"){
        #maintains genotype indexing
        return(paste0(s_info[[n]][1:sum(keep)], collapse = "/"))
      }
      if(grepl("SA", n)){
        return(paste0(s_info[[n]], collapse = ","))
      }
      #values that are for the other fields not including genotype
      val = s_info[[n]]
      #if there are elements for reference and alt, remove those that didn't pass filtering
      #assuming all ref passes the filtering criteria and if it didn't it was removed in line 49
      if(length(val) == length(keep)){
        return(paste0(val[keep], collapse=","))
      }
      #if only alternates are elements, remove "keep"s pointer to reference (1)
      if(length(val) < length(keep)){
        return(paste0(val[keep[-1]], collapse=","))
      }
    }) %>%
      unlist() %>% paste0(., collapse = ":")#collapse back to original format
    #if alts are removed to be only a single, return single, else return separated by comma
    tmp$ALT = ifelse(length(alt_allele[keep[-1]]) == 1,
                     alt_allele[keep[-1]],
                     paste0(alt_allele[keep[-1]], collapse = ","))
    if(FALSE %in% keep){
      tmp$FILTER = paste0(tmp$FILTER, ";adjusted")
    }
    #for element of INFO column
    tmp$INFO = lapply(seq(t), function(x_n){
      #get values associated with info element
      x = t[[x_n]]
      #if no inforamtion is available, just renturn name
      #Mostly noticed for STR
      if(is.null(x)){
        return(names(t)[x_n])
      }
      #if only single length for element keep single
      if(length(x) == 1){
        r = x
      }
      #if multiple length for element but less than length of "keep"
      #likely just including information for alternate allele, 
      #remove reference information and then filter based on "keep"
      if(length(x) == length(keep)-1){
        r=x[keep[-1]]
      }
      #if length of element is same as 'keep', has reference information
      #vector filter info element to 'keep'
      if(length(x) == length(keep)){
        r=x[keep]
      }
      #if there is only a single length of then just paste the name and element with equals between
      if(length(r) == 1){
        return(paste0(names(t)[x_n], "=", r))
      }
      #if length is greater than 1, collapse the elements to separated by comma then collapse to name with equals
      if(length(r) > 1){
        return(paste0(names(t)[x_n], "=", paste0(r, collapse = ",")))
      }
      
    }) %>%
      unlist() %>%
      paste0(., collapse = ";")#paste all info elements together with semi-colon
    #return fixed row if not NULL
    return(as.data.frame(tmp))
  }) %>%
    do.call(bind_rows, .) #collapse back to table
  
  vcf_file@fix = out[,1:8] %>% as.matrix()
  vcf_file@gt = out[,9:10] %>% as.matrix()
  f_name = paste0(f2, "AF04_F1R21_F2R11_10reads.vcf.gz")
  pass_adjusted = "##FILTER=<ID=adjusted,Description=\"Entry filtered to AF > 4 percent and at least 1 in F1R2 and F2R1\">"
  vcf_file@meta = append(vcf_file@meta, pass_adjusted, after = 1)
  write.vcf(vcf_file, file = paste0("Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads/", f_name))
  return(NULL)
}, mc.cores = 12)
