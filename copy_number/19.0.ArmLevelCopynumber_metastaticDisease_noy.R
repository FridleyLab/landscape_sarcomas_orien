rm(list=ls())
library(tidyverse)
library(parallel)
library(openxlsx)
library(data.table)
library(patchwork)

clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20240716_niceNames.xlsx") %>%
  filter(is.na(sarcoma), !is.na(somatic_file), primary_met == "Metastatic")
clinical$segment_file = str_split(clinical$somatic_file, pattern = "_") %>% 
  do.call(rbind, .) %>% data.frame() %>%
  rowwise() %>% mutate(F = paste0(X1, "_", X4, "_segments.txt")) %>% pull(F)
clin2 = clinical %>%
  filter(!is.na(somatic_file), #make sure that we have the somatic mutaiton file for the sample
         tumor_germline == "Tumor", #remove germline from clinical file
         is.na(sarcoma),
         !reviewer_remove) %>% #remove 'non-sarcoma' samples
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>% #remove everything after "." to get sample ID
  group_by(cdc_revision) %>% #group by the sarcoma histology Andrew collapsed
  mutate(new_collapsed = ifelse(n() < 5, "other", cdc_revision), #if there are less than 5 samples, make new collapsed "other"
         nice_name_collapsed = ifelse(n() < 5, "Other", cdc_nice_name)) %>% #with a nice name called "Other" for plotting
  ungroup()
#cosmic_genes = read.csv("Analysis_folders/Significantly_Mutated_Genes/Census_all_COSMIC.csv") #retrieve the COSMIC Tier 1 gene list
#created cnTable with the creating_cnTable.R script on the cluster
#cnTable = readRDS("Analysis_folders/CopyNumberVariation/cnTable_maftools_wNormal.rds") 


# Getting Centromere Locations --------------------------------------------
centromere = fread("data/hg38_centromereLocations_08may2023.tsv.gz", data.table = F)
centromere_simplified = centromere %>%
  group_by(chrom) %>%
  summarise(chromStart = min(chromStart),
            chromEnd = max(chromEnd))
segment_files = paste0("data/WES/somatic_CNV/", clinical$segment_file)
chromosome_order = fread(segment_files[1], data.table = F) %>% pull(chromosome) %>% unique()
#get chromosome arm
getArm = function(dat){
  dat[dat$end.pos < dat$chromStart,"Arm"] = "p"
  dat[dat$start.pos > dat$chromEnd,"Arm"] = "q"
  dat[is.na(dat$Arm), "Arm"] = "centromere intersect"
  return(dat)
}
#create big data of arm-level copy number changes
# out = mclapply(1:nrow(clin2), function(seg){
#   dat = fread(paste0("data/WES/somatic_CNV/",clin2$segment_file[seg]), data.table = F) %>%
#     select(chromosome, start.pos, end.pos, CNt) %>%
#     full_join(centromere_simplified, by = c("chromosome" = "chrom")) %>%
#     mutate(SegmentLength = end.pos - start.pos) %>%
#     filter(!is.na(CNt))
#   
#   dat2 = getArm(dat)
#   
#   dat3 = lapply(1:nrow(dat2), function(r){
#     wkn = dat2[r,]
#     if(wkn$Arm != "centromere intersect") return(wkn)
#     
#     if(wkn[,"start.pos"] < wkn[,"chromStart"]){
#       #if the overlap is on the p-arm
#       wkn$end.pos = wkn$chromStart-1
#       wkn$Arm = "p"
#     } else if(wkn[,"end.pos"] > wkn[,"chromEnd"]){
#       wkn$start.pos = wkn$chromEnd+1
#       wkn$Arm = "q"
#     } else {
#       return(NULL)
#     }
#     return(wkn)
#   }) %>% 
#     do.call(bind_rows, .)
#   
#   dat3 %>%
#     group_by(chromosome, Arm) %>%
#     mutate(ArmSegmentLength = sum(SegmentLength)) %>%
#     summarise(`Arm CNt` = sum(CNt * SegmentLength) / ArmSegmentLength,
#               SegStart = min(start.pos),
#               SegEnd = max(end.pos)) %>%
#     distinct() %>%
#     mutate(sample = clin2$wes[seg])
# }, mc.cores = 32) %>%
#   do.call(bind_rows, .)

#save data from cluster
# saveRDS(out, "Analysis_folders/CopyNumberVariation/segmentationFiles/data/CopyNumber_summarised-to-pq-arms_metastaticDisease.rds")
out = readRDS("Analysis_folders/CopyNumberVariation/segmentationFiles/data/CopyNumber_summarised-to-pq-arms_metastaticDisease.rds") %>%
  filter(chromosome != "chrY")

#making data
chromosome_order2 = paste(rep(chromosome_order, each = 2), c('p', 'q'), sep = " ")
plot_dat = out %>%
  mutate(cn = case_when(`Arm CNt` < 1 ~ "HOMDEL", 
                        `Arm CNt` >= 1 & `Arm CNt` < 2 ~ "HETLOSS",
                        `Arm CNt` == 2 ~ "NORMAL",
                        `Arm CNt` > 2 & `Arm CNt` <= 3 ~ "GAIN",
                        `Arm CNt` > 3 ~ "AMP"),
         cn = factor(cn, levels = c("HOMDEL", "HETLOSS", "NORMAL", "GAIN", "AMP"))) %>%
  unite("Chr_Arm", "chromosome", "Arm", sep = " ", remove = F) %>%
  mutate(Chr_Arm = factor(Chr_Arm, levels = chromosome_order2)) %>%
  left_join(clin2 %>% select(wes, nice_name_collapsed), by = c("sample" = "wes"))

# Cluster samples and histologies -----------------------------------------

#instead of using the amp-gain-etc, using the arm copyo number to determine 
between_samples_clustering = lapply(unique(clin2$nice_name_collapsed), function(x){
  tmp_dat = plot_dat[plot_dat$sample %in% clin2$wes[clin2$nice_name_collapsed == x],] %>%
    ungroup() %>%
    select(Chr_Arm, `Arm CNt`, sample) %>%
    spread("Chr_Arm", "Arm CNt") %>%
    column_to_rownames("sample") %>%
    as.matrix() %>% 
    dist(., method = "euclidean") %>%
    hclust(., method = "complete")
  tmp_dat$labels[tmp_dat$order]
})
names(between_samples_clustering) = unique(clin2$nice_name_collapsed)

#checking the between histology clustering orders now
between_histology_clustering = readRDS("Analysis_folders/CopyNumberVariation/segmentationFiles/data/histology_order.rds")

sample_order = c()
for(i in between_histology_clustering){
  sample_order = append(sample_order, between_samples_clustering[[i]])
}


sarcoma_segment_plot = plot_dat %>%
  mutate(chromosome = factor(chromosome, levels = chromosome_order),
         sample = factor(sample, levels = sample_order),
         nice_name_collapsed = factor(nice_name_collapsed, levels = between_histology_clustering)) %>%
  ggplot() +
  geom_segment(aes(x = SegStart, y = 0, xend = SegEnd, yend = 0, color = cn), linewidth = 3) + 
  theme_bw() +
  ylim(-0.25, 0.25) +
  ggh4x::facet_nested(nice_name_collapsed + sample ~ chromosome + Arm, scales = "free_x", space = "free_x") +
  #scale_color_gradientn(colors = color_fun(0:5)) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        #strip.text.y.right = element_blank(),
        strip.text.x.top = element_text(angle = 90),
        panel.spacing.y = unit(0, "lines")) +
  labs(title = "Chromosome Copy Number Change by Arm") +
  scale_color_manual(values = c("blue", "lightblue", "gray", "pink", "red")) +
  labs(x = "Genomic Location", y = "Samples")

pdf("Analysis_folders/CopyNumberVariation/segmentationFiles/results/All_Sarcoma_CN-by-arm_metastaticDisease_noy.pdf", width = 20, height = 30)
sarcoma_segment_plot
dev.off()

# whole genome amplification ----------------------------------------------

cn_data = readRDS("Analysis_folders/CopyNumberVariation/segmentationFiles/data/copy_number_level-by-length-noy.rds")
cn_data2 = cn_data %>%
  filter(Tumor_Sample_Barcode %in% clin2$Tumor_Sample_Barcode) %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(All_Segments_Length = sum(Length),
         cn_proportion = Length / All_Segments_Length) %>%
  filter(cn == "AMP") %>%
  mutate(WGA = ifelse(cn_proportion > 0.50, "Yes", "No")) %>%
  left_join(clin2)
table(cn_data2$WGA)
#plot it
pl = cn_data2 %>%
  mutate(sample = factor(wes, levels = sample_order),
         nice_name_collapsed = factor(nice_name_collapsed, levels = between_histology_clustering)) %>%
  ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 10, yend = 0, color = WGA), linewidth = 3) + 
  theme_bw() +
  ylim(-0.25, 0.25) +
  ggh4x::facet_nested(sample ~ ., scales = "free_x", space = "free_x") + #key for keeping the right order, for some reason doesn't like when the nice name is added but the order is the same as above without nice names
  #scale_color_gradientn(colors = color_fun(0:5)) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        #strip.text.y.right = element_blank(),
        strip.text.x.top = element_text(angle = 90),
        panel.spacing.y = unit(0, "lines")) +
  labs(title = "Greater than 75% of Segments with Amplification") +
  scale_color_manual(values = c("gray", "green")) +
  labs(x = "Arbitrary", y = "Samples")

pdf("Analysis_folders/CopyNumberVariation/segmentationFiles/results/All_Sarcoma_WGA_metastaticDisease_noy.pdf", width = 20, height = 30)
pl
dev.off()