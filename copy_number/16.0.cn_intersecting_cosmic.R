#Alex Soupir
#create heatmap for the copy number data
rm(list=ls())
library(tidyverse) #data processing
library(pbmcapply) #progress bar multicore apply
library(data.table) #data processing (faster than tidyverse but syntax is confusion)
library(openxlsx) #reading excel files of clinical data
library(ComplexHeatmap) #making heatmap

#import data
#the *_niceNames.xlsx file is the same is the one previouly used, just with the sacoma collapsed changed to Title format (first letter capitalized) and underscores removed in column called "nice_name_collapsed"
clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20240716_niceNames.xlsx")
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
cosmic_genes = read.csv("Analysis_folders/Significantly_Mutated_Genes/Census_all_COSMIC.csv") #retrieve the COSMIC Tier 1 gene list

#created cnTable with the creating_cnTable.R script on the cluster
cnTable = readRDS("Analysis_folders/CopyNumberVariation/cnTable_maftools.rds") 

#filter cnTable to have only COSMIC genes
#Much faster using data table rather than dplyr..
cnTable = cnTable %>% as.data.table(check.names = F)
#join clinical and the copy number table together
cnTable = merge(cnTable[genes %in% cosmic_genes$Gene.Symbol][,wes := gsub("\\_.*","", tsb)][,wes := gsub("^T", "", wes)], #mutate the tumor sample barcode ID from the copy number table to match the wes column in the clinical linkage to allow merging
                as.data.table(clin2 %>% #cleaned clinical data file
                                select(wes, nice_name_collapsed) %>% #select the wes column and nice names of histology subtypes collapsed by Andrew
                                group_by(nice_name_collapsed) %>% mutate(Samples = n())), by = "wes") #grouping by the collapsed subtype names, create a new column with number of samples
#different tables created that allow for checking frequency of copy number change by length of gene, as well as the heatmap
cn_cnhist = unique(cnTable[, freq := length(tsb) / Samples, by = list(genes, cn, nice_name_collapsed)][,list(genes, cn, nice_name_collapsed, length, freq)]) #calculate the frequency of copy number change for each gene/histology by each category of copy number change (up to 4 each gene)
cn_hist = unique(cnTable[, hist_freq := length(tsb) / Samples, by = list(genes, nice_name_collapsed)][,list(genes, nice_name_collapsed, hist_freq)]) #calculate the frequency of copy number change for each gene/histology
cnTable2 = merge(cn_cnhist, cn_hist) #merge above together

#breaking apart by histology
by_histology = lapply(unique(cnTable2$nice_name_collapsed), function(histology){ #for each histology subtype
  cnTable2[nice_name_collapsed == histology][order(-rank(freq))] #select said histology and sort the frequency in descending order
})
names(by_histology) = unique(cnTable2$nice_name_collapsed)#set names for this list

#plot each of the histology subtypes not in an x-y plot of sqrt length against frequency with a linear fit curve
plots = lapply(by_histology, function(dat){
  dat %>%
    ggplot() +
    geom_point(aes(x = sqrt(length), y = freq, color = cn), alpha = 0.25) +
    geom_smooth(aes(x = sqrt(length), y = freq), method = lm) +
    theme_bw() +
    facet_wrap(~cn, nrow = 2) + #create 4 panels, one for each class of copy number change
    labs(y = "Frequency", title = unique(dat$nice_name_collapsed)) +
    ylim(c(0, 1))
})

#saving plots
pdf("Analysis_folders/CopyNumberVariation/Figures/GeneCN_vs_GeneLength_COSMIC.pdf", width = 10, height = 10)
tmp = lapply(plots, print)
dev.off()


# Selecting Histology Specific Genes --------------------------------------
hist_genes = lapply(by_histology, function(dat){ #for each of the histolgy subtypes
  dat %>%
    filter(hist_freq > 0.75) #select only those genes that have copy number change (regardless of copy number category) in more 75% of samples
}) %>% do.call(bind_rows, .) %>% #bind all histology subtypes together
  pull(genes) %>% unique() #extract HUGO IDs and get the unique names


cosmic_table = cnTable[genes %in% hist_genes] #filter the big copy number table to include only those genes that were mutated in more 75% of samples
cosmic_mat_prep = cosmic_table %>%
  select(-tsb, -length, -Samples, -freq, -hist_freq) %>% #pull of the sample name, length, number of samples, per copy number category frequency and overall copy number frequency
  #merge it with all combinations of genes and samples to get those sample/gene combinations that may not exist for complete info
  full_join(expand.grid(wes = clin2$wes, 
                        genes = hist_genes) %>%
              full_join(clin2 %>% select(wes, nice_name_collapsed))) %>% #join with the histology subtype names for plotting
  spread("genes", "cn") %>% #convert data from long to wide with rows as samples, columns as genes mutated in at least 1 histology more than 75% samples, and values in each slot the copy number change for that sample/gene combo
  group_by(nice_name_collapsed) %>% #group by first column of histology subtypes
  mutate(num_samps = n(), #identify number of samples in each histology
         .before = nice_name_collapsed) %>%
  ungroup() %>% #removing grouping
  mutate(samp_muts = rowSums(!is.na(.[hist_genes])), #identify the number of overall samples that have some kind of copy number change from normal
         .before = nice_name_collapsed) %>%
  arrange(desc(num_samps), nice_name_collapsed, desc(samp_muts)) #sort the data first by most to least samples in histology, then alphabetical histology, then mutated sample rate
top_annot = cosmic_mat_prep[,1:4] #extract out annotation data in the first few columns (have to be removed for heat map matrix)
cosmic_mat = cosmic_mat_prep %>% 
  ungroup() %>% #remove grouping
  select(-nice_name_collapsed, -num_samps, -samp_muts) %>% #remove above mentioned columns except sample ID
  column_to_rownames("wes") %>% #make sample ID the row names
  as.matrix() %>% t() #convert to matrix and make genes rows and samples columns
#cosmic_mat[is.na(cosmic_mat)] = "None" #if missing copy number change (i.e. normal) replace 'none' category - didn't keep 
cosmic_mat2 = cosmic_mat
#this function converts the copy number categories to numbers and then is able to calculate distance
#with alphabetical categories there is no way to calculate difference in the hierarchical clustering so is necessary in this particular case 
#if wanting to sort the matrix and turn of clustering, not needed
dist_letters = function(x, y) {
  x = ifelse(x == "AMP", "A", ifelse(x == "GAIN", "B", ifelse(x == "None", "C", ifelse(x == "HETLOSS", "D", "E")))) #make sure to keep order of copy number change and replace with letters in said order
  y = ifelse(y == "AMP", "A", ifelse(y == "GAIN", "B", ifelse(y == "None", "C", ifelse(y == "HETLOSS", "D", "E")))) #same
  x = strtoi(charToRaw(paste(x, collapse = "")), base = 16) #convert letter ordering to number
  y = strtoi(charToRaw(paste(y, collapse = "")), base = 16) #same
  sqrt(sum((x - y)^2))#calculate difference in letter (number) distance
}
#get the genes in cosmic that remain in the filtered gene list of the sarcoma samples and convert the gene role if containing more than one role to just be "Multiple"
cosmic_cats = cosmic_genes %>%
  select(Gene.Symbol, Role.in.Cancer) %>%
  filter(Gene.Symbol %in% row.names(cosmic_mat)) %>%
  mutate(Role_simplified = ifelse(grepl(",", Role.in.Cancer), "Multiple", Role.in.Cancer))

#create an annotation for the columns that are the histology subtypes and turn off the legend - this is because when we group it'll throw the label on there anyway so no need for separate legend
col_anno = HeatmapAnnotation(df = data.frame(`Histology (Collapsed)` = top_annot$nice_name_collapsed), show_legend = FALSE)
#create an annotation for the roles the COSMIC gene plays in cancer and turn legend off - this is why there seems to be 4 groups, but multiple group has many colors within
row_anno = rowAnnotation(df = data.frame(Role.In.Cancer = cosmic_cats$Role.in.Cancer), show_legend = FALSE)
#construct the heatmap!
ht = Heatmap(cosmic_mat2, #heat map matrix with samples and columns and genes as rows filled with the copy number levels
             name = "Copy Number", #set legend title
             cluster_columns = T, #turn on the clustering for the columns, pretty sure its on by default
             clustering_distance_rows = dist_letters, #provide the clustering function for our ordinal category levels
             clustering_distance_columns = dist_letters, #same
             bottom_annotation = col_anno, #pass in the annotation for the the histology subtypes
             left_annotation = row_anno, #pass in the annotation for the COSMIC gene roles in cancer
             show_column_names = F, #Don't remember function but don't think in this case it does anything since column_split adds the titles back in
             col = c("AMP" = "red", "GAIN" = "pink3", "None" = "lightgray", "HETLOSS" = "lightblue4", "HOMDEL" = "blue"), #custom color scale for our copy number levels
             column_split = top_annot$nice_name_collapsed, #group our columns by histology subtype and cluster within
             row_split = cosmic_cats$Role_simplified, #group our rows by the gene role in cancer and cluster within
             column_title_rot = 45, #turn the names 45deg so that we can read them 
             row_title_rot = 45) #rotate names 45 degrees for style 
#save the plot
pdf("Analysis_folders/CopyNumberVariation/Figures/COSMIC_gt75-in-hist_heatmap.pdf", height = 20, width = 30)
ht
dev.off()