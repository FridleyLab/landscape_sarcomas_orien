rm(list=ls())
library(tidyverse)
library(openxlsx)
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
         nice_name_reassigned_collapsed = ifelse(n() < 5, "Other", cdc_nice_name),
         Tumor_ID = gsub("^T", "", Tumor_Sample_Barcode) %>%
           gsub("\\_t.*", "\\_t", .) %>%
           gsub("FT\\-", "FT\\.", .)) %>%
  ungroup()
#other parameters
hist_colors = readRDS("Analysis_folders/1.0.Histology_Reassignment/collapsed_histology_color-vec.rds")
#rename "other" to "Other"
names(hist_colors)[which(names(hist_colors) == "other")] = "Other"
#order of the histology subtypes
annotation_order = sort(table(clin2$nice_name_reassigned_collapsed), decreasing = T) %>% names()
#simulated annealing signature matric
SA = readRDS("Analysis_folders/Mutation_Signature/results/1170_SA-docomposition_86-cosmic-sigs.rds")
SA2 = SA %>% 
  rownames_to_column("rows") %>%
  data.frame(check.names = FALSE) %>% 
  column_to_rownames("rows") %>% 
  select(any_of(clin2$Tumor_ID)) %>% as.matrix()
#proposed signature details
# https://cancer.sanger.ac.uk/signatures/sbs/

#i think removing signatures that don't have any samples with a value will make the heatmap more to the point
sig_sums = rowSums(SA2)
SA3 = SA2[sig_sums!=0,]

# Heatmap -----------------------------------------------------------------

column_ha = HeatmapAnnotation(Histology = clin2$nice_name_reassigned_collapsed,
                              Primary_Met = clin2$primary_met,
                              col = list(Histology = hist_colors), show_legend = FALSE)
ht = Heatmap(SA3, name = "Scores", top_annotation = column_ha,
        show_column_names = FALSE, show_row_names = TRUE,
        column_split = factor(clin2$nice_name_reassigned_collapsed, levels = annotation_order), 
        cluster_columns = TRUE, 
        cluster_column_slices = FALSE, 
        column_title_rot = 45, 
        column_title_gp = gpar(fontsize = 10), 
        row_names_gp = gpar(fontsize = 10))
#names are long so increase right margin
ht2 = draw(ht, padding = unit(c(2,2,2,20), "mm"))

#save pdf
pdf("Analysis_folders/Mutation_Signature/results/signature_heatmap.pdf", height = 10, width = 15)
ht2
dev.off()