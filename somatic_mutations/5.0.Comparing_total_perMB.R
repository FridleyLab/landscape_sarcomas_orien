#libraries
library(tidyverse)
library(data.table)
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
clinical = read.xlsx("Analysis_folders/1.0.Histology_Reassignment/ClinicalLinkagewithFiles_20230731_niceNames.xlsx")
clin2 = clinical %>%
  filter(!is.na(somatic_file),
         tumor_germline == "Tumor",
         is.na(sarcoma)) %>%
  mutate(Tumor_Sample_Barcode = gsub("\\..*", "", somatic_file)) %>%
  group_by(changed_diagnosis_clean) %>%
  mutate(changed_diagnosis_collapsed = ifelse(n() < 5, "other", changed_diagnosis_clean),
         nice_name_reassigned_collapsed = ifelse(n() < 5, "other", nice_name_reassigned)) %>%
  ungroup()
#loop over the mafs and calculate the real mutations per MB using appropriate kit capture size
tmp = lapply(mafs, function(maf){
  laml = read.maf(paste0("Analysis_folders/FilteringVCFfiles/All_Sarcomas_MAFs/", maf))
  
  ##code to make the plots
  ##
  
  cohort = maf %>%
    gsub("pass_tumor_", "", .) %>%
    gsub("_maftools.maf", "", .)
  laml.mutload = tcgaCompare(maf = laml, cohortName = cohort, logscale = TRUE, capture_size = 39.671448)
  mut_dat = laml.mutload$mutation_burden_perSample %>%
    filter(cohort == !!cohort) %>%
    rowwise() %>%
    mutate(Capture_Region = ifelse(grep(Tumor_Sample_Barcode, kits, value = T) %>%
                                     grepl("IDT", .), IDT, NIM),
           total_perMB = total / Capture_Region)
  mut_dat %>%
    ungroup() %>%
    bind_rows(laml.mutload$mutation_burden_perSample %>%
                filter(cohort != !!cohort)) %>%
    write.csv(., paste0("Analysis_folders/FilteringVCFfiles/MutationPerMB10/", cohort, ".csv"), row.names = F)
})
#read in data and bind together
all_studies_total_perMB = lapply(list.files("Analysis_folders/FilteringVCFfiles/MutationPerMB10", full.names = T, pattern = "csv$"), 
                                 read.csv, check.names = F) %>%
  do.call(bind_rows, .) %>% distinct() #since all will have the TCGA, remove duplicates

#get the updated histology subtype assignments
sum(all_studies_total_perMB$Tumor_Sample_Barcode %in% clin2$Tumor_Sample_Barcode)
#result is 1168 because there are 2 samples with none
#all tcga samples have "TCGA" in the name which makes it easier
all_studies_total_perMB = clin2 %>%
  filter(Tumor_Sample_Barcode %in% all_studies_total_perMB$Tumor_Sample_Barcode) %>% 
  select(Tumor_Sample_Barcode, nice_name_reassigned_collapsed) %>% 
  rename("cohort" = 2) %>%
  full_join(all_studies_total_perMB %>%
              filter(!grepl("TCGA", Tumor_Sample_Barcode)) %>%
              select(-cohort)) %>%
  full_join(all_studies_total_perMB %>%
              filter(grepl("TCGA", Tumor_Sample_Barcode))) %>%
  relocate(cohort, .after = total)
  
#pull out the first few columns, remove column for kit
tcga.cohort = all_studies_total_perMB %>% 
  select(1:4)

#total mutations per each
raw_tmb = all_studies_total_perMB %>%
  right_join(clin2 %>%
               select(primary_met, Tumor_Sample_Barcode, nice_name_reassigned_collapsed)) %>%
  mutate(total_perMB = ifelse(is.na(total_perMB), 0, total_perMB),
         total = ifelse(is.na(total), 0, total))
mean(raw_tmb$total_perMB)
raw_tmb %>% group_by(primary_met) %>% summarise(`Mean TMB` = mean(total_perMB))
wilcox.test(total_perMB ~ primary_met, data = raw_tmb)

#i think these were part of the MAFtools `tcgaCompare` function but the test wasn't used for plotting i dont think
pt.test = pairwise.t.test(x = tcga.cohort$total_perMB, 
                          g = tcga.cohort$cohort, p.adjust.method = "fdr")
pt.test.pval = as.data.frame(pt.test$p.value)
data.table::setDT(x = pt.test.pval, keep.rownames = TRUE)
colnames(pt.test.pval)[1] = "Cohort"
pt.test.pval = data.table::melt(pt.test.pval, id.vars = "Cohort")
colnames(pt.test.pval) = c("Cohort1", "Cohort2", "Pval")

#picking back up with line 40
#pretty sure everything below is just making the data we have fig into the function `tcgaCompare`
tcga.cohort = as.data.table(tcga.cohort)
tcga.cohort[, `:=`(plot_total, total_perMB)]

decreasing = FALSE
tcga.cohort.med = tcga.cohort[, .(.N, median(plot_total)), 
                              cohort][order(V2, decreasing = decreasing)]
tcga.cohort$cohort = factor(x = tcga.cohort$cohort, levels = tcga.cohort.med$cohort)
colnames(tcga.cohort.med) = c("Cohort", "Cohort_Size", "Median_Mutations")
tcga.cohort$TCGA = ifelse(test = nchar(tcga.cohort$cohort %>% as.character()) > 4 | tcga.cohort$cohort == "RMS", yes = "Input", no = "TCGA")
tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
plot.dat = lapply(seq_len(length(tcga.cohort)), function(i) {
  x = tcga.cohort[[i]]
  x = data.table::data.table(rev(seq(i - 1, i, length.out = nrow(x))), 
                             x[order(plot_total, decreasing = TRUE), plot_total], 
                             x[, TCGA])
  x
})

names(plot.dat) = names(tcga.cohort)
logscale = T
if (logscale) {
  y_lims = range(log10(data.table::rbindlist(l = plot.dat)[V2 != 
                                                             0][, V2]))
} else {
  y_lims = range(data.table::rbindlist(l = plot.dat)[, 
                                                     V2])
}
y_min = floor(min(y_lims))
y_max = ceiling(max(y_lims))
y_lims = c(y_min, y_max)
y_at = pretty(y_lims)

pdf("Analysis_folders/TumorMutationBurden/AF04_F1R21_F2R11_10reads_wTCGApancancer_manuscript-figure.pdf", width = 15, height = 6)
plot(NA, NA, xlim = c(0, length(plot.dat)), ylim = y_lims, 
     axes = FALSE, xlab = NA, ylab = NA, frame.plot = TRUE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = grDevices::adjustcolor(col = "gray", alpha.f = 0.1))
col = c("gray70", "black")
bg_col = c("#EDF8B1", "#2C7FB8")
rect(xleft = seq(0, length(plot.dat) - 1, 1), ybottom = min(y_lims), 
     xright = seq(1, length(plot.dat), 1), ytop = y_max, 
     col = grDevices::adjustcolor(col = bg_col, alpha.f = 0.2), 
     border = NA)
abline(h = pretty(y_lims), lty = 2, col = "gray70")
lapply(seq_len(length(plot.dat)), function(i) {
  x = plot.dat[[i]]
  if (x[1, V3] == "Input") {
    if (logscale) {
      points(x$V1, log10(x$V2), pch = 16, cex = 0.4, 
             col = col[2])
    }
    else {
      points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[2])
    }
  }
  else {
    if (logscale) {
      points(x$V1, log10(x$V2), pch = 16, cex = 0.4, 
             col = col[1])
    }
    else {
      points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[1])
    }
  }
})
samp_sizes = lapply(plot.dat, nrow)
cohortFontSize = 0.8
axis(side = 1, at = seq(0.5, length(plot.dat) - 0.5, 1), 
     labels = names(plot.dat), las = 2, tick = FALSE, line = -0.8, 
     cex.axis = cohortFontSize)
axis(side = 3, at = seq(0.5, length(plot.dat) - 0.5, 1), 
     labels = unlist(samp_sizes), las = 2, tick = FALSE, 
     line = -0.8, font = 3, cex.axis = cohortFontSize)
tcga.cohort.med[, `:=`(Median_Mutations_log10, log10(Median_Mutations))]
if (decreasing) {
  sidePos = 2
  linePos = 2
} else {
  sidePos = 2
  linePos = 2
}
capture_size = "mine"
axisFontSize = 0.8
if (logscale) {
  if (is.null(capture_size)) {
    axis(side = 2, at = y_at, las = 2, line = -0.6, 
         tick = FALSE, labels = 10^(y_at), cex.axis = axisFontSize)
    mtext(text = "TMB", side = sidePos, line = linePos)
  } else {
    axis(side = 2, at = y_at, las = 2, line = -0.6, 
         tick = FALSE, labels = 10^(y_at), cex.axis = axisFontSize)
    mtext(text = "TMB (per MB)", side = sidePos, line = linePos)
  }
} else {
  if (is.null(capture_size)) {
    axis(side = 2, at = y_at, las = 2, line = -0.6, 
         tick = FALSE, cex.axis = axisFontSize)
    mtext(text = "TMB", side = sidePos, line = linePos)
  } else {
    axis(side = 2, at = y_at, las = 2, line = -0.6, 
         tick = FALSE, cex.axis = axisFontSize)
    mtext(text = "TMB (per MB)", side = sidePos, line = linePos)
  }
}
medianCol = "red"
if (logscale) {
  lapply(seq_len(nrow(tcga.cohort.med)), function(i) {
    segments(x0 = i - 1, x1 = i, y0 = tcga.cohort.med[i, 
                                                      Median_Mutations_log10], y1 = tcga.cohort.med[i, 
                                                                                                    Median_Mutations_log10], col = medianCol)
  })
} else {
  lapply(seq_len(nrow(tcga.cohort.med)), function(i) {
    segments(x0 = i - 1, x1 = i, y0 = tcga.cohort.med[i, 
                                                      Median_Mutations], y1 = tcga.cohort.med[i, Median_Mutations], 
             col = medianCol)
  })
}
tcga.cohort = data.table::rbindlist(l = tcga.cohort)
tcga.cohort[, `:=`(plot_total, NULL)]
tcga.cohort[, `:=`(TCGA, NULL)]
pt.test.pval = pt.test.pval[!is.na(Pval)][order(Pval, decreasing = FALSE)]
list(median_mutation_burden = tcga.cohort.med, mutation_burden_perSample = tcga.cohort, 
     pairwise_t_test = pt.test.pval)

dev.off()
