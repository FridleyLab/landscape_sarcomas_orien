#!/usr/bin/env bash

#from jamie
#1.
# You can exclude the --max_TCC_AF filter in vcf_mut_metrics.pl. It uses the old
# tumor-only targeted gene sequencing data, which isn’t as useful for whole
# exome. I usually keep it since we have it, but have been using it less and
# less.
#2.
# If you don’t have a GQ FORMAT field in your VCF, you should remove the --minGQ
# option in the vcf_mut_metrics.pl step (or change it as Aster’s GQ is likely
# different from ours.)
#3.
# vcf_mut_metrics.pl does require the “ANN” INFO field, which can be added by my
# annotation pipeline (details below). The 1000 Genomes allele frequency AF_1k
# INFO field (a useful filter for common germline snps) is also added by my
# annotation pipeline, but could be excluded.
#4.
# Here are my annotation pipeline scripts:
# /share/NGS/public/ngs_resource/scripts/MRC/wdl/singularity_slurm/*annotate*
#  “run_annotate_cromwell.sh” is the wrapper script to call Cromwell, a pipeline 
#  engine that runs WDL scripts. Annotate_single_vcf.wdl is the wdl script the 
#  includes all the steps and pointers to the singularity images. 
#  “annotate_single_vcf_inputs.json” is the config file that defines all the 
#  inputs. Copy these three files to your run location. Edit the JSON file to 
#  select the right references (probably GRCh38 with Gencode v30 is best) and 
#  delete the other reference blocks. Add the full path to your vcf input file 
#  at the top. Submit the pipeline with “sbatch run_annotate_cromwell.sh”. Once 
#  finished, the output file is in the Cromwell folder structure, which will be 
#  something like 
#      cromwell-executions/annotate_single_vcf/*/call_add_tcga/execution/*vcf.
#5.
# Once annotation is finished, use that file as input for the run_sigest.sh 
# script.



# create file of filenames for each subset
#run inside the mutational signature data folder
find ../../FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads/ -name "*.vcf.gz" > vcf_files.txt

#because they were just zipped and not bgzipped
cd ../../FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads/
#unzip
gunzip *.vcf.gz
#rezip
ls *.vcf | parallel -j 4 "bgzip -c {} > {}.gz"
#index
ls *.vcf.gz | parallel -j 4 "tabix -p vcf {} "

cd ../../Mutation_Signature/data

#merge samples together
bcftools merge --file-list "vcf_files.txt" -o "merge.1170_somatic_filtered.bcf"
#index the bcf
bcftools index "merge.1170_somatic_filtered.bcf"
#convert from bcf format to vcf format 
bcftools convert -O z -o merge.1170_somatic_filtered.vcf.gz merge.1170_somatic_filtered.bcf
#then run 2.0.annotate_wdl/run_annotate_cromwell.sh after editing the json in with the emrged VCF file

#output from annotation move to data folder 
# cp pgm/2.0.annotate_wdl/cromwell-executions/a48ed492-9ed0-4a7d-9d4d-823ecd99833c/call-add_tcga/execution/merge.1170_somatic_filtered.vcf.gz_ann_1k_cos_TCCnorm_ESP_CLNV_ExAC_TCGA.vcf data/

#some of the chromosomes aren't present in the reference, meaning that they are the contigs like 'chr7_KI270803v1_alt'
#have to remove
bcftools view merge.1170_somatic_filtered.vcf.gz_ann_1k_cos_TCCnorm_ESP_CLNV_ExAC_TCGA.vcf.gz -r $(paste -sd ',' target_regions/regions.txt) -o merge.1170_TCGA_annotated.vcf

#bgzip and tabix file
bgzip -c -@ 8 -l 9 merge.1170_TCGA_annotated.vcf > merge.1170_TCGA_annotated.vcf.gz
tabix -p vcf merge.1170_TCGA_annotated.vcf.gz

#need to separate into separate VCF files depending on the capture region that was used for the sample
#pgm/2.1.get-kit-sample-IDs.R

#split the 2.0 output into 2 files that are kit specific
bcftools view -S target_regions/WES_IDT_Sample-IDs.txt merge.1170_TCGA_annotated.vcf.gz > merge.IDT-samples_TCGA_annotated.vcf
bcftools view -S target_regions/WES_NIM_Sample-IDs.txt merge.1170_TCGA_annotated.vcf.gz > merge.NIM-samples_TCGA_annotated.vcf
#zip and index
bgzip -c -@ 8 -l 9 merge.IDT-samples_TCGA_annotated.vcf > merge.IDT-samples_TCGA_annotated.vcf.gz
bgzip -c -@ 8 -l 9 merge.NIM-samples_TCGA_annotated.vcf > merge.NIM-samples_TCGA_annotated.vcf.gz
tabix -p vcf merge.IDT-samples_TCGA_annotated.vcf.gz
tabix -p vcf merge.NIM-samples_TCGA_annotated.vcf.gz

