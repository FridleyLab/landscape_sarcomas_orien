#!/bin/bash

#SBATCH --job-name AnnoVar
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00-01:30:00
#SBATCH --mem 10G
#SBATCH --output /share/lab_fridley/AVATAR_sarcoma/Analysis_folders/FilteringVCFfiles/annotating_results.out
#SBATCH --array=1-1170

cd /share/lab_fridley/AVATAR_sarcoma/

file=Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads/$(ls Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads/ | grep "\\.vcf.gz$" | sed -n ${SLURM_ARRAY_TASK_ID}p) #grab actual file from array

#gunzip -k $file
vcf=$(sed "s/.gz//" <<< "$(basename $file)")
cp $file Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/$vcf.gz

name=$(sed "s/.vcf.gz//" <<< "$(basename $file)")

if $(zcat $file | grep -q "SentieonCommandLine")
then
	haplocommand=$(zcat $file | grep 'SentieonCommandLine')
	else
	echo -e "not containing SentieonCommandLine"
fi
if $(echo $haplocommand | grep -q "IDT_target")
then
	echo IDT
	target="Analysis_folders/VCFtoMAF/target_regions/hg38_IDT_targets_sorted_merged.bed"
	region_text="IDT" 
else
	echo False
	target="Analysis_folders/VCFtoMAF/target_regions/hg38_NIM_targets_sorted_merged.bed"	#set target name
	region_text="NIM"
fi

module load BEDTools/2.30.0-GCC-11.2.0
gunzip Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/$vcf.gz
bedtools intersect -b $target -a Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/$vcf -header > Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/${name}.${region_text}.vcf


module load BCFtools/1.12-GCC-10.2.0
./Analysis_folders/FilteringVCFfiles/SA_anntate_vcf_2maf_acs.sh -v Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/${name}.${region_text}.vcf -g hg38 -n $name.${region_text}

rm Analysis_folders/FilteringVCFfiles/somatic_vcfs_tumor_AF04_F1R21_F2R11_10reads_annotated/$vcf

