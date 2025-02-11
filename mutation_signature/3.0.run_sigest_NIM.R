
#!/bin/bash
#SBATCH --mem=12G
#SBATCH --time=72:00:00

# CHANGEME #
OUT=../results/merge.NIM-samples_TCGA_annotated # output file prefix
VCF=merge.NIM-samples_TCGA_annotated.vcf.gz  # bgzipped and tabix indexed
BED=target_regions/hg38_NIM_targets_sorted_merged_extended_sorted_merged.bed  # capture bed file
BEDNAME=NIM  # output name for capture
REF=target_regions/GCA_000001405.15_GRCh38_no_alt_add_chr.fa
mingq=3  # suggest 3 for MuTect/Strelka, 15 for GATK

#cp /share/NGS/public/ngs_resource/scripts/MRC/R/deconstructSigs_cosV3.R .
cwd=`pwd`

#module load bedtools
#module load tabix

singularity exec -C --pwd="$cwd" --bind /share/NGS/public:/share/NGS/public:ro,${cwd}:${cwd} \
/share/NGS/public/ngs_resource/scripts/MRC/singularity_images/tabix.sif \
perl -I /share/NGS/public/ngs_resource/scripts/MRC/ -I /share/NGS/public/apps/bam2mpg/lib/ \
/share/NGS/public/ngs_resource/scripts/MRC/perl/vcf_mut_metrics.pl \
--vcf $VCF \
--bed $BED \
--max_1Kgen 0.01 \
--max_TCC_AF 0.05 \
--extended_mfa $REF \
> ${OUT}.metrics.ext

/share/NGS/public/ngs_resource/scripts/MRC/perl/mut_sig2sigestimattion.pl \
${OUT}.metrics.ext \
> ${OUT}

#singularity exec -C --pwd="$cwd" --bind /share/NGS/public:/share/NGS/public:ro,/share/NGS/work/:/share/NGS/work/:ro,${cwd}:${cwd} \
#    /share/NGS/public/ngs_resource/scripts/MRC/singularity_images/bedtools.sif \
#        perl /share/NGS/public/ngs_resource/scripts/MRC/perl/bed2trinuc.pl \
#            --bed $BED \
#            --genome_fa $REF \
#            --out $BEDNAME \
#            > ${BEDNAME}.trinuc.txt
#gzip ${BEDNAME}.txt

#fails
# sinagularity exec -C --pwd="$cwd" -bind /share/NGS/public:/share/NGS/public:ro,${cwd}:${cwd} \
# /share/NGS/public/ngs_resource/scripts/MRC/singularity_images/signatureestimation_1.0.sif \
# R --vanilla <signature_estimation.R --args ${BEDNAME}.trinuc.txt ${OUT}



