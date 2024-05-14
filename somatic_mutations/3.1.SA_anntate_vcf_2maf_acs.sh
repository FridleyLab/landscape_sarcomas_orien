## Usage
usage ()
{
   cat << "EOF"

   -------------------------------------------------------------------------------------
      Usage : $0 -v <input_vcf> [ -g <genome_version> -d <OUT_DIR> -n <OUT_NAME> ]

        This function annotates a vcf and convert to maf

        -v <input_vcf>        vcf input
        -g <genome_version>   Reference genome: hg38 (default) or hg19
        -d <OUT_DIR>          Output directory. (default ./)
        -n <OUT_NAME>         Output name
   -------------------------------------------------------------------------------------

EOF
   exit
}

# defaults
GENOME=hg38

# if required variables are not passed directly (such as from qsub),
# Check command line arguments
if [ "$#" -lt 2 ]
   then
      echo VCF input is required.
      usage
else
   while [ "$1" != "" ]; do
      case $1 in
         -g )           shift
                        GENOME=$1
                        ;;
         -v )           shift
                        input_vcf=$1
                        ;;
         -d )           shift
                        OUTDIR=$1
                        ;;
         -n )           shift
                        OUT_NAME=$1
                        ;;
      esac
      shift
   done
fi

if [ -z "${input_vcf}" ]
then
   echo VCF input is required.
   usage
elif [ ! -f "${input_vcf}" ]
then
   echo Specified vcf does not exist.
   usage
fi

if(file $input_vcf | grep -q compressed);
then
   echo "File is compressed - unpacking"
   gunzip $input_vcf
   input_vcf=$(echo $input_vcf | sed 's/.gz//g')
else
   echo "File is uncompressed - continuing"
fi

echo $input_vcf

sleep 1

if [ -n "${OUTDIR}" ] 
then
   if [ -d "${OUTDIR}" ] 
   then
       cd ${OUTDIR}
   else
       echo Specified directory does not exist.
       usage
   fi
fi


## Validate genome
if [[ "${GENOME}" == "hg19" ]]
then
   # ref=/share/NGS/work/ref_seq/hs37d5/hs37d5.fa
   # ref_fai=/share/NGS/work/ref_seq/hs37d5/hs37d5.fa.fai
   # ref_dict=/share/NGS/work/ref_seq/hs37d5/hs37d5.dict
   # ref_bwt=/share/NGS/work/ref_seq/hs37d5/bwa_0.7.10/hs37d5.bwt
   # ref_sa=/share/NGS/work/ref_seq/hs37d5/bwa_0.7.10/hs37d5.sa
   # ref_prefix=/share/NGS/work/ref_seq/hs37d5/bwa_0.7.10/hs37d5
   # ref_ann=/share/NGS/work/ref_seq/hs37d5/bwa_0.7.10/hs37d5.ann
   # ref_pac=/share/NGS/work/ref_seq/hs37d5/bwa_0.7.10/hs37d5.pac
   # ref_amb=/share/NGS/work/ref_seq/hs37d5/bwa_0.7.10/hs37d5.amb
   # db_indel=/share/NGS/work/ref_seq/hs37d5/gatk_bundle_2.8/Mills_and_1000G_gold_standard.indels.b37.vcf
   # dbsnp=/share/NGS/work/ref_seq/hs37d5/gatk_bundle_2.8/dbsnp_138.b37.cleanHeader.vcf.gz
   # mutect1.cosmic=/share/NGS/work/COSMIC/CosmicCodingMuts_v68.vcf.gz
   # genotyperN.interval_bed=/share/NGS/public/ngs_resource/annotations/1000_gen/1000g_0.15MAF.refGene.uniq.bed
   # genotyperT.interval_bed=/share/NGS/public/ngs_resource/annotations/1000_gen/1000g_0.15MAF.refGene.uniq.bed
   # exclude_exp='INFO/1000g2015aug_all > 0.01 | INFO/ExAC_ALL > 0.01 | INFO/esp6500siv2_all > 0.01'
   # combine_exclude_exp='INFO/1000g2015aug_all > 0.01 | INFO/ExAC_ALL > 0.01 | INFO/esp6500siv2_all > 0.01 | QUAL < 2'
   humandb_annovar=/share/NGS/public/ngs_resource/annotations/annov_db/hg19_2019
   genome_version=hg19
   protocol=1000g2015aug_all,exac03,esp6500siv2_all,cosmic70,ensGene,refGene
   operation=f,f,f,f,g,g
elif [[  "${GENOME}" == "hg38" ]]
then
   # ref=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/download/GRCh38.d1.vd1.fa
   # ref_fai=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/download/GRCh38.d1.vd1.fa.fai
   # ref_dict=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/download/GRCh38.d1.vd1.dict
   # ref_bwt=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/index/bwa/GRCh38.d1.vd1.fa.bwt
   # ref_sa=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/index/bwa/GRCh38.d1.vd1.fa.sa
   # ref_prefix=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/index/bwa/GRCh38.d1.vd1.fa
   # ref_ann=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/index/bwa/GRCh38.d1.vd1.fa.ann
   # ref_pac=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/index/bwa/GRCh38.d1.vd1.fa.pac
   # ref_amb=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GRCh38.d1.vd1/index/bwa/GRCh38.d1.vd1.fa.amb
   # dbsnp=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GATK_bundle_hg38/dbsnp_138.hg38.vcf.gz
   # db_indel=/share/NGS/public/ngs_resource/genome_reference/GRCh38/GATK_bundle_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
   # mutect1.cosmic=/share/NGS/work/COSMIC/CosmicCodingMuts_v91_GRCh38_chr.vcf.gz
   # genotyperN.interval_bed=/share/NGS/public/ngs_resource/annotations/1000_gen/1000g_0.15MAF.refGene.uniq.hg38.bed
   # genotyperT.interval_bed=/share/NGS/public/ngs_resource/annotations/1000_gen/1000g_0.15MAF.refGene.uniq.hg38.bed
   # exclude_exp='INFO/1000g2015aug_all > 0.01 | INFO/ExAC_ALL > 0.01 | INFO/esp6500siv2_all > 0.01'
   # combine_exclude_exp='INFO/1000g2015aug_all > 0.01 | INFO/ExAC_ALL > 0.01 | INFO/esp6500siv2_all > 0.01 | QUAL < 2'
   humandb_annovar=/share/NGS/public/ngs_resource/genome_reference/GRCh38/humandb_annovar
   genome_version=hg38
   protocol=1000g2015aug_all,exac03,esp6500siv2_all,cosmic70,ensGene,refGene
   operation=f,f,f,f,g,g
else 
   echo "Error: Unknown genome version."
   usage
fi

## Derive output name from vcf
if [ -z "${OUT_NAME}" ]
then
   OUT_NAME=$(basename $input_vcf|cut -d. -f1)
fi

OUT_NAME=${OUT_NAME}.${genome_version}


## Validate vcf (simple)
if  [[ $(hostname -d) == "local"  ]]
then
   BCFTOOL=/share/apps/bcftools-1.9/bin/bcftools
elif [[ $(hostname -d) == "cm.cluster" ]]
then
   BCFTOOL=/app/eb/software/BCFtools/1.12-GCC-10.2.0/bin/bcftools
else
   BCFTOOL=bcftools
fi


TABLE_ANNOVAR=/share/NGS/public/apps/annovar/table_annovar.pl

bool=true
printf "%s\t" "sampleID" >$(echo $input_vcf | sed 's|\(.*\)/.*|\1|')/${OUT_NAME}.annotated.vcf

TMPDIR=$(echo $input_vcf | sed 's|\(.*\)/.*|\1|')/$OUT_NAME
mkdir -p ${TMPDIR}

for sample in  $(${BCFTOOL} query -l ${input_vcf} )
do
    #c1 only print AC > 1 variant #old INFO is removed 
    ${BCFTOOL} view -s ${sample} -o ${TMPDIR}/tmp-${sample}.vcf ${input_vcf} 

    ${BCFTOOL} annotate -x INFO -Oz -o ${TMPDIR}/${sample}.vcf.gz ${TMPDIR}/tmp-${sample}.vcf

    ${TABLE_ANNOVAR} \
        ${TMPDIR}/${sample}.vcf.gz \
        --vcfinput \
        --outfile ${TMPDIR}/${sample} \
        --buildver ${genome_version} \
        --otherinfo \
        --nastring "." \
        --remove \
        ${humandb_annovar} \
        --protocol ${protocol} \
        --operation ${operation}

    if [ "$bool" = true ]
    then    
        awk -v var="${TMPDIR}/${sample}" 'BEGIN{OFS="\t";getline;print}{print var, $0}' ${TMPDIR}/${sample}.${genome_version}_multianno.txt >>$(echo $input_vcf | sed 's|\(.*\)/.*|\1|')/${OUT_NAME}.annotated.vcf
        bool=false
    else
        awk -v var="${TMPDIR}/${sample}" 'BEGIN{OFS="\t";getline;}     {print var, $0}' ${TMPDIR}/${sample}.${genome_version}_multianno.txt >>$(echo $input_vcf | sed 's|\(.*\)/.*|\1|')/${OUT_NAME}.annotated.vcf
    fi 
    
    rm ${TMPDIR}/${sample}.*
done

#rm -rf ${TMPDIR}

gzip $(echo $input_vcf | sed 's|\(.*\)/.*|\1|')/${OUT_NAME}.annotated.vcf
rm ${input_vcf}

	
# #add all sample ID to maf file to include all samples in downstream analysis
# fieldNum=$( awk 'BEGIN{FS="\t"} NR==1 {print NF}' ${OUT_NAME}.tmp.vcf )
# fieldNum=$(($fieldNum-2))

# rm -f tmp-sample_NA.txt

# for sample in  $( ${BCFTOOL} query -l ${input_vcf} )
# do
#     printf "%s\t" ${sample}>>tmp-sample_NA.txt
#     for i in $( seq 1 $fieldNum )
#     do 
#         printf "%s\t" "NA" >> tmp-sample_NA.txt
#     done; 
#     printf "%s\n" "NA" >> tmp-sample_NA.txt
# done

# cat    ${OUT_NAME}.tmp.vcf tmp-sample_NA.txt > tmp-${OUT_NAME}.tmp.vcf

# mv tmp-${OUT_NAME}.tmp.vcf ${OUT_NAME}.tmp.vcf

rm $TMPDIR/tmp*
rmdir $TMPDIR


