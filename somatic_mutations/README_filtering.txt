This is to to document what has been done involving the somatic (and possibly other) VCF files to get to a point of proper mutation rates matching that of TCGA SARC.

Clinical linkage file that was locked on 20221103 in data/clinical was used to identify the samples that had somatic mutation files (out of the somatic mutation files that we have but aren't included in the clinical linkage file), removing those that are not sarcoma (based on Andrew's suggestion/mention).

The script at Analysis_folders/FilteringVCFfiles/FilteringVCF_AF04_FR1read_10supporting.R was run which was edited to filter to F1R2 + F2R1 requiring 5, 10, and 20 reads total. :
1. removed the germline sample from the VCF
2. filtered the results to only those that had a FILTER of "PASS"
3. removed alternate alleles that have less than a frequency of 4% (0.04, keeping those with 4% or more, also tested 5 and 10%)
4. removed alternate alleles with F1R2 or F2R1 with 0 supporting reads in either
5. filtered alleles that didn't meet the particular F1R2+F2R1 read numbers
6. added header entry for an "ADJUSTED" attribute in the header

Output from script was saved to Analysis_folders/FilteringVCFs/somatic_vcfs_tumor_AF04_F1R21_F2R11_*reads_annotated/*.vcf.

The 1170 VCF files of somatic mutations were then annotated with AnnoVar using Analysis_folders/FilteringVCFfiles/Annotate_VCFs_Annovar.slurm which pulls the filtered VCFs, intersects them with the capture region bed files per the kit that they used, then anntoates the results with Annovar.

To convert the annovar output to MAF, the script at Analysis_folders/FilteringVCFfiles/annovarToMaf.R was used. This script uses MAFtools built in function of reading annovar output right into a maf object which can than be saved/exported.

To plot the mutation burden of the top 10 sarcomas by sample number, Analysis_folders/Oncoplot/rpgm/Comparing_total_perMB.R was used and a PDF of the plot was saved to Analysis_folders/Oncoplot/MutationPerMB/Top10Sarcomas_mutationBurden.pdf and edited in Adobe Illustrator to make the text not falling off the page.

The number of mutations that were selected on Nov 21, 2022 by the whole group was the 10 reads in F1R2+F2R1
