module load bcftools-1.21
module load plink/1.9.0
#fst_calculation_script
cd /gpfs/data/user/yuvraj/selection/apoe/intersection/isec_output_AA_AMR
for group in AA DR TB IE IND; do
    for pop in AFR AMR EAS SAS EUR; do
bcftools merge 0000.vcf.gz 0001.vcf.gz -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz
plink --vcf merged.vcf.gz --make-bed --out merged_plink --allow-extra-chr
awk '{print $1, $2}' merged_plink.fam > sample_list.txt
bcftools query -l 0000.vcf.gz > pop1_samples.txt
bcftools query -l 0001.vcf.gz > pop2_samples.txt
bash cluster.sh
plink --bfile merged_plink --fst --within clusters.txt --out fst_results
