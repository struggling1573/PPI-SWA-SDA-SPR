#!/bin/bash

#SBATCH --job-name=whatshap
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# define your data
REF=Percocypris_pingi_hap2.genomic.fna

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate whatshap

[ -d pb/vcf_phased/ ] || mkdir -p pb/vcf_phased
cut -d$'\t' -f1 ${REF}.fai | parallel -j 8 'whatshap polyphase pb/vcf_filtered/{}.filtered.vcf pb/split/{}.rmdup.bam --ploidy 4 --reference Percocypris_pingi_hap2.genomic.fna -o pb/vcf_phased/{}.phased.vcf'

# MergeVcfs(combine all chromosomes)
gatk  --java-options -Xmx12G MergeVcfs $(for i in `cut -d$'\t' -f1 ${REF}.fai`; do echo "-I pb/vcf_phased/${i}.phased.vcf" ;done) -O variant.phased.vcf

bgzip variant.phased.vcf
tabix variant.phased.vcf.gz
#bcftools consensus -H 1 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype1.fasta&
#bcftools consensus -H 2 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype2.fasta&
#bcftools consensus -H 3 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype3.fasta&
#bcftools consensus -H 4 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype4.fasta&

