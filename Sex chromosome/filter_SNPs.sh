#!/bin/bash
#SBATCH --job-name=CombineGVCFs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=101G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err


ref=ref.fna
tmp=/tmp
cpu=2
GATK=/biosoft/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar
sample_vcf=`cut -f 1 sample_group.txt | sed 1d | awk '{print "--variant 05.gvcf_merge/"$i".vcf.gz"}' | sed ":a;N;s/\n/ /g;ta"`

# CombineGVCFs
java -XX:ParallelGCThreads=${cpu} -Xmx100g -Djava.io.tmpdir=${tmp} -jar ${GATK} CombineGVCFs -R $ref -O variant.CombineGVCFs.vcf $sample_vcf
# GenotypeGVCFs
java -XX:ParallelGCThreads=${cpu} -Xmx100g -Djava.io.tmpdir=${tmp} -jar ${GATK} GenotypeGVCFs -R $ref -V variant.CombineGVCFs.vcf -O variant.Genotype.vcf 
# VcftoolsFilter
java -XX:ParallelGCThreads=${cpu} -Xmx100g -Djava.io.tmpdir=${tmp} -jar ${GATK} VariantFiltration -V variant.Genotype.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O variants_filtered.vcf
# SelectSnp
java -XX:ParallelGCThreads=${cpu} -Xmx100g -Djava.io.tmpdir=${tmp} -jar ${GATK} SelectVariants -R $ref -select-type SNP -V variants_filtered.vcf -O snps_filtered.vcf
# SelectIndel
java -XX:ParallelGCThreads=${cpu} -Xmx100g -Djava.io.tmpdir=${tmp} -jar ${GATK} SelectVariants -R $ref -select-type INDEL -V variants_filtered.vcf -O indels_filtered.vcf
