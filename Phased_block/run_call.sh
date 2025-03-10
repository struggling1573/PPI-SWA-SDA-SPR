#!/bin/bash

#SBATCH --job-name=gatk
#SBATCH --nodelist=comput1
#SBATCH --cpus-per-task=20
#SBATCH --mem=101G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# define your data
REF=Percocypris_pingi_hap2.genomic.fna
IN_BAM=pacbioccs.fq.bam
TMPDIR=/tmp

# path to software
export PATH=/biosoft/gatk-4.2.6.1:$PATH
TMPDIR=/tmp
PICARD=/biosoft/picard.jar
GATK_jar=/biosoft/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar

############################################################################################################
###                                        variant calling(by chromosome)                                ###
############################################################################################################
# Build genome index
samtools faidx $REF
for i in `echo $REF | sed 's/\.fna//g'`;do gatk CreateSequenceDictionary -R ${i}.fna -O ${i}.dict;done

# split bam
[ -d pb/split ] || mkdir -p pb/split
[ -d pb/vcf/ ] || mkdir -p pb/vcf
samtools index $IN_BAM
cut -d$'\t' -f1 ${REF}.fai | parallel -j2 'samtools view -@ 10 -b pacbioccs.fq.bam {} > pb/split/{}.bam'
cut -d$'\t' -f1 ${REF}.fai | parallel -j2 'samtools index -@ 10 pb/split/{}.bam'
cut -d$'\t' -f1 ${REF}.fai | parallel -j5 'java -XX:ParallelGCThreads=4 -Xmx100g -jar /biosoft/picard.jar MarkDuplicates I=pb/split/{}.bam O=pb/split/{}.rmdup.bam M=pb/split/{}.m use_jdk_deflater=true use_jdk_inflater=true'
cut -d$'\t' -f1 ${REF}.fai | parallel -j2 'samtools index pb/split/{}.rmdup.bam'

############################################################################################################
###                                               gatk pipeline                                          ###
############################################################################################################
# HaplotypeCaller(#--emit-ref-confidence GVCF)
cut -d$'\t' -f1 ${REF}.fai | parallel -j 5 'gatk --java-options -Xmx12G HaplotypeCaller --native-pair-hmm-threads 8 --tmp-dir /tmp --reference Percocypris_pingi_hap2.genomic.fna --input pb/split/{}.rmdup.bam --output pb/vcf/{}.vcf --emit-ref-confidence GVCF --pcr-indel-model AGGRESSIVE --minimum-mapping-quality 60 --sample-ploidy 4 --intervals {}'

# GenotypeGVCFs
cut -d$'\t' -f1 ${REF}.fai | parallel -j 10 'gatk --java-options -Xmx12G GenotypeGVCFs --tmp-dir /tmp --reference Percocypris_pingi_hap2.genomic.fna --variant pb/vcf/{}.vcf --output pb/vcf/{}.genotype.vcf'

[ -d pb/vcf_filtered/ ] || mkdir -p pb/vcf_filtered
# VcftoolsFilter
cut -d$'\t' -f1 ${REF}.fai | parallel -j 10 'gatk --java-options -Xmx12G VariantFiltration --tmp-dir /tmp --reference Percocypris_pingi_hap2.genomic.fna -V pb/vcf/{}.genotype.vcf -O pb/vcf_filtered/{}.filtered.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter"'

############################################################################################################
###                                             whatshap phasing                                         ###
############################################################################################################
[ -d pb/vcf_phased/ ] || mkdir -p pb/vcf_phased
source /home1/miniconda/etc/profile.d/conda.sh
conda activate whatshap
cut -d$'\t' -f1 ${REF}.fai | parallel -j 10 'whatshap polyphase pb/vcf_filtered/{}.filtered.vcf pb/split/{}.rmdup.bam --ploidy 4 --reference Percocypris_pingi_hap2.genomic.fna -o pb/vcf_phased/{}.phased.vcf'

# MergeVcfs(combine all chromosomes)
gatk  --java-options -Xmx12G MergeVcfs $(for i in `cut -d$'\t' -f1 ${REF}.fai`; do echo "-I pb/vcf_phased/${i}.phased.vcf" ;done) -O variant.phased.vcf

bgzip variant.phased.vcf
tabix variant.phased.vcf.gz
#bcftools consensus -H 1 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype1.fasta&
#bcftools consensus -H 2 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype2.fasta&
#bcftools consensus -H 3 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype3.fasta&
#bcftools consensus -H 4 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype4.fasta&

