#!/bin/bash
#SBATCH --job-name=HapCaller
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

sample=$1
ref=ref.fna
index_dict=ref.dict

tmp=/tmp
cpu=2
ploid=4
GATK=/home/bio_soft/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar

# call SNP by chromosome
[ -d 03.ChrBed/ ] || mkdir -p 03.ChrBed
[ -d 04.gvcf_byChr/ ] || mkdir -p 04.gvcf_byChr
#cat ref.fna.fai | cut -f 1 |while read i ;do cat ref.fna.fai |awk '{print $1"\t"0"\t"$2}' |grep -w ${i}> 03.ChrBed/${i}.bed ;done
for chr in `cut -f 1 ref.fna.fai | cat`;do java -XX:ParallelGCThreads=${cpu} -Xmx40g -Djava.io.tmpdir=${tmp} -jar $GATK HaplotypeCaller --emit-ref-confidence GVCF -R $ref -O 04.gvcf_byChr/${sample}.${chr}.vcf.gz -I 02.markDup/${sample}.dedup.bam --sample-ploidy $ploid -L 03.ChrBed/${chr}.bed;done

# MergeVcfs
[ -d 05.gvcf_merge/ ] || mkdir -p 05.gvcf_merge
realpath 04.gvcf_byChr/${sample}.*gz > ${sample}.gvcf_list
java -XX:ParallelGCThreads=${cpu} -Xmx50g -Djava.io.tmpdir=${tmp} -jar $GATK MergeVcfs I=${sample}.gvcf_list O=05.gvcf_merge/${sample}.vcf.gz
