#!/bin/bash

#SBATCH --job-name=step01
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate rna_envs

### step01
i=$1
cpu=8

#hisat2-build -p 20 Percocypris_pingi.genomic.fna Percocypris_pingi
hisat2 -p ${cpu} --dta --phred33 --summary-file ${i}.stat --new-summary -x Percocypris_pingi -1 rna/${i}_R1.clean.fastq.gz -2 rna/${i}_R2.clean.fastq.gz | samtools view -@ 4 -bS -t Percocypris_pingi.genomic_v2.fna.fai | samtools sort -@ 4 -m 4000000000 -o ${i}.sorted.bam
samtools index -c ${i}.sorted.bam
stringtie ${i}.sorted.bam -p ${cpu} -o ${i}_temp_gtf -A ${i}.fpkm_tracking -G Percocypris_pingi.genomic.gtf
