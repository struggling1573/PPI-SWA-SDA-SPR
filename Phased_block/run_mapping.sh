#!/bin/bash

#SBATCH --job-name=mapping
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

ccs=m84208_240707_074535.fastq
minimap2 -t 20 -ax map-hifi -R @RG\\tID:m84208_240707_074535\\tPU:null\\tSM:m84208_240707_074535\\tPL:PACBIO\\tLB:hifi Percocypris_pingi_hap2.genomic.fna $ccs | samtools view -SbF 4 -@ 20 -b - | samtools sort -@ 20 -m 4G - > pacbioccs.fq.bam
