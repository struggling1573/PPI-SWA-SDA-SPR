#!/bin/bash

#SBATCH --job-name=gwas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate vcf2gwas 

vcf2gwas -v snps_filtered.vcf.gz -pf example.csv -ap -lmm -T 8
