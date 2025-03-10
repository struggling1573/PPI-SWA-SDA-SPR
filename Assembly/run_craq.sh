#!/bin/bash

#SBATCH --job-name=PPI_craq
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# define your data
REF=/path_to_data/PPI_total.fa
SMS=/path_to_data/PPI_ccs.fa
NGS_R1=/path_to_data/PPI_mu_3_R1.clean.fq.gz
NGS_R2=/path_to_data/PPI_mu_3_R2.clean.fq.gz
craq=/biosoft/CRAQ-v1.0.9-alpha/bin/craq

$craq -g $REF -sms $SMS -ngs $NGS_R1,$NGS_R2 --map map-hifi --plot --thread 20
