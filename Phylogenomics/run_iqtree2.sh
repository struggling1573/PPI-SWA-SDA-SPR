#!/bin/bash

#SBATCH --job-name=iqtree2
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate iqtree

iqtree2 -s allcds.fasta -T 8 -bb 1000
