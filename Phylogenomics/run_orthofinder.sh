#!/bin/bash

#SBATCH --job-name=orthofinder
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err


source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate orthofinder

orthofinder -f pep_rename -t 8 -a 8 -S blast -og
