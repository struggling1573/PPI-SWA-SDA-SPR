#!/bin/bash

#SBATCH --job-name=edta
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err


REF=/path_to_data/Percocypris_pingi.genomic.fna

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

EDTA.pl --genome $REF --species others --sensitive 0 --anno 1 --threads 16
