#!/bin/bash

#SBATCH --job-name=busco
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

ref=$1
prefix=$2
busco_db=/biosoft/busco/actinopterygii_odb10
source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate busco

busco -i $ref -o $prefix -m genome -l $busco_db -c 12 --offline
