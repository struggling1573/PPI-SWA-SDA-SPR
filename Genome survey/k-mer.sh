#!/bin/bash

#SBATCH --job-name=k-mer
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate genomescope
export PATH=/biosoft/KMC_v3.2.2/bin:$PATH

mkdir tmp_k21
kmc -k21 -t20 -m64 -ci1 -cs10000 @fastq db_k21 tmp_k21
kmc_tools transform db_k21 histogram kmcdb_k21.hist -cx10000
genomescope2 -i kmcdb_k21.hist -k 21 -p 4 -o ./ -n k21
