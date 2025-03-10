#!/bin/bash

#SBATCH --job-name=ploidy_survey
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate smudgeplot

L=100
U=500
echo $L $U

kmc_tools transform db_k21 -ci"$L" -cx"$U" dump -s db_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o db_L"$L"_U"$U" < db_L"$L"_U"$U".dump
smudgeplot.py plot -o k21 db_L"$L"_U"$U"_coverages.tsv
