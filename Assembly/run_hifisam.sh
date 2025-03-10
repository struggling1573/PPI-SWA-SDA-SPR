#!/bin/bash

#SBATCH --job-name=hifiasm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err


# define your data
ccs=/path_to_data/ccs_withHUST/PPI*.fa
ont=/path_to_data/ont_withHUST/PPI.pass.fq
hicR1=/path_to_data/PPI_hic*_R1.fastq
hicR2=/path_to_data/PPI_hic*_R2.fastq
pre=PPI

export PATH=/biosoft/hifiasm-0.19.8:$PATH

# run hifisam
hifiasm -o $pre -t 16 -l 2 $ccs --ul-rate 0.1 --ul $ont --h1 $hicR1 --h2 $hicR2

# --n-hap 3 for SWA and SPR 
#hifiasm -o $pre -t 20 -l 2 $ccs --ul-rate 0.1 --ul $ont --h1 $hicR1 --h2 $hicR2 --n-hap 3
