#!/bin/bash

#SBATCH --job-name=verkko
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
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

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate verkko

verkko -d $pre --hifi $ccs --nano $ont --hic1 $hicR1 --hic2 $hicR2 --threads 20 --local-memory 100
