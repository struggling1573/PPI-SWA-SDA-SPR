#!/bin/bash

#SBATCH --job-name=hifiasm-trio
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# define your data
ccs=/path_to_data/sandai_ccs/*fa
ont=/path_to_data/pass.fq.gz

export PATH=/bio_soft/hifiasm-0.24.0:$PATH

yak count -k31 -b37 -t16 -o pat.yak paternal.fq.gz
yak count -k31 -b37 -t16 -o mat.yak maternal.fq.gz
hifiasm -o h1 -t 16 -1 pat.yak -2 mat.yak $ccs --ul-rate 0.1 --ul $ont
