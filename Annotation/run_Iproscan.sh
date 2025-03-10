#!/bin/bash

#SBATCH --job-name=Iproscan
#SBATCH --nodelist=comput2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

infile=/path_to_data/Percocypris_pingi.genomic.longest.faa
prefix=Percocypris_pingi

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate interproscan
interproscan=/biosoft/interproscan-5.64-96.0/interproscan.sh
[ -d /tmp/tmp_haoyu/tmp_iproscan/ ] || mkdir -p /tmp/tmp_haoyu/tmp_iproscan

# run interproscan
$interproscan -i $infile -b ${prefix}.iprscan -goterms -iprlookup -pa -dp -cpu 20 --tempdir /tmp/tmp_haoyu/tmp_iproscan 
