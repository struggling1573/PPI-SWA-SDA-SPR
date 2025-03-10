#!/bin/bash

#SBATCH --job-name=kofamscan
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

infile=/path_to_data/Percocypris_pingi.genomic.longest.faa
prefix=Percocypris_pingi

ko_list=/path_to_database/kegg_annot/ko_list
profiles=/path_to_database/kegg_annot/profiles
source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan
[ -d /tmp/tmp_haoyu/tmp_kofamscan/ ] || mkdir -p /tmp/tmp_haoyu/tmp_kofamscan

# run exec_annotation
exec_annotation -k $ko_list -p $profiles -o ${prefix}.querry2KO --cpu 20 --format mapper -E 1e-5 $infile --tmp-dir /tmp/tmp_haoyu/tmp_kofamscan

