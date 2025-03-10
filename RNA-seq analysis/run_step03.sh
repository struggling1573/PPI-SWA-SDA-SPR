#!/bin/bash

#SBATCH --job-name=step03
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate rna_envs

[ -d index/ ] || mkdir -p index
rsem-prepare-reference -gtf predict/final.gtf Percocypris_pingi.genomic.fna --bowtie2 rsem_index/genome

# Quantification of gene sets
[ -d map_quantify/ ] || mkdir -p map_quantify
for i in `cat id`;do rsem-calculate-expression --phred33-quals -q --paired-end --strandedness none -p 8 --sort-bam-by-coordinate --bowtie2 --bowtie2-mismatch-rate 0.1 --bowtie2-k 200 --bowtie2-sensitivity-level sensitive rna/${i}_R1.clean.fastq.gz rna/${i}_R2.clean.fastq.gz rsem_index/genome map_quantify/${i};done

for i in `cat id`;do
    cut -f 1,7 map_quantify/${i}.genes.results > ${i}_genefpkm
    sed -i "s/FPKM/$i/g" ${i}_genefpkm
    cut -f 1,6 map_quantify/${i}.genes.results > ${i}_genetpm
    sed -i "s/TPM/$i/g" ${i}_genetpm
    cut -f 1,5 map_quantify/${i}.genes.results > ${i}_generc
done
