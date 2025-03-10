#!/bin/bash

#SBATCH --job-name=FunAnno
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

infile=/path_to_data/Percocypris_pingi.genomic.longest.faa
prefix=Percocypris_pingi

sprot_db=/path_to_database/uniprot_sprot_20231108
trembl_db=/path_to_database/uniprot_trembl_20231108
nr_db=/path_to_database/nr_20230728
[ -d ncbiOut/ ] || mkdir -p ncbiOut

## Swiss-Prot
[ -d /tmp/tmp_haoyu/tmp_sprot/ ] || mkdir -p /tmp/tmp_haoyu/tmp_sprot
diamond blastp --db $sprot_db --out ncbiOut/${prefix}_sprot.out --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles --query $infile --max-target-seqs 1 --evalue 1e-5 --sensitive --max-hsps 1 --threads 20 --block-size 20 --index-chunks 1 --tmpdir /tmp/tmp_haoyu/tmp_sprot
cat ncbiOut/${prefix}_sprot.out | awk -F '\t' '{print $1"\t"$17}' > ncbiOut/${prefix}_sprot.anno

## TrEMBL
[ -d /tmp/tmp_haoyu/tmp_trembl/ ] || mkdir -p /tmp/tmp_haoyu/tmp_trembl
diamond blastp --db $trembl_db --out ncbiOut/${prefix}_trembl.out --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles --query $infile --max-target-seqs 1 --evalue 1e-5 --sensitive --max-hsps 1 --threads 20 --block-size 20 --index-chunks 1 --tmpdir /tmp/tmp_haoyu/tmp_trembl
cat ncbiOut/${prefix}_trembl.out | awk -F '\t' '{print $1"\t"$17}' > ncbiOut/${prefix}_trembl.anno

## NR
[ -d /tmp/tmp_haoyu/tmp_nr/ ] || mkdir -p /tmp/tmp_haoyu/tmp_nr
diamond blastp --db $nr_db --out ncbiOut/${prefix}_nr.out --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles --query $infile --max-target-seqs 1 --evalue 1e-5 --sensitive --max-hsps 1 --threads 20 --block-size 20 --index-chunks 1 --tmpdir /tmp/tmp_haoyu/tmp_nr
cat ncbiOut/${prefix}_nr.out | sed 's/<>/\t/g' | awk -F '\t' '{print $1"\t"$17}' > ncbiOut/${prefix}_nr.anno

