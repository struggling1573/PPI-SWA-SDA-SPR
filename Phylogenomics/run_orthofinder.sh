#!/bin/bash

#SBATCH --job-name=orthofinder
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# Run orthofinder to cluster the protein sequences
orthofinder -f pep_rename -t 8 -a 8 -S blast -og

# Keep single-copy orthologous genes missing in some species
OTU=$1
A1=`csvtk stat -t Orthogroups.GeneCount.tsv | awk '{if(NR==2) print $2}'`
A2=`csvtk stat -t Orthogroups.GeneCount.tsv | awk '{if(NR==2) print $2}' | awk '{print $1-1}'`
cut -f 1-$A2 Orthogroups.GeneCount.tsv > tmp1.txt && awk -F'\t' -v a=$A2 '{for(i=2;i<=a;i++) if($i>1) next}1' tmp1.txt > tmp2.txt && awk 'NR==1{print > "tmp3.txt"; nextfile} {print > "tmp3.txt"}' tmp1.txt tmp2.txt 
awk -F'\t' -v a=$A2 '{sum=0; for(i=2;i<=a;i++) sum+=$i; print $0 "\t" sum}' tmp3.txt | awk -v a=$A1 -v b=$OTU '$a>=b' | cut -f1 > SingleCopyOrthologues_filtered_spe${OTU}.txt && rm tmp*
