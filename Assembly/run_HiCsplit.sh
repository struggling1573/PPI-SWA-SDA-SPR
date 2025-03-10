#!/bin/bash

#SBATCH --job-name=HiCsplit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err


# step01. run Juicer to get the .sam file
REF=PPI.hap2.fa 	# SWA.hap2.fa
juicer=/biosoft/juicer-1.6
$juicer/misc/generate_site_positions.py MboI $REF $REF
bwa index $REF
samtools faidx $REF
cut -f 1,2 ${REF}.fai > ${REF}.sizes
$juicer/CPU/juicer.sh -t 20 -D $juicer -d $PWD -y $PWD/${REF}_MboI.txt -z $PWD/$REF -p $PWD/${REF}.sizes -s MboI -S early


# step02. split hic reads accoding NM value of the alignment file
cut -f 1,2,3,6,12 hic.fastq.sam > hic.fastq.sam.info
cat hic.fastq.sam.info | sed 's/:/\t/g' | awk '{if($5=="NM") print $1"\t"$7}' > MaternalHiC_mapping.info
awk '{print $1"\t"$12}' hic.fastq.sam | sed '/@/d' | sed 's/:/\t/g' | awk '{if($2=="NM")print $1"\t"$4}' > MaternalHiC_mapping.out
python sum_by_first_column.py MaternalHiC_mapping.info MaternalHiC_mapping.info_sum
python sum_by_first_column.py PaternalHiC_mapping.info PaternalHiC_mapping.info_sum
python average_by_first_column.py MaternalHiC_mapping.info MaternalHiC_mapping.info_ava
python average_by_first_column.py PaternalHiC_mapping.info PaternalHiC_mapping.info_ava
less MaternalHiC_mapping.info_sum | sed 's/\"//g' | awk '{print $1"\t"$2"\t""Ma"}' | sed '/Group/d' > df_sum.txt
less PaternalHiC_mapping.info_sum | sed 's/\"//g' | awk '{print $1"\t"$2"\t""Pa"}' | sed '/Group/d' >> df_sum.txt
sort -k 1,1 -k 2n,2 df_sum.txt > df_sum_sorted.txt
python group_by_value.py df_ava_sorted.txt > df_ava_sorted.out

hic_R1=/path_to_data/hic_R1.fastq
hic_R2=/path_to_data/hic_R2.fastq
awk '{if($2=="Pa") print $1"/1"}' df_sum_sorted.out > paternal.reads_1 && seqtk subseq $hic_R1 paternal.reads_1 | gzip > paternal.reads_R1.fq.gz && rm paternal.reads_1 &
awk '{if($2=="Pa") print $1"/2"}' df_sum_sorted.out > paternal.reads_2 && seqtk subseq $hic_R2 paternal.reads_2 | gzip > paternal.reads_R2.fq.gz && rm paternal.reads_2 &
awk '{if($2=="Ma") print $1"/1"}' df_sum_sorted.out > maternal.reads_1 && seqtk subseq $hic_R1 maternal.reads_1 | gzip > maternal.reads_R1.fq.gz && rm maternal.reads_1 &
awk '{if($2=="Ma") print $1"/2"}' df_sum_sorted.out > maternal.reads_2 && seqtk subseq $hic_R2 maternal.reads_2 | gzip > maternal.reads_R2.fq.gz && rm maternal.reads_2 &
wait
echo 'Congratulation, all steps done!'
