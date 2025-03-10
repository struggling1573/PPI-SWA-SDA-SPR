#!/bin/bash
#SBATCH --job-name=haphic
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

draft_genome=/path_to_data/PPI_total.fa
HiC_R1=/path_to_data/PaternalHiC_R1.fq.gz
HiC_R2=/path_to_data/PaternalHiC_R2.fq.gz
CPU=12

# PATH to software
source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate haphic_py310
PATH=/biosoft/HapHiC:$PATH
PATH=/biosoft/HapHiC/utils:$PATH

# Align Hi-C data to the assembly
bwa index $draft_genome
bwa mem -t $CPU -5SP $draft_genome  $HiC_R1  $HiC_R2 | samblaster | samtools view - -@ $CPU -S -h -b -F 3340 -o HiC.bam
filter_bam HiC.bam 1 --nm 3 --threads $CPU | samtools view - -b -@ $CPU -o HiC.filtered.bam

# Run HapHiC scaffolding pipeline
haphic pipeline $draft_genome HiC.filtered.bam 49 --RE "GATC" --threads $CPU --remove_allelic_links 2
cd 04.build && bash 04.build/juicebox.sh
cat PPI_haphic_JBAT.liftover.agp | awk '{print $6"\t"$1}' > id_rename.txt
python fasta_transID.py id_rename.txt PPI_total.fa PPI_rename.fa
#python content_ReplaceBy_mappingID.py id_rename.txt merged_nodups.txt merged_nodups.rename
