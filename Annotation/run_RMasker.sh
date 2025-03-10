#!/bin/bash

#SBATCH --job-name=RMasker
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

REF=/path_to_data/Percocypris_pingi.genomic.fna
repeate_lib=/path_to_data/Percocypris_pingi.RM.lib
prefix=Percocypris_pingi

export PATH=/biosoft/RepeatMasker-4.1.2-p1:$PATH
export BLAST_USAGE_REPORT=false

RepeatMasker -pa 4 -gff -a -xsmall -gccalc -dir ${prefix}_RM.output -lib $repeate_lib $REF
