#!/bin/bash

#SBATCH --job-name=RModeler
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

REF=/path_to_data/Percocypris_pingi.genomic.fna
db=PPI_db
prefix=Percocypris_pingi
export PATH=/biosoft/RepeatModeler-2.0.3:$PATH
export BLAST_USAGE_REPORT=false

mkdir $db
BuildDatabase -name $db/$db $REF
RepeatModeler -database $db/$db -pa 4
less $db/${db}-families.fa |sed "s/rnd/${prefix}_rnd/g" | sed 's/ .*//g' >${prefix}.RModeler.fa
grep '>' $db/${db}-families.fa |sed "s/rnd/${prefix}_rnd/g" |sed 's/>//g' | sed 's/ (/\t(/g' >${prefix}.RModeler.anno
