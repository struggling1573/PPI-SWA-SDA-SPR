#!/bin/bash

#SBATCH --job-name=merqury
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

export PATH=/biosoft/meryl-1.4/bin:$PATH
export PATH=/biosoft/merqury-1.3:$PATH
export MERQURY=/biosoft/merqury-1.3

sh $MERQURY/trio/phase_block.sh hifiasm_trio maternal_compress.k30.hapmer.meryl paternal_compress.k30.hapmer.meryl out
sh $MERQURY/trio/hap_blob.sh maternal_compress.k30.hapmer.meryl paternal_compress.k30.hapmer.meryl SWA_total.fa PPI_total.fa out
sh $MERQURY/trio/hamming_error.sh out.hapmers.count SWA_total.fa PPI_total.fa
