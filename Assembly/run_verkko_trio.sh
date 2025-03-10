#!/bin/bash

#SBATCH --job-name=verkko-trio
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# define your data
paternal_ngs=/path_to_data/PPI_mu_3.total.fq
maternal_ngs=/path_to_data/SWA_mu_1.total.fq
child_ngs=/path_to_data/h2.total.fq
paternal_ont=/path_to_data/PPI_mu_3.ont.fq
maternal_ont=/path_to_data/SWA_mu_1.ont.fq
ccs=/path_to_data/sandai_ccs/*fa
ont=/path_to_data/pass.fq.gz

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate verkko
export PATH=/biosoft/meryl-1.4/bin:$PATH
export PATH=/biosoft/merqury-1.3:$PATH
export MERQURY=/biosoft/merqury-1.3

meryl count compress k=31 threads=20 memory=20 $maternal_ngs $maternal_ont output maternal_compress.k31.meryl
meryl count compress k=31 threads=20 memory=20 $paternal_ngs $paternal_ont output paternal_compress.k31.meryl
meryl count compress k=31 threads=20 memory=20 $child_ngs output child_compress.k31.meryl
$MERQURY/trio/hapmers.sh maternal_compress.k31.meryl paternal_compress.k31.meryl child_compress.k31.meryl

verkko -d verrkko_out --hifi $ccs --nano $ont --hap-kmers maternal_compress.k31.hapmer.meryl paternal_compress.k31.hapmer.meryl trio --threads 20 --local-memory 100
