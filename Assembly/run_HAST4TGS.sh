#!/bin/bash

#SBATCH --job-name=HAST4TGS
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
child_ccs=/path_to_data/ccs.fa 

build_unshared_kmers=/biosoft/HAST4TGS-master/00.build_unshare_kmers_by_meryl/build_unshared_kmers.sh
CLASSIFY_ONLY=/biosoft/HAST4TGS-master/01.classify_reads/CLASSIFY_ONLY.sh

$build_unshared_kmers --paternal $paternal_ngs --maternal $maternal_ngs --memory 100 --thread 20 --auto_bounds
$CLASSIFY_ONLY --paternal_mer paternal.unique.filter.mer --maternal_mer maternal.unique.filter.mer --offspring $child_ccs --thread 20

cat phasing.out | grep "ccs" | awk '{if($3=="haplotype0") print $2}' > id_ccs.haplotype0 &
cat phasing.out | grep "ccs" | awk '{if($3=="haplotype1") print $2}' > id_ccs.haplotype1 &
cat phasing.out | grep "ccs" | awk '{if($3=="ambiguous") print $2}' > id_ccs.ambiguous &
wait
seqtk subseq hifi_pass01.fq.gz id_ccs.haplotype0 | gzip > hifi_haplotype0.fq.gz &
seqtk subseq hifi_pass01.fq.gz id_ccs.haplotype1 | gzip > hifi_haplotype1.fq.gz &
seqtk subseq hifi_pass01.fq.gz id_ccs.ambiguous | gzip > hifi_ambiguous.fq.gz &
wait
seqkit stat hifi*.fq
