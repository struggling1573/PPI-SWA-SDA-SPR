#!/bin/bash

#SBATCH --job-name=LastzIdentity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err


# Define your data
ref_spe=hap1
query_spe=hap2

# Path to software
source /biosoft/miniconda/etc/profile.d/conda.sh
conda activate last
PATH=/biosoft/GenomeAlignmentTools_v1.0/kent/bin:$PATH

[ -d LastzIdentity/01.out_axt/ ] || mkdir -p LastzIdentity/01.out_axt
[ -d LastzIdentity/02.out_chain/ ] || mkdir -p LastzIdentity/02.out_chain
[ -d LastzIdentity/03.out_net/ ] || mkdir -p LastzIdentity/03.out_net
[ -d LastzIdentity/04.out_maf/ ] || mkdir -p LastzIdentity/04.out_maf
[ -d LastzIdentity/05.out_psl/ ] || mkdir -p LastzIdentity/05.out_psl
[ -d LastzIdentity/06.out_score/ ] || mkdir -p LastzIdentity/06.out_score

# Run lastz
seq -w 1 25 | parallel -j 10 "lastz $ref_spe/chr{}.fna $query_spe/chr{}.fna K=4500 L=3000 Y=15000 O=600 E=150 H=2000 T=2 --format=axt > LastzIdentity/01.out_axt/chr{}.axt"

# Get chain from axt
seq -w 1 25 | parallel -j 10 "axtChain LastzIdentity/01.out_axt/chr{}.axt $ref_spe/chr{}.2bit $query_spe/chr{}.2bit LastzIdentity/02.out_chain/chr{}.chain -minScore=5000 -linearGap=medium"

# Netting (Make alignment nets out of chains)
seq -w 1 25 | parallel -j 10 "chainNet LastzIdentity/02.out_chain/chr{}.chain -minSpace=1 $ref_spe/chr{}.sizes $query_spe/chr{}.sizes stdout /dev/null | netSyntenic stdin LastzIdentity/03.out_net/chr{}.net"

# Maffing (Convert net (and chain) to axt, to maf format)
seq -w 1 25 | parallel -j 10 "netToAxt LastzIdentity/03.out_net/chr{}.net LastzIdentity/02.out_chain/chr{}.chain $ref_spe/chr{}.2bit $query_spe/chr{}.2bit stdout | axtSort stdin LastzIdentity/04.out_maf/chr{}.axt"

seq -w 1 25 | parallel -j 10 "axtToMaf LastzIdentity/04.out_maf/chr{}.axt $ref_spe/chr{}.sizes $query_spe/chr{}.sizes LastzIdentity/04.out_maf/chr{}.maf -tPrefix=${ref_spe}. -qPrefix=${query_spe}.  "

# Convert maf to psl format
seq -w 1 25 | parallel -j 10 "mafToPsl ${ref_spe} ${query_spe} LastzIdentity/04.out_maf/chr{}.maf LastzIdentity/05.out_psl/chr{}.psl"

# Calculate web blat score from psl files
seq -w 1 25 | parallel -j 10 "pslScore LastzIdentity/05.out_psl/chr{}.psl > LastzIdentity/06.out_score/chr{}.score"

# filter
mkdir LastzIdentity/07.filter && cd LastzIdentity/07.filter
seq -w 1 25 | parallel -j 10 "awk '\$5>60 && \$6>55 && \$6<100 {a[\$1\"_\"\$2]=3} END { while (getline < \"../05.out_psl/chr{}.psl\") { if (a[\$14\"_\"\$16]==3) { print \$0 } } }' ../06.out_score/chr{}.score > chr{}.score.ide99.filt"
seq -w 1 25 | parallel -j 10 "python psl2sim.py chr{}.score.ide99.filt 50000 chr{}_query.txt chr{}_ref.txt"
cat chr*_query.txt > df_SDA_hap2_vs_hap1.iden_query.txt 
cat chr*_ref.txt > df_SDA_hap2_vs_hap1.iden_ref.txt 

# plot

