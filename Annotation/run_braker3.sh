#!/bin/bash

#SBATCH --job-name=braker3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

MaskedRef=/path_to_data/Percocypris_pingi.genomic.fna.masked
SpeciesID=Percocypris_pingi
HomoPep=/path_to_data/sp4_pep.purged_rename.faa
PATH2RNA=/path_to_data/rna_data
RNAID=pp_bl_1,pp_br_1,pp_fi_1,pp_gi_1,pp_gu_1,pp_he_1,pp_li_1,pp_mu_3,pp_mu_4,pp_sk_1,pp_sp_1,pp_tk_1
OutDIR=brakerout_rna
CPU=20

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate braker3
cp -r  /biosoft/miniconda3/envs/braker3/config .
export AUGUSTUS_CONFIG_PATH=$PWD/config
export PATH=/biosoft/GeneMark-ETP/bin/:$PATH
export PATH=/biosoft/GeneMark-ETP/tools/:$PATH
export PATH=/biosoft/GeneMark-ETP/bin/gmes/:$PATH
export PATH=/biosoft/ProtHint_v2.6.0/bin:$PATH

braker.pl --genome=$MaskedRef --species=$SpeciesID --prot_seq=$HomoPep --rnaseq_sets_ids=$RNAID --rnaseq_sets_dir=$PATH2RNA --workingdir=$OutDIR --gff3 --threads $CPU
