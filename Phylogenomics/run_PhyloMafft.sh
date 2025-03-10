#!/bin/bash

#SBATCH --job-name=PhyloMafft
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# PATH to software
source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate iqtree
export PATH=/biosoft/EMBOSS-6.6.0/bin:$PATH
export PATH=/biosoft/muscle_v5.1:$PATH

mkdir tmp_dir && cd tmp_dir
mkdir cds_single out_aln allcds_tmp

############################################################################################################
###                                        extract SingleCopyOrthologues                                 ###
############################################################################################################
ln -s ../../02.singleCopy/Orthogroups .
cat ../../01.cds_prepare/cds_rename/* >>species_all.cds
cat Orthogroups/SingleCopyOrthologues_filtered_spe10.txt | parallel -j 20 "grep {} Orthogroups/Orthogroups.txt | sed 's/ /\\n/g' | seqkit grep -f - species_all.cds >cds_single/{}.cds"

############################################################################################################
###                                                    MAFFT                                             ###
############################################################################################################
cat Orthogroups/SingleCopyOrthologues_filtered_spe10.txt | parallel -j 20 'mafft --maxiterate 1000 --localpair cds_single/{}.cds >out_aln/{}.cds.aln'
cat Orthogroups/SingleCopyOrthologues_filtered_spe10.txt | parallel -j 20 'trimal -in out_aln/{}.cds.aln -out out_aln/{}.cds.aln.filter -gt 0.1'
cat Orthogroups/SingleCopyOrthologues_filtered_spe10.txt | parallel -j 20 'seqkit seq out_aln/{}.cds.aln.filter -w 0 > out_aln/{}.new.cds'

############################################################################################################
###                                           allcds Concatenation                                       ###
############################################################################################################
cp -r out_aln/*.new.cds allcds_tmp/

sed -i 's/>Acro_fasc_gene_.*/>Acro_fasc_gene/g' allcds_tmp/*
sed -i 's/>Acro_wenc_gene_.*/>Acro_wenc_gene/g' allcds_tmp/*
sed -i 's/>Aspi_lati_gene_.*/>Aspi_lati_gene/g' allcds_tmp/*
sed -i 's/>Cara_aura_subA_.*/>Cara_aura_subA/g' allcds_tmp/*
sed -i 's/>Cara_aura_subB_.*/>Cara_aura_subB/g' allcds_tmp/*
sed -i 's/>Cypr_carp_subA_.*/>Cypr_carp_subA/g' allcds_tmp/*
sed -i 's/>Cypr_carp_subB_.*/>Cypr_carp_subB/g' allcds_tmp/*
sed -i 's/>Dani_reri_gene_.*/>Dani_reri_gene/g' allcds_tmp/*
sed -i 's/>Gymn_eckl_eckl_.*/>Gymn_eckl_eckl/g' allcds_tmp/*
sed -i 's/>Gymn_prze_hap1_.*/>Gymn_prze_hap1/g' allcds_tmp/*
sed -i 's/>Luci_capi_subA_.*/>Luci_capi_subA/g' allcds_tmp/*
sed -i 's/>Luci_capi_subB_.*/>Luci_capi_subB/g' allcds_tmp/*
sed -i 's/>Onyc_macr_gene_.*/>Onyc_macr_gene/g' allcds_tmp/*
sed -i 's/>Onyc_raru_hap1_.*/>Onyc_raru_hap1/g' allcds_tmp/*
sed -i 's/>Oxyg_stew_gene_.*/>Oxyg_stew_gene/g' allcds_tmp/*
sed -i 's/>Perc_ping_hap1_.*/>Perc_ping_hap1/g' allcds_tmp/*
sed -i 's/>Perc_ping_hap2_.*/>Perc_ping_hap2/g' allcds_tmp/*
sed -i 's/>Poro_huan_gene_.*/>Poro_huan_gene/g' allcds_tmp/*
sed -i 's/>Proc_raba_subA_.*/>Proc_raba_subA/g' allcds_tmp/*
sed -i 's/>Proc_raba_subB_.*/>Proc_raba_subB/g' allcds_tmp/*
sed -i 's/>Punt_tetr_gene_.*/>Punt_tetr_gene/g' allcds_tmp/*
sed -i 's/>Scap_acan_gene_.*/>Scap_acan_gene/g' allcds_tmp/*
sed -i 's/>Schi_mala_gene_.*/>Schi_mala_gene/g' allcds_tmp/*
sed -i 's/>Schi_pylz_gene_.*/>Schi_pylz_gene/g' allcds_tmp/*
sed -i 's/>Schi_youn_geno_.*/>Schi_youn_geno/g' allcds_tmp/*
sed -i 's/>Schi_davi_hap1_.*/>Schi_davi_hap1/g' allcds_tmp/*
sed -i 's/>Schi_davi_hap2_.*/>Schi_davi_hap2/g' allcds_tmp/*
sed -i 's/>Schi_lant_gene_.*/>Schi_lant_gene/g' allcds_tmp/*
sed -i 's/>Schi_ocon_gene_.*/>Schi_ocon_gene/g' allcds_tmp/*
sed -i 's/>Schi_pren_hap1_.*/>Schi_pren_hap1/g' allcds_tmp/*
sed -i 's/>Schi_pren_hap2_.*/>Schi_pren_hap2/g' allcds_tmp/*
sed -i 's/>Schi_pren_hap3_.*/>Schi_pren_hap3/g' allcds_tmp/*
sed -i 's/>Schi_wang_hap1_.*/>Schi_wang_hap1/g' allcds_tmp/*
sed -i 's/>Schi_wang_hap2_.*/>Schi_wang_hap2/g' allcds_tmp/*
sed -i 's/>Schi_wang_hap3_.*/>Schi_wang_hap3/g' allcds_tmp/*
sed -i 's/>Spin_sine_subA_.*/>Spin_sine_subA/g' allcds_tmp/*
sed -i 's/>Spin_sine_subB_.*/>Spin_sine_subB/g' allcds_tmp/*

ls allcds_tmp/* >allcds.list
catsequences allcds.list
mv allseqs.fas ../allcds.fasta
mv allseqs.partitions.txt ../allcds.partitions.txt

cd .. && rm -rf tmp_dir
