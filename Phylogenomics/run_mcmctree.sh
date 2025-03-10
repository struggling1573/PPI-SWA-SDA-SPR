#!/bin/bash

#SBATCH --job-name=mcmctree
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# define yourdata
baseml_time=0.90
in_aln=`realpath allcds.fasta`
ref_tree=`realpath allcds.fasta.treefile.reroot`
iqtree_log=`realpath allcds.fasta.log`
time_tree=`realpath time_lable.tree`

cp -r /biosoft/pipeline_mcmctree/* .
source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate iqtree

############################################################################################################
###                                                step01                                                ###
############################################################################################################
cd step01
nw_topology -I $ref_tree >baseml.tree
sp_num=`nw_labels baseml.tree | wc -l`
cp $in_aln baseml.fas
# modify you treefile
sed -i '1i '"$sp_num"' 1' baseml.tree
export baseml_time
sed -i "s/;/'@"$baseml_time"';/g" baseml.tree
# run baseml
baseml baseml.ctl
# calculate rgene_gamma
Substitution_mlb=`grep -A 2 Substitution mlb |grep '+-' |awk '{print $1}'`
rgene_gamma=`echo "scale=8;1/$Substitution_mlb" | bc`
export rgene_gamma
# calculate Gamma_shape
export ref_tree
export iqtree_log
Gamma_shape=`less $iqtree_log | grep "alpha" | tail -n 1 | awk '{print $4}'`
export Gamma_shape

############################################################################################################
###                                                step02                                                ###
############################################################################################################
cd ../step02
cp ../step01/baseml.fas all.fa
cp $time_tree all.tre
# modify mcmctree.ctl
sed -i "s/0.28/$baseml_time/g" mcmctree.ctl
sed -i "s/3.021/$Gamma_shape/g" mcmctree.ctl
sed -i "s/2.59224504/$rgene_gamma/g" mcmctree.ctl
# run mcmctree
mcmctree mcmctree.ctl
cd ..

############################################################################################################
###                                                step03                                                ###
############################################################################################################
mkdir step03 && cd step03
# modify mcmctree.ctl
ln -s ../step02/all.* .
ln -s ../step02/out.BV in.BV
cp ../step02/mcmctree.ctl .
sed -i 's/usedata = 3/usedata = 2/g' mcmctree.ctl
# run mcmctree
mcmctree mcmctree.ctl
cd ..

############################################################################################################
###                                                step04                                                ###
############################################################################################################
mkdir step04 && cd step04
# modify mcmctree.ctl
ln -s ../step02/all.* .
ln -s ../step02/out.BV in.BV
cp ../step02/mcmctree.ctl .
sed -i 's/usedata = 3/usedata = 2/g' mcmctree.ctl
# run mcmctree
mcmctree mcmctree.ctl
