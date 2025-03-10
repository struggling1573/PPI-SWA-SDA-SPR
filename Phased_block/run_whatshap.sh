#!/bin/bash

#SBATCH --job-name=whatshap
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

# define your data
REF=Percocypris_pingi_hap2.genomic.fna

source /biosoft/miniconda3/etc/profile.d/conda.sh
conda activate whatshap

[ -d pb/vcf_phased/ ] || mkdir -p pb/vcf_phased
cut -d$'\t' -f1 ${REF}.fai | parallel -j 8 'whatshap polyphase pb/vcf_filtered/{}.filtered.vcf pb/split/{}.rmdup.bam --ploidy 4 --reference Percocypris_pingi_hap2.genomic.fna -o pb/vcf_phased/{}.phased.vcf'

# MergeVcfs(combine all chromosomes)
gatk  --java-options -Xmx12G MergeVcfs $(for i in `cut -d$'\t' -f1 ${REF}.fai`; do echo "-I pb/vcf_phased/${i}.phased.vcf" ;done) -O variant.phased.vcf

bgzip variant.phased.vcf
tabix variant.phased.vcf.gz
#bcftools consensus -H 1 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype1.fasta&
#bcftools consensus -H 2 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype2.fasta&
#bcftools consensus -H 3 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype3.fasta&
#bcftools consensus -H 4 -f Percocypris_pingi_hap2.genomic.fna variant.phased.vcf.gz > haplotype4.fasta&

seq -w 1 25 | parallel -j 10 "bgzip -c hap2_chr{}.phased.vcf >hap2_chr{}.phased.vcf.gz"
seq -w 1 25 | parallel -j 10 "tabix hap2_chr{}.phased.vcf.gz"

# Extract the largest phased block from the phased vcf
mkdir vcf_phased_LargestBlock && cd vcf_phased_LargestBlock && ln -s ../vcf_phased/* .
for i in `seq -w 1 25`;do
  PS=`grep -v '^#' hap2_chr${i}.phased.vcf | grep 'PS' | cut -d':' -f11 | sort | uniq -cd | sort -k1nr | head -1 | awk '{print $2}'`
  grep -E "^#|$PS$" hap2_chr${i}.phased.vcf > hap2_chr${i}.phased.largestBlock.vcf
done

seq -w 1 25 | parallel -j 10 "bgzip -c hap2_chr{}.phased.largestBlock.vcf >hap2_chr{}.phased.largestBlock.vcf.gz"
seq -w 1 25 | parallel -j 10 "tabix hap2_chr{}.phased.largestBlock.vcf.gz"

# Add phased tags to the bam file
mkdir bam_haplotag_largestBlock/
seq -w 1 25 | parallel -j 10 "whatshap haplotag --ploidy 4 --reference ../Percocypris_pingi_hap2.genomic.fna vcf_phased_LargestBlock/hap2_chr{}.phased.largestBlock.vcf.gz split/hap2_chr{}.rmdup.bam -o bam_haplotag_largestBlock/hap2_chr{}.haplotag.bam"

# Reads of the corresponding hap were extracted according to the label
mkdir largestBlock2fa && cd largestBlock2fa/ && ln -s ../bam_haplotag_largestBlock/* .
for i in `seq -w 1 25`;do samtools view hap2_chr${i}.haplotag.bam |  grep 'HP:i:1' | awk '{print ">"$1"\n"$10}' >chr${i}_hap1.fasta && samtools view hap2_chr${i}.haplotag.bam |  grep 'HP:i:2' | awk '{print ">"$1"\n"$10}' >chr${i}_hap2.fasta && samtools view hap2_chr${i}.haplotag.bam |  grep 'HP:i:3' | awk '{print ">"$1"\n"$10}' >chr${i}_hap3.fasta && samtools view hap2_chr${i}.haplotag.bam |  grep 'HP:i:4' | awk '{print ">"$1"\n"$10}' >chr${i}_hap4.fasta;done

# Run LAST alignment with Onychostoma macrolepis as the reference
#seq -w 1 25 | parallel -j 8 "lastdb -P40 -uNEAR -cR11 OMAchr{} OMAchr{}.fna"
seq -w 1 25 | parallel -j 8 "lastal -P8 -m100 -E0.05 OMAsplit/OMAchr{} PPIchr{}_hap1.fasta | last-split -fMAF+ > chr{}_hap1.maf"
seq -w 1 25 | parallel -j 8 "maf-swap chr{}_hap1.maf | last-split | maf-swap | maf-sort > chr{}_hap1.filter.maf"
seq -w 1 25 | parallel -j 8 "perl /home/user/haoyu/01.script/maf.rename.species.S.pl chr{}_hap1.filter.maf OMAchr{} PPIchr{}_hap1 chr{}_hap1.rename.maf"
seq -w 1 25 | parallel -j 8 "lastal -P8 -m100 -E0.05 OMAsplit/OMAchr{} PPIchr{}_hap2.fasta | last-split -fMAF+ > chr{}_hap2.maf"
seq -w 1 25 | parallel -j 8 "maf-swap chr{}_hap2.maf | last-split | maf-swap | maf-sort > chr{}_hap2.filter.maf"
seq -w 1 25 | parallel -j 8 "perl /home/user/haoyu/01.script/maf.rename.species.S.pl chr{}_hap2.filter.maf OMAchr{} PPIchr{}_hap2 chr{}_hap2.rename.maf"
seq -w 1 25 | parallel -j 8 "lastal -P8 -m100 -E0.05 OMAsplit/OMAchr{} PPIchr{}_hap3.fasta | last-split -fMAF+ > chr{}_hap3.maf"
seq -w 1 25 | parallel -j 8 "maf-swap chr{}_hap3.maf | last-split | maf-swap | maf-sort > chr{}_hap3.filter.maf"
seq -w 1 25 | parallel -j 8 "perl /home/user/haoyu/01.script/maf.rename.species.S.pl chr{}_hap3.filter.maf OMAchr{} PPIchr{}_hap3 chr{}_hap3.rename.maf"
seq -w 1 25 | parallel -j 8 "lastal -P8 -m100 -E0.05 OMAsplit/OMAchr{} PPIchr{}_hap4.fasta | last-split -fMAF+ > chr{}_hap4.maf"
seq -w 1 25 | parallel -j 8 "maf-swap chr{}_hap4.maf | last-split | maf-swap | maf-sort > chr{}_hap4.filter.maf"
seq -w 1 25 | parallel -j 8 "perl /home/user/haoyu/01.script/maf.rename.species.S.pl chr{}_hap4.filter.maf OMAchr{} PPIchr{}_hap4 chr{}_hap4.rename.maf"
for i in `seq -w 1 25`;do python mafConcatenated_byMultiz.py -o chr${i}.maf chr${i}_hap1.rename.maf chr${i}_hap2.rename.maf chr${i}_hap3.rename.maf chr${i}_hap4.rename.maf;done
for i in `seq -w 1 25`;do python maf_ReOrder_v2.py chr${i}.maf PPIchr${i}_hap1,PPIchr${i}_hap2,PPIchr${i}_hap3,PPIchr${i}_hap4,OMAchr${i} > chr${i}.trans.maf;done
for i in `seq -w 1 25`;do python Maf2ConcatFasta.py PPIchr${i}_hap1,PPIchr${i}_hap2,PPIchr${i}_hap3,PPIchr${i}_hap4,OMAchr${i} < chr${i}.trans.maf > chr${i}.trans.maf.fasta; done

# Construct phylogenetic trees from multiple sequence alignment files
for i in `seq -w 1 25`;do iqtree2 -s chr${i}.trans.maf.fasta -T 20 -bb 1000 -o OMAchr${i} ;done
for i in `seq -w 1 25`;do nw_reroot chr${i}.trans.maf.fasta.treefile OMAchr${i} > chr${i}.trans.maf.fasta.treefile.reroot;done
