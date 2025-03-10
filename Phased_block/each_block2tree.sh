# Get block position
for i in `seq -w 1 25`;do grep -v '^#' ckhap2_chr${i}.phased.vcf | grep 'PS' | cut -d':' -f11 | sort | uniq -cd | sort -k1nr | awk '{print $2}' > ckhap2_chr${i}.phased.block_start;done
for i in `seq -w 1 25`;do
        for start_id in `cat ckhap2_chr${i}.phased.block_start`;do
                end_id=`grep -E "$start_id$" ckhap2_chr${i}.phased.vcf | awk '{print $1"\t"$2"\t"$10}' | sed 's/:/\t/g' | awk '$8<=a' a=$start_id | tail -n 1 | awk '{print $2}'`
                echo -e "$start_id\t$end_id" >> ckhap2_chr${i}.phased.blockpos;
        done;
done

# Extract the phased vcf file corresponding to each block in a chromosome
[ -d result_block/ ] || mkdir -p result_block
num=`cat hap2_chr19.phased.blockpos | wc -l`

for i in `seq 1 $num`;do
        start_id=`cat hap2_chr19.phased.blockpos | head -n $i | tail -n 1 | awk '{print $1}'`
        end_id=`cat hap2_chr19.phased.blockpos | head -n $i | tail -n 1 | awk '{print $2}'`
        grep -E "^#" hap2_chr19.phased.vcf >> result_block/chr19_${start_id}_${end_id}.vcf
        grep -E "$start_id$" hap2_chr19.phased.vcf | awk -v a=$end_id '$2<=a' >> result_block/chr19_${start_id}_${end_id}.vcf
done

# Get haplotype sequences from the phased vcf files
REF=/path_to_data/Percocypris_pingi_hap2.genomic.fna
InBam=/path_to_data/hap2_chr19.rmdup.bam
ls result_block/*vcf | sed 's/\.vcf//g' | sed 's/\//\t/g' | cut -f2 | parallel -j 20 "bgzip -c result_block/{}.vcf > result_block/{}.vcf.gz"
ls result_block/*vcf | sed 's/\.vcf//g' | sed 's/\//\t/g' | cut -f2 | parallel -j 20 "tabix result_block/{}.vcf.gz"
i=$1
whatshap haplotag --ploidy 4 --reference $REF result_block/${i}.vcf.gz $InBam -o result_block/${i}.haplotag.bam
samtools view result_block/${i}.haplotag.bam |  grep 'HP:i:1' | awk '{print ">"$1"\n"$10}' > result_block/${i}_hap1.fasta&
samtools view result_block/${i}.haplotag.bam |  grep 'HP:i:2' | awk '{print ">"$1"\n"$10}' > result_block/${i}_hap2.fasta&
samtools view result_block/${i}.haplotag.bam |  grep 'HP:i:3' | awk '{print ">"$1"\n"$10}' > result_block/${i}_hap3.fasta&
samtools view result_block/${i}.haplotag.bam |  grep 'HP:i:4' | awk '{print ">"$1"\n"$10}' > result_block/${i}_hap4.fasta&
wait
rm result_block/${i}.haplotag.bam;
#ls result_block/*vcf | sed 's/\.vcf//g' | sed 's/\//\t/g' | cut -f2 | parallel -j 4 "bash 02.get_fasta.sh {}"

# Filter fasta
cd result_block
ll *fasta | awk '{print $9"\t"$5}' > ../fasta_count.txt
cd ..
cat fasta_count.txt | awk '$2=="0"' | sed 's/_hap.*//g' | sort | uniq | grep -v -f - fasta_count.txt > fasta_count.filter.txt
mkdir block2fa
for i in `cat fasta_count.filter.txt | cut -f 1`;do cp result_block/$i block2fa/;done

# Run LAST alignment with Onychostoma macrolepis as the reference
i=$1
[ -d WGA/$i/ ] || mkdir -p WGA/$i/
lastal -P8 -m100 -E0.05 00.ppha2_for_Omac/Omac_chr01 result_block_filter/${i}_hap1.fasta | last-split -fMAF+ >WGA/$i/${i}_hap1.maf
maf-swap WGA/$i/${i}_hap1.maf | last-split | maf-swap | maf-sort > WGA/$i/${i}_hap1.filter.maf
perl maf.rename.species.S.pl WGA/$i/${i}_hap1.filter.maf Omac pp_hap1 WGA/$i/${i}_hap1.rename.maf
lastal -P8 -m100 -E0.05 00.ppha2_for_Omac/Omac_chr01 result_block_filter/${i}_hap2.fasta | last-split -fMAF+ >WGA/$i/${i}_hap2.maf
maf-swap WGA/$i/${i}_hap2.maf | last-split | maf-swap | maf-sort > WGA/$i/${i}_hap2.filter.maf
perl maf.rename.species.S.pl WGA/$i/${i}_hap2.filter.maf Omac pp_hap2 WGA/$i/${i}_hap2.rename.maf
lastal -P8 -m100 -E0.05 00.ppha2_for_Omac/Omac_chr01 result_block_filter/${i}_hap3.fasta | last-split -fMAF+ >WGA/$i/${i}_hap3.maf
maf-swap WGA/$i/${i}_hap3.maf | last-split | maf-swap | maf-sort > WGA/$i/${i}_hap3.filter.maf
perl maf.rename.species.S.pl WGA/$i/${i}_hap3.filter.maf Omac pp_hap3 WGA/$i/${i}_hap3.rename.maf
lastal -P8 -m100 -E0.05 00.ppha2_for_Omac/Omac_chr01 result_block_filter/${i}_hap4.fasta | last-split -fMAF+ >WGA/$i/${i}_hap4.maf
maf-swap WGA/$i/${i}_hap4.maf | last-split | maf-swap | maf-sort > WGA/$i/${i}_hap4.filter.maf
perl maf.rename.species.S.pl WGA/$i/${i}_hap4.filter.maf Omac pp_hap4 WGA/$i/${i}_hap4.rename.maf
multiz M=1 WGA/$i/${i}_hap1.rename.maf WGA/$i/${i}_hap2.rename.maf 0 U1_${i} U2_${i} > WGA/$i/${i}.tmp1.maf
multiz M=1 WGA/$i/${i}.tmp1.maf WGA/$i/${i}_hap3.rename.maf 0 U1_${i} U2_${i} > WGA/$i/${i}.tmp2.maf
multiz M=1 WGA/$i/${i}.tmp2.maf WGA/$i/${i}_hap4.rename.maf 0 U1_${i} U2_${i} > WGA/$i/${i}.final.maf
python Maf2ConcatFasta.py Omac,pp_hap1,pp_hap2,pp_hap3,pp_hap4 < WGA/$i/${i}.final.maf > WGA/$i/${i}.final.fasta
rm U*
rm WGA/$i/${i}.tmp*

# Construct phylogenetic trees from multiple sequence alignment files
ls 00.data/*fasta | sed 's/\.final\.fasta//g' | sed 's/\//\t/g' | cut -f2 | parallel -j 20 'trimal -in 00.data/{}.final.fasta -out 00.data/{}.filter.aln'
ls 00.data/*fasta | sed 's/\.final\.fasta//g' | sed 's/\//\t/g' | cut -f2 | parallel -j 1 'iqtree2 -s 00.data/{}.filter.aln -T 20 -bb 1000'
ls 00.data/*fasta | sed 's/\.final\.fasta//g' | sed 's/\//\t/g' | cut -f2 | parallel -j 20 'nw_reroot 00.data/{}.filter.aln.treefile Omac > 00.data/{}.reroot.tree'
#ls 00.data/*fasta | sed 's/\.final\.fasta//g' | sed 's/\//\t/g' | cut -f2 | parallel -j 20 'nw_topology -I 00.data/{}.reroot.tree > 00.data/{}.topology.tree'
