# Merge the GFF file
ls *_temp_gtf > mergelist.txt
stringtie --merge -p 4 -o merge.gtf mergelist.txt
mkdir merge predict
sed 's/gene_name ".*"//g' Percocypris_pingi.genomic.gtf > ref_new.gtf
cuffcompare -r ref_new.gtf -e 100 -C -s Percocypris_pingi.genomic.fna -o compared.stringtie.gtf merge.gtf
awk '{if($3 ~ /^(u|i|o|j)$/) print $1"\t"$2"\t"$3"\t"$4"\t"$5}' compared.stringtie.gtf.merge.gtf.tmap > gene_transcript.list
python get_new_transcript.py -i merge.gtf -fa Percocypris_pingi.genomic.fna -l gene_transcript.list -o new_transcript
conda activate agat
agat_convert_sp_gff2gtf.pl --gff new_transcript.gtf -o new_transcript.tmp1
agat_convert_sp_gxf2gxf.pl -g new_transcript.tmp1 -o new_transcript.tmp2
python gff_addCDS.py -i new_transcript.tmp2 -o new_transcript.tmp3
cat new_transcript.tmp3 | awk -F '\t' -v OFS=',' '{$1=$1;print}' | sort -t, -k1,1 -k4,5n | awk -F ',' -v OFS='\t' '{$1=$1;print}' > new_transcript.tmp4
grep -v "#" new_transcript.tmp4 | cut -f1 | sort | uniq | awk '{print $1"\t"$1}' > change.bed
python gff_GeneRname.py -g new_transcript.tmp4 -c change.bed -a 10 -o new_transcript.tmp5 -s newG -m_gene gene_mapping.txt -m_mRNA mRNA_mapping.txt
cp new_transcript.tmp5 new_transcript_final.gff && rm *tmp* change.bed
gffread -T new_transcript_final.gff > new_transcript_final.gtf
gffread -g ../Percocypris_pingi.genomic.fna new_transcript_final.gff -x new_transcript_final.cds -y new_transcript_final.faa

# Get longest transcript sequence of the predicted new transcript
gffio view -L new_transcript_final.gff > new_transcript.longest.gff    
gffread -g ../Percocypris_pingi.genomic.fna new_transcript.longest.gff -x new_transcript.longest.cds -y new_transcript.longest.faa

# Get the final GFF file
python /biosoft/CPC2_standalone-1.0.1/bin/CPC2.py -i new_transcript.fa -o cpc2.predict.result  
python get_final_gtf.py -i cpc2.predict.result -o . -r_gtf Percocypris_pingi.genomic.gtf -n_gtf new_transcript.gtf -l gene_transcript.list
agat_convert_sp_gff2gtf.pl --gff final.gtf -o final.tmp1
agat_convert_sp_gxf2gxf.pl -g final.tmp1 -o final.tmp2
python gff_addCDS.py -i final.tmp2 -o final.gff3
rm final.tmp*



