ls *_temp_gtf > mergelist.txt
stringtie --merge -p 4 -o merge.gtf mergelist.txt

mkdir merge predict
sed 's/gene_name ".*"//g' Percocypris_pingi.genomic.gtf > ref_new.gtf
cuffcompare -r ref_new.gtf -e 100 -C -s Percocypris_pingi.genomic.fna -o compared.stringtie.gtf merge.gtf

awk '{if($3 ~ /^(u|i|o|j)$/) print $1"\t"$2"\t"$3"\t"$4"\t"$5}' compared.stringtie.gtf.merge.gtf.tmap > gene_transcript.list
python get_new_transcript.py -i merge.gtf -fa Percocypris_pingi.genomic.fna -l gene_transcript.list -o new_transcript
python /biosoft/CPC2_standalone-1.0.1/bin/CPC2.py -i new_transcript.fa -o cpc2.predict.result  
python get_final_gtf.py -i cpc2.predict.result -o . -r_gtf Percocypris_pingi.genomic.gtf -n_gtf new_transcript.gtf -l gene_transcript.list
