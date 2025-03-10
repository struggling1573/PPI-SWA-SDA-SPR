paste fe_pool.win50k.stat ma_pool.win50k.stat | cut -f 1,2,3,7,15 | sed '/#/d' | awk '{print $1"\t"$2"\t"$3"\t"$4/$5}' | sed '1i Chr\tStart\tEnd\tCoverage' > df_Coverage.txt
paste fe_pool.win50k.stat ma_pool.win50k.stat | cut -f 1,2,3,8,16 | sed '/#/d' | awk '{print $1"\t"$2"\t"$3"\t"$4/$5}' | sed '1i Chr\tStart\tEnd\tMeanDepth' > df_MeanDepth.txt

cat df_Coverage.txt | sed 1d | sed 's/hap2_chr0//g' | sed 's/hap2_chr//g' | awk '{print "rs"NR"\t"$1"\t"$3"\t"$4}' | sed '1i SNP\tCHR\tBP\tP' > df_Coverage_trans.txt
cat df_MeanDepth.txt | sed 1d | sed 's/hap2_chr0//g' | sed 's/hap2_chr//g' | awk '{print "rs"NR"\t"$1"\t"$3"\t"$4}' | sed '1i SNP\tCHR\tBP\tP' > df_MeanDepth_trans.txt
