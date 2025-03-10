vcftools --vcf snps_filtered.vcf --weir-fst-pop id_female.txt --weir-fst-pop id_male.txt --out P_fema_male --fst-window-size 100000 --fst-window-step 10000
mv P_fema_male.windowed.weir.fst P_fema_male.win100k10k.weir.fst

#less P_fema_male.win100k10k.weir.fst | sed 's/chr//g' | sed 's/hap2_//g'| sed '/F/d' | sed '/superscaf/d' > P_fema_male.win100k10k_plot.txt
sed "1d" P_fema_male.win100k10k.weir.fst | awk '{if($5<0)print $1"\t"$2"\t0";else print $1"\t"$2"\t"$5}' | sed 's/chr//g' | sed '/F/d' | sed '/superscaf/d' > P_fema_male.win100k10k_plot.txt
Rscript manhattan_SNPs.r P_fema_male_plot.txt P_fema_male

awk '{print $1"\trs"NR"\t"$2"\t"$3}' P_fema_male.win100k10k_plot.txt | sed '1i CHR\tSNP\tBP\tP' > df_fst.xls
