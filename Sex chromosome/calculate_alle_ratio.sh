vcftools --vcf snps_filtered.vcf --keep id_female.txt --freq --out female_allele_frequencies
vcftools --vcf snps_filtered.vcf --keep id_male.txt --freq --out male_allele_frequencies

paste female_allele_frequencies.frq male_allele_frequencies.frq | cut -f 1,2,5,11 | sed 1d | sed '/nan/d' | sed 's/:/\t/g' | cut -f 1,2,4,6 | awk '$3!="" && $4!=""' | sed '1i CHR\tPOS\tfemale\tmale' > df.txt
