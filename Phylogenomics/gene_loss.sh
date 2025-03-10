# Count the number of missing genes in hap1
cat Orthogroups.GeneCount.tsv | awk '$2>0 && $2<100 && $3=="0" && $4<100 && $4>0' | cut -f 1 | grep -f - Orthogroups.txt | sed 's/ /\n/g' | grep -v 'OG\|PPI' | grep "chr" | sed 's/G.*//g' | sort | uniq -c | awk '{print $2"\t"$1}'
# Count the number of missing genes in hap2
cat Orthogroups.GeneCount.tsv | awk '$2>0 && $2<100 && $3>0 && $3<100 && $4=="0"' | cut -f 1 | grep -f - Orthogroups.txt | sed 's/ /\n/g' | grep -v 'OG\|PPI' | grep "chr" | sed 's/G.*//g' | sort | uniq -c | awk '{print $2"\t"$1}'

#cat Orthogroups.GeneCount.tsv | awk '$2>0 && $2<100 && $3=="0" && $4<100 && $4>0 && $5<100 && $5>0' | cut -f 1 | grep -f - Orthogroups.txt | sed 's/ /\n/g' | grep -v 'OG\|SWA' | grep "chr" | sed 's/G.*//g' | sort | uniq -c | awk '{print $2"\t"$1}'
#cat Orthogroups.GeneCount.tsv | awk '$2>0 && $2<100 && $3>0 && $3<100 && $4=="0" && $5=="0"' | cut -f 1 | grep -f - Orthogroups.txt | sed 's/ /\n/g' | grep -v 'OG\|SWA' | grep "chr" | sed 's/G.*//g' | sort | uniq -c | awk '{print $2"\t"$1}'
