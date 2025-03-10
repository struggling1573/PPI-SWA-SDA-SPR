#!/bin/bash

# Usage: for i in `seq 1 25`;do bash AlleCount_ploid6.sh snp_gene_num.txt $i >> tmp;done && sed -i '1i CHR\tAlleNumber\tHexasomic\tTetrasomic\tDisomic\tNoAlle\tOther' tmp && mv tmp snp_gene_num.count

input=$1
chr_name=$2

tmp1=`awk -v a=$chr_name 'NR==a' $input  | cut -f 1`
tmp2=`awk -v a=$chr_name 'NR==a' $input | sed 's/\t/\n/g' | grep -v '/' | grep -v '-' | awk '{sum += $1};END {print sum}'`
tmp3=`awk -v a=$chr_name 'NR==a' $input | sed 's/\t/\n/g' | grep -A 1 '0/0/0/0/0/1\|0/1/1/1/1/1' | grep -v '/' | grep -v '-' | awk '{sum += $1};END {print sum}'`
tmp4=`awk -v a=$chr_name 'NR==a' $input | sed 's/\t/\n/g' | grep -A 1 '0/0/0/0/1/1\|0/0/1/1/1/1' | grep -v '/' | grep -v '-' | awk '{sum += $1};END {print sum}'`
tmp5=`awk -v a=$chr_name 'NR==a' $input | sed 's/\t/\n/g' | grep -A 1 '0/0/0/1/1/1\|1/1/1/0/0/0' | grep -v '/' | grep -v '-' | awk '{sum += $1};END {print sum}'`
tmp6=`awk -v a=$chr_name 'NR==a' $input | sed 's/\t/\n/g' | grep -A 1 '0/0/0/0/0/0\|1/1/1/1/1/1' | grep -v '/' | grep -v '-' | awk '{sum += $1};END {print sum}'`
tmp7=`expr $tmp2 - $tmp3 - $tmp4 - $tmp5 - $tmp6`
echo -e "$tmp1\t$tmp2\t$tmp3\t$tmp4\t$tmp5\t$tmp6\t$tmp7"
