#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wanghyx666@163.com
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

cpu=10
REF=ref.fna

bwa index $REF
samtools faidx $REF
export PATH=/biosoft/gatk-4.2.6.1:$PATH
ref=`echo $REF | sed 's/.fna//g'`
gatk CreateSequenceDictionary -R ${ref}.fna -O ${ref}.dict

[ -d 01.bwa_out/ ] || mkdir -p 01.bwa_out
[ -d 02.markDup/ ] || mkdir -p 02.markDup

for sample in `cat id`;do
        bwa mem -t $cpu -R $(echo "@RG\tID:$sample\tSM:$sample\tLB:$sample"_"$sample\tPL:ILLUMINA") $REF reseq/${sample}_R1.clean.fastq.gz reseq/${sample}_R2.clean.fastq.gz  |  samtools sort -@ $cpu -O BAM -o 01.bwa_out/${sample}.sorted.bam  -
        pandepth -i 01.bwa_out/${sample}.sorted.bam -o ${sample} -t $cpu
        pandepth -i 01.bwa_out/${sample}.sorted.bam -o ${sample} -t $cpu -w 50000
        mv ${sample}.chr.stat.gz 01.bwa_out/
        mv ${sample}.win.stat.gz 01.bwa_out/${sample}.win50k.stat.gz
        picard SortSam -I 01.bwa_out/${sample}.sorted.bam -O 02.markDup/${sample}.sorted.bam -SO coordinate -VALIDATION_STRINGENCY LENIENT --TMP_DIR /tmp
        picard MarkDuplicates -I 02.markDup/${sample}.sorted.bam -O 02.markDup/${sample}.dedup.bam -M 02.markDup/${sample}.m --TMP_DIR /tmp
        samtools index 02.markDup/${sample}.dedup.bam
done
