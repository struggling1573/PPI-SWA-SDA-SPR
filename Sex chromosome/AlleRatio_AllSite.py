import argparse
import pysam
from collections import defaultdict

def read_samples(group_file):
    with open(group_file, 'r') as f:
        return [line.strip() for line in f]

def count_alleles(vcf, samples_group1, samples_group2, outfile, max_missing_ratio):
    vcf_in = pysam.VariantFile(vcf)
    num_samples_group1 = len(samples_group1)
    num_samples_group2 = len(samples_group2)
    total_samples = num_samples_group1 + num_samples_group2

    with open(outfile, 'w') as outf:
        outf.write("chr\tpos")
        for i in range(10):
            outf.write(f"\tgroup1_{i}")
        for i in range(10):
            outf.write(f"\tgroup2_{i}")
        outf.write("\n")

        for record in vcf_in.fetch():
            num_alleles = max([len(gt) for gt in record.samples.values() if gt is not None])
            counts_group1 = defaultdict(int)
            counts_group2 = defaultdict(int)
            total_group1 = 0
            total_group2 = 0
            num_missing_group1 = 0
            num_missing_group2 = 0

            for sample in record.samples:
                gt = record.samples[sample].get('GT')

                if gt is None or gt == (None, None):
                    if sample in samples_group1:
                        num_missing_group1 += 1
                    elif sample in samples_group2:
                        num_missing_group2 += 1
                    continue

                if sample in samples_group1:
                    total_group1 += len([allele for allele in gt if allele is not None])
                    for allele in gt:
                        if allele is not None:
                            counts_group1[allele] += 1
                elif sample in samples_group2:
                    total_group2 += len([allele for allele in gt if allele is not None])
                    for allele in gt:
                        if allele is not None:
                            counts_group2[allele] += 1

            if (num_missing_group1 / num_samples_group1) > (1 - max_missing_ratio) or \
               (num_missing_group2 / num_samples_group2) > (1 - max_missing_ratio):
                continue

            freq_group1 = {allele: round(count / total_group1, 4) if total_group1 > 0 else 0 for allele, count in counts_group1.items()}
            freq_group2 = {allele: round(count / total_group2, 4) if total_group2 > 0 else 0 for allele, count in counts_group2.items()}
            outf.write("{}\t{}".format(record.chrom, record.pos))
            for i in range(num_alleles):
                outf.write(f"\t{freq_group1.get(i, 0)}")
            for i in range(num_alleles):
                outf.write(f"\t{freq_group2.get(i, 0)}")
            outf.write("\n")

def parse_args():
    author_info = """
Author: Hangyu Wang
Date: Nov 25 2024
Unit: Southwest University
Contact: wanghyx666@163.com
"""
    parser = argparse.ArgumentParser(
        description="Calculate allele frequencies for two groups in a VCF file, adapting to various ploidies.",
        epilog=author_info,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True, help="Path to the input VCF file.")
    parser.add_argument("--group1", type=str, required=True, help="File containing sample names for group1.")
    parser.add_argument("--group2", type=str, required=True, help="File containing sample names for group2.")
    parser.add_argument("--out", type=str, required=True, help="Path to the output file.")
    parser.add_argument("--misratio", type=float, default=0.75, help="Maximum allowed missing ratio (default: 0.75).")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    samples_group1 = read_samples(args.group1)
    samples_group2 = read_samples(args.group2)
    count_alleles(args.vcf, samples_group1, samples_group2, args.out, args.misratio)
  
