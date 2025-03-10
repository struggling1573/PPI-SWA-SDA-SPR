import sys
import re
import gzip
import argparse


def parse_vcf_and_write_output(vcf_file, num_file_prefix, ratio_file_prefix, sample_name=None):
    if vcf_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt' 
    else:
        open_func = open
        mode = 'r'

    chromosomes_data = {}
    with open_func(vcf_file, mode) as vcf:
        sample_header = None
        for line in vcf:
            if line.startswith('#CHROM'):
                sample_header = line.strip().split('\t')[9:]
                break
        if sample_header is None:
            print("Error: No sample header found in VCF file.")
            return
        for line in vcf:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom = parts[0]
                if sample_name:
                    try:
                        sample_index = sample_header.index(sample_name)
                        genotype = parts[sample_index + 9]
                    except ValueError:
                        print(f"Error: Sample name '{sample_name}' not found in the sample columns.")
                        return
                else:
                    genotype = parts[9]
                main_genotype = re.match(r'([^\:]+)', genotype).group(1)
                if chrom not in chromosomes_data:
                    chromosomes_data[chrom] = {}
                if main_genotype not in chromosomes_data[chrom]:
                    chromosomes_data[chrom][main_genotype] = 0
                chromosomes_data[chrom][main_genotype] += 1

    num_file = f"{sample_name}.{num_file_prefix}" if sample_name else num_file_prefix
    ratio_file = f"{sample_name}.{ratio_file_prefix}" if sample_name else ratio_file_prefix
    with open(num_file, 'w') as num_out, open(ratio_file, 'w') as ratio_out:
        for chrom, genotypes in chromosomes_data.items():
            num_out.write(f"{chrom}")
            ratio_out.write(f"{chrom}")
            total = sum(genotypes.values())
            for genotype, count in genotypes.items():
                num_out.write(f"\t{genotype}\t{count}")
                ratio_out.write(f"\t{genotype}\t{count / total if total > 0 else 0}")
            num_out.write('\n')
            ratio_out.write('\n')


def main():
    author_info = """
Author: Hangyu Wang
Date: Nov 25 2024
Unit: Southwest University
Contact: wanghyx666@163.com
"""
    epilog = f"""
Usage example:
    python geneAlle_byCHR_v2.py input.vcf[.gz] [sampleID]
{author_info}
"""
    parser = argparse.ArgumentParser(
        description="This script is used to calculate the genotype frequencies on each chromosome in a VCF file. It reads a VCF file (compressed or uncompressed), optionally filters by a specified sample name, and then outputs the genotype counts and ratios for each chromosome to separate files.",
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input_file", help="Path to the input VCF file, which can be compressed (.gz) or uncompressed.")
    parser.add_argument("sample_name", nargs='?', default=None, help="Optional parameter. Specify the name of the sample to be analyzed.")
    args = parser.parse_args()
    vcf_filename = args.input_file
    gene_num_file = "gene_num.txt"
    gene_ratio_file = "gene_ratio.txt"
    sample_name = args.sample_name

    parse_vcf_and_write_output(vcf_filename, gene_num_file, gene_ratio_file, sample_name)

if __name__ == "__main__":
    main()
