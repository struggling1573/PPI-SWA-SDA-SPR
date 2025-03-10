import pysam
import pandas as pd
import argparse

def count_genotypes_in_windows(vcf_file, genotype_file, sample_name, output_file, window_size=50000):
    result = []
    vcf = pysam.VariantFile(vcf_file)
    with open(genotype_file, 'r') as f:
        genotypes_to_count = [tuple(map(int, line.strip().split('/'))) for line in f.readlines()]
    for contig in vcf.header.contigs:
        chr_name = contig
        contig_length = vcf.header.contigs[contig].length
        start_pos = 0
        while start_pos < contig_length:
            end_pos = start_pos + window_size
            if end_pos > contig_length:
                end_pos = contig_length
            genotype_counts = {genotype: 0 for genotype in genotypes_to_count}
            for record in vcf.fetch(chr_name, start_pos, end_pos):
                if sample_name in record.samples:
                    genotype = record.samples[sample_name]['GT']
                    if genotype in genotype_counts:
                        genotype_counts[genotype] += 1
            result.append([chr_name, start_pos, end_pos] + [genotype_counts[genotype] for genotype in genotypes_to_count])
            start_pos += window_size
    vcf.close()
    output_df = pd.DataFrame(result, columns=['chr', 'start pos', 'end pos'] + [str(genotype) for genotype in genotypes_to_count])
    output_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    author_info = """
Author: Hangyu Wang
Date: Nov 25 2024
Unit: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script counts the frequencies of specified genotypes in sliding windows of a VCF file. It reads a VCF file, a list of genotypes to count, and a sample name, then counts the occurrences of each genotype within each window and outputs the results to a tab - separated file.",
        epilog=author_info,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file.')
    parser.add_argument('--genetype', type=str, required=True, help='Path to the file containing genotypes to count. Each line of this file should represent a genotype in the format of integers separated by a slash, e.g., 0/1.')
    parser.add_argument('--sample', type=str, required=True, help='Name of the sample in the VCF file to count genotypes for.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file. The output will be a tab - separated file.')
    parser.add_argument('--window_size', type=int, default=50000, help='Size of the window in base pairs. The default value is 50000.')

    args = parser.parse_args()

    count_genotypes_in_windows(args.vcf, args.genetype, args.sample, args.output, args.window_size)
