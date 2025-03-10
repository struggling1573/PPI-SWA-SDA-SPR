import argparse
import pysam
import os

def rename_chromosomes_in_bam(mapping_file, input_bam_file, output_bam_file):
    chromosome_map = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            old_chrom, new_chrom = line.strip().split()
            chromosome_map[old_chrom] = new_chrom

    input_bam = pysam.AlignmentFile(input_bam_file, 'rb')
    output_bam_header = input_bam.header.copy()
    sq_entries = []

    for sq in output_bam_header['SQ']:
        old_chrom = sq['SN']
        new_chrom = chromosome_map.get(old_chrom, old_chrom)
        sq_entries.append({'SN': new_chrom, 'LN': sq['LN']})
    new_header_dict = output_bam_header.to_dict()
    new_header_dict['SQ'] = sq_entries
    with pysam.AlignmentFile('temp.bam', 'wb', header=pysam.AlignmentHeader.from_dict(new_header_dict)) as temp_bam:
    temp_bam = pysam.AlignmentFile('temp.bam', 'rb')
    correct_header = temp_bam.header
    temp_bam.close()
    os.remove('temp.bam')
    output_bam = pysam.AlignmentFile(output_bam_file, 'wb', header=correct_header)

    for read in input_bam:
        output_bam.write(read)
    input_bam.close()
    output_bam.close()


if __name__ == "__main__":
    author_info = """
Author: Haoyu Wang
Date: Aug 24  2023
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script renames chromosomes in a BAM file according to a mapping file. It reads the mapping file to create a chromosome name mapping dictionary, then updates the @SQ lines in the BAM header with the new chromosome names. Finally, it writes the alignment records from the input BAM file to the output BAM file with the updated header.",
        epilog=author_info,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('mapping_file', type=str, help='Path to the mapping file. Each line in this file should contain two fields separated by whitespace: the old chromosome name and the new chromosome name.')
    parser.add_argument('input_bam_file', type=str, help='Path to the input BAM file.')
    parser.add_argument('output_bam_file', type=str, help='Path to the output BAM file where the BAM file with renamed chromosomes will be written.')
    args = parser.parse_args()

    rename_chromosomes_in_bam(args.mapping_file, args.input_bam_file, args.output_bam_file)
