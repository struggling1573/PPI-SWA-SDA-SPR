import argparse
import pysam

def parse_args():
    author_info = """
Author: Hangyu Wang
Date: Nov 25 2024
Unit: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script filters a VCF file according to the genotypes of male and female samples. It keeps the records where all male samples have genotypes (0, 1) or (1, 0), and all female samples have the genotype (0, 0).",
        epilog=author_info,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True, help="Path to the input VCF file.")
    parser.add_argument("--male", type=str, required=True, help="File containing male sample names, one name per line.")
    parser.add_argument("--female", type=str, required=True, help="File containing female sample names, one name per line.")
    parser.add_argument("--out", type=str, required=True, help="Path to the output VCF file.")
    return parser.parse_args()

def read_samples(group_file):
    with open(group_file, 'r') as f:
        return [line.strip() for line in f]

def filter_vcf(vcf_file, male_samples, female_samples, out_file):
    vcf_in = pysam.VariantFile(vcf_file)
    vcf_out = pysam.VariantFile(out_file, 'w', header=vcf_in.header)
    for record in vcf_in.fetch():
        male_condition_met = True
        female_condition_met = True
        for sample in male_samples:
            gt = record.samples[sample].get('GT')
            if gt not in [(0, 1), (1, 0)]:
                male_condition_met = False
                break
        for sample in female_samples:
            gt = record.samples[sample].get('GT')
            if gt not in [(0, 0)]:
                female_condition_met = False
                break
        if male_condition_met and female_condition_met:
            vcf_out.write(record)
    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    args = parse_args()
    male_samples = read_samples(args.male)
    female_samples = read_samples(args.female)
    filter_vcf(args.vcf, male_samples, female_samples, args.out)

