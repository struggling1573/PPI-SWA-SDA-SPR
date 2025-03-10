import argparse

def parse_args():
    author_info = """
Author: Hangyu Wang
Date: Nov 25 2024
Unit: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script defines two genotype frequency intervals (Interval A and Interval B) and extracts rows from the input file where the genotype frequencies of group1 and group2 meet specific conditions. Specifically, for at least one of the four genotype frequencies (A, T, G, C), the frequency in group1 should be greater than that in group2, and the frequency in group1 should fall within Interval A while the frequency in group2 should fall within Interval B.",
        epilog=author_info,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", type=str, required=True, help="Path to the input file which contains genotype frequency data.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output file where the filtered rows will be saved.")
    return parser.parse_args()


def filter_rows(input_file, output_file):
    interval_A = (0.9, 1.0) 
    interval_B = (0.7, 0.8) 

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = next(infile) 
        outfile.write(header) 
        for line in infile:
            parts = line.strip().split('\t')
            group1_values = parts[2:6]  
            group2_values = parts[6:10] 
            keep_row = False  
            for i in range(4):
                group1_val = float(group1_values[i])
                group2_val = float(group2_values[i])
                if group1_val > group2_val and interval_A[0] < group1_val < interval_A[1] and interval_B[0] < group2_val < interval_B[1]:
                    keep_row = True
                    break  
            if keep_row:
                outfile.write(line)

if __name__ == "__main__":
    args = parse_args()
    filter_rows(args.input, args.output)

