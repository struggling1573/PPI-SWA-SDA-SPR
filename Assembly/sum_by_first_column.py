import argparse

def sum_by_first_column(input_file, output_file):
    current_key = None
    current_sum = 0
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split('\t')
            key = parts[0]
            value = int(parts[1])
            if key != current_key:
                if current_key is not None:
                    f_out.write(f"{current_key}\t{current_sum}\n")
                current_key = key
                current_sum = value
            else:
                current_sum += value
        if current_key is not None:
            f_out.write(f"{current_key}\t{current_sum}\n")


if __name__ == "__main__":
    author_info = """
Author: Haoyu Wang
Date: Aug 13  2023
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script reads a tab-separated input file, sums the values in the second column grouped by the values in the first column, and writes the results to an output file. Each line of the input file should contain two columns separated by a tab: the first column is a key, and the second column is an integer value.",
        epilog=author_info,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('input_file', type=str, help='Path to the input file. The file should be a tab-separated text file with two columns: the first column is a key, and the second column is an integer value.')
    parser.add_argument('output_file', type=str, help='Path to the output file where the results will be written.')
    args = parser.parse_args()

    sum_by_first_column(args.input_file, args.output_file)
