import argparse
import pandas as pd


def replace_column(input_file, output_file, mapping_file, column_index):
    mapping_df = pd.read_csv(mapping_file, sep='\t', header=0)
    mapping = dict(zip(mapping_df.iloc[:, 0], mapping_df.iloc[:, 1]))
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.iloc[:, column_index - 1] = df.iloc[:, column_index - 1].astype(str)

    def replace_ids(x):
        for original_id, renamed_id in mapping.items():
            x = x.replace(original_id, renamed_id)
        return x

    df.iloc[:, column_index - 1] = df.iloc[:, column_index - 1].apply(replace_ids)
    df.to_csv(output_file, sep='\t', header=None, index=False)


def main():
    author_info = """
Author: Haoyu Wang
Date: May 23  2024
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="Replace the specified column in the input file according to the correspondence file. If there is no matching item in the correspondence file, the information in the original file remains unchanged.",
        epilog=author_info,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Path to the input file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file')
    parser.add_argument('-m', '--mapping', required=True, help='Path to the correspondence file')
    parser.add_argument('-c', '--column', type=int, required=True, help='Index of the column to be replaced (starting from 1)')
    args = parser.parse_args()

    replace_column(args.input, args.output, args.mapping, args.column)


if __name__ == "__main__":
    main()