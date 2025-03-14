import sys
import glob
import pandas as pd
import argparse

def main(compare_txt, file_suffix, output_file):
    group = {}
    group_key = []
    all_samples = []
    processed_samples = set()
    with open(compare_txt) as f:
        for line in f:
            tmp = line.strip().split(':')
            g1, g2 = tmp[0].split('_vs_')
            group[g1] = tmp[1].split(',')
            group[g2] = tmp[2].split(',')
            if g1 not in group_key:
                group_key.append(g1)
            if g2 not in group_key:
                group_key.append(g2)
            all_samples.extend(tmp[1].split(','))
            all_samples.extend(tmp[2].split(','))

    order_df = []
    all_df = glob.glob(f'*_{file_suffix}')
    for sample in all_samples:
        if sample in processed_samples:
            continue
        file_path = f'{sample}_{file_suffix}'
        if file_path in all_df:
            try:
                df = pd.read_csv(file_path, sep='\t')
                df = df[['gene_id'] + [col for col in df.columns if col != 'gene_id']]
                df = df.rename(columns={col: sample if col != 'gene_id' else 'gene_id' for col in df.columns})
                order_df.append(df)
                print(f"Successfully read file {file_path}")
                processed_samples.add(sample)
            except Exception as e:
                print(f"Error reading file {file_path}: {e}")
        else:
            print(f"File {file_path} not found.")

    if not order_df:
        print("No valid data files found. Please check your input.")
        return

    new_df = order_df[0]
    for n in range(1, len(order_df)):
        new_df = pd.merge(new_df, order_df[n], on='gene_id', how='outer')

    columns = ['gene_id'] + all_samples
    valid_columns = [col for col in columns if col in new_df.columns]
    new_df = new_df[valid_columns]

    new_df.to_csv(output_file, index=0, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge gene expression files based on a comparison group file.",
                                     epilog="Author: Haoyu Wang\n"
                                            "Date: Dec 30 2024\n"
                                            "Affiliation: Southwest University\n"
                                            "Contact: wanghyx666@163.com\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("compare_txt", help="Path to the compare group text file (e.g., compare_group.txt).")
    parser.add_argument("file_suffix", choices=['generc', 'genefpkm', 'genetpm'],
                        help="Suffix of the gene expression files. Must be one of generc, genefpkm or genetpm.")
    parser.add_argument("output_file", help="Output file name for the merged gene expression matrix.")
    args = parser.parse_args()

    main(args.compare_txt, args.file_suffix, args.output_file)
