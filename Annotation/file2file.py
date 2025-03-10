import argparse


def combine_files(file1_path, file1_column, file2_path, output_path):
    data1_dict = {}
    with open(file1_path, 'r') as file1:
        for line in file1:
            columns1 = line.strip().split('\t')
            data1_dict.setdefault(columns1[file1_column - 1], []).append(columns1)

    data2_dict = {}
    with open(file2_path, 'r') as file2:
        for line in file2:
            columns2 = line.strip().split('\t')
            for key in columns2:
                data2_dict.setdefault(key, []).append(columns2)

    with open(output_path, 'w') as output:
        for key, rows1 in data1_dict.items():
            if key in data2_dict:
                for row1 in rows1:
                    for row2 in data2_dict[key]:
                        combined_row = row1 + row2
                        output.write('\t'.join(combined_row) + '\n')
            else:
                for row1 in rows1:
                    combined_row = row1 + ['NA'] * len(data2_dict.get(next(iter(data2_dict.keys()), []), []))
                    output.write('\t'.join(combined_row) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='''
    This script is used to merge two tab-separated files. It searches for values in a specified column of the first file and looks for matches in the second file.
    When a match is found, the corresponding rows from both files are combined and written to the output file. If no match is found in the second file, "NA" is added.
    ''',
        epilog="Author: Haoyu Wang\n"
               "Date: Feb 02  2023\n"
               "Affiliation: Southwest University\n"
               "Contact: wanghyx666@163.com\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('infile1', type=str, help='Path to the first input file.')
    parser.add_argument('infile2', type=str, help='Path to the second input file.')
    parser.add_argument('file1_column', type=int, default=1, help='Column number in the first file used for querying (default value is 1).')
    parser.add_argument('out', type=str, default='output.txt', help='Path to the output file (default value is output.txt).')

    args = parser.parse_args()
    combine_files(args.infile1, args.file1_column, args.infile2, args.out)