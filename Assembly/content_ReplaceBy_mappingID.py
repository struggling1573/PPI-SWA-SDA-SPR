import argparse

def replace_columns_with_mapping(mapping_file, data_file, output_file):
    with open(mapping_file, 'r') as mf:
        mapping = {}
        for line in mf:
            original_id, new_id = line.strip().split()
            mapping[original_id] = new_id

    with open(output_file, 'w') as of:
        with open(data_file, 'r') as df:
            for line in df:
                fields = line.strip().split()
                replaced_fields = [mapping.get(field, field) for field in fields]
                of.write('\t'.join(replaced_fields) + '\n') 


if __name__ == "__main__":
    author_info = """
Author: Haoyu Wang
Date: Aug 24  2023
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script replaces columns in a data file based on a mapping file. It reads the mapping file to create a mapping dictionary, then replaces each field in the data file with its corresponding new ID from the mapping if available. If no mapping is found for a field, the original field is kept unchanged.",
        epilog=author_info,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('mapping_file', type=str, help='Path to the mapping file. Each line in this file should contain two fields separated by tabs or spaces: the original ID and the new ID.')
    parser.add_argument('input_file', type=str, help='Path to the input data file. Fields in this file should be separated by tabs or spaces.')
    parser.add_argument('output_file', type=str, help='Path to the output file where the data with replaced columns will be written.')
    args = parser.parse_args()

    replace_columns_with_mapping(args.mapping_file, args.input_file, args.output_file)
