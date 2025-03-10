import argparse

def process_file(input_file):
    current_key = None
    ma_value = None
    pa_value = None
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            key = parts[0]
            value = float(parts[1])
            tag = parts[2]
            if key != current_key:
                if current_key is not None:
                    if ma_value is not None and pa_value is None:
                        print(f"{current_key}\tMa")
                    elif ma_value is None and pa_value is not None:
                        print(f"{current_key}\tPa")
                    elif ma_value < pa_value:
                        print(f"{current_key}\tMa")
                    elif ma_value > pa_value:
                        print(f"{current_key}\tPa")
                    elif ma_value == pa_value:
                        print(f"{current_key}\tMa")
                        print(f"{current_key}\tPa")
                current_key = key
                ma_value = value if tag == 'Ma' else None
                pa_value = value if tag == 'Pa' else None
            else:
                if tag == 'Ma':
                    ma_value = value
                elif tag == 'Pa':
                    pa_value = value
    if current_key is not None:
        if ma_value is not None and pa_value is None:
            print(f"{current_key}\tMa")
        elif ma_value is None and pa_value is not None:
            print(f"{current_key}\tPa")
        elif ma_value < pa_value:
            print(f"{current_key}\tMa")
        elif ma_value > pa_value:
            print(f"{current_key}\tPa")
        elif ma_value == pa_value:
            print(f"{current_key}\tMa")
            print(f"{current_key}\tPa")


if __name__ == "__main__":
    author_info = """
Author: Haoyu Wang
Date: Aug 13  2023
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script processes an input file where each line contains a key, a numerical value, and a tag (either 'Ma' or 'Pa') separated by tabs. It groups the data by key and compares the values associated with 'Ma' and 'Pa' tags for each key. Based on the comparison results, it prints the key along with the appropriate tag ('Ma' or 'Pa') or both if the values are equal.",
        epilog=author_info,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('input_file', type=str, help='Path to the input file. The file should be a tab-separated text file where each line contains a key, a numerical value, and a tag (either "Ma" or "Pa").')
    args = parser.parse_args()

    process_file(args.input_file)
