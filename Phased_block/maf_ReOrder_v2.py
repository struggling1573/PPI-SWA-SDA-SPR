import sys
import argparse


def main():
    author_info = """
Author: Hangyu Wang
Date: Dec 28 2024
Unit: Southwest University
Contact: wanghyx666@163.com
"""
    epilog = f"""
Usage example:
    python script.py input.maf order1,order2,...
{author_info}
"""
    parser = argparse.ArgumentParser(
        description="This tool is used to sort the blocks in a MAF file according to the specified order. It reads an input MAF file and a sorting order string, then sorts the 's' lines within each block of the MAF file based on the given order.",
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input_file", help="Path to the input MAF file.")
    parser.add_argument("order_str", help="Sorting order, separated by commas, e.g., order1,order2,...")

    args = parser.parse_args()

    input_file = args.input_file
    order_str = args.order_str
    order_list = order_str.split(',')

    with open(input_file, 'r') as f:
        header = []
        line = f.readline()
        while line and line.startswith('#'):
            header.append(line)
            line = f.readline()
        print(''.join(header), end='')

        blocks = []
        current_block = []

        while line:
            stripped = line.strip()
            if not stripped: 
                line = f.readline()
                continue

            if stripped.startswith('a'):
                if current_block:
                    blocks.append(current_block)
                current_block = [stripped]
            elif stripped.startswith('s'):
                if current_block:
                    current_block.append(stripped)
            else:  
                if current_block:
                    current_block.append(stripped)
                    print(line, end='')
            line = f.readline()

        if current_block:
            blocks.append(current_block)
        for i, block in enumerate(blocks):
            a_line = block[0]
            other_lines = block[1:]
            s_lines = []
            non_s_lines = []
            for line in other_lines:
                if line.startswith('s'):
                    s_lines.append(line)
                else:
                    non_s_lines.append(line)
            try:
                sorted_s = sorted(
                    s_lines,
                    key=lambda x: order_list.index(x.split()[1].split('.')[0])
                )
            except ValueError as e:
                print(f"Error: An unspecified prefix was found (original line: {line})", file=sys.stderr)
                print(f"Current order list: {order_list}", file=sys.stderr)
                sys.exit(1)
            new_block = [a_line] + sorted_s + non_s_lines
            if i == 0:
                print('\n'.join(new_block))
            else:
                print('\n' + '\n'.join(new_block))
        print()

if __name__ == '__main__':
    main()
