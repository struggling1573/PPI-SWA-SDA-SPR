import argparse
import os
import re
import glob


def parse_busco_file(file_path):
    filename = os.path.basename(file_path)
    if 'result_' not in filename:
        return None
    species = filename.split('result_')[1].split('.')[0]
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    result_line = None
    for line in lines:
        if line.startswith('C:'):
            result_line = line
            break
    if not result_line:
        print(f"Warning: Result line not found in {file_path}")
        return None
    c_match = re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\]', result_line)
    f_match = re.search(r'F:([\d.]+)%', result_line)
    m_match = re.search(r'M:([\d.]+)%', result_line)
    n_match = re.search(r'n:(\d+)', result_line)

    if not all([c_match, f_match, m_match, n_match]):
        print(f"Warning: Failed to parse percentages in {file_path}")
        return None
    try:
        percent_C = "{:.2f}%".format(float(c_match.group(1)))
        percent_S = "{:.2f}%".format(float(c_match.group(2)))
        percent_D = "{:.2f}%".format(float(c_match.group(3)))
        percent_F = "{:.2f}%".format(float(f_match.group(1)))
        percent_M = "{:.2f}%".format(float(m_match.group(1)))
        count_total = int(n_match.group(1))
    except ValueError:
        print(f"Warning: Invalid numeric format in {file_path}")
        return None
        
    counts = {'C': None, 'S': None, 'D': None, 'F': None, 'M': None}
    patterns = {
        'C': r'^\s*(\d+)\s+Complete BUSCOs \(C\)',
        'S': r'^\s*(\d+)\s+Complete and single-copy BUSCOs \(S\)',
        'D': r'^\s*(\d+)\s+Complete and duplicated BUSCOs \(D\)',
        'F': r'^\s*(\d+)\s+Fragmented BUSCOs \(F\)',
        'M': r'^\s*(\d+)\s+Missing BUSCOs \(M\)'
    }

    for line in lines:
        for key in patterns:
            if counts[key] is None:
                match = re.match(patterns[key], line)
                if match:
                    try:
                        counts[key] = int(match.group(1))
                    except ValueError:
                        pass

    if None in counts.values():
        print(f"Warning: Missing counts in {file_path}")
        return None

    return {
        'species': species,
        'percent_C': percent_C,
        'percent_S': percent_S,
        'percent_D': percent_D,
        'percent_F': percent_F,
        'percent_M': percent_M,
        'count_C': counts['C'],
        'count_S': counts['S'],
        'count_D': counts['D'],
        'count_F': counts['F'],
        'count_M': counts['M'],
        'count_total': count_total
    }


def main():
    author_info = """
Author: Haoyu Wang
Date: Dec 30  2024
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description='Parse BUSCO results into a summary table',
        epilog=author_info
    )
    parser.add_argument('-datapath', required=True, help='Directory containing BUSCO result files')
    parser.add_argument('-output', default='busco_summary.tsv', help='Output TSV file name')
    args = parser.parse_args()

    file_pattern = os.path.join(args.datapath, 'short_summary.specific.actinopterygii_odb10.result_*.txt')
    files = glob.glob(file_pattern)

    results = []
    for file in files:
        data = parse_busco_file(file)
        if data:
            results.append(data)

    if not results:
        print("No valid data to write.")
        return

    headers = [
        'species', 'percent_C', 'percent_S', 'percent_D', 'percent_F', 'percent_M',
        'count_C', 'count_S', 'count_D', 'count_F', 'count_M', 'count_total'
    ]

    with open(args.output, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for res in results:
            row = [
                res['species'],
                res['percent_C'],
                res['percent_S'],
                res['percent_D'],
                res['percent_F'],
                res['percent_M'],
                str(res['count_C']),
                str(res['count_S']),
                str(res['count_D']),
                str(res['count_F']),
                str(res['count_M']),
                str(res['count_total'])
            ]
            f.write('\t'.join(row) + '\n')

    print(f"Successfully wrote {len(results)} records to {args.output}")


if __name__ == '__main__':
    main()
    
