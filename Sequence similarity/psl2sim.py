import sys
import math
from collections import defaultdict
import gzip
import argparse

def psl_calc_milli_bad(size_mul, q_end, q_start, t_end, t_start, q_num_insert, t_num_insert, matches, rep_matches, mis_matches, is_mrna):
    q_ali_size = q_end - q_start
    t_ali_size = t_end - t_start
    ali_size = min(q_ali_size, t_ali_size)
    if ali_size <= 0:
        return (0, 0)

    size_dif = q_ali_size - t_ali_size
    if size_dif < 0:
        if is_mrna:
            size_dif = 0
        else:
            size_dif = -size_dif

    insert_factor = q_num_insert
    if not is_mrna:
        insert_factor += t_num_insert

    total = matches + rep_matches + mis_matches
    if total == 0:
        return (0, 0)

    round_away_from_zero = 3 * math.log(1 + size_dif)
    rounded = round(round_away_from_zero)

    mis2 = mis_matches + insert_factor + rounded
    total2 = total
    return (mis2, total2)

def process_file(input_file, win_size, output_file, use_query=True):
    mis_all = defaultdict(int)
    total_all = defaultdict(int)
    chrom_all = defaultdict(str)

    opener = gzip.open if input_file.endswith('.gz') else open
    mode = 'rt' if input_file.endswith('.gz') else 'r'

    with opener(input_file, mode) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 17:
                continue

            matches = int(fields[0])
            mis_matches = int(fields[1])
            rep_matches = int(fields[2])
            q_num_insert = int(fields[4])
            t_num_insert = int(fields[6])
            q_start = int(fields[11])
            q_end = int(fields[12])
            t_start = int(fields[15])
            t_end = int(fields[16])

            if use_query:
                current_start = q_start
                chrom = fields[9]  # qName
            else:
                current_start = t_start
                chrom = fields[13]  # tName

            mis2, total2 = psl_calc_milli_bad(
                1, q_end, q_start, t_end, t_start,
                q_num_insert, t_num_insert,
                matches, rep_matches, mis_matches, 1
            )

            w = current_start // win_size
            if not chrom_all[w]:
                chrom_all[w] = chrom

            mis_all[w] += mis2
            total_all[w] += total2

    with open(output_file, 'w') as f_out:
        for w in sorted(mis_all.keys()):
            start = w * win_size
            chrom = chrom_all.get(w, '')
            mis = mis_all[w]
            total = total_all[w]

            if total == 0:
                percent = 100.0
            else:
                percent = 100.0 - (mis / total) * 100.0

            f_out.write(f"{w}\t{start}\t{chrom}\t{mis}\t{total}\t{percent:.2f}\n")

def main():
    author_info = """
Author: Haoyu Wang
Date: Feb 03 2024
Affiliation: Southwest University
Contact: wanghyx666@163.com
"""
    description = 'Calculate sequence similarity in sliding windows for PSL alignments.'
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=author_info
    )
    parser.add_argument('input_file', help='Input PSL file (supports .gz compression)')
    parser.add_argument('window_size', type=int, help='Window size in base pairs')
    parser.add_argument('output_query', help='Output file using query coordinates (columns 14-15)')
    parser.add_argument('output_ref', help='Output file using reference coordinates (columns 10-11)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    process_file(args.input_file, args.window_size, args.output_query, use_query=True)
    process_file(args.input_file, args.window_size, args.output_ref, use_query=False)

if __name__ == "__main__":
    main()