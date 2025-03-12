#!/usr/bin/env python3
import argparse
import re
import gzip
import sys
from collections import defaultdict

def determine_file_type(filename):
    if filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(filename, 'rt') as f:
        first_char = f.read(1)
    if first_char == '>':
        return 'fa'
    elif first_char == '@':
        return 'fq'
    else:
        raise ValueError(f"Unknown file format for {filename}")

def read_fasta(file_obj):
    seq_id = None
    seq = []
    for line in file_obj:
        line = line.strip()
        if line.startswith('>'):
            if seq_id is not None:
                yield (seq_id, ''.join(seq).upper())
            seq_id = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if seq_id is not None:
        yield (seq_id, ''.join(seq).upper())

def read_fastq(file_obj):
    while True:
        header = file_obj.readline().strip()
        if not header:
            break
        seq = file_obj.readline().strip().upper()
        file_obj.readline()  # +
        file_obj.readline()  # quality
        yield (header[1:], seq)

def process_sequence(seq, stats, scaf_cutoff, contig_cutoff, check_n):
    seq_len = len(seq)

    # Scaffold statistics
    if seq_len >= scaf_cutoff:
        stats['scaffold']['lengths'].append(seq_len)
        stats['scaffold']['total_length'] += seq_len
        stats['scaffold']['gc_count'] += seq.count('G') + seq.count('C')
        stats['scaffold']['count'] += 1
        if seq_len >= 2000:
            stats['scaffold']['count_2k'] += 1

        # Split into contigs
        contigs = re.split('N+', seq)
        for contig in contigs:
            contig_len = len(contig)
            if contig_len >= contig_cutoff:
                stats['contig']['lengths'].append(contig_len)
                stats['contig']['total_length'] += contig_len
                stats['contig']['gc_count'] += contig.count('G') + contig.count('C')
                stats['contig']['count'] += 1
                if contig_len >= 2000:
                    stats['contig']['count_2k'] += 1

    # N statistics
    if check_n:
        n_segments = re.findall('N+', seq)
        for n in n_segments:
            n_len = len(n)
            if stats['N']['min'] is None or n_len < stats['N']['min']:
                stats['N']['min'] = n_len
            if stats['N']['max'] is None or n_len > stats['N']['max']:
                stats['N']['max'] = n_len
            stats['N']['total'] += n_len

def calculate_n50(sorted_lengths, total_length):
    n_stats = {}
    if not sorted_lengths or total_length == 0:
        return n_stats

    cumulative = 0
    for idx, length in enumerate(sorted_lengths, 1):
        cumulative += length
        for n in range(1, 10):
            threshold = total_length * n / 10
            if cumulative >= threshold and str(n) not in n_stats:
                n_stats[str(n)] = {'length': length, 'number': idx}
    return n_stats

def format_number(num, thousand_sep=True):
    if thousand_sep:
        return "{:,}".format(num)
    return str(num)

def main():
    author_info = """
Author: Hangyu Wang
Date: Dec 08 2023
Unit: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(description='Genome assembly statistics calculator', epilog=author_info, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_files', nargs='*', help='Input FASTA/Q files')
    parser.add_argument('--inlist', help='File containing list of input files')
    parser.add_argument('--scacut', type=int, default=100, help='Scaffold length cutoff (default: 100)')
    parser.add_argument('--cutoff', type=int, help='Contig length cutoff (default: fa=100, fq=1)')
    parser.add_argument('--Ncheck', choices=['y', 'n'], default='y', help='Check N statistics (default: y)')
    parser.add_argument('--showGC', choices=['y', 'n'], default='y', help='Show GC content (default: y)')
    parser.add_argument('--thousp', choices=['y', 'n'], default='y', help='Use thousand separators (default: y)')
    parser.add_argument('--justiC', choices=['y', 'n'], default='y', help='Column justification (default: y)')
    parser.add_argument('--strout', help='Output file for tab-separated stats')
    parser.add_argument('--sample', help='Sample name for output')
    args = parser.parse_args()

    # Process input files
    input_files = args.input_files
    if args.inlist:
        with open(args.inlist) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    input_files.append(line)

    if not input_files:
        print("Error: No input files specified", file=sys.stderr)
        sys.exit(1)

    # Initialize statistics
    stats = {
        'scaffold': {
            'lengths': [],
            'total_length': 0,
            'gc_count': 0,
            'count': 0,
            'count_2k': 0,
            'max_length': 0,
        },
        'contig': {
            'lengths': [],
            'total_length': 0,
            'gc_count': 0,
            'count': 0,
            'count_2k': 0,
            'max_length': 0,
        },
        'N': {
            'total': 0,
            'min': None,
            'max': None,
        },
        'total_sequences': 0,
    }

    # Process each file
    for filename in input_files:
        try:
            file_type = determine_file_type(filename)
        except ValueError as e:
            print(e, file=sys.stderr)
            continue

        opener = gzip.open if filename.endswith('.gz') else open
        with opener(filename, 'rt') as f:
            if file_type == 'fa':
                reader = read_fasta(f)
                contig_cutoff = args.cutoff if args.cutoff is not None else 100
            elif file_type == 'fq':
                reader = read_fastq(f)
                contig_cutoff = args.cutoff if args.cutoff is not None else 1
            for _, seq in reader:
                stats['total_sequences'] += 1
                process_sequence(
                    seq,
                    stats,
                    args.scacut,
                    contig_cutoff,
                    args.Ncheck == 'y'
                )

    # Post-process statistics
    # Sort lengths and calculate max lengths
    for stype in ['scaffold', 'contig']:
        if stats[stype]['lengths']:
            stats[stype]['lengths'].sort(reverse=True)
            stats[stype]['max_length'] = stats[stype]['lengths'][0]
        else:
            stats[stype]['max_length'] = 0

    # Calculate N50-N90
    scaffold_n50 = calculate_n50(
        stats['scaffold']['lengths'],
        stats['scaffold']['total_length']
    )
    contig_n50 = calculate_n50(
        stats['contig']['lengths'],
        stats['contig']['total_length']
    )

    # Prepare output data
    output_data = {
        'scaffold': {
            'total_length': stats['scaffold']['total_length'],
            'max_length': stats['scaffold']['max_length'],
            'n_values': scaffold_n50,
            'gc_rate': stats['scaffold']['gc_count'] / stats['scaffold']['total_length'] if stats['scaffold']['total_length'] else 0,
            'count': stats['scaffold']['count'],
            'count_2k': stats['scaffold']['count_2k'],
        },
        'contig': {
            'total_length': stats['contig']['total_length'],
            'max_length': stats['contig']['max_length'],
            'n_values': contig_n50,
            'gc_rate': stats['contig']['gc_count'] / stats['contig']['total_length'] if stats['contig']['total_length'] else 0,
            'count': stats['contig']['count'],
            'count_2k': stats['contig']['count_2k'],
        },
        'N': stats['N'],
        'total_sequences': stats['total_sequences'],
    }

    # Generate output
    thousand_sep = args.thousp == 'y'

    # Simple table output
    if args.strout:
        with open(args.strout, 'w') as f:
            f.write("Sample name\tTotal number\tTotal number(>2kb)\tTotal bases(bp)\tMax length(bp)\tMean length(bp)\tN50(bp)\tN90(bp)\tGC content(%)\n")
            sample = args.sample if args.sample else 'sample'
            mean_len = stats['scaffold']['total_length'] // stats['total_sequences'] if stats['total_sequences'] else 0
            gc_percent = output_data['scaffold']['gc_rate'] * 100
            f.write(f"{sample}\t"
                    f"{format_number(stats['total_sequences'], thousand_sep)}\t"
                    f"{format_number(output_data['scaffold']['count_2k'], thousand_sep)}\t"
                    f"{format_number(output_data['scaffold']['total_length'], thousand_sep)}\t"
                    f"{format_number(output_data['scaffold']['max_length'], thousand_sep)}\t"
                    f"{format_number(mean_len, thousand_sep)}\t"
                    f"{format_number(scaffold_n50.get('5', {}).get('length', 0), thousand_sep)}\t"
                    f"{format_number(scaffold_n50.get('9', {}).get('length', 0), thousand_sep)}\t"
                    f"{gc_percent:.2f}\n")

    # Detailed formatted output
    if args.justiC == 'y':
        print("=" * 70)
        print(" " * 26 + "scaffold" + " " * 19 + "contig")
        print("{:<16}{:>15}{:>10}{:>15}{:>10}".format("", "length(bp)", "number", "length(bp)", "number"))

        # Max length
        print("{:<16}{:>15}{:>10}{:>15}{:>10}".format(
            "max_len",
            format_number(output_data['scaffold']['max_length'], thousand_sep),
            "",
            format_number(output_data['contig']['max_length'], thousand_sep),
            ""
        ))

        # N values
        for n in range(1, 10):
            key = str(n)
            scaf_n = scaffold_n50.get(key, {'length': 0, 'number': 0})
            cont_n = contig_n50.get(key, {'length': 0, 'number': 0})
            print("{:<16}{:>15}{:>10}{:>15}{:>10}".format(
                f"N{n}0",
                format_number(scaf_n['length'], thousand_sep),
                format_number(scaf_n['number'], thousand_sep),
                format_number(cont_n['length'], thousand_sep),
                format_number(cont_n['number'], thousand_sep)
            ))

        # Total length
        print("{:<16}{:>15}{:>10}{:>15}{:>10}".format(
            "Total_length",
            format_number(output_data['scaffold']['total_length'], thousand_sep),
            "",
            format_number(output_data['contig']['total_length'], thousand_sep),
            ""
        ))

        # Count >= cutoff
        print("{:<16}{:>15}{:>10}{:>15}{:>10}".format(
            f"number>={args.scacut}bp",
            "",
            format_number(output_data['scaffold']['count'], thousand_sep),
            "",
            format_number(output_data['contig']['count'], thousand_sep)
        ))

        # Count >= 2000
        print("{:<16}{:>15}{:>10}{:>15}{:>10}".format(
            "number>=2,000bp",
            "",
            format_number(output_data['scaffold']['count_2k'], thousand_sep),
            "",
            format_number(output_data['contig']['count_2k'], thousand_sep)
        ))

        # GC rate
        if args.showGC == 'y':
            print("=" * 70)
            print("{:<16}{:>15}{:>10}{:>15}{:>10}".format(
                "GC_rate",
                "",
                f"{output_data['scaffold']['gc_rate']:.3f}",
                "",
                f"{output_data['contig']['gc_rate']:.3f}"
            ))

        # N statistics
        if args.Ncheck == 'y':
            print("=" * 70)
            print("Total N bases: {}\tMin N: {}\tMax N: {}".format(
                format_number(output_data['N']['total'], thousand_sep),
                format_number(output_data['N']['min'] if output_data['N']['min'] is not None else 0, thousand_sep),
                format_number(output_data['N']['max'] if output_data['N']['max'] is not None else 0, thousand_sep)
            ))
        print("=" * 70)

if __name__ == '__main__':
    main()
    
