import sys

USAGE = f"""
Usage: python rename_fasta.py mapping.txt sequences.fasta output.fasta

This script is used to rename sequence IDs in a FASTA file based on a mapping file.

Arguments:
  mapping.txt       Path to the mapping file, which should contain two columns separated by tabs.
                    The first column is the original ID, and the second column is the new ID.
  sequences.fasta   Path to the input FASTA file.
  output.fasta      Path to the output FASTA file with renamed sequence IDs.

To show this help message, use the -h option:
  python rename_fasta.py -h

Author: Haoyu Wang
Date: Mar 13  2025
Affiliation: Southwest University
Contact: wanghyx666@163.com
"""

if '-h' in sys.argv:
    print(USAGE)
    sys.exit(0)
    
if len(sys.argv) < 4:
    print(USAGE)
    sys.exit(1)

input_mapping = sys.argv[1]
input_fasta = sys.argv[2]
output_fasta = sys.argv[3]

ID_mapping = {}
try:
    with open(input_mapping, 'r', encoding='utf-8') as f1:
        for line in f1:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                ID_mapping[parts[0]] = parts[1]
except Exception as e:
    print(f"Error reading mapping file: {e}")
    sys.exit(1)

try:
    with open(input_fasta, 'r', encoding='utf-8') as f2, open(output_fasta, 'w', encoding='utf-8') as f3:
        for line in f2:
            if line.startswith('>'):
                fasta_id = line.strip()[1:]
                new_id = ID_mapping.get(fasta_id, fasta_id)
                f3.write('>{}\n'.format(new_id))
            else:
                f3.write(line)
except Exception as e:
    print(f"Error processing FASTA file: {e}")
    sys.exit(1)
