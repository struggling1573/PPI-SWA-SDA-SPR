import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_sequences(id_list_file, fasta_file):
    utglist = []
    with open(id_list_file, 'r') as id_file:
        for line in id_file:
            utglist.append(line.strip())

    d = {}
    seqs = SeqIO.parse(fasta_file, "fasta")
    for record in seqs:
        d[str(record.id)] = str(record.seq)

    for id in utglist:
        if id in d:
            print(f">{id}")
            print(d[id])
        else:
            print(f"Warning: ID {id} not found in the FASTA file.")


if __name__ == "__main__":
    author_info = """
Author: Haoyu Wang
Date: Aug 24  2023
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script extracts sequences from a FASTA file based on a list of sequence IDs. It reads a file containing a list of sequence IDs and a FASTA file, then prints the sequences corresponding to the provided IDs.",
        epilog=author_info,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('id_list_file', type=str, help='Path to the file containing a list of sequence IDs, one ID per line.')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file from which sequences will be extracted.')
    args = parser.parse_args()

    extract_sequences(args.id_list_file, args.fasta_file)
