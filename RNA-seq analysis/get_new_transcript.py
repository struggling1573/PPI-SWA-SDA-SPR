import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re

def fun_one(genome):
    return str(Seq(genome).reverse_complement())

def fun_two(gtf_file, fa_file, list_file, out_fa, out_gtf):
    df_list = pd.read_csv(list_file, sep='\t', header=None)
    tid_list = set(df_list[4].tolist())
    fasta = SeqIO.to_dict(SeqIO.parse(fa_file, 'fasta'))
    transcript_regex = re.compile(r'(.*)transcript_id "(.*?)";')
    exon_regex = re.compile(r'(.*)transcript_id "(.*)"; exon_number')

    append_fa = set()
    with open(out_fa, 'w') as f1, open(out_gtf, 'w') as f2:
        for chunk in pd.read_csv(gtf_file, sep='\t', comment='#', header=None, chunksize=10000):
            chunk.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
            transcript_chunk = chunk[chunk['feature'] == 'transcript']
            for _, row in transcript_chunk.iterrows():
                match = transcript_regex.search(row['attribute'])
                if match:
                    tid = match.group(2)
                    if tid in tid_list:
                        seq = str(fasta[row['seqname']].seq)
                        f1.write(f'>{tid}\n')
                        sequence = seq[row['start'] - 1:row['end']]
                        if row['strand'] == '+':
                            f1.write(f'{sequence}\n')
                        else:
                            f1.write(f'{fun_one(sequence)}\n')
                        append_fa.add(tid)
                        f2.write('\t'.join(row.astype(str).tolist()) + '\n')
            exon_chunk = chunk[chunk['feature'] == 'exon']
            for _, row in exon_chunk.iterrows():
                match = exon_regex.search(row['attribute'])
                if match:
                    tid = match.group(2)
                    if tid in append_fa:
                        f2.write('\t'.join(row.astype(str).tolist()) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script is used to process genomic data. It extracts specific transcripts and their corresponding sequences and exon information based on the provided GTF file, reference FASTA file, and a list of transcript IDs. The script can handle both forward and reverse-stranded transcripts and generate output files in FASTA and GTF formats.",
                                     epilog="Author: Haoyu Wang\n"
                                            "Date: May 23  2024\n"
                                            "Affiliation: Southwest University\n"
                                            "Contact: wanghyx666@163.com\n",
                                            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='merge gtf', dest='i', type=str)
    parser.add_argument('-o', '--output', help='output prefix', dest='o', default='new_transcript', type=str)
    parser.add_argument('-fa', '--ref_fa', help='ref fasta', required=True, dest='fa', type=str)
    parser.add_argument('-l', '--trans_list', help='new transcript list', required=True, dest='l', type=str)
    args = parser.parse_args()

    fun_two(args.i, args.fa, args.l, args.o + ".fa", args.o + ".gtf")
