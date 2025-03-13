import sys
from Bio import SeqIO

USAGE = f"""
Usage: python %s genbank out.CDS.fasta out.pep.fasta prefix

This script is used to extract CDS and protein sequences from a GenBank file.

Arguments:
  genbank         Path to the input GenBank file.
  out.CDS.fasta   Path to the output file for CDS sequences in FASTA format.
  out.pep.fasta   Path to the output file for protein sequences in FASTA format.
  prefix          Prefix to be added to the gene names in the output files.

Example:
  python %s LY_10A.gb LY_10A.gb.CDS.fasta LY_10A.gb.pep.fasta LY_10A

To show this help message, use the -h option:
  python %s -h

Author: Hangyu Wang
Date: Feb 05 2024
Unit: Southwest University
Contact: wanghyx666@163.com
""" % (sys.argv[0], sys.argv[0], sys.argv[0])

if '-h' in sys.argv:
    print(USAGE)
    sys.exit()

if len(sys.argv) != 5:
    print(USAGE)
    sys.exit()

record = SeqIO.read(sys.argv[1], "genbank")
fw_CDS = open(sys.argv[2], "w")
fw_pep = open(sys.argv[3], "w")
prefix = sys.argv[4]
for feature in record.features:
    if feature.type == 'CDS':
        gene_name = feature.qualifiers['gene'][0]
        cds_seq = feature.extract(record.seq)
        pep_seq = feature.qualifiers['translation'][0]
        fw_CDS.write(f">{prefix} {gene_name}\n{cds_seq}\n")
        fw_pep.write(f">{prefix} {gene_name}\n{pep_seq}\n")
      
fw_CDS.close()
fw_pep.close()
