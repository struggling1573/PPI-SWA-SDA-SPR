import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

USAGE = f"""
Usage: python %s genbank out.CDS.fasta out.pep.fasta

This script is used to extract CDS and protein sequences from a GenBank file using /product as primary ID.

Arguments:
  genbank         Path to the input GenBank file.
  out.CDS.fasta   Path to the output file for CDS sequences in FASTA format.
  out.pep.fasta   Path to the output file for protein sequences in FASTA format.

Example:
  python %s LY_10A.gb LY_10A.gb.CDS.fasta LY_10A.gb.pep.fasta

To show this help message, use the -h option:
  python %s -h

Author: Haoyu Wang
Date: Mar 13  2025
Affiliation: Southwest University
Contact: wanghyx666@163.com
""" % (sys.argv[0], sys.argv[0], sys.argv[0])

def main():
    if '-h' in sys.argv:
        print(USAGE)
        sys.exit(0)

    if len(sys.argv) != 4:
        print(USAGE)
        sys.exit(1)

    input_file = sys.argv[1]
    cds_output = sys.argv[2]
    protein_output = sys.argv[3]

    cds_records = []
    protein_records = []

    for record in SeqIO.parse(input_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                identifiers = [
                    feature.qualifiers.get('product', [None])[0],
                    feature.qualifiers.get('locus_tag', [None])[0],
                    feature.qualifiers.get('gene', [None])[0],
                    feature.qualifiers.get('protein_id', [None])[0]
                ]
                feature_id = next(
                    (str(i) for i in identifiers if i),
                    f"{record.id}_CDS_{feature.location.start}_{feature.location.end}"
                )

                cds_seq = feature.location.extract(record.seq)
                cds_records.append(SeqRecord(cds_seq, id=feature_id, description=""))

                if 'translation' in feature.qualifiers:
                    protein_seq = Seq(feature.qualifiers['translation'][0])
                    protein_records.append(SeqRecord(protein_seq, id=feature_id, description=""))

    SeqIO.write(cds_records, cds_output, "fasta")
    SeqIO.write(protein_records, protein_output, "fasta")

if __name__ == "__main__":
    main()
