import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(
        description='Extract CDS and protein sequences from GenBank file using /product as primary ID',
        epilog='Author: Haoyu Wang\nDate: Dec 30 2024\nAffiliation: Southwest University\nContact: wanghyx666@163.com',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Input GenBank file path')
    parser.add_argument('-c', '--cds', required=True, help='Output FASTA file for CDS sequences')
    parser.add_argument('-p', '--protein', required=True, help='Output FASTA file for protein sequences')
    args = parser.parse_args()
    
    cds_records = []
    protein_records = []
    for record in SeqIO.parse(args.input, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                identifiers = [
                    feature.qualifiers.get('product', [None])[0],  # 第一优先级
                    feature.qualifiers.get('locus_tag', [None])[0],
                    feature.qualifiers.get('gene', [None])[0],
                    feature.qualifiers.get('protein_id', [None])[0]
                ]
                feature_id = next(
                    (str(i) for i in identifiers if i),
                    f"{record.id}_CDS_{feature.location.start}_{feature.location.end}"
                )
                cds_seq = feature.location.extract(record.seq)
                cds_record = SeqRecord(cds_seq, id=feature_id, description="")
                cds_records.append(cds_record)

                if 'translation' in feature.qualifiers:
                    protein_seq = Seq(feature.qualifiers['translation'][0])
                    protein_record = SeqRecord(protein_seq, id=feature_id, description="")
                    protein_records.append(protein_record)
    SeqIO.write(cds_records, args.cds, "fasta")
    SeqIO.write(protein_records, args.protein, "fasta")

if __name__ == "__main__":
    main()
