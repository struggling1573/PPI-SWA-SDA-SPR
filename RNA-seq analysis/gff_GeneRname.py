import re
import sys
import argparse
import logging
import gffutils
from collections import defaultdict


def get_actual_transcript_types(gff_file):
    possible_transcript_types = ["transcript", "mRNA"]
    actual_transcript_types = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                feature_type = fields[2]
                if feature_type in possible_transcript_types and feature_type not in actual_transcript_types:
                    actual_transcript_types.append(feature_type)
    return actual_transcript_types


def rename(args):
    seqid2name = dict()
    for line in open(args.change, 'r'):
        tem = line.strip().split()
        seqid2name[tem[0]] = tem[1]

    db = gffutils.create_db(args.gff, ':memory:', force=True, keep_order=True, merge_strategy="create_unique",
                            sort_attribute_values=True)
    mRNA_children = ("five_prime_UTR", "three_prime_UTR", "CDS", "exon")
    idmap = {
        "CDS": "cds",
        "exon": "exon",
        "five_prime_UTR": "5utrp",
        "three_prime_UTR": "3utrp"
    }

    actual_transcript_types = get_actual_transcript_types(args.gff)

    f_out = open(args.output, 'w')
    gene_mapping = {}
    transcript_mapping = {}

    seqid = None
    for gene in db.features_of_type("gene", order_by=('seqid', 'start', 'end')):
        if gene.seqid != seqid:
            genenum = 0
        genenum += args.addnum
        seqid = gene.seqid

        original_gene_id = gene.id
        genename = '{0}{1}{2:06}'.format(seqid2name[seqid], args.separator, genenum)
        gene_mapping[original_gene_id] = genename
        f_out.write('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={geneid}\n'.format(
            seqid=gene.seqid, source=gene.source, featuretype=gene.featuretype, start=gene.start, end=gene.end,
            score=gene.score, strand=gene.strand, frame=gene.frame, geneid=genename))

        for transcript_featuretype in actual_transcript_types:
            for t, transcript in enumerate(
                    db.children(gene, featuretype=transcript_featuretype, order_by=('seqid', 'start', 'end'))):
                original_transcript_id = transcript.id
                transcript_num = t + 1
                transcript_id = '{genename}.mRNA{num}'.format(genename=genename, num=transcript_num)
                transcript_mapping[original_transcript_id] = transcript_id
                f_out.write('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={transcript_id};Parent={geneid}\n'.format(
                    seqid=transcript.seqid, source=transcript.source, featuretype=transcript.featuretype,
                    start=transcript.start, end=transcript.end, score=transcript.score, strand=transcript.strand,
                    frame=transcript.frame, transcript_id=transcript_id, geneid=genename))

                numdict = defaultdict(int)
                for child in db.children(transcript, featuretype=mRNA_children, order_by=("start", 'end')):
                    numdict[child.featuretype] += 1
                    child_id = '{genename}.{childid}{num}'.format(genename=genename,
                                                                  childid=idmap[child.featuretype],
                                                                  num=numdict[child.featuretype])
                    f_out.write('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={child_id};Parent={transcript_id}\n'.format(
                        seqid=child.seqid, source=child.source, featuretype=child.featuretype, start=child.start,
                        end=child.end, score=child.score, strand=child.strand, frame=child.frame, child_id=child_id,
                        transcript_id=transcript_id))

    f_out.close()

    if args.m_gene:
        with open(args.m_gene, 'w') as gene_mapping_file:
            gene_mapping_file.write("Original_Gene_ID\tRenamed_Gene_ID\n")
            for original_id, renamed_id in gene_mapping.items():
                gene_mapping_file.write(f"{original_id}\t{renamed_id}\n")

    if args.m_mRNA:
        with open(args.m_mRNA, 'w') as transcript_mapping_file:
            transcript_mapping_file.write("Original_Transcript_ID\tRenamed_Transcript_ID\n")
            for original_id, renamed_id in transcript_mapping.items():
                transcript_mapping_file.write(f"{original_id}\t{renamed_id}\n")


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="rename gff3 file"
    )

    parser.add_argument('-g', '--gff', required=True, help='gff3 file')
    parser.add_argument('-c', '--change', required=True, help='a file, correspondence between sequence name and gene name prefix')
    parser.add_argument('-a', '--addnum', type=int, default=1, help='diff in gene number, such as if addnum = 10, xx1G000010, xx1G000020')
    parser.add_argument('-o', '--output', required=True, help='Output GFF file name') 
    parser.add_argument('-s', '--separator', default='G', help='Separator used in gene ID, e.g., G or newG')
    parser.add_argument('-m_gene', help='Output file for mapping of original and renamed gene IDs')
    parser.add_argument('-m_mRNA', help='Output file for mapping of original and renamed mRNA/transcript IDs')

    author_info = """
Author: Haoyu Wang
Date: May 23  2024
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser.epilog = author_info

    args = parser.parse_args()

    rename(args)


if __name__ == "__main__":
    main()
