import sys, os, re
import argparse

def read_ref_gtf(ref_gtf):
    gene_trans = {}
    trans_exon = {}
    with open(ref_gtf, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            lines = line.strip().split("\t")
            if lines[2] == "transcript":
                gid = re.search(r'gene_id "(.*?)"', lines[8]).group(1)
                tid = re.search(r'transcript_id "(.*?)"', lines[8]).group(1)
                if gid not in gene_trans:
                    gene_trans[gid] = {}
                gene_trans[gid][tid] = line.strip()
                trans_exon[tid] = []
            elif lines[2] == "exon":
                tid = re.search(r'transcript_id "(.*?)"', lines[8]).group(1)
                trans_exon[tid].append(line.strip())

    return gene_trans, trans_exon

def get_trans(mydict):
    t_length = {}
    for k in mydict:
        lines = mydict[k].split("\t")
        t_length[k] = int(lines[4]) - int(lines[3]) + 1
    longest_trans_id = get_longest_trans(t_length)
    longest_info = mydict[longest_trans_id].strip().split("\t")
    return (longest_trans_id, longest_info[0], longest_info[1], longest_info[3], longest_info[4], longest_info[5],
            longest_info[6], ",".join(t_length.keys()))


def sum_exon(exons):
    length = 0
    for line in exons:
        info = line.strip().split("\t")
        length += (int(info[4]) - int(info[3]) + 1)
    return length


def get_longest_trans(length):
    for key, value in length.items():
        if value == max(length.values()):
            return key


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is designed to process and analyze genomic data related to gene transcripts. It selects coding transcripts from CPC (Coding Potential Calculator) results, extracts relevant information from cuffcompare results, and then combines and modifies GTF (Gene Transfer Format) files. The script aims to generate a final GTF file that includes necessary information about genes, transcripts, and exons.",
        epilog="Output files:\n"
               "  - final.gtf: The final processed GTF file that includes necessary information about genes,\n"
               "               transcripts, and exons. Each gene entry contains additional attributes\n"
               "               such as gene length, number of transcripts, and length of the longest transcript.\n"
               "  - {final_gtf_name}_merged.gtf: A merged GTF file that combines the reference GTF and the new\n"
               "                                transcript GTF without considering the CPC results. It retains\n"
               "                                all transcript information for comprehensive analysis.\n"
               "  - {final_gtf_name}_coding.gtf: A GTF file containing only the coding transcripts filtered by\n"
               "                                the CPC results. It extracts coding transcripts from the new\n"
               "                                transcript GTF and updates the gene IDs according to the new\n"
               "                                transcript list.\n"
               "  - 02final.gtf: A filtered version of the final GTF file, where gene lines are removed. It\n"
               "                 only contains transcript and exon information, suitable for analyses that\n"
               "                 only require this level of detail.\n"
               "Note: {final_gtf_name} is the base name of the final GTF file specified by the -o option.\n"
               "\n\n"
               "Author: Haoyu Wang\n"
               "Date: May 23  2024\n"
               "Affiliation: Southwest University\n"
               "Contact: wanghyx666@163.com\n",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--input', help='cpc result', dest='cpc', type=str)
    parser.add_argument('-o', '--output', help='final gtf', dest='final_gtf', default='final.gtf', type=str)
    parser.add_argument('-n_gtf', '--new_gtf', help='new transcript gtf', required=True, dest='n_gtf', type=str)
    parser.add_argument('-r_gtf', '--ref_gtf', help='ref gtf', required=True, dest='r_gtf', type=str)
    parser.add_argument('-l', '--trans_list', help='new transcript list', required=True, dest='trans_list', type=str)
    args = parser.parse_args()

    coding_list = []
    if args.cpc:
        with open(args.cpc) as f:
            f.readline()
            for line in f:
                tmp = line.strip().split('\t')
                if tmp[-1] == 'coding':
                    coding_list.append(tmp[0])

    novel_gene = {}
    novel_isoform = {}
    if args.cpc:
        with open(args.trans_list) as f:
            for line in f:
                tmp = line.strip().split('\t')
                if tmp[-1] in coding_list:
                    if tmp[0] == '-':
                        novel_gene[tmp[-1]] = tmp[-2]
                    else:
                        novel_isoform[tmp[-1]] = tmp[0]

    with open(args.r_gtf) as f, open(args.n_gtf) as f1, open(args.final_gtf, 'w') as f2:
        f2.write(f.read().strip() + '\n')
        for line in f1:
            tmp = line.strip().split('\t')

            tid = re.search(r'transcript_id "(.*?)"', tmp[8]).group(1)
            if args.cpc and tid in novel_gene.keys():
                f2.write(line)
            elif args.cpc and tid in novel_isoform.keys():
                tmp[8] = 'gene_id "%s"; transcript_id "%s";' % (novel_isoform[tid], tid)
                f2.write('\t'.join(tmp) + '\n')
            elif not args.cpc:
                f2.write(line)

    merged_gtf_path = os.path.splitext(args.final_gtf)[0] + "_merged.gtf"
    with open(args.r_gtf) as f_ref, open(args.n_gtf) as f_new, open(merged_gtf_path, 'w') as f_merged:
        f_merged.write(f_ref.read().strip() + '\n')
        f_merged.write(f_new.read().strip())

    if args.cpc:
        coding_transcript = []
        with open("%s/cpc2.predict.result.txt" % os.path.dirname(args.final_gtf), "r") as f2:
            for line in f2:
                lines = line.strip().split("\t")
                if lines[7] == "coding":
                    coding_transcript.append(lines[0])
        print("coding transcript num:%s" % len(coding_transcript))

        new_transcript_gtf = {}
        append_list = []
        with open(args.n_gtf, "r") as f1:
            for line in f1:
                if line.startswith("#"):
                    continue
                lines = line.split("\t")
                if lines[2] == "transcript" and lines[6] != ".":
                    match = re.search(r'transcript_id "(.*?)"', lines[8])
                    new_transcript_gtf[match.group(1)] = [line.strip()]
                    append_list.append(match.group(1))
                elif lines[2] == "exon" and lines[6] != ".":
                    match = re.search(r'transcript_id "(.*)"', lines[8])
                    if match.group(1) in append_list:
                        new_transcript_gtf[match.group(1)].append(line.strip())

        coding_gtf_path = os.path.splitext(args.final_gtf)[0] + "_coding.gtf"
        with open(args.trans_list) as f4, open(coding_gtf_path, "w") as f3:
            for line in f4:
                lines = line.strip().split("\t")
                tid = lines[4]
                if tid in coding_transcript and tid in new_transcript_gtf:
                    old_gid = lines[3]
                    if lines[0] != "-":
                        new_gid = lines[0]
                    else:
                        new_gid = lines[3]
                        print(new_gid)
                    txt = "\n".join(new_transcript_gtf[tid])
                    old_txt = "gene_id \"%s\"" % old_gid
                    new_txt = "gene_id \"%s\"" % new_gid
                    txt = txt.replace(old_txt, new_txt)
                    f3.write(txt)
                    f3.write("\n")

        os.system("cat %s %s > %s/tmp.gtf" % (
            args.r_gtf, coding_gtf_path, os.path.dirname(args.final_gtf)))
    else:
        os.system("cat %s %s > %s/tmp.gtf" % (
            args.r_gtf, args.n_gtf, os.path.dirname(args.final_gtf)))

    gene_trans, trans_exon = read_ref_gtf("%s/tmp.gtf" % os.path.dirname(args.final_gtf))
