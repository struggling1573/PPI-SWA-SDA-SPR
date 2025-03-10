import argparse
import pysam

def main():
    author_info = """
Author: Haoyu Wang
Date: Oct 04  2023
Affiliation: Southwest University
Contact: wanghyx666@163.com
    """
    parser = argparse.ArgumentParser(
        description="This script extracts reads from a BAM file corresponding to a list of contigs provided in a file. It creates a new BAM file containing only the reads mapped to the specified contigs.",
        epilog=author_info,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('contig_list_file', type=str, help='Path to the file containing a list of contig names, one per line.')
    parser.add_argument('bam_file', type=str, help='Path to the input BAM file.')
    args = parser.parse_args()
    contig_list_file = args.contig_list_file
    bam_file = args.bam_file
    try:
        with open(contig_list_file, 'r') as file:
            contigs = [line.strip() for line in file]
    except FileNotFoundError:
        print(f"Error: Contig list file '{contig_list_file}' not found.")
        return
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
    except FileNotFoundError:
        print(f"Error: BAM file '{bam_file}' not found.")
        return
    except OSError as e:
        print(f"Error opening BAM file: {e}")
        return
    outfname = f"{contig_list_file.split('.')[0]}.sub.bam"
    outf = pysam.AlignmentFile(outfname, "wb", template=samfile)
    for contig in contigs:
        try:
            for read in samfile.fetch(contig):
                outf.write(read)
        except ValueError:
            print(f"Warning: Contig '{contig}' not found in BAM file.")
            continue
    outf.close()
    samfile.close()
    print(f"Data extraction complete. Output file created: {outfname}")


if __name__ == "__main__":
    main()