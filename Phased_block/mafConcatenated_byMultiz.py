import argparse
import subprocess
import os
import sys


def merge_maf_files(output_file, input_files):
    temp_files = []
    try:
        current = input_files[0]
        for i in range(1, len(input_files)):
            tmp_file = f"tmp{i}.maf" if i < len(input_files) - 1 else output_file
            temp_files.append(tmp_file)
            next_file = input_files[i]
            cmd_str = f"multiz M=1 {current} {next_file} 0 U1 U2 > {tmp_file}"
            print(f"\n🚀 Executing step {i}: {cmd_str}")
            cmd_list = ["multiz", "M=1", current, next_file, "0", "U1", "U2"]
            with open(tmp_file, 'w') as f:
                proc = subprocess.run(
                    cmd_list,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True
                )
            if proc.returncode != 0:
                error_msg = f"Command execution failed (code {proc.returncode}):\n"
                error_msg += f"Error message:\n{proc.stderr.strip()}"
                raise RuntimeError(error_msg)
            current = tmp_file
        print(f"\n✅ Final output: {os.path.abspath(output_file)}")

    finally:
        # Clean up temporary files (keep the final output)
        if temp_files:
            print("\n🧹 Cleaning up intermediate files...")
            for f in temp_files[:-1]:  # Don't delete the final output file
                if os.path.exists(f):
                    os.remove(f)
                    print(f"Deleted: {f}")

def main():
    author_info = """
Author: Hangyu Wang
Date: Dec 28 2024
Unit: Southwest University
Contact: wanghyx666@163.com
"""
    epilog = f"Usage example:\n" \
             "  python merge_maf.py -o merged.maf hap1.maf hap2.maf hap3.maf\n" \
             "Output example:\n" \
             "  🚀 Executing step 1: multiz M=1 hap1.maf hap2.maf 0 U1 U2 > tmp1.maf\n" \
             "  🚀 Executing step 2: multiz M=1 tmp1.maf hap3.maf 0 U1 U2 > merged.maf\n" \
             f"{author_info}"

    parser = argparse.ArgumentParser(
        description="MAF file merging tool (displays the full command). This tool merges multiple MAF files step - by - step using the 'multiz' command and shows the actual commands executed for debugging purposes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog
    )
    parser.add_argument("-o", "--output", required=True, help="Path to the output file.")
    parser.add_argument("input_files", nargs="+", help="List of input MAF files (in the order of merging).")
    args = parser.parse_args()
    if len(args.input_files) < 2:
        print("✖ Error: At least 2 input files are required.")
        sys.exit(1)
    for f in args.input_files:
        if not os.path.isfile(f):
            print(f"✖ Error: Input file does not exist - {f}")
            sys.exit(1)
    try:
        merge_maf_files(args.output, args.input_files)
    except Exception as e:
        print(f"\n❌ Merging failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
