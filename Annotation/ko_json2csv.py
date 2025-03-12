import argparse
import json
import re
import pandas as pd


def parse_kegg_json(input_file):
    with open(input_file) as f:
        ko_map_data = json.load(f)

    data = []
    for level1 in ko_map_data["children"]:
        m = re.match(r"(\S+)\s+([\S\w\s]+)", level1["name"])
        level1_pathway_id = m.group(1).strip()
        level1_pathway_name = m.group(2).strip()
        for level2 in level1["children"]:
            m = re.match(r"(\S+)\s+([\S\w\s]+)", level2["name"])
            level2_pathway_id = m.group(1).strip()
            level2_pathway_name = m.group(2).strip()
            for level3 in level2["children"]:
                m = re.match(r"(\S+)\s+([^\[]*)", level3["name"])
                level3_pathway_id = m.group(1).strip()
                level3_pathway_name = m.group(2).strip()
                if "children" in level3:
                    for ko in level3["children"]:
                        m = re.match(r"(\S+)\s+(\S+);\s+([^\[]+)\s*(\[EC:\S+(?:\s+[^\[\]]+)*\])*", ko["name"])
                        if m:
                            ko_id = m.group(1).strip()
                            ko_name = m.group(2).strip()
                            ko_des = m.group(3).strip()
                            ec = m.group(4) if m.group(4) else "-"
                            data.append([
                                level1_pathway_id, level1_pathway_name,
                                level2_pathway_id, level2_pathway_name,
                                level3_pathway_id, level3_pathway_name,
                                ko_id, ko_name, ko_des, ec
                            ])
    return data


def write_to_output(data, output_file):
    df = pd.DataFrame(data, columns=[
        "level1_pathway_id", "level1_pathway_name",
        "level2_pathway_id", "level2_pathway_name",
        "level3_pathway_id", "level3_pathway_name",
        "KO", "KO_name", "KO_des", "ec"
    ])
    df = df.drop_duplicates()
    df.to_csv(output_file, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Parse KEGG pathway JSON data and output a tab-separated text file.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Author: Hangyu Wang\nDate: Jan 04 2024\nUnit: Southwest University\nContact: wanghyx666@163.com"
    )
    parser.add_argument('input_file', type=str, help='Path to the input KEGG JSON file (e.g., ko00001.json)')
    parser.add_argument('output_file', type=str, help='Path to the output tab-separated text file')
    args = parser.parse_args()

    data = parse_kegg_json(args.input_file)
    write_to_output(data, args.output_file)
    
