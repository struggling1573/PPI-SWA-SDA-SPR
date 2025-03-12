import argparse
import obonet
import pandas as pd


def parse_go_basic_obo(input_file):
    go_data = []
    graph = obonet.read_obo(input_file)
    for node, data in graph.nodes(data=True):
        go_id = node
        go_name = ' '.join(word.capitalize() for word in data.get('name', '').split())
        go_class = data.get('namespace', '').replace('_', ' ')
        go_data.append([go_id, go_name, go_class])
    return go_data


def write_to_output(go_data, output_file):
    df = pd.DataFrame(go_data, columns=['GO ID', 'GO Name', 'GO Class'])
    df.to_csv(output_file, sep='\t', na_rep='nan', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Parse the go-basic.obo file and write the results to the output file.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Author: Hangyu Wang\nDate: Jan 04 2024\nUnit: Southwest University\nContact: wanghyx666@163.com"
    )
    parser.add_argument('input_file', type=str, help='Path to the input go-basic.obo file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    args = parser.parse_args()

    go_data = parse_go_basic_obo(args.input_file)
    write_to_output(go_data, args.output_file)
