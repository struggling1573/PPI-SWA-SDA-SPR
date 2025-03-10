import re

def parse_go_basic_obo(input_file):
    go_ids = []
    go_names = []
    go_classes = []
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('id:'):
                go_id = line.strip().split(' ')[1]
                go_ids.append(go_id)
            elif line.startswith('name:'):
                go_name = line.strip().replace('name:', '').strip()
                go_name = ' '.join([word.capitalize() for word in go_name.split(' ')])
                go_names.append(go_name)
            elif line.startswith('namespace:'):
                go_class = line.strip().split(' ')[1]
                go_class = go_class.replace('_', ' ')
                go_classes.append(go_class)

    return go_ids, go_names, go_classes


def write_to_output(go_ids, go_names, go_classes, output_file):
    with open(output_file, 'w', encoding='utf-8') as f:
        for go_id, go_name, go_class in zip(go_ids, go_names, go_classes):
            f.write(f"{go_id}\t{go_name}\t{go_class}\n")


if __name__ == "__main__":
    input_file = "go-basic.obo"
    output_file = "GO.library"
    go_ids, go_names, go_classes = parse_go_basic_obo(input_file)
    write_to_output(go_ids, go_names, go_classes, output_file)
