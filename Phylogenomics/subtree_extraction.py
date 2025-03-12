from ete3 import Tree
import argparse


def has_bootstrap_info(tree):
    for node in tree.traverse():
        if hasattr(node, 'support') and node.support is not None:
            return True
    return False


def has_dist_info(tree):
    for node in tree.traverse():
        if hasattr(node, 'dist') and node.dist != 0:
            return True
    return False


def extract_subtree(tree_file, species_file, keep_bootstrap=True, keep_dist=True):
    try:
        tree = Tree(tree_file, format=1)

        with open(species_file, 'r') as f:
            # 过滤掉空行
            species_list = [line.strip().split()[0] for line in f.readlines() if line.strip()]

        existing_species = []
        non_existing_species = []
        for species in species_list:
            node = tree.search_nodes(name=species)
            if node:
                existing_species.append(species)
            else:
                non_existing_species.append(species)

        if non_existing_species:
            for species in non_existing_species:
                print(f"Note: {species} is not in the input tree!")

        target_nodes = []
        for species in existing_species:
            node = tree.search_nodes(name=species)
            if node:
                target_nodes.append(node[0])

        subtree = tree.get_common_ancestor(target_nodes)
        subtree.prune(existing_species)

        if has_bootstrap_info(subtree) and not keep_bootstrap:
            for node in subtree.traverse():
                if hasattr(node, 'support'):
                    node.support = None

        if has_dist_info(subtree) and not keep_dist:
            for node in subtree.traverse():
                if hasattr(node, 'dist'):
                    node.dist = 0

        return subtree

    except FileNotFoundError as e:
        print(f"Error: File {e.filename} not found.")
    except Exception as e:
        print(f"An unknown error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract a subtree corresponding to a subset of species from a rooted phylogenetic tree, and optionally keep or remove bootstrap values and genetic distances.',
        epilog="Author: Haoyu Wang\n"
               "Date: Dec 14  2023\n"
               "Affiliation: Southwest University\n"
               "Contact: wanghyx666@163.com\n",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('input_tree', type=str, help='Path to the input phylogenetic tree file (in Newick format)')
    parser.add_argument('species_file', type=str, help='Path to the file containing the species IDs to be extracted. The first column should be the species ID.')
    parser.add_argument('output_tree', type=str, help='Path to the output subtree file')
    parser.add_argument('--bootstrap', type=str, default='y', choices=['y', 'n'],
                        help='Whether to keep the bootstrap values. Default is y (keep).')
    parser.add_argument('--dist', type=str, default='y', choices=['y', 'n'],
                        help='Whether to keep the genetic distances. Default is y (keep).')

    args = parser.parse_args()

    keep_bootstrap = args.bootstrap.lower() == 'y'
    keep_dist = args.dist.lower() == 'y'

    subtree = extract_subtree(args.input_tree, args.species_file, keep_bootstrap, keep_dist)
    if subtree:
        try:
            with open(args.output_tree, 'w') as f:
                f.write(subtree.write(format=1))
            print(f"The subtree has been successfully saved to {args.output_tree}")
        except Exception as e:
            print(f"Error occurred while saving the subtree file: {e}")
    
