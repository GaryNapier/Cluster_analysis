import sys
import ete3
import argparse
from tqdm import tqdm

class tree:
    def __init__(self,treefile):
        self.treefile = treefile
        self.t = ete3.Tree(treefile)

    def distance_from_root(self):
        dists = {}
        for leaf in tqdm(self.t.get_leaves()):
            dists[leaf.name] = self.t.get_distance(leaf)
        return dists

    def extract_distance_matrix(self):
        dists = {}
        for i,leafi in enumerate(self.t.get_leaves()):
            dists[leafi.name] = {}
            for j,leafj in tqdm(enumerate(self.t.get_leaves())):
                if i>=j: continue
                dists[leafi.name][leafj.name] = leafi.get_distance(leafj)
        return dists




def main_extract_leaf_root_dists(args):
    t = tree(args.tree)
    dists = t.distance_from_root()
    for sample in dists:
        print(f"{sample}\t{dists[sample]}")

def main_extract_paired_dists(args):
    t = tree(args.tree)
    print(t.extract_distance_matrix())

parser = argparse.ArgumentParser(description='Tree analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('root_dist', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--tree',required=True)
parser_sub.set_defaults(func=main_extract_leaf_root_dists)

parser_sub = subparsers.add_parser('paired_dist', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--tree',required=True)
parser_sub.set_defaults(func=main_extract_paired_dists)


args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
