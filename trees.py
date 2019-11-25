import sys
import ete3
import argparse
from tqdm import tqdm
from statistics import median


def extract_distance_matrix(t):
    dists = {}
    for i,leafi in enumerate(t.get_leaves()):
        dists[leafi.name] = {}
        for j,leafj in tqdm(enumerate(t.get_leaves())):
            if i>=j: continue
            dists[leafi.name][leafj.name] = leafi.get_distance(leafj)
    return dists


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
        return extract_distance_matrix(self)




def main_extract_leaf_root_dists(args):
    t = tree(args.tree)
    dists = t.distance_from_root()
    for sample in dists:
        print(f"{sample}\t{dists[sample]}")

def main_extract_paired_dists(args):
    t = tree(args.tree)
    print(t.extract_distance_matrix())

def main_subsample(args):
    t = tree(args.tree)
    dists = t.distance_from_root()
    median_dist = median(list(dists.values()))
    dist_cutoff = args.arg1*median_dist
    for node in t.t.traverse():
        farthest_dist = node.get_farthest_leaf()[1]
        if farthest_dist > dist_cutoff: continue
        dists = extract_distance_matrix(node)
        print(dists)


parser = argparse.ArgumentParser(description='Tree analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('root_dist', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--tree',required=True)
parser_sub.set_defaults(func=main_extract_leaf_root_dists)

parser_sub = subparsers.add_parser('paired_dist', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--tree',required=True)
parser_sub.set_defaults(func=main_extract_paired_dists)

parser_sub = subparsers.add_parser('subsample', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--tree',required=True)
parser_sub.add_argument('--arg1',default=0.1,type=float,required=True)
parser_sub.set_defaults(func=main_subsample)


args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
