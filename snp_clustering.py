import csv
from networkx.algorithms.components.connected import connected_components
import networkx
import argparse
import json
import os
import os.path
import sys
import subprocess
import random
rand_generator = random.SystemRandom()
import statistics
import subprocess as sp
# version = 1.0

# ----------------
# Universal setup/variables
# ----------------
round_place = 3
path_command = "find -name "

def samples2itol(l, append_samps_only=False):
    arr = []
    text = """DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL    samples
COLOR    #ff0000

LEGEND_TITLE    samples
LEGEND_SHAPES    1
LEGEND_COLORS    black
LEGEND_LABELS    Sample

DATA
"""
    if append_samps_only is False:
        arr.append(text)
        for x in l:
            arr.append("%s\tblack" % x.rstrip())
        return "\n".join(arr)
    else:
        for x in l:
            arr.append("%s\tblack" % x.rstrip())
        return "\n".join(arr)


def get_random_file(prefix=None, extension=None):
    randint = rand_generator.randint(1, 999999)
    if prefix:
        if extension:
            return "%s.%s%s" % (prefix, randint, extension)
        else:
            return "%s.%s.txt" % (prefix, randint)
    else:
        if extension:
            return "%s.tmp%s" % (randint, extension)
        else:
            return "%s.tmp.txt" % (randint)


def filecheck(filename):
    """
    Check if file is there and quit if it isn't
    """
    if not os.path.isfile(filename):
        print("Can't find %s" % filename)
        exit(1)
    else:
        return filename


def nofile(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isfile(filename):
        return True
    else:
        return False


def index_vcf(vcf_file, threads=4, overwrite=False):
    """
    Indexing a bam file
    """
    cmd = "bcftools index --threads %s -f %s" % (threads, vcf_file)
    if filecheck(vcf_file):
        if nofile(vcf_file+".csi"):
            run_cmd(cmd)
        elif (os.path.getmtime(vcf_file+".csi") < os.path.getmtime(vcf_file)) or overwrite:
            run_cmd(cmd)


def run_cmd(cmd, verbose=1, target=None):
    """
    Wrapper to run a command using subprocess with 3 levels of verbosity and
    automatic exiting if command failed
    """
    cmd = "set -u pipefail; " + cmd
    if verbose == 2:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stdout = open("/dev/stdout", "w")
        stderr = open("/dev/stderr", "w")
    elif verbose == 1:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stdout = open("/dev/null", "w")
        stderr = open("/dev/null", "w")
    else:
        stdout = open("/dev/null", "w")
        stderr = open("/dev/null", "w")

    res = subprocess.call(cmd, shell=True, stderr=stderr, stdout=stdout)
    stderr.close()
    if res != 0:
        print("Command Failed! Please Check!")
        exit(1)


class vcf:
    def __init__(self, filename, prefix=None, threads=4, keepfile=None):
        self.samples = []
        self.filename = filename
        self.threads = threads
        self.keepfile = keepfile
        if prefix is None:
            if filename[-4:] == ".bcf":
                self.prefix = filename[:-4]
            elif filename[-5:] == ".gbcf":
                self.prefix = filename[:-5]
            elif filename[-7:] == ".vcf.gz":
                self.prefix = filename[:-7]
            elif filename[-4:] == ".vcf":
                self.prefix = filename[:-4]
            else:
                self.prefix = filename
        else:
            self.prefix = prefix
        self.prefix = self.prefix
        self.temp_file = get_random_file()
        index_vcf(filename, self.threads)
        cmd = "bcftools query -l %(filename)s > %(temp_file)s" % vars(self)
        run_cmd(cmd)
        for l in open(self.temp_file):
            self.samples.append(l.rstrip())
        os.remove(self.temp_file)
        self.vcf = "%s.vcf" % self.prefix


    def get_plink_dist(self):
        tmpfile = get_random_file()
        # keep_cmd = " --keep %s " % keepfile if keepfile else  ""
        keep_cmd = " --keep %s " % self.keepfile if self.keepfile else  ""
        cmd = "plink --vcf %s %s --distance square --allow-extra-chr --out %s --double-id" % ( self.filename, keep_cmd, tmpfile)
        run_cmd(cmd)
        O = open("%s.dist" % (self.prefix), "w")
        dists = []
        for l in open("%s.dist" % tmpfile):
            row = [float(d)/2 for d in l.rstrip().split()]
            O.write("%s\n" % "\t".join([str(x) for x in row]))
            dists.append(row)
        O.close()
        run_cmd("rm %s*" % tmpfile)
        return dists

    def get_clusters(self, cutoff=10, remove_singletons=False):
        dists = self.get_plink_dist()
        edges = []
        tmp_node_set = set()
        for i in range(len(dists)):
            for j in range(len(dists)):
                if j >= i:
                    continue
                if dists[i][j] < cutoff:
                    edge = {
                        "source": self.samples[i], "target": self.samples[j], "snps": dists[i][j]}
                    tmp_node_set.add(self.samples[i])
                    tmp_node_set.add(self.samples[j])
                    edges.append(edge)
        nodes = [{"id": s} for s in tmp_node_set] if remove_singletons else [
            {"id": s} for s in self.samples]
        graph = {"nodes": nodes, "edges": edges}
        json.dump(graph, open("%s.ters.json" % self.prefix, "w"))
        return graph


class transmission_graph:
    def __init__(self, filename):
        tmp = json.load(open(filename))
        self.json_graph = tmp
        self.graph = networkx.Graph()
        self.graph.add_nodes_from([x["id"] for x in tmp["nodes"]])
        for edge in tmp["edges"]:
            self.graph.add_edges_from([(edge["source"], edge["target"])])
        self.clusters = sorted(list(connected_components(self.graph)), key=lambda x: len(x), reverse=True)

    def add_meta_data(self, csvfile):
        for row in csv.DictReader(open(csvfile)):
            found = False
            for i, node in enumerate(self.json_graph["nodes"]):
                if node["id"] == row["id"]:
                    found = True
                    break
            if found:
                for column in set(row.keys())-set(["id"]):
                    if column not in self.json_graph["nodes"][i]:
                        self.json_graph["nodes"][i][column] = row[column]

    def extract_clusters(self):
        print(self.clusters)


# def main_stats(args):
#     graph = transmission_graph(args.graph)
#     print("Num clusters: %s" % len(graph.clusters))
#     print("Mean clusters size: %s" %
#           (sum([len(x) for x in graph.clusters])/len(graph.clusters)))

def main_stats(args):
    graph = args.graph
    filecheck(graph)
    graph = transmission_graph(graph)
    min_clust_size = args.cluster_minimum
    stats_output_file = args.stats_output_file
    len_clusts = []
    for clust in graph.clusters:
        if len(clust) >= min_clust_size:
            len_clusts.append(len(clust))

    # Save and tidy variables. Max and stdev are fussy so use try/except
    n_samps = sum(len_clusts)
    n_clusts = len(len_clusts)
    max_clust_sz = max(len_clusts or [0])
    min_clust_sz = min(len_clusts or [0])
    try:
        mean_clust_sz = round((sum(len_clusts) / len(len_clusts)), round_place)
    except:
        mean_clust_sz = 0

    try:
        sd_clust_sz = round(statistics.stdev(len_clusts), round_place)
    except:
        sd_clust_sz = 0

    # Put in dictionary
    stats_dict = {"N samples in clusters": n_samps,
    "Number of clusters": n_clusts,
    "Max clust size": max_clust_sz,
    "Min clust size": min_clust_sz,
    "Mean clusters size": mean_clust_sz,
    "s.d. clusters": sd_clust_sz}

    print(stats_dict)

    # Write dictionary to file:
    if nofile(stats_output_file):
        with open(stats_output_file, 'w') as f:
            writer = csv.DictWriter(f, fieldnames = list(stats_dict.keys()), delimiter = '\t')
            writer.writeheader()
            writer.writerows([stats_dict])
    else:
        with open(stats_output_file, 'a+') as f:
            writer = csv.DictWriter(f, fieldnames = list(stats_dict.keys()), delimiter = '\t')
            writer.writerows([stats_dict])

    # Save samples to itol file
    itol_file = args.itol_file
    for cluster in graph.clusters:
        if len(cluster)>=min_clust_size:
            if nofile(itol_file):
                open(itol_file,"w").write(samples2itol(list(cluster)))
            else:
                open(itol_file,"a+").write(samples2itol(list(cluster)))

def main_add_meta(args):
    graph = transmission_graph(args.graph)
    graph.add_meta_data(args.csv)
    json.dump(graph.json_graph, open(args.out, "w"))


def main_view(args):
    graph = transmission_graph(args.graph)
    graph.extract_clusters()


def main_downsample(args):
    graph = transmission_graph(args.graph)
    for i,cluster in enumerate(graph.clusters):
        if len(cluster)>args.min_size:
            open(f"{args.graph}.{i}.cluster.itol.txt","w").write(samples2itol(list(cluster)))


def main_get_largest_cluster(args):
    graph = transmission_graph(args.graph)
    sys.stdout.write("%s\n" % ("\n".join(graph.clusters[0])))


def main_vcf2clusters(args):
    # vcf_file = sp.getoutput(path_command+args.vcf)
    # keepfile = sp.getoutput(path_command+args.keepfile)
    vcf_file = args.vcf
    keepfile = args.keepfile
    filecheck(vcf_file)
    filecheck(keepfile)
    vcf_class = vcf(vcf_file,keepfile=keepfile)
    vcf_class.get_clusters(cutoff=args.cutoff)


parser = argparse.ArgumentParser(
    description='Tree analysis', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser(
    'vcf2clusters', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--vcf', required=True)
parser_sub.add_argument('--cutoff', type=int, default=10)
parser_sub.add_argument('--keepfile', type=str)
parser_sub.set_defaults(func=main_vcf2clusters)

parser_sub = subparsers.add_parser(
    'stats', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('graph')
parser_sub.add_argument("--out", help = "output file name", dest = "stats_output_file", type = str)
parser_sub.add_argument("--itol_file", help = "name of itol file output", dest = "itol_file", type = str)
parser_sub.add_argument("--clust_min", help="minimum number of samples in cluster", dest="cluster_minimum", type=int, default=1)
parser_sub.set_defaults(func=main_stats)

parser_sub = subparsers.add_parser(
    'annotate', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('graph')
parser_sub.add_argument('csv')
parser_sub.add_argument('out')
parser_sub.set_defaults(func=main_add_meta)

parser_sub = subparsers.add_parser(
    'view', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('graph')
parser_sub.set_defaults(func=main_view)

parser_sub = subparsers.add_parser(
    'downsample', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('graph')
parser_sub.add_argument('--min-size',type=int,default=50)
parser_sub.set_defaults(func=main_downsample)

parser_sub = subparsers.add_parser(
    'largest_cluster', help='Calculate stats', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('graph')
parser_sub.set_defaults(func=main_get_largest_cluster)

args = parser.parse_args()
if vars(args) == {}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
