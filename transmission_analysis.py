#! /usr/bin/env python
import subprocess
import argparse
import numpy as np
import pandas as pd
import csv
import os
import sys
import subprocess as sp
import ete3
from ete3 import Tree, ClusterTree, TreeStyle, NodeStyle, AttrFace, ProfileFace, TextFace, faces
from ete3.treeview.faces import add_face_to_node

def run_cmd(cmd,verbose=1,target=None):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	if target and filecheck(target): return True
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/stdout","w")
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")
	else:
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
	stderr.close()
	if res!=0:
		print("Command Failed! Please Check!")
		exit(1)

def read_tree(tree_file):
	with open(tree_file, mode='r') as newick_file:
		t = Tree(newick_file.read())
	# Midpoint root the tree
	# Calculate the midpoint node
	md_pt = t.get_midpoint_outgroup()
	# and set it as tree outgroup
	t.set_outgroup(md_pt)
	# Write midpoint rooted as newick file
	with open(tree_file+".rooted", mode='w') as newick_file:
		newick_file.write(t.write())
	return t

def run(args):

    # ------------------
    # Setup/variables
    # ------------------
    path_command = "find -name "
    lineage = args.lineage
    dist_file = sp.getoutput(path_command+args.dist_file)
    mds_file = sp.getoutput(path_command+args.mds_file)
    tree_file = sp.getoutput(path_command+args.tree_file)
    lineage_file = sp.getoutput(path_command+args.lineage_file)

    transmission_threshold = args.transmission_threshold
    bootstrap_cutoff = args.bootstrap_cutoff
    node_cutoff = args.node_cutoff

    # -------------------------
    # Pull transmission samples
    # -------------------------
    # mds file to get sample names
    with open(mds_file, mode='r') as mds_file:
    	mds = list(csv.reader(mds_file, delimiter = "\t"))

    # Get sample names from first column of mds file
    samples = list()
    for line in mds:
        line = ' '.join(line) # Collapse spaces
        line = [i.strip() for i in line.split(' ')] # Make into individual strings
        line = [i for i in line if i] # Remove empty strings(formerly spaces)
        samples.append(line[0]) # Save first entry (col)
    samples = samples[1:] # Remove column header

    # Read in distance file
    dist = pd.read_csv(dist_file, sep='\t', names=samples)
    dist.index = samples

    transmission_samples = []
    for row_ind, row in enumerate(dist.index.values):
      for col_ind, col in enumerate(dist.columns.values):
        if col >= row: # Lower triangle
          continue
        #if subdist_within_table.loc[row, col] <= outliers_within_stat[0] or subdist_within_table.loc[row, col] >= outliers_within_stat[1]:
        if dist.loc[row, col] <= transmission_threshold:
          transmission_samples.append([row, col])

    # Get unique
    transmission_samples = list(set([item for sublist in transmission_samples for item in sublist]))

    # print out (use "> file.txt" to save)
    for i in transmission_samples:
        print(i)

    # ----------
    # Get stats
    # ----------
    # Read in list of samples in lineage/clade
    lineage_samples =  pd.read_csv(lineage_file, sep='\t', names = None, header = None)
    # Read in metadata for drug resistance
    meta = pd.read_csv("all_lins/csv/tb10k.meta.csv", sep=',', na_filter = False)

    # Join transmission samples to drug resistance data
    trans_dr = pd.merge(pd.DataFrame(transmission_samples, columns = ['id']), meta[['id', 'drtype']], on='id')
    # Count drug-resistant samples
    drug_resistance_pivot = pd.DataFrame(pd.pivot_table(trans_dr,index=["drtype"], values = ["drtype"], aggfunc = [len])).reset_index()

    print(drug_resistance_pivot['drtype' == 'Drug-resistant'])

    # Put together table

    # trans_df = pd.DataFrame({"Lin/cluster": lineage,
    # "n in lineage/cluster": len(lineage_samples),
    # "n samples in transmission": len(transmission_samples),
    # "n drug-resistant": })

    # print(df_test)


    #
    # len(transmission_samples)/len(lineage_samples)

    # # --------------
    # # Tree analysis
    # # --------------
    # # Traverse tree and get samples belonging to clade greater than 'node_cutoff'
    # t = read_tree(tree_file)
    # # Traverse tree to get min cluster size
    # leaf_list = []
    # for node in t.traverse("preorder"):
    #     # Get leaves
    #     if len(node) > node_cutoff and node.support >= bootstrap_cutoff:
    #         leaf_list.append([n.name for n in node.get_leaves()])

    # # -------------
    # # Put together
    # # -------------
    # final_transmission_samps = []
    # for trans_samp_i in transmission_samples:
    #     for leaf_list_i in leaf_list:
    #         # print(trans_samp_i)
    #         # print("------------")
    #         # print(leaf_list_i)
    #         # print("------------ ")
    #         # print(trans_samp_i in leaf_list_i)
    #         if trans_samp_i in leaf_list_i:
    #             final_transmission_samps.append(trans_samp_i)

    # # Save unique
    # final_transmission_samps = list(set(final_transmission_samps))
    # # print out (use "> file.txt" to save)
    # for i in final_transmission_samps:
    #     print(i)

def main():
    parser=argparse.ArgumentParser(description="Find tree clusters")
    parser.add_argument("-l", help = "lineage", dest = "lineage", type = str, required = True)
    parser.add_argument("-dm",help="distance matrix file", dest="dist_file", type = str, required=True)
    parser.add_argument("-mds",help="mds file", dest="mds_file", type=str, required=True)
    parser.add_argument("-tf",help="tree file", dest="tree_file", type = str)
    parser.add_argument("-lf", help = "lineage samples file", dest = "lineage_file", type = str, required = True)
    parser.add_argument("-tt",help="transmission_threshold", dest="transmission_threshold", type = int, default=10)
    parser.add_argument("-ct", help="minimum cluster size", dest = "node_cutoff", type = int, default = 5)
    parser.add_argument("-bs", help = "minimum bootstrap", dest = "bootstrap_cutoff", type = int, default = 75)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
	main()

#datetime.datetime(2020, 2, 13, 19, 34, 34, 740336)
