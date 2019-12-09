#! /usr/bin/env python
import subprocess
import argparse
import numpy as np
import ete3
import csv
from skbio import DistanceMatrix
from skbio.stats.distance import anosim
import pandas as pd
from statsmodels import robust
from itertools import chain
import sys
from numpy import mean
from scipy import stats
from numpy import var
from math import sqrt
import seaborn as sns
import matplotlib.pyplot as plt
from ete3 import Tree, ClusterTree, TreeStyle, NodeStyle, AttrFace, ProfileFace, TextFace, faces
from ete3.treeview.faces import add_face_to_node
from itertools import chain
import scipy
from scipy.stats import iqr

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

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=10)
        faces.add_face_to_node(N, node, 0, position="aligned")

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

def get_leaves_in_clusts(tree_file, node_cutoff, bootstrap_cutoff):

	# Read tree file
	t = read_tree(tree_file)
	# Parse tree and get lists of samples for each clade ('leaves')
	# Clades have to be > 50 bootstrap an min of 20 samples (by default)
	leaves = list()
	for node in t.traverse("preorder"):
		# Do some analysis on node 
		if len(node) > node_cutoff and node.support >= bootstrap_cutoff:
			# print(node) 
			# print(node.name)
			# print(node.support)
			leaves.append([n.name for n in node.get_leaves()])
	return leaves

def do_anosim(leaf_list, all_other_samps, dist, test_stat, perm=99):
	# http://scikit-bio.org/docs/0.5.5/generated/skbio.stats.distance.anosim.html#skbio.stats.distance.anosim
	# Convert pandas distance matrix to DistanceMatrix object
	dist = DistanceMatrix(dist)
	# Need to create two groups for anosim - samples in leaf list and all other samples 
	grouping = [["Cluster"]*len(leaf_list), ["Everything_else"]*len(all_other_samps)]
	grouping = list(chain.from_iterable(grouping)) # Unlist groups into one list
	# Do anosim
	anosim_res = anosim(dist, grouping, permutations = perm)
	return anosim_res
	# return clust_samps_list

def make_itol(samples, file_nm, lineage):
	col = ['#ff0000']*len(samples)    
	df = pd.DataFrame(list(zip(samples, col)))  
	df.to_csv(file_nm, sep="\t", quoting=csv.QUOTE_NONE, index = False, header = False)   

	run_cmd("cp itol_templates/dataset_color_strip_template.txt lin_%s/results/%s" % (lineage, file_nm))
	run_cmd("cat %s >> lin_%s/results/%s" % (file_nm, lineage, file_nm) )

def do_clustering(tree_file, matrix_file, samples):
	# Load matrix file. Must be in same format as here: https://github.com/etetoolkit/ete/blob/master/examples/clustering/diauxic.array
	# i.e. the sample names are the first column and the column headers are the first row. 
	# The samples and column headers are therefore PART of the data, not just indices, so have to do this convoluted operation:

	# Read in the matrix without headers (assuming there are no headers), but add headers as samples
	matrix = pd.read_csv(matrix_file, sep='\t', header = None, names=samples) 
	# Add samples as row indices/row names
	matrix.index = samples
	# Write this to a temporary file as a tab-separated text file
	matrix.to_csv("matrix_file_tmp.txt", sep = "\t")  
	# Need a string of text to now fill the first position (because the headers and row names are now part of the data)
	# Needs to be '#NAMES'
	# Headers and row names should now be part of data
	run_cmd("sed -i '1s/^/#NAMES/' matrix_file_tmp.txt")

	# Read tree
	# tree = read_tree(tree_file)
	# Loads tree and array
	t = ClusterTree(tree_file, "matrix_file_tmp.txt")
	# Remove temp file
	run_cmd("rm matrix_file_tmp.txt")	
	# nodes are linked to the array table
	array = t.arraytable
	# Calculates some stats on the matrix. Needed to establish the colour gradients.
	matrix_dist = [i for r in range(len(array.matrix)) for i in array.matrix[r] if np.isfinite(i)]
	matrix_max = np.max(matrix_dist)
	matrix_min = np.min(matrix_dist)
	matrix_avg = matrix_min+((matrix_max-matrix_min)/2)
	# Creates a profile face that will represent node's profile as as heatmap
	profileFace  = ProfileFace(matrix_max, matrix_min, matrix_avg, 200, 14, "heatmap")
	cbarsFace = ProfileFace(matrix_max,matrix_min,matrix_avg,200,70,"cbars")
	nameFace = AttrFace("name", fsize=8)
	# Creates my own layout function that uses previous faces
	def mylayout(node):
	    # If node is a leaf
	    if node.is_leaf():
	        # And a line profile
	        add_face_to_node(profileFace, node, 0, aligned=True)
	        node.img_style["size"]=0
	        add_face_to_node(nameFace, node, 1, aligned=True)

	    # If node is internal
	    else:
	        # If silhouette is good, creates a green bubble
	        if node.silhouette>0:
	            validationFace = TextFace("Silh=%0.2f" %node.silhouette,
	                                      "Verdana", 10, "#056600")
	            node.img_style["fgcolor"]="#056600"
	        # Otherwise, use red bubbles
	        else:
	            validationFace = TextFace("Silh=%0.2f" %node.silhouette,
	                                      "Verdana", 10, "#940000")
	            node.img_style["fgcolor"]="#940000"

	        # Sets node size proportional to the silhouette value.
	        node.img_style["shape"]="sphere"
	        if node.silhouette<=1 and node.silhouette>=-1:
	            node.img_style["size"]= 15+int((abs(node.silhouette)*10)**2)

	        # If node is very internal, draw also a bar diagram
	        # with the average expression of the partition
	        add_face_to_node(validationFace, node, 0)
	        if len(node)>100:
	            add_face_to_node(cbarsFace, node, 1)

	# Use my layout to visualize the tree
	ts = TreeStyle()
	ts.layout_fn = mylayout
	t.show(tree_style=ts)
	return ts

# Function for effect size https://machinelearningmastery.com/effect-size-measures-in-python/
def cohend(d1, d2):
	# calculate the size of samples
	n1, n2 = len(d1), len(d2)
	# calculate the variance of the samples
	s1, s2 = var(d1, ddof=1), var(d2, ddof=1)
	# calculate the pooled standard deviation
	s = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
	# calculate the means of the samples
	u1, u2 = mean(d1), mean(d2)
	# calculate the effect size
	return (u1 - u2) / s


def run(args):
	lineage = args.input # these match the "dest": dest="input"
	# tree_file = "lin_%s/results/lin_%s.iqtree.treefile" % (lineage, lineage)
	tree_file = "lin_%s/results/lin_%s.iqtree.treefile" % (lineage, lineage)
	mds_file = "lin_%s/results/lin_%s_results.filt.pca.mds" % (lineage, lineage)
	dist_file = "lin_%s/results/lin_%s_results.filt.dist" % (lineage, lineage)
	node_cutoff = args.node_cutoff 
	bootstrap_cutoff = args.bootstrap_cutoff
	# "R" value - 0.75 < R < 1 - highly different; 0.5 < R < 0.75 -  different; 0.25 < R < 0.5 -  different with some overlap; 0.1 < R < 0.25 - high overlap; R < 0.1 - similar
	# See https://www.researchgate.net/post/Which_R_value_is_considered_to_show_a_strong_difference_between_the_groups_in_ANOSIM
	test_stat = args.test_stat 
	do_t_test = args.do_t_test
	show_plots = args.show_plots
	do_anosim_analysis = args.do_anosim_analysis
	metadata = args.metadata
	do_clustering_analysis = args.do_clustering_analysis
	if metadata == "":
		do_clustering_analysis = 0


	off_diag = -1 # Have to specify off-diagonal with '-1'. 

	if test_stat < -1 or test_stat > 1:
		sys.exit("Please choose R value (-r) between -1 and 1")

	# Load files:

	# mds file to get samples
	with open(mds_file, mode='r') as mds_file:
		mds = list(csv.reader(mds_file, delimiter = "\t"))

	# Get samples from first column of mds file
	samples = list()
	for line in mds:
		line = ' '.join(line) # Collapse spaces
		line = [i.strip() for i in line.split(' ')] # Make into individual strings
		line = [i for i in line if i] # Remove empty strings(formerly spaces)
		samples.append(line[0]) # Save first entry (col)
	
	samples = samples[1:] # Remove column header

	# Read in lineages table
	lineages_master = pd.read_csv("all_lins/txt/lineages_master.txt", sep='\t') 
	# Get unique sublineages
	sublins = lineages_master.loc[lineages_master["first_dec"] == lineage, ].sub_lineage.unique()
	# Subset lins table to only include those in the tree/from mds file
	lineages_master = lineages_master[lineages_master["sample"].isin(samples)]

	# Read in distance matrix setting and set header and index as sample IDs. 
	dist = pd.read_csv(dist_file, sep='\t', names=samples) 
	dist.index = samples

	# Get tree 
	t = read_tree(tree_file)

	# Traverse tree, subset distance matrix according to leaf names and do tests on distance matrix

	# Setup
	clust_samps_list = []
	clust_num = 0 # Loop counter
	ts = TreeStyle() # Tree display setup
	ts.mode = "c" # draw tree in circular mode
	ts.show_leaf_name = False 
	ts.show_branch_length = True
	ts.show_branch_support = True
	ts.title.add_face(TextFace("Lineage %s" % lineage, fsize=len(t)/4), column=0) # Add title to tree
	ts.layout_fn = layout
	colours = sns.hls_palette(len(sublins), l = 0.9, s = 0.7).as_hex()
	

	# Traverse tree loop 
	for node in t.traverse("preorder"):

		# Get leaves
		if len(node) > node_cutoff and node.support >= bootstrap_cutoff:
			leaf_list = [n.name for n in node.get_leaves()]
			print("-------------")
			print(leaf_list)
			print("-------------")
		else:
			continue

		# Get samples that are NOT in clade of interest
		all_other_samps = list(set(samples).difference(set(leaf_list)))

		# # Subset raw data
		subdist_within_table = dist.loc[leaf_list, leaf_list]
		subdist_within_lt = np.tril(subdist_within_table, k=off_diag) # Get lower triangle
		subdist_within_lt = subdist_within_lt[np.tril_indices(len(subdist_within_lt), k=off_diag)] # Actually get lower triangle. Pathetic.

		subdist_btwn_table = dist.loc[leaf_list, all_other_samps] # Save for later
		subdist_btwn = np.array(dist.loc[leaf_list, all_other_samps])
		subdist_btwn = subdist_btwn.reshape(subdist_btwn.shape[0]*subdist_btwn.shape[1], ) # Turn to vector


		# Do t-tests
		if do_t_test == 1:

			t_test = stats.ttest_ind(subdist_within_lt, subdist_btwn, equal_var=False)
			effect_sz = cohend(subdist_within_lt, subdist_btwn)

			if t_test.pvalue <= 0.05:
				clust_num += 1

				# Stats stuff: 

				within_dict = {"Within":"Within", "n_samps":len(leaf_list), "n_dist":len(subdist_within_lt), "Min": np.min(subdist_within_lt), "Max":np.max(subdist_within_lt), \
				 "Sum":np.sum(subdist_within_lt), "Mean": np.mean(subdist_within_lt), "Median": np.median(subdist_within_lt), "sd": np.std(subdist_within_lt, ddof=1), "mad": robust.mad(subdist_within_lt)}
				btwn_dict = {"Between":"Between", "n_samps":len(all_other_samps), "n_dist":len(subdist_btwn), "Min": np.min(subdist_btwn), "Max":np.max(subdist_btwn), "Sum":np.sum(subdist_btwn), \
				 "Mean": np.mean(subdist_btwn), "median": np.median(subdist_btwn), "sd": np.std(subdist_btwn, ddof=1), "mad": robust.mad(subdist_btwn)}

				# Get outliers
				q1_within = np.percentile(subdist_within_lt, 25)
				q3_within = np.percentile(subdist_within_lt, 75)
				iqr = abs(q1_within - q3_within) # interquartile range (Q1 - Q3)
				outliers_within_stat = [q1_within - (1.5*iqr), q3_within+(1.5*iqr)]

				q1_btwn = np.percentile(subdist_btwn, 25)
				q3_btwn = np.percentile(subdist_btwn, 75)
				iqr = abs(q1_btwn - q3_btwn)
				outliers_btwn_stat = [q1_btwn - (1.5*iqr), q3_btwn+(1.5*iqr)]

				print("Outliers within: ", outliers_within_stat)
				print("Outliers between: ", outliers_btwn_stat)

				# print(subdist_within_table)

				outliers_within = []
				for row_ind, row in enumerate(subdist_within_table.index.values):
					for col_ind, col in enumerate(subdist_within_table.columns.values):
						if col >= row: # Lower triangle
							continue
						if subdist_within_table.loc[row, col] <= outliers_within_stat[0] or subdist_within_table.loc[row, col] >= outliers_within_stat[1]:
							outliers_within.append([row, col])

				outliers_within = list(set(list(chain.from_iterable(outliers_within)))) # Unlist and get unique
				print("Outliers within, clust %s: " % clust_num, outliers_within)

				outliers_btwn = []
				for row_ind, row in enumerate(subdist_btwn_table.index.values):
					for col_ind, col in enumerate(subdist_btwn_table.columns.values):
						if subdist_btwn_table.loc[row, col] <= outliers_btwn_stat[0] or subdist_btwn_table.loc[row, col] >= outliers_btwn_stat[1]:
							outliers_btwn.append([row, col])

				outliers_btwn = list(set(list(chain.from_iterable(outliers_btwn))))
				print("Outliers between, clust %s: " % clust_num, outliers_btwn)				


				print("Clust num: ", clust_num)
				print(t_test)
				print("Effect size (d):", effect_sz)
				print(within_dict)
				print(btwn_dict)
				print("Diff: ", abs(within_dict["Sum"] - btwn_dict["Sum"]))

				# Do boxplots
				if show_plots == 1:
					bxplt = sns.boxplot(data=[subdist_within_lt, subdist_btwn]) #.set_ylim([0, np.max(dist)])
					bxplt.set_title("Clust"+str(clust_num)) # set(xlabel = ["Within", "Between"]) #, x = "Group", y = "Distance", order = ["Within", "Between"])
					bxplt.set_ylim([-100, dist.max().max()+0.1*dist.max().max()])	
					# fig = img.get_figure()
					# fig.savefig()
					# fig.clf()
					plt.show()

				# Tree stuff:

				# l_width = abs(t_test.statistic) * 0.05 # Set thickness of line based on test stat degree
				l_width = abs(effect_sz) * 10 # Set thickness of line based on effect size
				colour = "#ff0000"

				style = NodeStyle()
				style["fgcolor"] = "#0f0f0f"
				style["size"] = 0
				style["vt_line_color"] = colour
				style["hz_line_color"] = colour
				# style["vt_line_width"] = l_width
				style["hz_line_width"] = l_width
				style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
				style["hz_line_type"] = 0
				node.set_style(style)
				clust_num_annotation = TextFace("\n"+str(clust_num))
				node.add_face(clust_num_annotation, column=1, position = "branch-bottom")

				# make_itol(leaf_list, "clust_"+str(clust_num)+".txt", lineage)

			# else:

				# style = NodeStyle()
				# style["fgcolor"] = "#0f0f0f"
				# style["size"] = 0
				# style["vt_line_color"] = "blue"
				# style["hz_line_color"] = "blue"
				# style["vt_line_width"] = 1
				# style["hz_line_width"] = 1
				# style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
				# style["hz_line_type"] = 0
				# node.set_style(style)			

		# Do anosim:
		if do_anosim_analysis:
			anosim_res = do_anosim(leaf_list, all_other_samps, dist, test_stat)
			# Save samples 
			# if < 0.05 p-val and > 0.5 # R test stat ("Different groups" - https://www.researchgate.net/post/Which_R_value_is_considered_to_show_a_strong_difference_between_the_groups_in_ANOSIM)
			if anosim_res.loc["p-value"] < 0.05 and anosim_res.loc["test statistic"] >= test_stat:
				clust_samps_list.append(leaf_list)
				print("-------------")
				print(anosim_res)

	for i in clust_samps_list:
		print(i)
		print("-------------")		

	# Clustering
	# Nb, as long as the read_tree() function has run on the tree_file, then the .rooted file should be available. 
	# get_leaves_in_clusts() runs the read_tree() function
	if do_clustering_analysis == 1:
		do_clustering(tree_file+".rooted", metadata, samples)

	# Show tree in GUI window:
	# Add colour to last common ancestor of each sublin
	for i, sublin in enumerate(sublins):
		samps_sublin = lineages_master.loc[lineages_master["sub_lineage"] == sublin, "sample"] # Subset samples by sublineage
		nst = NodeStyle()
		nst["bgcolor"] = colours[i]
		t.get_common_ancestor(list(samps_sublin)).set_style(nst)

	t.show(tree_style=ts)


def main():
	parser=argparse.ArgumentParser(description="Find tree clusters")
	parser.add_argument("-l",help="Lineage to analyse" ,dest="input", type=str, required=True)
	# parser.add_argument("-out",help="output filename" ,dest="output", type=str, required=True)
	parser.add_argument("-s",help="Min number of samples in clusters" ,dest="node_cutoff", type=int, default=20)
	parser.add_argument("-b",help="Min bootstrap for clusters" ,dest="bootstrap_cutoff", type=int, default=51)
	parser.add_argument("-r",help="R value cutoff (-1<0<1) for anosim", dest="test_stat", type=float, default=0.5)
	parser.add_argument("-t",help="0 = don't do t-test, 1 = do t-test, default = 0", dest="do_t_test", type=int, default=0)
	parser.add_argument("-p", help="0 = don't show plots, 1 = show plots, default = 0", dest="show_plots", type=int, default=0)
	parser.add_argument("-a", help="0 = don't do anosim analysis, 1 = do anosim analysis", dest="do_anosim_analysis", type=int, default=0)
	parser.add_argument("-c", help="0 = don't do cluster analysis, 1 = do clustering analysis", dest="do_clustering_analysis", type=int, default=0)
	parser.add_argument("-m", help="Metadata file for clustering analysis", dest="metadata", type=str, default="")
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()


# test_samps = ['ERR046997', 'ERR072039', 'ERR1034655', 'ERR1034661']
	# all_samps = test_samps + ['SRR058419', 'SRR058415', 'SRR057733', 'SRR058370']
	# test_samps = ['SRR058417', 'SRR057731', 'SRR058373', 'SRR058116', 'SRR058399']
	# all_samps  =['SRR058417', 'SRR057731', 'SRR058373', 'SRR058116', 'SRR058399', 'SRR058081', 'SRR058375', 'SRR057619', 'SRR058372', 'SRR057770', 'SRR058420']
	# test_samps = ['ERR046997', 'ERR072039', 'ERR2517277', 'ERR551657', 'ERR551574', 'ERR1034655', 'ERR2516176', 'ERR1034661']
	# all_samps = ['ERR046997', 'ERR072039', 'ERR2517277', 'ERR551657', 'ERR551574', 'ERR1034655', 'ERR2516176', 'ERR1034661', 'ERR2512852', 'ERR2517447', 'ERR2510830', 'ERR553390', 'ERR2514189', 'SRR2100564', 'ERR2517392', 'ERR2517428']
