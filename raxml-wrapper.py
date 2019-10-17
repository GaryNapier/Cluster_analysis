import sys
import argparse
import subprocess
import random


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

def main(args):
	run_cmd("./raxml-ng --model GTR+G --msa %s --parse" % (args.msa))
	fine_grain_threads = None
	for l in open("%s.raxml.log" % args.msa):
		row = l.rstrip().split()
		if "Recommended number of threads" in l:
			fine_grain_threads = int(row[8])
	sys.stderr.write("Using %s threads for fine grain parallelisation\n" % fine_grain_threads)
	if fine_grain_threads<args.threads:
		course_grain_threads = int(args.threads//fine_grain_threads)
	else:
		course_grain_threads = 1
		fine_grain_threads = args.threads
	bs_trees_per_process = int(args.bs_trees//course_grain_threads)
	sys.stderr.write("Performing %s trees per bootstrap process\n" % bs_trees_per_process)
# 20 / 4 - 5
# 100/5 = 20
	with open("%s.parallel.sh" % args.msa,"w") as O:
		for i in range(int(args.starting_trees/2)):
			O.write("./raxml-ng --model GTR+G --msa %s.raxml.rba --search --threads %s --tree rand{1} --prefix %s.search.%s --seed %s\n" % (args.msa,fine_grain_threads,args.msa,i,random.randint(1,10000)))
		for i in range(int(args.starting_trees/2),args.starting_trees):
			O.write("./raxml-ng --model GTR+G --msa %s.raxml.rba --search --threads %s --tree pars{1} --prefix %s.search.%s --seed %s\n" % (args.msa,fine_grain_threads,args.msa,i,random.randint(1,10000)))
		for i in range(course_grain_threads):
			O.write("./raxml-ng --model GTR+G  --msa %s.raxml.rba --bootstrap --seed %s --bs-trees %s --prefix %s.bootstraps.%s --threads %s\n" % (args.msa,random.randint(1,10000),bs_trees_per_process, args.msa,i,fine_grain_threads))
		if args.bs_trees%bs_trees_per_process>0:
			O.write("./raxml-ng --bootstrap --msa %s.raxml.rba --seed %s --bs-trees %s --prefix %s.bootstraps.%s --threads %s\n" % (args.msa,random.randint(1,10000),(args.bs_trees%bs_trees_per_process), args.msa,(i+1),fine_grain_threads))
	run_cmd("cat %s.parallel.sh | parallel -j %s" % (args.msa,course_grain_threads))
	best_tree_lik = None
	best_tree = None
	for i in range(args.starting_trees):
		for l in open("%s.search.%s.raxml.log" % (args.msa,i)):
			if "Final LogLikelihood" in l:
				row = l.strip().split()
				lik = float(row[2])
				if best_tree_lik==None:
					best_tree_lik = lik
					best_tree = i
				elif lik>best_tree_lik:
					best_tree_lik = lik
					best_tree = i

	print("Best tree search: %s" % best_tree)
	run_cmd("cp %s.search.%s.raxml.bestTree %s.raxml.bestTree" % (args.msa,i,args.msa))
	run_cmd("cat %(msa)s.bootstraps.*.raxml.bootstraps > %(msa)s.bootstraps" % vars(args))
	run_cmd("./raxml-ng --support --tree %(msa)s.raxml.bestTree --bs-trees %(msa)s.bootstraps --prefix %(msa)s --threads 1 " % vars(args))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--msa',help='Prefix for files',required=True)
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.add_argument('--starting-trees',default=24, type=int, help='Number of starting trees')
parser.add_argument('--bs-trees',default=100, type=int, help='Number of starting trees')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
