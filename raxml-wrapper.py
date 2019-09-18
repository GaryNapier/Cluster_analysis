import sys
import argparse
import subprocess

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
	run_cmd("raxml-ng --model GTR+G --msa %s --parse" % (args.msa))
	fine_grain_threads = None
	for l in open("%s.raxml.log" % args.msa):
		row = l.rstrip().split()
		if "Recommended number of threads" in l:
			fine_grain_threads = int(row[8])
	sys.stderr.write("Using %s threads for fine grain parallelisation\n" % fine_grain_threads)
	if fine_grain_threads<args.threads:
		course_grain_threads = args.threads//fine_grain_threads
	else:
		course_grain_threads = 1
		fine_grain_threads = args.threads
	with open("%s.parallel.sh" % args.msa,"w") as O:
		for i in range(int(args.starting_trees/2)):
			O.write("raxml-ng --model GTR+G --msa %s --search --threads %s --tree rand{1} --prefix %s.search%s --seed $RANDOM\n" % (args.msa,fine_grain_threads,args.msa,i))
		for i in range(int(args.starting_trees/2),args.starting_trees):
			O.write("raxml-ng --model GTR+G --msa %s --search --threads %s --tree pars{1} --prefix %s.search%s --seed $RANDOM\n" % (args.msa,fine_grain_threads,args.msa,i))
	run_cmd("cat %s.parallel.sh | parallel -j %s" % (args.msa,course_grain_threads))
	run_cmd("grep Final %(msa)s*search*.log > %(msa)s.final_likelihoods.log" % vars(args))
	best_tree_lik = None
	best_tree = None
	for i in range(args.starting_trees):
		for l in open("%s.search%s.raxml.log" % (args.msa,i)):
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
	run_cmd("cp %s.search%s.raxml.bestTree %s.raxml.bestTree" % (args.msa,i,args.msa))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--msa',help='Prefix for files',required=True)
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.add_argument('--starting-trees',default=24, type=int, help='Number of starting trees')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
