import sys
import subprocess
import argparse
import random
import os
rand_generator = random.SystemRandom()

# version = 1.0

def write_set_GT_script():
    open("setGT.py","w").write("""#! /usr/bin/env python
import sys
from tqdm import tqdm
import argparse
import re
def main(args):
    ad_cutoff = args.fraction
    for l in tqdm(sys.stdin):
        something_changed = False
        row = l.strip().replace("|","/").split()
        if l[0]=="#":
            sys.stdout.write(l)
            continue
        alleles = [row[3]]+row[4].split(",")
        if len(alleles)>9: continue

        if "*" in row[4]:
            something_changed = True
            for i in range(9,len(row)):
                if row[i][0]=="." or row[i][2]==".": continue
                if alleles[int(row[i][0])]=="*" or alleles[int(row[i][2])]=="*":
                    tmp = list(row[i])
                    tmp[0] = "."
                    tmp[2] = "."
                    row[i]="".join(tmp)

        uniq_mixed_genotypes = set([x for x in re.findall("[0-9]/[0-9]",l) if x[0]!=x[2]])
        if len(uniq_mixed_genotypes)>1:
            for i in range(9,len(row)):
                if row[i][:3] not in uniq_mixed_genotypes: continue
                fmt = row[i].split(":")
                ad = [int(x) for x in fmt[1].split(",")]
                total_ad = sum(ad)
                if total_ad==0:continue
                adf = [ad[j]/total_ad for j in range(len(ad))]
                if max(adf)>=ad_cutoff:
                    new_gt = adf.index(max(adf))
                    fmt[0] = f"{new_gt}/{new_gt}"
                    something_changed = True
                    row[i] = ":".join(fmt)
        if something_changed:
            sys.stdout.write("\t".join(row)+"\n")
        else:
            sys.stdout.write(l)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fraction',default=0.7,type=float,help='Fraction of coverage to assign major')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)""")


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


def log(msg, ext=False):
    sys.stderr.write("\n"+str(msg)+"\n")
    if ext:
        exit(1)


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


def nofile(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isfile(filename):
        return True
    else:
        return False


def nofolder(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isdir(filename):
        return True
    else:
        return False


class vcf_class:
    def __init__(self, filename, threads=4):
        self.samples = []
        self.filename = filename
        self.threads = threads
        self.prefix = filename[:-7]
        if nofile(filename+".csi"):
            run_cmd("bcftools index  %(filename)s" % vars(self))
        self.temp_file = get_random_file()
        run_cmd("bcftools query -l %(filename)s > %(temp_file)s" % vars(self))
        for l in open(self.temp_file):
            self.samples.append(l.rstrip())
        os.remove(self.temp_file)

    def vcf_to_fasta(self, ref_file, threads=4, chunk_size=50000):
        self.ref_file = ref_file
        self.chunk_size = chunk_size
        self.cmd_split_chr = "bedtools makewindows -g %(ref_file)s.fai -w %(chunk_size)s -s  %(chunk_size)s | awk '{print $1\":\"$2\"-\"$3}'" % vars(
            self)
        self.tmp_file = "%s.tmp.txt" % self.prefix
        self.threads = threads
        cmd = "%(cmd_split_chr)s | parallel -j %(threads)s \"bcftools view  %(filename)s -r {} -Ou | bcftools query -f '%%POS[\\t%%IUPACGT]\\n' | sed 's/\*[\/|]\*/\.\/\./g' |  datamash transpose > %(prefix)s.{}.tmp.txt\"" % vars(
            self)
        run_cmd(cmd)
        cmd = "paste `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'` > %(tmp_file)s" % vars(
            self)
        run_cmd(cmd)
        cmd = "rm `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'`" % vars(
            self)
        run_cmd(cmd)
        with open(self.prefix+".snps.fa", "w") as O:
            for i, l in enumerate(open(self.tmp_file)):
                row = l.rstrip().split()
                if i == 0:
                    continue
                s = self.samples[i-1]
                seq = "".join(row).replace("./.", "N")
                O.write(">%s\n%s\n" % (s, seq))
        run_cmd("rm %s" % self.tmp_file)

    def vcf_to_matrix(self,):
        self.matrix_file = self.prefix+".mat"
        self.binary_matrix_file = self.prefix+".mat.bin"
        O = open(self.matrix_file, "w").write(
            "chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
        run_cmd(
            "bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./N/g' >> %(matrix_file)s" % vars(self))
        O = open(self.binary_matrix_file, "w").write(
            "chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
        run_cmd(
            "bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%GT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./N/g' | sed 's/1\/1/1/g' | sed 's/0\/0/0/g' >> %(binary_matrix_file)s" % vars(self))

    def get_plink_dist(self, pca=True, mds=True):
        self.tempfile = get_random_file(extension=".vcf")
        run_cmd("bcftools view %(filename)s > %(tempfile)s" % vars(self))
        run_cmd("plink --vcf %(tempfile)s --distance square --double-id --allow-extra-chr --out %(prefix)s.temp" % vars(self))
        O = open("%(prefix)s.dist" % vars(self), "w")
        dists = []
        for l in open("%(prefix)s.temp.dist" % vars(self)):
            row = [float(d)/2 for d in l.rstrip().split()]
            O.write("%s\n" % "\t".join([str(x) for x in row]))
            dists.append(row)
        O.close()
        if pca:
            run_cmd(
                "plink --vcf %(tempfile)s --pca --double-id --allow-extra-chr --out %(prefix)s.pca" % vars(self))
        if mds:
            run_cmd("plink --vcf %(tempfile)s --mds-plot 10 eigendecomp --cluster --double-id --allow-extra-chr --out %(prefix)s.pca" % vars(self))

        run_cmd(
            "rm %(tempfile)s* %(prefix)s.temp* %(prefix)s.pca.log %(prefix)s.pca.nosex" % vars(self))
        return dists


def main(args):
    write_set_GT_script()
    FAILED_SAMPLES = open("%s.failed_samples.log" % args.prefix, "w")
    params = {"threads": args.threads, "prefix": args.prefix, "ref": args.ref}
    params["map_file"] = "%s.map" % (args.prefix)
    if args.redo:
        samples = [x.rstrip() for x in open(params["map_file"]).readlines()]
    else:
        with open(params["map_file"], "w") as O:

            # Set up list to hold sample names
            samples = []

            # Loop through sample-file and do (1) append samples to list, (2) write sample to map file and (3) check for VCF index
            for line in open(args.sample_file):
                sample = line.rstrip()
                if args.ignore_missing and nofile("%s/%s%s" % (args.vcf_dir, sample, args.vcf_extension)):
                    continue
                if args.no_validate or args.redo:
                    pass
                else:
                    if not os.path.isfile(f"{args.vcf_dir}/{sample}{args.vcf_extension}.validated"):
                        FAILED_SAMPLES.write(sample+"\n")
                        continue

                samples.append(sample)
                O.write("%s\t%s/%s%s\n" %
                        (sample, args.vcf_dir, sample, args.vcf_extension))
                if nofile("%s/%s%s.tbi" % (args.vcf_dir, sample, args.vcf_extension)):
                    run_cmd("bcftools index --tbi %s/%s%s" %
                            (args.vcf_dir, sample, args.vcf_extension))
    stages = {"dbimport": 1, "genotype": 2,
              "filtering": 3, "fasta": 4, "matrix": 5, "pca": 6}
    # Create .dict file (GATK fasta index) has been created for the reference
    if nofile("%s.dict" % args.ref.replace(".fasta", "").replace(".fa", "")):
        run_cmd("gatk CreateSequenceDictionary -R %(ref)s" % params)
    # Create .fai file (SAMtools fasta index) has been created for the reference
    if nofile("%s.fai" % args.ref.replace(".fasta", "").replace(".fa", "")):
        run_cmd("samtools faidx %(ref)s" % params)
    if nofolder("%(prefix)s_genomics_db" % params) or stages[args.redo] <= 1:
        run_cmd("gatk GenomicsDBImport --genomicsdb-workspace-path %(prefix)s_genomics_db -L Chromosome --sample-name-map %(map_file)s --reader-threads %(threads)s --batch-size 500" % params, verbose=2)
    if nofile("%(prefix)s.raw.vcf.gz" % params) or stages[args.redo] <= 2:
        run_cmd("gatk --java-options \"-Xmx40g\" GenotypeGVCFs -R %(ref)s -V gendb://%(prefix)s_genomics_db -O %(prefix)s.raw.vcf.gz" % params, verbose=2)
    if nofile("%(prefix)s.filt.vcf.gz" % params) or stages[args.redo] <= 3:
        run_cmd("bcftools view -V indels %(prefix)s.raw.vcf.gz | python setGT.py | bcftools view -a | bcftools filter -e 'GT=\"het\"' -S . | bcftools view -i 'F_PASS(GT!=\"mis\")>0.9' | bcftools view -c 1 | bcftools view -a -Oz -o %(prefix)s.filt.vcf.gz" % params)
    vcf = vcf_class("%s.filt.vcf.gz" % (args.prefix))
    if nofile("%(prefix)s.snps.fa" % vars(vcf)) or stages[args.redo] <= 4:
        vcf.vcf_to_fasta(args.ref)
    if nofile("%(prefix)s.mat" % vars(vcf)) or stages[args.redo] <= 5:
        vcf.vcf_to_matrix()
    if nofile("%(prefix)s.pca.eigenvec" % vars(vcf)) or stages[args.redo] <= 6:
        vcf.get_plink_dist()


parser = argparse.ArgumentParser(
    description='TBProfiler pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--sample-file', help='sample file', required=True)
parser.add_argument('--prefix', help='Prefix for files', required=True)
parser.add_argument('--ref', help='reference file', required=True)
parser.add_argument('--vcf-dir', default="./vcf/",
                    type=str, help='VCF firectory')
parser.add_argument('--vcf-extension', default=".gatk.vcf.gz",
                    type=str, help='VCF extension')
parser.add_argument('--threads', default=4, type=int,
                    help='Number of threads for parallel operations')
parser.add_argument('--ignore-missing', action="store_true",
                    help='If this option is set, missing samples are ignored')
parser.add_argument('--redo', type=str,
                    choices=["dbimport", "genotype", "filtering", "fasta", "matrix", "pca"])
parser.add_argument('--no-validate', action="store_true",)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
