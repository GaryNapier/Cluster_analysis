import csv
import argparse
import seaborn as sns

def main(args):
	meta = {}

	for row in csv.DictReader(open(args.csv)):
		test = False
		for lin in args.lineages.split(","):
			if lin in row["sublineage"]:
				test = True
		if test:
			meta[row["id"]] = row["sublineage"]

	values = sorted(list(set(meta.values())))
	nvalues = len(values)
	colours = {val:c for val,c in zip(values,sns.color_palette("Set2",nvalues).as_hex())}
	itol_text = """DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL\tLineage
COLOR\t#ff0000

LEGEND_TITLE\tLineage
LEGEND_SHAPES\t%s
LEGEND_COLORS\t%s
LEGEND_LABELS\t%s

DATA
""" % (
	"\t".join(["1" for _ in values]),
	"\t".join([c for c in colours.values()]),
	"\t".join([l for l in values])

)
	for s in meta:
		itol_text = itol_text + "%s\t%s\n" % (s,colours[meta[s]])

	print(itol_text)


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--csv',type=str,help='CSV with lineages',required=True)
parser.add_argument('--lineages',type=str,help='lineages you want to select',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
