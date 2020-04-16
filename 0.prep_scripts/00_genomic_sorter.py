
import argparse, pysam, gzip
parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True, help="index of genome version")
parser.add_argument("--window", type=int, required=True, help="the length of division")
parser.add_argument("--output", type=str, required=True, help="sorted data")

args = parser.parse_args()

in_file = open(args.input)
bin_size = args.window
out_file = open(args.output, "w")

for line in in_file.readlines():
	indi = line.rstrip().split("\t")
	chrom = indi[0]
	gene_size = int(indi[1])
	for i in range(1,int(gene_size)):
		gap = int(gene_size) - bin_size*i
		if gap >= 0:
			start = 1 + int(bin_size*(i-1))
			end = int(bin_size*i)
			out_file.write(chrom + "\t" + str(start) + "\t" + str(end) + "\n")
		else:
			start_resi = 1 + int(bin_size*(i-1))
			end_resi = int(gene_size)
			out_file.write(chrom + "\t" + str(start_resi) + "\t" + str(end_resi)+ "\n")
			break
