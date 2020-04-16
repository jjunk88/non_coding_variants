
import numpy as np
import pandas as pd
import argparse
import math
from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu
from toolz import unique

parser = argparse.ArgumentParser()
parser.add_argument("--mut_metrics", type=str, required=True, help="your mutation list")
parser.add_argument("--exp_metrics", type=str, required=True, help="total expression array")
parser.add_argument("--ref", type=str, required=True, help="coordinates for genome")
parser.add_argument("--bin_size", type=int, required=True, help="window size for selecting target genes")
parser.add_argument("--output", type=str, required=True, help="name of output")

args = parser.parse_args()

dat_mut = open(args.mut_metrics)
dat_exp = open(args.exp_metrics)
out_file = open(args.output,'w')

window = args.bin_size

exp_array = pd.read_csv(dat_exp,sep="\t")
samples_list = []
for col in exp_array.columns:
	samples_list.append(col)
samples_list = samples_list[1:]

def average(list):
	return (sum(list) / len(list))


#out_file.write("chr" + "\t" +  "start" + "\t" +  "end" + "\t" +  "target" + "\t"+ "log2_fold_change" + "\t" + "p_value" + "\t" +  "samples_ID" + "\n")

for line in dat_mut:
	indi = line.rstrip().split("\t")
	chromosome = indi[0]
	start = indi[1]
	start_bin = int(start) - window
	end = indi[2]
	end_bin = int(end) + window
	samples_MT = indi[3]
	samples_MT_ind = samples_MT.split(",")
	samples_WT = list(set(samples_list) - set(samples_MT_ind))
	samples_WT_srt = ",".join(samples_WT)
	samples_WT_ind = samples_WT_srt.split(",")
	dat_ref = open(args.ref)
	for line_ref in dat_ref:
		ind = line_ref.rstrip().split("\t")
		ref_chr = ind[0]
		ref_start = ind[1]
		ref_end = ind[2]
		gene_list = ind[3]
		gene_list_sep = gene_list.split(",")
		if chromosome == ref_chr:
			if int(start_bin) < int(ref_start) and int(ref_end) < int(end_bin):
				for i in range(len(gene_list.split(","))):
					target = gene_list_sep[i]
					is_target = exp_array[exp_array['gene_symbol'] == target]
					if len(is_target) == 1:
						energ_MT = []
						energ_WT = []
						for k in range(1,227):
							samples_array = is_target.columns[k]
							for j in range(len(samples_MT.split(","))):
								sample_each_MT = samples_MT_ind[j]
								if samples_array == sample_each_MT:
									energ_MT.append([is_target.iloc[0,k],1])
								else:
									continue
							for l in range(len(samples_WT_srt.split(","))):
								sample_each_WT = samples_WT_ind[l]
								if samples_array == sample_each_WT:
									energ_WT.append([is_target.iloc[0,k],0])
								else:
									continue
						df_MT = map(list,unique(map(tuple, energ_MT)))
						df_WT = map(list,unique(map(tuple, energ_WT)))
						all_list = df_WT+df_MT
						res = np.array(all_list)
						group1 = res[:, 1] == 0
						group1 = res[group1][:, 0]
						group2 = res[:, 1] == 1
						group2 = res[group2][:, 0]
						print target, group1, group2
						fold_change = math.log(average(group2)/average(group1),2)
						u, p_value = mannwhitneyu(group1, group2)
						out_file.write(chromosome + "\t" +  str(start) + "\t" +  str(end) + "\t" +  target + "\t"+ str(fold_change) + "\t" + str(p_value) + "\t" +  samples_MT + "\n")
			else:
				continue
		else:
			continue
