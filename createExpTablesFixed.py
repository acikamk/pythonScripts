import numpy as np
import pandas as pd
import math
from collections import OrderedDict as Odict
import random
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import itertools
import pdb

'''
Script takes input vector simulated to steadystate and ids vector,
and than creates X experiment tables where all n-log(i) parameters 
are set to fix while the rest log(i) are used for optimization
The fixed parameters are taken from the init_file where they are run to
steady state.
Cost function is fixed.
Exp tables are generated. 
number of parameters = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2028, 3000]

Inputs
	- init vectors file init_vector_Hasenauer2_cellDiv_fold.txt
	- default experiment table precompiled with simtools
	- model, identifiers.py?
Outputs
	- experiment tables?
	- config files
		to change exp_table info 
		to change duration_sec info
		to change report_each_sec
		to change job id info
'''

def main(input_vectors_file, exp_table_file):
	rand = True
	f = open(input_vectors_file)
	data = f.readlines()
	par_ids = eval([l for l in data if l.startswith('PAR_IDS')][0].split('=')[-1])
	var_ids = eval([l for l in data if l.startswith('VAR_PAR_IDS')][0].split("=")[-1])
	vector = eval([l for l in data if l.startswith('RND_VECTOR')][0].split('=')[-1])
	var_idx = eval([l for l in data if l.startswith('VAR_PAR_IDX')][0].split('=')[-1])
	rnd_vect = [j for i, j in enumerate(vector) if i in var_idx] 
	pdb.set_trace()	

	assert len(var_ids) == len(rnd_vect), "Something wrong with input data!"

	exp_table = pd.read_table(exp_table_file, index_col=0, low_memory=False)

	nparam = [1, 2]

	for i in list(reversed(nparam)):
		tmp_dict = Odict()
		print i
		if rand:
			var = random.sample(var_ids, i)
			fix = [k for k in var_ids if k not in var]
		else:
			var = var_ids[0:i]
			fix = var_ids[i:]
		
		for j in fix:
			if not j in tmp_dict:
				tmp_dict[j] = [rnd_vect[var_ids.index(j)] for m in xrange(0, len(exp_table.columns))]
		
		# pdb.set_trace()	
	 	tmp_df = pd.DataFrame.from_dict(tmp_dict, orient='index', dtype='float64')	
	 	tmp_df.columns = exp_table.columns
	 	res = pd.concat([exp_table, tmp_df])
	 	res.index.rename('ID', inplace=True) 
	 	print len(tmp_dict), exp_table.shape, tmp_df.shape, res.shape
	 	name = ('.').join(exp_table_file.split('.')[:-1])+"_5_"+str(i)+".csv"
	 	print name
		res.to_csv(name, sep = '\t')

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
