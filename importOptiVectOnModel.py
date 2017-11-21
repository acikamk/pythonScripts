import sys
import pandas as pd
import numpy as np
import pdb
import pickle
import re
import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join


# def sort_input(opti_params, par_ids_mod, vectors, n):
# 	output = sorted(opti_params, key=lambda x: par_ids_mod.index(x))

# 	for vect in vectors:
# 		vector_modi = [None]*len(opti_params)
# 		for i in xrange(len(vector[16:])):
# 			vector_modi[i] = vector[opti_params.index(output[i])]
# 		vector[16:]=vector_modi

# 	return output, vectors

# def fill_input(opti_params, par_ids_mod, vectors):
# 	missing_ids = [par_ids_mod.index(i) for i in par_ids_mod \
# 		if i not in opti_params]

	# if len(opti_params) == len(par_ids_mod): 
	# 	if opti_params != par_ids_mod:
	# 		print "Input vector is not correctly sorted..."
	# 		[opti_params, vectors] = sort_input(opti_params, par_ids_mod, vectors)
	# else:		
	# 	print "Input vector is not correct size..."
	# 	[opti_params, vectors] = fill_input(opti_params, par_ids_mod, vectors)


def main(input_file, id_file_org, id_file_mod):

	assert os.path.isfile(input_file), \
		 "Input file" + input_file + " doesn't exists!"
	assert os.path.isfile(id_file_org), \
		 "Input file" + id_file_org + " doesn't exists!"
	assert os.path.isfile(id_file_mod), \
		 "Input file" + id_file_mod + " doesn't exists!"

	progress_file = input_file[:-4]+'_progress_trace_wTrl.txt'

	with open(input_file, 'r') as f:
		data=[l.strip().split('\t') for l in f]

	with open(id_file_org, 'r') as f:
		org = Odict()
		execfile(id_file_org, Odict(), org)

	with open(id_file_mod, 'r') as f:
		mod = Odict()
		execfile(id_file_mod, Odict(), mod)

	par_ids = org['par_ids']
	par_ids_mod = {k:v for k,v in mod['par_ids'].items() if \
		not 'RPKM2protein' in k \
		and not 'GeneSpecificScaling' in k \
		and not 'DrugTargetInteraction' in k}

	pdb.set_trace()
	# optimized parameters data 
	n = sum([1 for i in data[0] if 'kCD' in i ])
	k_ids = data[0][0:n]
	model_ids = data[0][n:]
	vectors = [[float(i) for i in j] for j in data[1:]]

	if 'kCD0' in k_ids:
		kCDs = ['kCD'+str(i) for i in xrange(n)]
	else:
		kCDs = ['kCD'+str(i+1) for i in xrange(n)]

	index = kCDs + par_ids_mod.keys()
	columns = [i for i in xrange(len(vectors))]
	matrix = np.empty((len(index),len(columns),))
	matrix.fill(np.nan)
	df_res = pd.DataFrame(matrix, index=index, columns=columns)

	pdb.set_trace()


	# list of pybios_ids ordered as in the
	# input optimized vector
	pybios_ids = []
	for k in model_ids:
		pybios_ids.append(\
			par_ids.keys()[zip(*par_ids.values())[0].index(int(k[1:])-1)])

	missing_ids = {k: float(val[1]) for k, val in par_ids_mod.items() \
		if k not in par_ids.keys()}
	s_missing = pd.Series(missing_ids)	
	df_missing=pd.concat([s_missing] * len(columns), axis=1)
	pdb.set_trace()

	items = []
	keys = k_ids + pybios_ids
	for i, k in enumerate(keys):
		vals = [i for i in zip(*vectors)[i]]
		items.append((k,vals))

	df_in = pd.DataFrame.from_items(items, columns=columns, orient='index')
	# pdb.set_trace()
	df_res.update(df_in)
	df_res.fillna(df_missing, inplace=True)
	# df_res.dropna(inplace=True)
	pdb.set_trace()

	with open(progress_file, "w") as prog_file:
		prog_file.write("VARIABLE PAR IDS:\t{}\n\n".format(df_res.index.tolist()))
		for column in df_res:
			prog_file.write("VECTOR:\t{}\n\n".format(df_res[column].tolist()))


	return


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],sys.argv[3])
