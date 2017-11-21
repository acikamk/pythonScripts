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

'''
translates sbml file with ids vector and optimzied vectors to simtools optimized vector
'''

def main(input_file, id_file_org, mapping_file):

	assert os.path.isfile(input_file), \
		 "Input file" + input_file + " doesn't exists!"
	assert os.path.isfile(id_file_org), \
		 "Input file" + id_file_org + " doesn't exists!"
	assert os.path.isfile(mapping_file), \
		 "Input file" + mapping_file + " doesn't exists!"

	progress_file = input_file[:-4] + '_progress_trace.txt'
	transl = True

	with open(input_file, 'r') as f:
		vectors = []
		for line in f:
			if 'VECTOR' in line:
				sbml_ids=eval(line.split('=')[1])
			else:
				
				vectors.append(eval(line.split('=')[1]))
	
	with open(id_file_org, 'r') as f:
		org = Odict()
		execfile(id_file_org, Odict(), org)

	par_ids = org['par_ids']

	map_df = pd.read_csv(mapping_file, sep = '\t', index_col = 1)
	map_dict = x=map_df['ExperimentTable_id'].to_dict()

	# optimized parameters data 
	k_ids = [k for k in sbml_ids if k.startswith('k')]
	pdb.set_trace()
	pybios_ids = [map_dict[i.split('_reaction')[0]] for i in sbml_ids if not i.startswith('k') ]

	index = pybios_ids + k_ids

	columns = [i for i in xrange(len(vectors))]
	matrix = np.empty((len(index),len(columns)))
	matrix.fill(np.nan)
	df_res = pd.DataFrame(matrix, index=index, columns=columns)

	pdb.set_trace()

	items = []
	
	assert len(index) == len(vectors[0]), "Problems with param_ids vector!"
	for i, k in enumerate(index):
		vals = [i for i in zip(*vectors)[i]]
		items.append((k, vals))

	df_in = pd.DataFrame.from_items(items, columns=columns, orient='index')
	

	df_res.update(df_in)
	pdb.set_trace()
	if df_res.isnull().values.any():
		print "NaN values exists!"
		pdb.set_trace()

	with open(progress_file, "w") as prog_file:
		prog_file.write("VARIABLE PAR IDS:\t{}\n\n".format(df_res.index.tolist()))
		for column in df_res:
			prog_file.write("VECTOR:\t{}\n\n".format(df_res[column].tolist()))

	# with open(progress_file, "w") as prog_file:
	# prog_file.write("VARIABLE PAR IDS:\t{}\n\n".format(index))
	# for vect in vectors:
	# 	prog_file.write("VECTOR:\t{}\n\n".format(vect))

	return


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
