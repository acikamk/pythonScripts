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
Translates pybios like ids file with 1 offset from Fabian
 to simtools optimized vector
 Inpuit:
 	- input_file: the file provided with S,P,K naminng with 1 offset
 	- id_file_or: identifiers.py file

 	optional: transl 
 		- parameter that tells whether the Translocation paramters from 
 		optimizaiton should be included in the resulting progress_trace file 

 NOTE: kCDX od kX can be used as weight parameters so this might needs
 to be addapted 
'''

def main(input_file, id_file_org):

	transl = True


	assert os.path.isfile(input_file), \
		 "Input file" + input_file + " doesn't exists!"
	assert os.path.isfile(id_file_org), \
		 "Input file" + id_file_org + " doesn't exists!"

	
	if transl:
		progress_file = input_file[:-4]+'_progress_trace_wTrl.txt'
	else: 
		progress_file = input_file[:-4]+'_progress_trace.txt'

	# read in the data from Fabian
	with open(input_file, 'r') as f:
		data=[l.strip().split('\t') for l in f]

	# read in the identifiers.py file
	with open(id_file_org, 'r') as f:
		org = Odict()
		execfile(id_file_org, Odict(), org)

	par_ids = org['par_ids']

	# handle the optimized parameters data 
	n = sum([1 for i in data[0] if 'kCD' in i ])
	k_ids = data[0][0:n]
	model_ids = data[0][n:]
	vectors = [[float(i) for i in j] for j in data[1:]] #!!!!!

	# check if kCDs where used
	if 'kCD0' in k_ids:
		kCDs = ['kCD'+str(i) for i in xrange(n)]
	else:
		kCDs = ['kCD'+str(i+1) for i in xrange(n)]

	pybios_ids = []
	for k in model_ids:
		pybios_ids.append(\
			par_ids.keys()[zip(*par_ids.values())[0].index(int(k[1:])-1)])

	# pdb.set_trace()
	index = kCDs + pybios_ids
	if not transl: 
		index = kCDs + [i for i in pybios_ids if "DrugTranslocation" not in i]

	# create empty dataframe with correct sorted parameter ids
	columns = [i for i in xrange(len(vectors))]
	matrix = np.empty((len(index),len(columns),))
	matrix.fill(np.nan)
	df_res = pd.DataFrame(matrix, index=index, columns=columns)

	
	# create dataframe with mapped optimized values to pybios_ids
	items = []
	keys = k_ids + pybios_ids
	assert len(keys) == len(vectors[0]), "Problems with param_ids vector!"
	for i, k in enumerate(keys):
		if not transl and "DrugTranslocation" in k:
				continue
		else:
			vals = [i for i in zip(*vectors)[i]]
			items.append((k, vals))

	df_in = pd.DataFrame.from_items(items, columns=columns, orient='index')
	
	# fill in the empty df
	df_res.update(df_in)
	# pdb.set_trace()
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
	try:
		main(sys.argv[1], sys.argv[2])
	except:
		print 'Translates pybios like ids file with 1 offset from Fabian\n\
		to simtools optimized vector \n\
		Inpuit: \n\
			\t- input_file: the file provided with S,P,K naminng with 1 offset\n\
			\t- id_file_or: identifiers.py file\n\
		\t 	optional: transl \n \
		\t\t- parameter that tells whether the Translocation paramters from \n\
		 		optimizaiton should be included in the resulting progress_trace file \n\n\
		NOTE: kCDX od kX can be used as weight parameters so this might needs to be addapted\n'