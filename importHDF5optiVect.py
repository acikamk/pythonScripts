import sys
import pandas as pd
import numpy as np
import pdb
import operator
import re
import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
import h5py
import warnings

'''
Script that maps hdf5 optimized vectors file
to pybios ids and progress_trace.txt file

Input: 
	- model identifiers.py file
	- model xx.ID_mapping_parameters.txt file that
		maps pybios to sbml parameters ids
	- hdf5 file with optimized vectors
	- hdf5 file with model data, i.e. parameters list
	- taketop10 - False to exaluate all, else True for the top 10
Output:
	- hdf5_progress_trace.txt file

****TODO: Handling of input parameters****
'''

mapping_ids_file = '/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/Speedy_v3_r403445.ID_mapping_parameters.txt'
identifiers_file = '/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/identifiers.py'

file_h5_opti = "/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/model-Speedy_v3_r403445_v2-srv23ib.712174.0-multistarts.h5"
file_h5_names = "/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/Speedy_v3_r403445_v1_2.h5"
taketop10 = False

progress_file = '/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/hdf5_progress_trace.txt'

# evaluate the indentifiers.py file
with open(identifiers_file, 'r') as f:
	org = Odict()
	execfile(identifiers_file, Odict(), org)

par_ids = org['par_ids']
par_ids = sorted(par_ids.items(), key=operator.itemgetter(1))
par_ids = Odict(par_ids)


# read data from the hdf5 files

hdf5_opti = h5py.File(file_h5_opti, 'r')

# get cost values
cost = np.array(hdf5_opti["/finalCost"][0])
cost[cost == 0.0] = np.nan

# get optmized vectors
opti_vectors = hdf5_opti["/finalParameters"]

# take the best 10 vectors only
if taketop10:
	indexes = np.argsort(cost)[0:10]
	opti_vectors = opti_vectors[:,np.sort(indexes)]
	cost = cost[np.sort(indexes)]

# read in the hdf5 names
hdf5_names = h5py.File(file_h5_names, 'r')
opti_names = np.array(hdf5_names["/parameters/parameterNames"])
 
# read in the mapping dict of ids
sbml2pybios_df = pd.read_table(mapping_ids_file, sep='\t', index_col=1)
sbml2pybios = sbml2pybios_df.to_dict()['ExperimentTable_id']

# map the hdf5 names to pybios
# check if needs to be adapted to some changes in sbml definitions
hdf52pybios = {i:sbml2pybios[j] for i in opti_names \
	for j in sbml2pybios.keys() if (i[:i.rindex("_")]==\
									j+"_reaction" or i==j)}
# generate hdf5 dict of names and opti values
dict_opti = dict(zip(opti_names, opti_vectors))	
#pdb.set_trace()

# generate dict of pybios names and opti values, and add cost value
pybios2opti = {v1:dict_opti[k1] for k1,v1 in hdf52pybios.iteritems()}
pybios2opti['fitness'] = cost

# pybios_ids = [i for i in par_ids.keys() if i in pybios2opti.keys()]
# add kcds values and set correct pybios_ids order
# however they are not properly handled/imported form the optimized vectors
n_kCDs = len([i for i in pybios2opti.keys() if i.startswith('k')])
kCDs = ['kCD'+str(i) for i in range(0, n_kCDs)]
pybios_order = ['fitness'] + kCDs + [i for i in par_ids.keys() \
			if i in pybios2opti.keys()] 

# convert to dataframe and trasfrom them to normal scale
df_res = pd.DataFrame.from_dict(pybios2opti, orient='index')
# pdb.set_trace()
df_res.dropna(axis=1, how='all', inplace=True)
df_res = df_res.reindex(pybios_order)
df_res.loc[kCDs] = np.power(10, df_res.loc[kCDs])

# save to to file
with open(progress_file, "w") as prog_file:
		prog_file.write("VARIABLE PAR IDS:\t{}\n\n".\
				format(df_res.index.tolist()[1:]))
		for column in df_res:
			prog_file.write("FITNESS:\t{}\n".\
					format(df_res.loc['fitness',column]))
			prog_file.write("VECTOR:\t{}\n\n".\
					format(df_res.loc[kCDs[0]:,column].tolist()))


