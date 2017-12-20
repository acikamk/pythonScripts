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

def getClustersFromCorr(candids_corr, corr_thresh):
# from sklearn.cluster import KMeans
# kmZ=KMeans(n_clusters=20).fit(z)
# kmZ.cluster_centers_
	ini_clusters = {}
	final_clusters = {}
	pdb.set_trace()
	n_samps=candids_corr.shape[0]
	for i in range(1,n_samps):
		col=candids_corr[:,i]
		#col(i)=0;
		ind=np.where(col>corr_thresh)
		ini_clusters[i]=ind

	pdb.set_trace()
	np_nodes=range(1,n_samps)
	# nodes not yet parsed
	cls_num=0
	pdb.set_trace()
	while np_nodes:
		crt_node=np_nodes[0]
		# put in a queue the current node's neighbors
		queue=ini_clusters[crt_node][ini_clusters[crt_node] != crt_node]
		final_clusters[cls_num] = [crt_node, queue];
		pdb.set_trace()
		while queue:
			node_to_parse=queue[-1]
			del queue[-1]
			new_nb=set(ini_clusters[node_to_parse])-set(final_clusters[cls_num])
			queue=[new_nb, queue]
			final_clusters[cls_num]=[final_clusters[cls_num], new_nb]
		np_nodes=set(np_nodes) - set(final_clusters[cls_num])
		cls_num=cls_num + 1
	pdb.set_trace()
	count=1
	for i in range(1,len(final_clusters)):
		if len(final_clusters[i]) > 1:
			out_clusters[count]=final_clusters[i]
			count=count+1
	return out_clusters


mapping_ids_file = '/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/Speedy_v3_r403445.ID_mapping_parameters.txt'
identifiers_file = '/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/identifiers.py'


with open(identifiers_file, 'r') as f:
	org = Odict()
	execfile(identifiers_file, Odict(), org)

par_ids = org['par_ids']
par_ids = sorted(par_ids.items(), key=operator.itemgetter(1))
par_ids = Odict(par_ids)

file_h5_opti = "/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/model-Speedy_v3_r403445_v2-srv23ib.712174.0-multistarts.h5"
file_h5_names = "/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/Speedy_v3_r403445_v1_2.h5"

hdf5_opti = h5py.File(file_h5_opti, 'r')

cost = np.array(hdf5_opti["/finalCost"][0])
cost[cost == 0.0] = np.nan

opti_vectors = hdf5_opti["/finalParameters"]

progress_file = '/home/H0001-1/a.kovachev/simtools/data/models/Speedy_v3_r403445/hdf5_progress_trace.txt'


hdf5_names = h5py.File(file_h5_names, 'r')
opti_names = np.array(hdf5_names["/parameters/parameterNames"])
 
sbml2pybios_df = pd.read_table(mapping_ids_file, sep='\t', index_col=1)
sbml2pybios = sbml2pybios_df.to_dict()['ExperimentTable_id']

hdf52pybios = {i:sbml2pybios[j] for i in opti_names for j in sbml2pybios.keys() if (i[:i.rindex("_")]==j+"_reaction" or i==j)}

dict_opti = dict(zip(opti_names, opti_vectors))	
# pdb.set_trace()

pybios2opti = {v1:dict_opti[k1] for k1,v1 in hdf52pybios.iteritems()}
pybios2opti['fitness'] = cost

# pybios_ids = [i for i in par_ids.keys() if i in pybios2opti.keys()]
n_kCDs = len([i for i in pybios2opti.keys() if i.startswith('k')])
kCDs = ['kCD'+str(i) for i in range(0, n_kCDs)]
pybios_order = ['fitness'] + kCDs + [i for i in par_ids.keys() if i in pybios2opti.keys()] 

df_res = pd.DataFrame.from_dict(pybios2opti, orient='index')
df_res.dropna(axis=1, inplace=True)
df_res = df_res.reindex(pybios_order)
df_res.loc[kCDs] = np.power(10, df_res.loc[kCDs])
# pdb.set_trace()
z=opti_vectors[:,~np.all(np.isnan(opti_vectors), axis=0)]
x=np.corrcoef(z, rowvar=False)
thr = 0.9
getClustersFromCorr(x,thr)

from sklearn.cluster import MeanShift
mpp = MeanShift().fit(np.transpose(x))
mpp.labels_
len(mpp.cluster_centers_)

# print to file
# with open(progress_file, "w") as prog_file:
# 		prog_file.write("VARIABLE PAR IDS:\t{}\n\n".format(df_res.index.tolist()))
# 		for column in df_res:
# 			prog_file.write("FITNESS:\t{}\n".format(df_res.loc['fitness',column]))
# 			prog_file.write("VECTOR:\t{}\n\n".format(df_res.loc[kCDs[0]:,column].tolist()))


