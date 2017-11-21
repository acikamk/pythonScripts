import collections 
from collections import OrderedDict as Odict
from collections import deque 
import os
import sys
import argparse, ast
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
from libsbml import *
import re
import pandas as pd
import numpy as np
import pdb

prot2trans = {'L-04':'L-03',
			'L-06':'L-05',
			'L-02':'L-01',
			'L-08':'L-07',
			'L-10':'L-09',
			'L-12':'L-11',
			'L-14':'L-13',
			'L-16':'L-15',
			'L-18':'L-17'}


def updateCostFunc(cost_file, data, sample_dict, tc_pairs):

	f_out = open(cost_file[:cost_file.rindex('.')] +'_proteomix_v3.txt', 'w')
	
	with open(cost_file) as f:
		for line in f:
			if line.startswith('ID:') and 'ENSMUSG' in line:
				# pdb.set_trace()
				gene_id = line.split()[0].split(':')[1].split('_')[0].strip()
				spec = line.split()[0].split(':')[1].split('_')[1].strip()
				data_t = [i for i in data.keys() if gene_id in i and 'Conc' in i]
				data_c = [i for i in data.keys() if gene_id in i and 'Conc' not in i]
				# print bio_samp_mean_key
				if not data_t:
					print "ID not found for Gene: " + gene_id
					continue
				sample_dict_r = {v:k for k,v in sample_dict.iteritems()}
				f_out.write(line)
				# pdb.set_trace()
				for tumor in data_t:
					cs = sample_dict[pairs[sample_dict_r['-'.join(tumor.split('-')[1:9])]]]
					# pdb.set_trace()
					control = [c for c in data_c if cs in c][0]
					fold_change = data[tumor]/data[control]
					# print fold_change	
					prot_c = control[control.index('L'):control.index('L')+4]
					prot_t = tumor[tumor.index('L'):tumor.index('L')+4]
					# pdb.set_trace()
					t = tumor.replace(gene_id,'').replace(prot_t, prot2trans[prot_c]).replace('__','_')
					c = control.replace(gene_id,'').replace(prot_c, prot2trans[prot_c]).replace('_','')
					f_out.write("{}/{}\t{}\n".format(t, c, fold_change))

			# pdb.set_trace()

def parseTCFile(tc_file):

	tc_df = pd.read_table(tc_file, sep = '\t')
	tc_df.dropna(how='all', inplace=True) 
	controls = set(tc_df['Control'].tolist())
	treatments = set(tc_df['Treatment'].tolist())
	drugs = set(tc_df['Drug'].tolist())
	gf = set(tc_df['Growthfactor'].tolist())
	pairs_df =  tc_df[['Treatment', 'Control']] 
	# pdb.set_trace()
	pairs = {x[0]:x[1] for x in pairs_df.values}
	# 
	return tc_df, controls, treatments, pairs

def parseProteomixFile(p_files, controls, treatments, tc_df):

	Conc = {'MK2206':'1000.0nM', 'Wortmannin':'5000.0nM'}
	bio_samp = Odict()
	sample_dict = Odict()
	for file in p_files:
		with open(file) as f:
			df = pd.read_table(f, sep = '\t', low_memory = False)
			genes =  set(df['Genes'].tolist())
		for i in df.index:
			#C-13-29-N-01-L-08
			sample_id = df.iloc[i]['Sample ID']
			biosource = '-'.join(sample_id.split('-')[1:4])
			replicate = '-'.join(sample_id.split('-')[4])
			sample = '-'.join(sample_id.split('-')[5:])
			gene_id = df.iloc[i]['Genes'].split(';')[0]
			# print gene_id
			if '13' in biosource:
				replicates = '010203'
			if '11' in biosource:
				replicates = '040506'
			sample_name = '_'.join(['-'.join(['TUMOR', biosource, 'GEM-00-C', sample, 'Mean', replicates ]), gene_id])
			#TUMOR-11-28-N-GEM-00-C-L-09-Mean-010203	
			#TUMOR-11-28-N-GEM-00-C-L-09-Mean-010203_MK-2206_Conc1000.0nM	
			#TUMOR-11-28-N-GEM-00-C-L-09-Mean-010203_Wortmannin_Conc5000.0nM
			# pdb.set_trace()
			if sample_id in controls:
				sample_dict[sample_id] = '-'.join([biosource, 'GEM-00-C', sample])
			elif sample_id in treatments:
				# pdb.set_trace()
				drug = tc_df.loc[tc_df['Treatment'] == sample_id]['Drug'].values[0]
				# gf = tc_df.loc[tc_df['Treatment'] == sample_id]['Growthfactor'].values[0]
				sample_name = sample_name + '_' + drug + '_Conc' + Conc[drug]
				sample_name = sample_name.replace('MK2206', 'MK-2206')
				sample_dict[sample_id] = '-'.join([biosource, 'GEM-00-C', sample])
			else:
				print 'Sample name not present in tc_table!'
				pdb.set_trace()
			# pdb.set_trace()
			if not sample_name in bio_samp:
				bio_samp[sample_name] = []
			bio_samp[sample_name].append(df.iloc[i]['Absquant.fmol'])
	
	
	bio_samp_mean = {k: np.mean(v) for k, v in bio_samp.iteritems()}
	bio_samp_median = {k: np.median(v) for k, v in bio_samp.iteritems()}
	return bio_samp_mean, bio_samp_median, sample_dict 




tc_df, controls, treatments, pairs = parseTCFile(sys.argv[1])
# pdb.set_trace()
bio_samp_mean, bio_samp_median, sample_dict = parseProteomixFile(sys.argv[3:], controls, treatments, tc_df)
updateCostFunc(sys.argv[2], bio_samp_mean, sample_dict, pairs)