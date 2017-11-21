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


def main(input_vectors_file, progress_trace, last):
	last = int(last)
	if not (last == -1 or last == 0 or last == -2):
		last = -1
	f = open(input_vectors_file)
	data = f.readlines()
	par_ids = eval([l for l in data if l.startswith('PAR_IDS')][0].split('=')[-1])
	var_ids = eval([l for l in data if l.startswith('VAR_PAR_IDS')][0].split("=")[-1])
	vector = eval([l for l in data if l.startswith('RND_VECTOR')][0].split('=')[-1])
	var_idx = eval([l for l in data if l.startswith('VAR_PAR_IDX')][0].split('=')[-1])
	rnd_vect = [j for i, j in enumerate(vector) if i in var_idx] 
	assert len(var_ids) == len(rnd_vect), "Something wrong with input data!"
	f.close()
	# 	
	f = open(progress_trace)
	data = f.readlines()
	var_par_ids = eval([l for l in data if l.startswith('VARIABLE PAR IDS:')][0].split('\t')[-1])
	if last !=-2:
		all_vect = [eval([l for l in data if l.startswith('VECTOR:')][last].split('\t')[-1])]
	else:
		all_vect = [eval(i.split('\t')[-1]) for i in [l for l in data if l.startswith('VECTOR:')]]
	f.close()
	# pdb.set_trace()	
	# corr_matrix = []
	
	
	for v in all_vect:
		c=0
		a = []
		b = []
		for var in var_par_ids:
			if var in var_ids:
				# print rnd_vect[var_ids.index(var)], math.pow(10, last_vector[var_par_ids.index(var)])
				rel_err = abs(1-math.pow(10, v[var_par_ids.index(var)])/rnd_vect[var_ids.index(var)])
				# print abs(1-math.pow(10, last_vector[var_par_ids.index(var)])/rnd_vect[var_ids.index(var)])
				if rel_err < 0.03:
					print var
					# print abs(1-math.pow(10, last_vector[var_par_ids.index(var)])/rnd_vect[var_ids.index(var)])
					# print rnd_vect[var_ids.index(var)], math.pow(10, last_vector[var_par_ids.index(var)])
					c+=1
				a.append(rnd_vect[var_ids.index(var)])
				b.append(math.pow(10, v[var_par_ids.index(var)]))
	# pdb.set_trace()
		if not len(a)<4:
			print c
			print np.corrcoef(a, b)
		else:
			print [(i, j, abs(1-i/j)) for (i,j) in zip(a, b)]  


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
