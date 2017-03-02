import sys
import pandas
import numpy
import argparse
from argparse import RawTextHelpFormatter
import collections 
from collections import OrderedDict as Odict
sys.path.append("/home/H0001-1/a.kovachev/simtools/src")
import runopt
from modelstuff import PyBiosModel
import time
import matplotlib.pyplot as plt 
"""
Script that runs a model to steady state for each given samples
from experiment table independently and then creates 
cost function file cosisting of cost function for each species. 
Later on the created file sould be used for optimization to 
steady states.

NOTE: The experiment table should be output from simtools, i.e. 
all species/parameters should be present in order do be correctly 
assigned.

Input:
 - input folder - Folder where the identifiers.py and exp_table.csv are present

Output:
 - cost_func_conc.txt - Cost function file 

 ### TO DO: Decide how to handle in/out params
"""

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def main(inputs_folder):
	try:
		inputs_folder = "/home/H0001-1/a.kovachev/simtools/data/models/Notch_Signalling_r377624/"
		model_id = inputs_folder[inputs_folder[:-1].rindex("/")+1:-1]
		identifiers_file = inputs_folder + "identifiers.py"
		exp_table_file = inputs_folder + "exp_table.csv"

		cost_func =  inputs_folder + "cost_func_empty.txt"
		sys.path.append(inputs_folder)

		f = open(identifiers_file)
		identifiers_code = compile(f.read(), identifiers_file, 'exec')
		f.close()
		eval(identifiers_code)
		diff_species = locals()['diff_var_ids']
		pyid = {j[0]: k for k, j in diff_species.iteritems()}
	except:	
		print bcolors.WARNING \
		 + "Problems with analysis of the identifiers.py file!"\
		 + bcolors.ENDC
		return

	exp_table_full = pandas.read_table(exp_table_file, sep='\t', index_col=0)

	conc = Odict()
	for sampleId in exp_table_full.columns:
		print bcolors.OKGREEN \
		 + "Running simulations for: " + sampleId
		 + bcolors.ENDC

		exp_table_sample = exp_table_full[sampleId].to_frame()
		model = PyBiosModel(model_id, identifiers_file, exp_table_sample, cost_func)
		conc[sampleId] = run_simulation( model, sampleId)

	create_cost_file(pyid, conc, conc.keys())	

def create_cost_file(pyid, conc_dict, samples):
	weight = 1
	f  = open("cost_func_conc.txt", "w")
	for key, val in pyid.iteritems():
		f.write("ID:cost_fct_{}\tkCD{}*{}\n".format(key,key,val))
		for s in samples:
			if any(conc_dict[s]):
				f.write("{}\t{}\t{}\n".format(s,conc_dict[s][key],weight))
	f.close()

def run_simulation(model, sampleId):

	import modeltools
	diff_quot = 0
	sim_time = 1e4	
	step = sim_time/100
	time_vect = numpy.arange(0, float(sim_time)+1 , step, dtype=numpy.float) 
	bounds = (-1,1)
	mp = modeltools.ModelProcessor()

	varKParVect = numpy.power(10, bounds[0] \
			+ numpy.random.rand(len(model.variable_indexK_arr)) 
			*(bounds[1] - bounds[0]) ) 
	fixed_Ki = numpy.setdiff1d(xrange(model.dimK), 
			model.variable_indexK_arr)
	K = numpy.zeros(model.dimK)
	if len(model.variable_indexK_arr) > 0:
		K[model.variable_indexK_arr] = varKParVect

	sample = model.samples[sampleId]
	if len(fixed_Ki) > 0:
		K[fixed_Ki] = sample.K[fixed_Ki] 
	res1 = mp.simulate(sample.S0, sim_time, sample.F, K)
	if not res1.success:
		print bcolors.FAIL\
		 + "Simulation failed with flag {} for sample {}"\
		 .format(res1.exitcode, sampleId)
		 + bcolors.ENDC
		return []
	else:
		return res1.finalState 

	return

if __name__ == "__main__":
    main(sys.argv[1])
