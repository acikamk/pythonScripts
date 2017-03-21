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

def main(inputs_folder, sufx):
	try:
		# inputs_folder = "/home/H0001-1/a.kovachev/simtools/data/models/Notch_Signalling_r377624/"
		model_id = inputs_folder[inputs_folder[:-1].rindex("/")+1:-1]
		identifiers_file = inputs_folder + "identifiers.py"
		exp_table_file = inputs_folder + "exp_table_body.csv"

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
	# model stuff
	bounds = (-1,1)
	model = PyBiosModel(model_id, identifiers_file, exp_table_full, cost_func)
	varKParVect = numpy.power(10, bounds[0] \
	  + numpy.random.rand(len(model.variable_indexK_arr)) 
	  *(bounds[1] - bounds[0]) ) 
	fixed_Ki = numpy.setdiff1d(xrange(model.dimK), 
			model.variable_indexK_arr)
	K = numpy.zeros(model.dimK)
	if len(model.variable_indexK_arr) > 0:
		K[model.variable_indexK_arr] = varKParVect

	var_par_ids_list = model.get_var_par_ids_list(with_kcd=False)
	par_ids_list = [pId for pId, par in model.parameters.items() if par.arrayId == 'K']
	f = open("./data/init_vector_" + sufx +".txt", "w")
	f.write("PAR_IDS={}\n\n".format(par_ids_list))
	f.write("VAR_PAR_IDX={}\n\n".format(model.variable_indexK_arr.tolist()))
	f.write("VAR_PAR_IDS={}\n\n".format(var_par_ids_list))
	f.write("RND_VECTOR={}\n\n".format(K.tolist()))
	# exp_table_sample = exp_table_full[sampleId].to_frame()
	for sampleId in exp_table_full.columns:
		print bcolors.OKGREEN \
		 + "Running simulations for: " + sampleId\
		 + bcolors.ENDC

		
		sample = model.samples[sampleId]
		if len(fixed_Ki) > 0:
			K[fixed_Ki] = sample.K[fixed_Ki] 
			# print K
		conc[sampleId] = run_simulation(sample, sampleId, K, f)
	f.close()
	create_cost_file(pyid, conc, conc.keys(), sufx)	

def create_cost_file(pyid, conc_dict, samples, sufx):
	weight = 1
	f = open("./data/cost_func_"+sufx+".txt", "w")
	try: 
		print bcolors.OKGREEN \
	 	 + "Saving to file " + f.name\
	 	 + bcolors.ENDC
		for key, val in pyid.iteritems():
			f.write("ID:cost_fct_{}\tkCD{}*{}\n".format(key,key,val))
			for s in samples:
				if any(conc_dict[s]):
					f.write("{}\t{}\t{}\n".format(s,conc_dict[s][key],weight))
	except:
		print bcolors.ERROR \
	 	 + "Problem with saving to file " + f.name\
	 	 + bcolors.ENDC
	f.close()
	 	 
def run_simulation(sample, sampleId, K, f):

	import modeltools
	sim_time = 1e6	
	step = sim_time/100
	mp = modeltools.ModelProcessor()

	res = mp.simulate(sample.S0, sim_time, sample.F, K)

	# import pdb; pdb.set_trace()
	
	f.write("CONC_{}={}\n\n".format(sampleId, res.finalState.tolist()))
	if not res.success:
		print bcolors.FAIL\
		 + "Simulation failed with flag {} for sample {}"\
		 .format(res.exitcode, sampleId)\
		 + bcolors.ENDC
		return []
	else:
		return res.finalState 

	return

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
