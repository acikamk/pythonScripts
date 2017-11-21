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
import os
import sys
import pdb
import re
"""
Script that runs a model to steady state for each given samples
from experiment table independently and then creates 
cost function file cosisting of cost function for each species. 
Later on the created file sould be used for optimization to 
steady states.

NOTE: The experiment table should be output from simtools, i.e. 
all species/parameters should be present in order do be correctly 
assigned.
Also the program expects that the naming difference between
control/tumor and treatment is in the sign "_"
Ex:TUMOR-HCT20152-thyroid-NF-01-01
   TUMOR-HCT20152-thyroid-NF-01-01_20percKO_PSENEN

   CONTROL-152-CB-P-TF-01-04-PRI-01-cellline-01-01	
   TUMOR-152-CB-P-TF-01-04-PRI-01-cellline-01-01	
   TUMOR-152-CB-P-TF-01-04-PRI-01-cellline-01-01_Cetuximab_Conc0.003nM	

Input:
 - input folder - Folder where the identifiers.py and exp_table.csv are present
 - sufx - Suffix name to add to the output files
 - fold (optional) - whether the cost function to use fold changes 

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

def main():
	parser = argparse.ArgumentParser(description='Process input data.')
	parser.add_argument('inputs_folder', type=existingDir,
	                    help='Folder where the identifiers.py and exp_table.csv are present')
	parser.add_argument('sufix', type=str,
	                    help='Suffix name to add to the output files')
	parser.add_argument('--fold', type=bool, default=False,
	                    help='(optional) - whether the cost function to use fold changes')
	parser.add_argument('--cost_func', type=str, default=None,
	                    help='(optional) - cost function to be used: ex. "MAML1Complexes	kCD0*[48]+kCD1*[70]+kCD2*[72]+kCD3*[74]+0.00001"')
	parser.set_defaults(func=process_data)
	args = parser.parse_args()
	args.func(args)

	return

def process_data(args):
	# try:
		# MAML1Complexes	kCD0*[48]+kCD1*[70]+kCD2*[72]+kCD3*[74]+0.00001
		# CellDivisionFunction	(k0*[131]+k1*[1293]+k2*[130]+k3*[326]+k4*[902]+k5*[115]+k6*[325]+k7*[323]+k8*[829]+k9*[127]+k10*[919]+k11*[129])/(k12*[18]+k13*[19]+k14*[20]+k15*[21])
		# inputs_folder = "/home/H0001-1/a.kovachev/simtools/data/models/Notch_Signalling_r377624/"
		
	inputs_folder = str(args.inputs_folder)
	model_id = inputs_folder[inputs_folder[:-1].rindex("/")+1:-1]
	identifiers_file = existingFilePath(args.inputs_folder + "/identifiers.py")
	exp_table_file = existingFilePath(args.inputs_folder + "/exp_table_short.csv")
	cost_func_file =  existingFilePath(args.inputs_folder + "/cost_func_empty.txt")
	if args.cost_func:
		assert (len(str(args.cost_func).split())==2), "Incorect cost function format!"

	# except:	
	# 	print bcolors.WARNING \
	# 	 + "Problems with parsing files!"\
	# 	 + bcolors.ENDC
	# 	return
	sys.path.append(args.inputs_folder)
	# 
	try:
		f = open(identifiers_file)
		identifiers_code = compile(f.read(), identifiers_file, 'exec')
		f.close()
		eval(identifiers_code)
		diff_species = locals()['diff_var_ids']
		
		if not args.cost_func:
			pyid = {j[0]: k for k, j in diff_species.iteritems()}
		else:
			cost_function = args.cost_func.split()[1]
			if 'kCD' in cost_function:
				cost_function_noK = re.sub('kCD(\d+)\*', '', cost_function)	
			else:
				cost_function_noK = re.sub('k(\d+)\*', '', cost_function)			
			pybios_ids = re.findall(r"\[\d+\]", cost_function_noK)
			pyid =  {k: j[0] for k, j in diff_species.iteritems() if k in pybios_ids}
	except:	
		print bcolors.WARNING \
		 + "Problems with parsing identifiers.py file!"\
		 + bcolors.ENDC
		return

	exp_table_full = pandas.read_table(exp_table_file, sep='\t', index_col=0)

	conc = Odict()
	# model stuff
	bounds = (-1,1)
	model = PyBiosModel(model_id, identifiers_file, exp_table_full, cost_func_file)
	varKParVect = numpy.power(10, bounds[0] \
	  + numpy.random.rand(len(model.variable_indexK_arr)) 
	  *(bounds[1] - bounds[0]) ) 
	fixed_Ki = numpy.setdiff1d(xrange(model.dimK), 
			model.variable_indexK_arr)
	K = numpy.zeros(model.dimK)
	# import pdb; pdb.set_trace()
	if len(model.variable_indexK_arr) > 0:
		K[model.variable_indexK_arr] = varKParVect

	var_par_ids_list = model.get_var_par_ids_list(with_kcd=False)
	par_ids_list = [pId for pId, par in model.parameters.items() if par.arrayId == 'K']
	f = open("./data/init_vector_" + args.sufix +".txt", "w")
	f.write("PAR_IDS={}\n\n".format(par_ids_list))
	f.write("VAR_PAR_IDX={}\n\n".format(model.variable_indexK_arr.tolist()))
	f.write("VAR_PAR_IDS={}\n\n".format(var_par_ids_list))
	f.write("RND_VECTOR={}\n\n".format(K.tolist()))
	# exp_table_sample = exp_table_full[sampleId].to_frame()
	for sampleId in exp_table_full.columns:
		if "CONTROL" in sampleId:
			continue
		print bcolors.OKGREEN \
		 + "Running simulations for: " + sampleId\
		 + bcolors.ENDC

		
		sample = model.samples[sampleId]
		if len(fixed_Ki) > 0:
			K[fixed_Ki] = sample.K[fixed_Ki] 
		conc[sampleId] = run_simulation(sample, sampleId, K, f)
	f.close()
	weight = 1
	f = open("./data/cost_func_"+ args.sufix + ".txt", "w")
	# try: 
	print bcolors.OKGREEN \
 	 + "Saving to file " + f.name\
 	 + bcolors.ENDC
 	# import pdb; pdb.set_trace()
	if not args.cost_func:
		create_cost_file_species(pyid, conc, f, weight, args)	
	else:
		create_cost_file_formula(pyid, conc, f, weight, args, cost_function_noK)
		# import pdb; pdb.set_trace()
	# except:
	# 	print bcolors.FAIL \
	#  	 + "Problem with saving to file " + f.name\
	#  	 + bcolors.ENDC
	f.close()
	return

def create_cost_file_formula(pyid, conc_dict, f, weight, args, cost_function_noK):
	
	f.write("ID:{}\t{}\n".format(args.cost_func.split()[0], args.cost_func.split()[1]))
	for sample in conc_dict.keys():
		try:
			# import pdb; pdb.set_trace()
			if not args.fold:
				value = eval(reduce(lambda x, y: x.replace(y, str(conc_dict[sample][pyid[y]])), pyid.keys(), cost_function_noK))
				# print value
				# import pdb; pdb.set_trace()
				f.write("{}\t{}\t{}\n".format(sample,value,weight))
				# import pdb; pdb.set_trace()
			elif "_" in sample:
				control = sample.split("_")[0]
			
				value_t = eval(reduce(lambda x, y: x.replace(y, str(conc_dict[sample][pyid[y]])), pyid.keys(), cost_function_noK))
				value_c = eval(reduce(lambda x, y: x.replace(y, str(conc_dict[control][pyid[y]])), pyid.keys(), cost_function_noK))
				value =value_t/value_c
				f.write("{}/{}\t{}\t{}\n".format(sample,control,value,weight))
		except:
			print ("Problem with sample {}".format(sample))

	return

def create_cost_file_species(pyid, conc_dict, f, weight, args):

	for key, val in pyid.iteritems():
		f.write("ID:cost_fct_{}\tkCD{}*{}+0.00001\n".format(key,key,val))
		for sample in conc_dict.keys():
			try:
				# import pdb; pdb.set_trace()
				if not args.fold:
					f.write("{}\t{}\t{}\n".format(sample,conc_dict[sample][key],weight))
				elif "_" in sample:
					control = sample.split("_")[0]
					if conc_dict[control][key] == 0.0:
						conc_dict[control][key] = 1e9
					else:	
						value = conc_dict[sample][key]/conc_dict[control][key]
					f.write("{}/{}\t{}\t{}\n".format(sample,control,value,weight))
			except:
				print ("Problem with sample {}".format(sample))

	 	 
def run_simulation(sample, sampleId, K, f):

	import modeltools
	sim_time = 1e9	
	mp = modeltools.ModelProcessor()

	res = mp.simulate(sample.S0, sim_time, sample.F, K)

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

def existingDir(dir_path):
	dir_name = os.path.abspath(os.path.expanduser(dir_path))
	assert os.path.isdir(dir_name), "Cannot find dir: %s" % dir_path
	return dir_name

def existingFilePath(file_path):
	file_path = os.path.abspath(os.path.expanduser(file_path))
	assert os.path.isfile(file_path), "File %s was not found" % file_path
	return file_path

if __name__ == "__main__":
    main()
