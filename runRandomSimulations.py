import numpy as np
import datetime
import re
import os
import sys
import math
import collections 
from collections import OrderedDict as Odict
import subprocess
import shutil
import pdb
import pandas as pd

"""
Script that runs forward simulations for a given number of random vectors or 
for given optimized vector by increasing each parameter at a time for a given
percentage. 

Inputs:
model_path = path to the pybios model
simtools_dir = path to simtools
config_path = path to temporary config file (local verson) in this folder is tmp_cfg.py
tmp_progress_trace_dir = path to the directory where temporary 
						 progress_trace is stored  
numb_vect = 20
opti_flag = False
perc = 30

Output:
simlated_vector files stored in the res_dir given in the tmp.cfg or simtools/data/results dir

*****TODO: handle the input parameters from console if necessary!******

"""

model_path = "/project/V0001-1/modcell/users/a.kovachev/Model_Optimization_Mai2017/Model_Mai2017_sim_thomas_270318_val_best3/"
simtools_dir = "/home/H0001-1/a.kovachev/simtools/"
config_path = "/home/H0001-1/a.kovachev/simtools/data/results/config_files/"
tmp_progress_trace_dir = "/project/V0001-1/modcell/users/a.kovachev/Model_Optimization_Mai2017/"
numb_vect = 10
opti_flag = True
perc = -90
# needed to generate the random vector

def test_param_sensitivity(model_path, numb_vect, opti_flag, perc):

	model_input = model_path+"inputs/"
	model_output = model_path+"outputs/"
	sim_vect_file = model_output+"simulated_vectors_file.txt"
	sys.path.append(model_input)

	import modeltools 
	mp = modeltools.ModelProcessor(maxStepNum=100000)
	total_results = Odict()
	
	_, _, ids_vector, opti_vector, _ = parse_simulation_file(sim_vect_file)
	n_param = len([i for i in ids_vector if "{" in i])
	n_kCD = len([i for i in ids_vector if i.startswith("k")])
	
	if opti_flag:
		n_repeats = 1 
		i_size = len(opti_vector[-1])
		isrand =  'opti_'
	else:
		# n_param = 4072
		# n_kCD = 16
		i_size = n_param + n_kCD
		n_repeats = numb_vect
		isrand = 'rand_'
	tmp_opti = {}
	# pdb.set_trace()	
	for r in xrange(0, n_repeats):
		# tmp_progress_trace_rand0.txt
		org_rand_sim_file = tmp_progress_trace_dir + "tmp_progress_trace_" + isrand + str(r) +".txt"
		tmp_opti[r] = []
		if opti_flag:
			for i in range(i_size):
			 	tmp_opti_vect = list(opti_vector[-1]) # clean way to copy a list
			 	tmp_opti_vect[i] = tmp_opti_vect[i] + tmp_opti_vect[i]*perc/100.0 
			 	# pdb.set_trace()
				tmp_opti[r].append(tmp_opti_vect)
		else:
			# use the next 2 lines to generate new random vector
			# tmp_opti_vect = np.concatenate((np.random.rand(n_kCD), np.random.uniform(-1,1,n_param))) 
			# tmp_opti[r].append(tmp_opti_vect)
			# use the next 2 lines to use existing random vector
			_, _, _, tmp_opti_vect, _  = parse_simulation_file(org_rand_sim_file)
			tmp_opti[r].append(tmp_opti_vect[0])
			# pdb.set_trace()
			for i in xrange(i_size):
				# pdb.set_trace()
				tmp_mod = list(tmp_opti[r][0]) 
				tmp_mod[i] = tmp_mod[i] + tmp_mod[i]*perc/100.0 
				tmp_opti[r].append(tmp_mod)
			# pdb.set_trace()	
	# pdb.set_trace()	
	for k, vect in tmp_opti.iteritems():
		sufx = '_rand_'+str(k) if not opti_flag else '_opti_' + str(k)
		if perc>0:
			direct = 'incr_'
		else:
			direct = 'decr_'
		tmp_opti_dir = tmp_progress_trace_dir + 'Model_Mai2017_sim_' + str(abs(perc)) + direct + isrand + str(k)+'_best3/'
		tmp_prog_file = tmp_progress_trace_file(tmp_progress_trace_dir, 
												ids_vector, 
												vect,
												sufx+direct)
		# set the random progres trace file in the config file
		config_file = modify_config_file(config_path, 
										["VectorsFile","JobId"],
										[tmp_prog_file, tmp_opti_dir])
		# pdb.set_trace()
		simtools_exe = simtools_dir + "simcontrol.sh"
		subprocess.call([simtools_exe, 'startjob', config_file])
		# set the config file back to "VectorsFile"

		# print "Finished simulation at {}".\
		# format(datetime.datetime.strftime(datetime.datetime.now(), "%H:%M %d-%m-%Y"))
		# pdb.set_trace()
		res_dir, res_file = get_config_data(config_file)


		_, _, ids_vector2, _, fitness2 = parse_simulation_file(res_file)

		if not opti_flag:
			ids_vector2 = ['orig_fitness'] + ids_vector2

		change = [el/fitness2[0] for el in fitness2]

		res = pd.DataFrame(data={'IDS': ids_vector2, 'FITNESS': fitness2, 'CHANGE': change})
		res.sort_values(by=['CHANGE'], inplace=True)
		res.to_csv(res_dir+'results_'+str(abs(perc)) + direct + isrand + str(k), sep='\t', index=False)

		config_file = modify_config_file(config_path, 
										[tmp_prog_file, tmp_opti_dir], 
										["VectorsFile","JobId"])	
		shutil.move(tmp_prog_file, res_dir)				
	return 

	
	# constr_dict, states_dict, tmp1, tmp2 = parse_simulation_file(res_file)
	# total_results[i] = Odict()
	# total_results[i]["constraints"] = constr_dict
	# total_results[i]["states"] = states_dict
	# tmp_sim_file = tmp_progress_trace_dir + "rand_vectors_file_" + ".txt"

	# shutil.move(res_file, tmp_sim_file)
		
	# shutil.rmtree(res_dir)
 	

def parse_simulation_file(simulated_vectors_file):

	constr_dict = Odict()
	states_dict = Odict()
	ids_vector = []
	opti_vectors = []
	fitness = []

	with open(simulated_vectors_file) as file:
		for line in file:
			if "CONSTRAINTS:" in line:
				constr_dict_raw = eval(line.split('CONSTRAINTS:')[-1])
				for k, dd in constr_dict_raw.items():   
					# check if the experiments where fold or abs
					if len(k) == 2: 
						constr_dict['/'.join(k)] = dd
					elif len(k) == 1:
						constr_dict[''.join(k)] = dd
			elif "VARIABLE PAR IDS:" in line:
				ids_vector = eval(line[line.index("["):])
			elif "FITNESS:" in line:
				fitness.append(eval(line[line.index(":")+1:]))
			elif "VECTOR:" in line:
				opti_vectors.append(eval(line[line.index(":")+1:]))	
			elif "STATES:" in line: 
				key = line[line.index(":")+1:line.rindex(":")]
				value = eval(line[line.rindex(":")+1:])
				states_dict[key] = value
	
	return constr_dict, states_dict, ids_vector, opti_vectors, fitness




def tmp_progress_trace_file(tmp_progress_trace_dir, ids_vector, opti_vector, sufx=""):
	
	tmp_prog_file = tmp_progress_trace_dir + "tmp_progress_trace_"+ sufx +".txt"
	tmp_prog = open(tmp_prog_file, "w")

	tmp_prog.write("VARIABLE PAR IDS:\t" + str(ids_vector) + "\n")
	for i, vector in enumerate(opti_vector):
		tmp_prog.write("VECTOR:\t" + str(list(vector)) + "\n")

	tmp_prog.close()

	return tmp_prog.name

def modify_config_file(cfg_path, search, text):

	cfg_file= cfg_path + "tmp_cfg.py"   
	cfg = open(cfg_file, "r")
	cfg_txt = cfg.read()

	rep = {}
	if len(search)!=len(text):
		print "Search and replace text not of equal length. Please check!"
		pdb.set_trace()
	for (k,v) in zip(search, text):
		rep[k] = v
	# pdb.set_trace()
	# Replace in text
	rep = dict((re.escape(k), v) for k, v in rep.iteritems())
	pattern = re.compile("|".join(rep.keys()))

	tmp_cfg_txt = pattern.sub(lambda m: rep[re.escape(m.group(0))], cfg_txt)
	cfg = open(cfg_file, "w")
	cfg.write(tmp_cfg_txt)

	return cfg.name 

def get_config_data(config_file):
	res_dir = ""
	model_id = ""
	job_id = ""
	with open(config_file) as file:
		for l in file:
			if l.startswith("model_id"):
				model_id =  l.split("=")[1].strip().replace('"',"")
			if l.startswith("job_id"):
				job_id = l.split("=")[1].strip().replace('"',"")
			if l.startswith("results_dir"):
				res_dir = l.split("=")[1].strip().replace('"',"")
	if res_dir:
		# res_dir =  res_dir + "/" + model_id + "/" + job_id + "/"
		res_dir = job_id
		res_file = res_dir + "outputs/simulated_vectors_file.txt"
	else:
		res_dir = simtools_dir + "data/results/" + model_id + "/" + job_id + "/"
		res_file = res_dir + "outputs/simulated_vectors_file.txt"

	return res_dir, res_file

test_param_sensitivity(model_path, 
					   numb_vect, 
					   opti_flag, 
					   perc)