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

"""
Script that runs forward simulations for a given number of random vectors or 
for given optimized vector by increasing each parameter at a time for a given
percentage. 

Inputs:
model_path = path to the pybios model
simtools_dir = path to simtools
config_path = path to temporary config file (local verson)
tmp_progress_trace_dir = path to the directory where temporary 
						 progress_trace is stored  
numb_vect = 20
opti_flag = False
perc = 30

Output:
simlated_vector files stored in the res_dir given in the tmp.config or simtools

"""

model_path = "/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/OncoTrack_sim_23_12_16/"
simtools_dir = "/home/H0001-1/a.kovachev/simtools/"
config_path = "/home/H0001-1/a.kovachev/simtools/data/results/config_files/"
tmp_progress_trace_dir = "/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/"
numb_vect = 20
opti_flag = False
perc = 30

def  test_param_sensitivity(model_path, numb_vect, opti_flag, perc):

	model_input = model_path+"inputs/"
	model_output = model_path+"outputs/"
	sim_vect_file = model_output+"simulated_vectors_file.txt"
	sys.path.append(model_input)
	import modeltools 
	mp = modeltools.ModelProcessor(maxStepNum=100000)
	total_results = Odict()

	tmp1, tmp2, ids_vector, opti_vector = parse_simulation_file(sim_vect_file)

	if opti_flag: 
		i_size = len(opti_vector)
	else:
		i_size = numb_vect

	for i in range(i_size):
		if opti_flag:
			 tmp_opti_vect = opti_vector
			 tmp_opti_vect[i] = tmp_opti_vect[i] + tmp_opti_vect[i] * perc/100 #  change + to - for decrease
		else:
			tmp_opti_vect = np.random.uniform(-1,1,17875) 

		tmp_prog_file = tmp_progress_trace_file(tmp_progress_trace_dir, 
												ids_vector, 
												tmp_opti_vect)
		# set the random progres trace file in the config file
		config_file = modify_config_file(config_path, 
										"VectorsFile",
										 tmp_prog_file)

		simtools_exe = simtools_dir + "simcontrol.sh"
		subprocess.call([simtools_exe, 'startjob', config_file])

		print "Finished simulation for id: {} at {}".\
			format(i, datetime.datetime.strftime(datetime.datetime.now(), "%H:%M %d-%m-%Y"))
		res_dir, res_file = get_config_data(config_file)

		constr_dict, states_dict, tmp1, tmp2 = parse_simulation_file(res_file)
		total_results[i] = Odict()
		total_results[i]["constraints"] = constr_dict
		total_results[i]["states"] = states_dict
		tmp_sim_file = tmp_progress_trace_dir + "rand_vectors_file_4id_" + str(i) + ".txt"

		shutil.move(res_file, tmp_sim_file)
		
		shutil.rmtree(res_dir)
		# set the config file bacj to "VectorsFile"
		config_file = modify_config_file(config_path, 
										tmp_prog_file, 
										"VectorsFile")					
	return 

def parse_simulation_file(simulated_vectors_file):

	constr_dict = Odict()
	states_dict = Odict()
	ids_vector = []
	opti_vector = []

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
			elif "VECTOR:" in line:
				opti_vector = eval(line[line.index(":")+1:])		
			elif "STATES:" in line: 
				key = line[line.index(":")+1:line.rindex(":")]
				value = eval(line[line.rindex(":")+1:])
				states_dict[key] = value
	
	return constr_dict, states_dict, ids_vector, opti_vector




def tmp_progress_trace_file(tmp_progress_trace_dir, ids_vector, opti_vector):
	
	tmp_prog_file = tmp_progress_trace_dir + "tmp_progress_trace.txt"
	tmp_prog = open(tmp_prog_file, "w")

	tmp_prog.write("VARIABLE PAR IDS:\t" + str(ids_vector) + "\n")
	tmp_prog.write("VECTOR:\t" + str(list(opti_vector)) + "\n")

	tmp_prog.close()

	return tmp_prog.name

def modify_config_file(cfg_path, search, text):

	cfg_file= cfg_path + "tmp_cfg.py"   
	cfg = open(cfg_file, "r")
	cfg_txt = cfg.read()

	rep = {}
	rep[search] = text

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
		res_dir =  res_dir + "/" + model_id + "/" + job_id + "/"
		res_file = res_dir + "outputs/simulated_vectors_file.txt"
	else:
		res_dir = simtools_dir + "data/results/" + model_id + "/" + job_id + "/"
		res_file = res_dir + "outputs/simulated_vectors_file.txt"

	return res_dir, res_file

test_param_sensitivity(model_path, 
					   numb_vect, 
					   opti_flag, 
					   perc)