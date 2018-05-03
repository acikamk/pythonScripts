import numpy as np
import datetime
import re
import os
import sys
import math
import matplotlib.pyplot as plt
import collections 
from collections import OrderedDict as Odict
import subprocess
import shutil
from os import listdir
from os.path import isfile, join

"""
Script that runs anaysis on the random forward simulations results of the 
runRandomSimulations.py sript.

Inputs:
expected_values: dict of the respective cost function and the exp. values
results_path = path to the random analysis results
opti_sim_file = path to one optimized vector simulation results
spec_id = id of the cost function species
plot = if True then the plotsimresults.py script will be used


optional:
# path to the plot script and config file
# ** test.cfg must be addapted 
plot_exe = "/home/H0001-1/a.kovachev/simtools/src/plotsimresults.py"
config_file = "/home/H0001-1/a.kovachev/simtools/data/results/config_files/test.cfg"

numb_vect = 20
opti_flag = False
perc = 30

Output:
simlated_vector files stored in the res_dir given in the tmp.config or simtools

"""

expected_values = {"TUMOR-327-MB-P-TF-01-03-PRI-01-cellline-01-01": 100.0,           
			"TUMOR-327-MB-P-TF-01-03-PRI-01-cellline-01-01_Afatinib_Conc3.048nM": 106.78529526, 
			"TUMOR-327-MB-P-TF-01-03-PRI-01-cellline-01-01_Afatinib_Conc60000.0nM": 0.1955416504, 
			"TUMOR-327-MB-P-TF-01-03-PRI-01-cellline-01-01_Afatinib_Conc6666.67nM": 2.1314039891, 
			"TUMOR-327-MB-P-TF-01-03-PRI-01-cellline-01-01_Afatinib_Conc740.741nM": 8.955807587, 
			"TUMOR-327-MB-P-TF-01-03-PRI-01-cellline-01-01_Afatinib_Conc82.305nM": 22.252639812, 
			"TUMOR-364-CB-M-MF-01-04-PRI-01-cellline-01-01": 100.0, 
			"TUMOR-364-CB-M-MF-01-04-PRI-01-cellline-01-01_Afatinib_Conc3.048nM": 97.662889518, 
			"TUMOR-364-CB-M-MF-01-04-PRI-01-cellline-01-01_Afatinib_Conc60000.0nM": 0.1922298665, 
			"TUMOR-364-CB-M-MF-01-04-PRI-01-cellline-01-01_Afatinib_Conc6666.67nM": 37.950222582, 
			"TUMOR-364-CB-M-MF-01-04-PRI-01-cellline-01-01_Afatinib_Conc740.741nM": 46.489275596, 
			"TUMOR-364-CB-M-MF-01-04-PRI-01-cellline-01-01_Afatinib_Conc82.305nM": 68.494536624,
			}

results_path = "/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/random_vectors_analysis/"
opti_sim_file = "/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/OncoTrack_sim_23_12_16/outputs/simulated_vectors_file.txt"

plot_exe = "/home/H0001-1/a.kovachev/simtools/src/plotsimresults.py"
config_file = "/home/H0001-1/a.kovachev/simtools/data/results/config_files/test.cfg"
spec_id = 3520
plot = False

# additional functions

# extract simulated file data
def parse_simulation_file(simulated_vectors_file):

	constr_dict = Odict()
	states_dict = Odict()
	ids_vector = []
	opti_vector = []
	fitness = ""

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
				fitness = eval(line[line.index(":")+1:])
			elif "VECTOR:" in line:
				opti_vector = eval(line[line.index(":")+1:])		
			elif "STATES:" in line: 
				key = line[line.index(":")+1:line.rindex(":")]
				value = eval(line[line.rindex(":")+1:])
				states_dict[key] = value
	
	return constr_dict, states_dict, fitness, ids_vector, opti_vector

# run analysis for the total random resutls and opti simulated file
"""
total_results - dict of the reuslts of the random sim
tumors_dict - results from the expected vals and optimized simulation file
fitness_list
"""
def run_statistical_analysis(total_results, tumors_dict, fitness_list):

	for tumor in tumors_dict.keys():

		print "Plotting results for sample: " + tumor
		fig, ax = plt.subplots()
		# the histogram of the data
		num_bins = range(len(tumors_dict[tumor]["states"]))
		bins = ax.bar(num_bins, tumors_dict[tumor]["states"])
		plt.setp(bins[0], 'facecolor', 'r')
		# ax.set_ylim(10e-6, 10)
		plt.yscale('log', nonposy='clip')
		fig.tight_layout()
		plt.draw()
		plt.savefig(os.path.join(output_path, tumor + "_states.png"), 
							dpi=200, bbox_inches='tight')
		plt.close(fig)
		fig, ax = plt.subplots()
		# the histogram of the data
		num_bins = range(len(tumors_dict[tumor]["constraints"]))
		bins = ax.bar(num_bins, tumors_dict[tumor]["constraints"])
		plt.setp(bins[0], 'facecolor', 'r')
		plt.setp(bins[1], 'facecolor', 'g')
		# ax.set_ylim(10e-6, 10)
		plt.yscale('log', nonposy='clip')
		fig.tight_layout()
		plt.draw()
		plt.savefig(os.path.join(output_path, tumor + "_constr.png"), 
							dpi=200, bbox_inches='tight')
		plt.close(fig)

	print "Plotting results for fitness"
	fig, ax = plt.subplots()
	# the histogram of the data
	num_bins = range(len(fitness_list))
	bins = ax.bar(num_bins, fitness_list)
	plt.setp(bins[0], 'facecolor', 'r')
	plt.yscale('log', nonposy='clip')
	# Tweak spacing to prevent clipping of ylabel
	fig.tight_layout()
	plt.draw()
	plt.savefig(os.path.join(output_path, "fitness_stats.png"), 
						dpi=200, bbox_inches='tight')
	plt.close(fig)
	return

# start analysis

output_path = results_path + "outputs/"
sim_file = output_path + "simulated_vectors_file.txt"

files = [f for f in listdir(output_path) if isfile(join(output_path, f))]
total_results = Odict()
tumors_dict = Odict()
fitness_list = []

# get all otpimized sim vector results and set them as first
constr_dict, states_dict, fitness, tmp1, opti_vector = parse_simulation_file(opti_sim_file)
for tumor in expected_values.keys():
	tumors_dict[tumor] = Odict()
	tumors_dict[tumor]["constraints"] = []
	tumors_dict[tumor]["states"] = []
	tumors_dict[tumor]["constraints"].append(expected_values[tumor])
	tumors_dict[tumor]["constraints"].append(constr_dict[tumor]["CellProliferation"])
	tumors_dict[tumor]["states"].append(states_dict[tumor][spec_id])
fitness_list.append(float(fitness))

print "Avg val of sim optimized vector: " + str(np.mean(np.array(opti_vector)))
print "Standard deviation of sim optimized vector: " + str(np.std(np.array(opti_vector)))

sim_opti_vector = opti_vector
sum_rand = np.zeros(len(opti_vector))

# get all random simulated data results
for f in files:
	tmp_file = output_path + f

	if plot:
		shutil.move(tmp_file, sim_file)
	constr_dict, states_dict, fitness, tmp1, opti_vector = parse_simulation_file(tmp_file)

	total_results[f] = Odict()
	total_results[f]["constraints"] = constr_dict
	total_results[f]["states"] = states_dict
	total_results[f]["fitness"] = float(fitness)
	total_results[f]["opti_vector"] = opti_vector
	sum_rand +=np.array(opti_vector)

	if plot:
		subprocess.call(["python", plot_exe, "-c" + config_file])
		shutil.move(sim_file, tmp_file)
		shutil.move(output_path + "analyzed_simulation_results.txt", output_path + "an_res_"+f[11:-4]+".txt")
		shutil.move(output_path + "cost_func_evaluation.txt", output_path + "cost_res_"+f[11:-4]+".txt")
		shutil.move(output_path + "Function_id_[3274].png",  output_path + "cost_plot_"+f[11:-4]+".png")
		shutil.move(output_path + "MAX-001:P[S62]-MYC-001 [comple_[complex in nucleus]_id_3274.png", output_path + "myc_plot_"+f[11:-4]+".png")

	fitness_list.append(float(fitness))

	for tumor in states_dict.keys():
		if not tumor in tumors_dict.keys():
			tumors_dict[tumor] = Odict()
			tumors_dict[tumor]["constraints"] = []
			tumors_dict[tumor]["states"] = []
		tumors_dict[tumor]["constraints"].append(constr_dict[tumor]["CellProliferation"])
		tumors_dict[tumor]["states"].append(states_dict[tumor][spec_id])

	print "Avg val of rand optimized vector: " + str(np.mean(np.array(opti_vector)))
	print "Standard deviation of sim optimized vector: " + str(np.std(np.array(opti_vector)))

avg_rand = sum_rand/len(files)
st_rand = np.std()


plt.vlines(range(50),[-1], avg_rand[:50])
plt.vlines(np.arange(50)+0.2, [-1], sim_opti_vector[:50], colors='r')
plt.draw()
plt.savefig(os.path.join(output_path, "rand_vals_stats.png"), 
						dpi=200, bbox_inches='tight')

run_statistical_analysis(total_results, constr_dict, states_dict, tumors_dict, fitness_list)
	
