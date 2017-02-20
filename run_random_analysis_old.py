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


TABLEAU64 = [(0,0,0),(158,0,142),(1,0,103),(213,255,0),(255,0,86),(14,76,161),(255,229,2),
	(0,95,57),(0,255,0),(149,0,58),(255,147,126),(164,36,0),(0,21,68),(145,208,203),(98,14,0),
	(107,104,130),(0,0,255),(0,125,181),(106,130,108),(0,174,126),(194,140,159),(190,153,112),
	(0,143,156),(95,173,78),(255,0,0),(255,0,246),(255,2,157),(104,61,59),(255,116,163),
	(150,138,232),(152,255,82),(167,87,64),(1,255,254),(255,238,232),(254,137,0),(189,198,255)
	,(1,208,255),(187,136,0),(117,68,177),(165,255,210),(255,166,254),(119,77,0),(122,71,130)
	,(38,52,0),(0,71,84),(67,0,44),(181,0,255),(255,177,103),(255,219,102),(144,251,146),
	(126,45,210),(189,211,147),(229,111,254),(222,255,116),(0,255,120),(0,155,255),(0,100,1),
	(0,118,255),(133,169,0),(0,185,23),(120,130,49),(0,255,198),(255,110,65),(232,94,190)]

for i in range(len(TABLEAU64)):
    r, g, b = TABLEAU64[i]
    TABLEAU64[i] = (r/255., g/255., b/255.)

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
output_path = results_path + "outputs/"
sim_file = output_path + "simulated_vectors_file.txt"
opti_sim_file = "/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/OncoTrack_sim_23_12_16/outputs/simulated_vectors_file.txt"
plot_exe = "/home/H0001-1/a.kovachev/simtools/src/plotsimresults.py"
config_file = "/home/H0001-1/a.kovachev/simtools/data/results/config_files/test.cfg"
myc_id = 3520

files = [f for f in listdir(output_path) if isfile(join(output_path, f))]
total_results = Odict()
tumors_dict = Odict()
fitness_list = []

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

def run_statistical_analysis(total_results, constr_dict, states_dict, tumors_dict, fitness_list):

	# for tumor in tumors_dict.keys():

	# 	print "Plotting results for sample: " + tumor
	# 	fig, ax = plt.subplots()
	# 	# the histogram of the data
	# 	num_bins = range(len(tumors_dict[tumor]["states"]))
	# 	# import pdb; pdb.set_trace()
	# 	bins = ax.bar(num_bins, tumors_dict[tumor]["states"])
	# 	plt.setp(bins[0], 'facecolor', 'r')
	# 	# ax.set_ylim(10e-6, 10)
	# 	plt.yscale('log', nonposy='clip')
	# 	# Tweak spacing to prevent clipping of ylabel
	# 	fig.tight_layout()
	# 	plt.draw()
	# 	plt.savefig(os.path.join(output_path, tumor + "_states.png"), 
	# 						dpi=200, bbox_inches='tight')
	# 	plt.close(fig)
	# 	fig, ax = plt.subplots()
	# 	# the histogram of the data
	# 	num_bins = range(len(tumors_dict[tumor]["constraints"]))
	# 	# import pdb; pdb.set_trace()
	# 	bins = ax.bar(num_bins, tumors_dict[tumor]["constraints"])
	# 	plt.setp(bins[0], 'facecolor', 'r')
	# 	plt.setp(bins[1], 'facecolor', 'g')
	# 	# ax.set_ylim(10e-6, 10)
	# 	plt.yscale('log', nonposy='clip')
	# 	# Tweak spacing to prevent clipping of ylabel
	# 	fig.tight_layout()
	# 	plt.draw()
	# 	plt.savefig(os.path.join(output_path, tumor + "_constr.png"), 
	# 						dpi=200, bbox_inches='tight')
	# 	plt.close(fig)

	# print "Plotting results for fitness"
	# fig, ax = plt.subplots()
	# # the histogram of the data
	# num_bins = range(len(fitness_list))
	# bins = ax.bar(num_bins, fitness_list)
	# plt.setp(bins[0], 'facecolor', 'r')
	# plt.yscale('log', nonposy='clip')
	# # Tweak spacing to prevent clipping of ylabel
	# fig.tight_layout()
	# plt.draw()
	# plt.savefig(os.path.join(output_path, "fitness_stats.png"), 
	# 					dpi=200, bbox_inches='tight')
	# plt.close(fig)
	# return

constr_dict, states_dict, fitness, tmp1, opti_vector = parse_simulation_file(opti_sim_file)
for tumor in expected_values.keys():
	tumors_dict[tumor] = Odict()
	tumors_dict[tumor]["constraints"] = []
	tumors_dict[tumor]["states"] = []
	tumors_dict[tumor]["constraints"].append(expected_values[tumor])
	tumors_dict[tumor]["constraints"].append(constr_dict[tumor]["CellProliferation"])
	tumors_dict[tumor]["states"].append(states_dict[tumor][myc_id])
fitness_list.append(float(fitness))

print "Avg val of sim optimized vector: " + str(np.mean(np.array(opti_vector)))
print "Standard deviation of sim optimized vector: " + str(np.std(np.array(opti_vector)))
sim_opti_vector = opti_vector

sum_rand = np.zeros(len(opti_vector))
for f in files:
	tmp_file = output_path + f
	# shutil.move(tmp_file, sim_file)
	constr_dict, states_dict, fitness, tmp1, opti_vector = parse_simulation_file(tmp_file)

	total_results[f] = Odict()
	total_results[f]["constraints"] = constr_dict
	total_results[f]["states"] = states_dict
	total_results[f]["fitness"] = float(fitness)
	total_results[f]["opti_vector"] = opti_vector
	sum_rand +=np.array(opti_vector)

	fitness_list.append(float(fitness))

	for tumor in states_dict.keys():
		if not tumor in tumors_dict.keys():
			tumors_dict[tumor] = Odict()
			tumors_dict[tumor]["constraints"] = []
			tumors_dict[tumor]["states"] = []
		tumors_dict[tumor]["constraints"].append(constr_dict[tumor]["CellProliferation"])
		tumors_dict[tumor]["states"].append(states_dict[tumor][myc_id])

	print "Avg val of rand optimized vector: " + str(np.mean(np.array(opti_vector)))
	print "Standard deviation of sim optimized vector: " + str(np.std(np.array(opti_vector)))

avg_rand = sum_rand/len(files)
st_rand = np.std()
# import pdb; pdb.set_trace()

plt.vlines(range(50),[-1], avg_rand[:50])
plt.vlines(np.arange(50)+0.2, [-1], sim_opti_vector[:50], colors='r')

plt.show()
# import pdb; pdb.set_trace()
# run_statistical_analysis(total_results, constr_dict, states_dict, tumors_dict, fitness_list)
	

	# for f in files:
	# tmp_file = output_path + f
	# shutil.move(tmp_file, sim_file)
	# constr_dict, states_dict, fitness, tmp1, tmp2 = parse_simulation_file(tmp_file)

	# total_results[f] = Odict()
	# total_results[f]["constraints"] = constr_dict
	# total_results[f]["states"] = states_dict
	# total_results[f]["fitness"] = fitness
	# subprocess.call(["python", plot_exe, "-c" + config_file])
	# shutil.move(sim_file, tmp_file)
	# shutil.move(output_path + "analyzed_simulation_results.txt", output_path + "an_res_"+f[11:-4]+".txt")
	# shutil.move(output_path + "cost_func_evaluation.txt", output_path + "cost_res_"+f[11:-4]+".txt")
	# shutil.move(output_path + "Function_id_[3274].png",  output_path + "cost_plot_"+f[11:-4]+".png")
	# shutil.move(output_path + "MAX-001:P[S62]-MYC-001 [comple_[complex in nucleus]_id_3274.png", output_path + "myc_plot_"+f[11:-4]+".png")


# import pdb; pdb.set_trace()

