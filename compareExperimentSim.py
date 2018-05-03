import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
"""
Script that compares the final steady state concentrations of the species 
from optimized vector simulation with Project Manager and simtools.
Input: 
study_path = path to the project manager study 
simulation_path = path to the simtools forward simulaton

optional:
identifiers_file = path to the pybios model indentifiers.py file

output:
stout output
The absolute and relative tolerances of the total species concentrations 
and the sample vise species concetrations
"""

study_path = "/project/V0001-1/modcell/modcell_data/project_manager/5_TESTTEST/2016.12.06_09.34.29_499_Notch_Model_Testing/Studies/2605_Notch_simulation_test" 
#"/project/V0001-1/modcell/modcell_data/project_manager/5_TESTTEST/2017.03.03_13.37.29_554_Speedy_r378088_tests/Studies/2332_Speedy_r378088_sim/"
simulation_path = "/project/V0001-1/modcell/users/a.kovachev/Notch_Signalling_r377624/Project_Manager_Tests/Notch_Debug_t_f_sim_test4/"
#"/project/V0001-1/modcell/users/a.kovachev/Speedy_v1_r378088/Speedy_v1_r378088_sim_02_02_17_rel"

tabfiles_path = study_path + "/AnalysisResults/Tabfiles/"
mapping_table = study_path + "/MappingTable.txt"
simulated_file = simulation_path + "/outputs/simulated_vectors_file.txt"
f = open(mapping_table)
info = f.readlines()
f.close()
tumor_files = {"Tabfile_" + i.split()[0] + ".txt" : i.split()[1] \
				for i in info if "TUMOR" in i}

conc_dict = Odict()
conc_dict_samples = Odict()

# analyse tabfiles and create sample : conc dict and species : conc dict
for file in tumor_files.keys():
	# try:

	f =open(tabfiles_path + file)
	info = f.readlines()
	conc_dict_samples[tumor_files[file]] = \
		[float(i.split()[-1]) for i in info]
	keys = [i.split('\t')[0] for i in info]
	vals = [float(i.split('\t')[-1]) for i in info]
	# import pdb; pdb.set_trace()
	if not conc_dict:
		for k in keys:
			conc_dict[k] = []
	[conc_dict[k].append(v) for k, v in zip(keys,vals)]
	f.close()
	# except:
	# 		print "Error parsing " + file

# analyse simulated_vector_file.txt and create state : conc dict
states_dict = {}
f = open(simulated_file)
info = f.readlines()
states_dict = {i[i.index(":")+1:i.rindex(":")]: eval(i[i.rindex(":")+1:])\
 				for i in info if "STATES" in i}
f.close()

abs_tol = Odict()
rel_tol = Odict()
tab_total_sum = sum([j for i in conc_dict_samples.values() for j in i])
stat_total_sum = sum([j for i in states_dict.values() for j in i])  
# import pdb; pdb.set_trace()
print "Tabfiles sum: " + str(tab_total_sum) + "\tSimtools sum:" + str(stat_total_sum) + "\n"
print "TOTAL SUM ABSOLUTE TOLERANCE: {0:.5f}\tRELATIVE TOLERANCE: {1:.5f}".\
		format(abs(tab_total_sum-stat_total_sum),
			   abs(1-tab_total_sum/stat_total_sum))
try:
	tab_value = max(len(k) for k in states_dict.keys())
	print "\n{1:{0}}\tabs tol\t\t\treal tol".format(tab_value, 'Sample Name')
	for k in states_dict.keys():		
		abs_tol[k] = abs(sum(states_dict[k])-sum(conc_dict_samples[k]))
		rel_tol[k] = abs(1-sum(states_dict[k])/sum(conc_dict_samples[k]))
		print "{1:{0}}\t{2:5.5f}\t\t{3:5.5f}".\
				format(tab_value, k, abs_tol[k], rel_tol[k])

except:
	print "SAMPLE NAMES IN SIMULATED FILE AND MAPPING TABLE NOT CONSISNTENT!"

#### additional code for macthing tabfiles to indetifiers.py
#### extract species data from identifiers.py

# identifiers_file = "/project/V0001-1/modcell/modcell_data/project_manager/5_TESTTEST/2016.12.06_09.34.29_499_Notch_Model_Testing/Model/identifiers.py"
# f = open(identifiers_file)
# identifiers_code = compile(f.read(), identifiers_file, 'exec')
# f.close()
# eval(identifiers_code)
# diff_species = locals()['diff_var_ids']
# diff_sp_dict =  {j[1]: k for k,j in diff_species.iteritems()}   

# # create pybios_id : conc dict 
# spec_conc_dict = {v1: conc_dict[k] for k, v1 in diff_sp_dict.iteritems()}
