import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
import compiler
import re
import random

def write_results(key, value, f, absolute, controls):

	if absolute:	
		f.write(key + '\t' + str(value)+ '\n')
	else:
		if key in controls.keys():
			controls[key] = value
		else:
			if any(control in key for control in controls.keys()):	 
				tumor = eval(reduce(lambda x, y: x.replace(y, str(conc_dict_samples[key][y])), conc_dict_samples[key], cost_function_noK))
				c, control = [(c, controls[c]) for c in controls.keys() if key.startswith(c)][0]
				value = tumor/control
				f.write(key + '/' + c  + '\t' + str(value)+ '\n')
			else:
				print "ERROR: key"  + "control" + "not found!\n"
	return


study_path = "/project/V0001-1/modcell/modcell_data/project_manager/7_OncoTrack/2017.03.20_15.32.00_573_Model_Test_for_Optimization/" 

tabfiles_path = study_path + "Studies/2467_Testing_Model_Level1/AnalysisResults/Tabfiles/"
mapping_table = study_path + "Studies/2467_Testing_Model_Level1/MappingTable.txt"
identifiers_file = study_path + "Model/identifiers.py"
n_sim = 3
cost_function = "(k0*[205]+k1*[710]+k2*[1180]+k3*[1658]+k4*[1638]+k5*[207]+k6*[1599])/(k7*[1028]+k8*[1025]+k9*[1074]+k10*[1027]+k11*[1026]+k12*[1077])"
cost_function_name = "CellDivisionFunction"
output_file = "cost_func.txt"
absolute = False
training = True

cost_function_noK = re.sub('k(\d+)\*', '', cost_function)
pybios_ids = re.findall(r"\[\d+\]", cost_function_noK)
nominator = cost_function.split('/')[0][1:-1]
denominator = cost_function.split('/')[1][1:-1]

# parse indentifiers.py file
f = open(identifiers_file)
identifiers_code = compile(f.read(), identifiers_file, 'exec')
f.close()
eval(identifiers_code)
diff_species = locals()['diff_var_ids']
diff_sp_dict =  {j[1]: k for k,j in diff_species.iteritems()}  

# parse mapping table
f = open(mapping_table)
info = f.readlines()
f.close()
tumor_files = {"Tabfile_" + i.split()[0] + ".txt" : i.split()[1] \
				for i in info if "TUMOR" in i}

# analyse tabfiles and create conc_dict_samples 
# with samples names as keys and new
# dictionary with species names as key and 
# average concentration as values

conc_dict_samples = Odict()
controls = Odict()

for file in tumor_files.keys():
	try:
		if not tumor_files[file] in conc_dict_samples.keys():
			conc_dict_samples[tumor_files[file]] = Odict()
			temp_dict = {}
			if "CONTROL" in tumor_files[file] or not "_" in tumor_files[file]:
				controls[tumor_files[file]] = None
		f=open(tabfiles_path + file)
		for line in f:
			temp_dict[diff_sp_dict[' '.join(line.split()[0:-n_sim])]]= sum(float(x) for x in line.split()[-n_sim:])/n_sim
			conc_dict_samples[tumor_files[file]] = temp_dict
		f.close()
	except:
			print "Error parsing " + file
			f.close()

# import pdb; pdb.set_trace()
if training:
	bucket_size = len(controls.keys())
	training_size = int(round(bucket_size*80/100))
	training_list = random.sample(controls.keys(), training_size)
	print sorted(training_list)
	testing_list = list(set(controls.keys())-set(training_list))
	print sorted(testing_list)

import pdb; pdb.set_trace()

# print results to file
if training:
	f_train = open(output_file[:-4]+"_train.txt", 'w') 
	f_test = open(output_file[:-4]+"_test.txt", 'w') 
	f_train.write("ID:" + cost_function_name + '\t' + cost_function + '\n')
	f_test.write("ID:" + cost_function_name + '\t' + cost_function + '\n')
else:
	f = open(output_file, 'w')
	f.write("ID:" + cost_function_name + '\t' + cost_function + '\n')

# import pdb; pdb.set_trace()
for key in sorted(conc_dict_samples.keys()):
	value = eval(reduce(lambda x, y: x.replace(y, str(conc_dict_samples[key][y])), conc_dict_samples[key], cost_function_noK))
	if training:
		if key.split("_")[0] in training_list:
			write_results(key, value, f_train, absolute, controls)
		elif key.split("_")[0] in testing_list:
			write_results(key, value, f_test, absolute, controls)
		else:
			print "Key Error :)! " +key +"\n"
	else:
		write_results(key, value, f, absolute, controls, conc_dict_samples)

if training:
	f_train.close() 
	f_test.close()
else:
	f.close()



# formula = re.compile(r'\b(' + '|'.join(conc_dict_samples[key].keys()) + r')\b')
#  # re.compile('|'.join(d.keys()))
# result = pattern.sub(lambda x: conc_dict_samples[key][x.group()], cost_function_noK)

# def get_value(formula, S, F):
#     return eval(formula)
# def evaluateFormula(formula)
# 	for 


# conc_dict_samples[tumor_files[file]] = \
# 	[[float(x) for x in i.split()[-n_sim:]] for i in info]
# keys = [' '.join(i.split()[0:-n_sim]) for i in info]
# vals =[sum([float(x) for x in i.split()[-n_sim:]])/n_sim for i in info]
# if not conc_dict:
# 	for k in keys:
# 		conc_dict[k] = []
# [conc_dict[k].append(v) for k, v in zip(keys,vals)]