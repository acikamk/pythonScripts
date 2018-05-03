import pandas as pd
from collections import OrderedDict as Odict
import collections
import sys 
import pdb
'''
Compares two outputs from forward simulation and 
input experiment tables for different runs.
Checks for additional parameters and their values.

'''

def get_progTrace(file):
	opti_vect = []
	with open(file) as f:	
		for line in f:
			if line.startswith("VARIABLE PAR IDS:"):
				ids_vect = eval(line[line.index("["):])
			if line.startswith("VECTOR:"):
				vect = eval(line[line.index("["):])
				assert len(ids_vect) == len(vect) 
				opti_vect.append(vect)
	return ids_vect, opti_vect

def main(model_path_1, model_path_2):

	prog_1_file = model_path_1 +"/outputs/simulated_vectors_file.txt"
	exp_1_file  = model_path_1 + "/inputs/exp_table.csv"
	prog_2_file = model_path_2 +"/outputs/simulated_vectors_file.txt"
	exp_2_file  = model_path_2 + "/inputs/exp_table.csv"

	exp_1 = pd.read_table(exp_1_file, index_col=0, low_memory=False)
	exp_2 = pd.read_table(exp_2_file, index_col=0, low_memory=False)

	ids_vect_1, opti_vect_1 = get_progTrace(prog_1_file)
	ids_vect_2, opti_vect_2 = get_progTrace(prog_2_file)

	# check model sizes
	if len(ids_vect_1) > len(ids_vect_2):
		prim = ids_vect_1 
		prim_v = opti_vect_1
		sec = ids_vect_2 
		sec_v = opti_vect_2
		print model_path_1 + " has larger size vector."
	elif len(ids_vect_1) < len(ids_vect_2):
		prim = ids_vect_2
		prim_v = opti_vect_2
		sec = ids_vect_1  
		sec_v = opti_vect_1
		print model_path_2 + " has larger size vector."
	else:
		prim = ids_vect_1 
		prim_v = opti_vect_1
		sec = ids_vect_2 
		sec_v = opti_vect_2
		print "Models have same size. " + model_path_1 + " was sellected as primary."
	vals_dict = Odict()
	for k in prim:
		if not k in vals_dict:
			vals_dict[k] = []
		# pdb.set_trace()
		for i in range(len(prim_v)):
			if k in sec:
				vals_dict[k].append((prim_v[i][prim.index(k)], sec_v[i][sec.index(k)]))
				if abs(prim_v[i][prim.index(k)]-sec_v[i][sec.index(k)]) > 0.00001:
					print "Problem found with parameters!"
					pdb.set_trace()
			else:
				vals_dict[k].append((prim_v[i][prim.index(k)], ""))
				print k + " not in the secondary vector." 
				if prim_v[i][prim.index(k)] != 1.0:
					print "Problem found with drugTransl!"
					pdb.set_trace()		
	# check experiment tables
	if exp_1.columns.tolist() != exp_2.columns.tolist():
		print "Exp tables columns are not equal..."
		pdb.set_trace()	

	if exp_1.index.tolist() != exp_2.index.tolist():
		l1 = [i for i in exp_1.index.tolist() if i not in exp_2.index.tolist()]
		l2 = [i for i in exp_2.index.tolist() if i not in exp_1.index.tolist()]
		print "Additionaly in exp_1.index "+ l1
		print "Additionaly in exp_2.index "+ l2
		# pdb.set_trace()				
	for el in l1:
		print [ i for i in exp_1.loc[el].values if i!=1.0]

	for el in l2:
		print [ i for i in exp_2.loc[el].values if i!=1.0]

	pdb.set_trace()

if __name__ == "__main__":
	try:
		main(sys.argv[1], sys.argv[2])
	except:
		print "This function takes two input arguments, paths to simulated model folders."
		print "Usage ex.: py compare2inVecotors.py path_to_simulated_model1 path_to_simulated_model2"