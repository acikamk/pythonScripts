import pandas
import sys
import re
import os
import pdb

"""
Script that converts experiment table format from simtools to AMICI.
It modifies the IDs column to SBML_IDs. 
NOTE: Does not 100% garanties that the ids will match.

Input:
 - Simtools experiment table file
 - Pybios indentifiers.py file

Output:
 - exp_tab_name_sbmlIDS.csv - Experiment table wiht sbml_ids

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


def convert2sbmlIDS(exp_tb, cost_f, iden_py):
	compartments = {'cellular component]' : '1',
					'extracellular]': '2',
					'golgi apparatus]': '11',
					'endoplasmic reticulum]':'12',
					'peroxisome]': '13',
					'plasma membrane]': '3',
					'intracellular]':'4',
					'cytoplasm]':'5',
					'nucleus]':'6',
					'nucleolus]':'7',
					'ribosome]':'8',
					'mitochondrion]':'9',
					'lysosome]':'10'}
                       
	f = open(iden_py)
	try:
		identifiers_code = compile(f.read(), iden_py, 'exec')
		f.close()
		eval(identifiers_code)
		# pdb.set_trace()
		diff_sp_dict = locals()['diff_var_ids']
		fixed_sp_dict = locals()['fixed_var_ids']
		pars_sp_dict = locals()['par_ids']
		
		
		sbml_dict_diff =  {k:"SP_"+ str(j[3]) + "_" + \
			compartments[j[1].split(" in ")[1]] for k,j in diff_sp_dict.iteritems()}   
		sbml_dict_fixed =  {k: "SP_"+ str(j[3]) + "_" + \
			compartments[j[1].split(" in ")[1]] for k,j in fixed_sp_dict.iteritems()} 

		unique_parids = sorted(set([k[1:-1].split('|')[0] for k in pars_sp_dict.keys()]))
		sbml_dict_pars = {k: "par_" + k[1:-1].split('|')[0] + "_" + \
			k[1:-1].split('|')[-1] + "_reaction_" + k[1:-1].split('|')[0]\
			# k[1:-1].split('|')[-1] + "_reaction_" + str(unique_parids.index(k[1:-1].split('|')[0])+1)\
			for k,j in pars_sp_dict.iteritems()}    
		
	except:
		print  bcolors.FAIL \
		+ "The identifiers.py file has incompatible structure!" \
	 	+ bcolors.ENDC
		return

	try:
		# import pdb; pdb.set_trace()
		if isinstance(exp_tb, pandas.DataFrame):
			exp_df = exp_tb
		else:
			exp_df = pandas.read_table(exp_tb, sep='\t', index_col=0)
	except:
		print bcolors.FAIL \
		+ "Inproper file format!" \
		+ bcolors.ENDC
		return

	as_list = exp_df.index.tolist()	
	for i, el in enumerate(as_list):
		if el in sbml_dict_fixed.keys():
			as_list[i] = sbml_dict_fixed[el] 
		elif el in sbml_dict_diff.keys():
			as_list[i] = sbml_dict_diff[el]
		elif el in sbml_dict_pars.keys(): 
			as_list[i] = sbml_dict_pars[el]
		else:
			print "Unknown array ID: '{}'".format(el)

	exp_df.index = as_list
	exp_df.index.rename('ID', inplace=True) 
	exp_file = exp_tb[:exp_tb.rindex(".")] + "_sbmlIDS.csv"
	exp_df.to_csv(exp_file, sep='\t')
	print bcolors.OKGREEN \
	+ "Conversion of experiment table suceeded!" \
	+ bcolors.ENDC
	
	cost_file = open((os.path.splitext(cost_f)[0])+"_sbmlIDS.txt", "w")
	
	try:
		with open(cost_f) as f:
			for l in f:
				if l.startswith('ID'):
					cost_id = l.split()[0]
					cost_form = l.split()[1]
					cost_sbml = re.sub('\[\d+\]', lambda x: sbml_dict_diff[x.group()], cost_form)
					# import pdb; pdb.set_trace()
					cost_file.write(cost_id + "\t" + cost_sbml + "\n")
				else:
					cost_file.write(l)
		print bcolors.OKGREEN \
		+ "Conversion of cost function suceeded!" \
		+ bcolors.ENDC
		# import pdb; pdb.set_trace()
	except:
		print bcolors.FAIL \
		+ "The cost function file was not parsed correctly!" \
	 	+ bcolors.ENDC
		return
	return

if __name__ == "__main__":
	try:
		convert2sbmlIDS(sys.argv[1], sys.argv[2], sys.argv[3])
	except:
		print bcolors.FAIL \
		+ "Please insert the input files in correct order:" \
		+ "experiment_table_file, cost_function_file, identifiers.py_file." \
	 	+ bcolors.ENDC
		
	