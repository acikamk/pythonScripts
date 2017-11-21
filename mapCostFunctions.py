import pandas
from collections import OrderedDict as Odict
import collections
import sys 
import pdb
'''
Function that modifies/removes/adds parameters
to specific index places in optimized vectors
Inputs:
 - model_path =  path to the model structure
 - funct = remove, add or modify parameters
 - val =  which value to bi assign for add
 - ids = the pybios_ids model ids to bi handled
Output:
 - File
			+ 1.0										 +kCD1	      kCD5 		1.0          kCD6       kCD11      kCD8		    kCD9       kCD10       kCD7 		 kCD12      kCD13      kCD14      kCD15
 (kCD0*[131]+kCD1*[2227]+kCD2*[130]+kCD3*[902]+kCD4*[326]+kCD5*[1293]+kCD6*[115]+kCD7*[2229]+kCD8*[325]+kCD9*[129]+kCD10*[919]+kCD11*[127]+kCD12*[829]+kCD13*[323])/(kCD14*[18]+kCD15*[19]+kCD16*[20]+kCD17*[21])
 (kCD0*[131]+kCD1*[1293]+kCD2*[130]+kCD3*[902]+kCD4*[326]+kCD5*[115]+kCD6*[325]+kCD7*[323]+kCD8*[919]+kCD9*[127]+kCD10*[829]+kCD11*[129])/(kCD12*[18]+kCD13*[19]+kCD14*[20]+kCD15*[21])
			+kCD5*										 +kCD6        kCD8      kCD13      kCD19	  kCD11      kCD12       kCD9		   kCD14      kCD15      kCD16      kCD17
'''

def get_formula(file):
	with open(file) as f:
		for line in f:
			if line.startswith("ID"):
				cost_str = line.strip().split('\t')[-1]
				cost_str = cost_str.replace("(","").replace(")","")
				break
	return cost_str

def get_costDict(cost_str):
	param_dict = Odict()
	for nominator in cost_str.split("/"):
		for pair in nominator.split("+"):
			if pair.split("*")[1] in param_dict.keys():
				raise ValueError("Error! Same species present in cost multiple times!")
			# pdb.set_trace()
			param_dict[pair.split("*")[1]] = pair.split("*")[0]
	return param_dict

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
	'''
	modify model_1 progress_trace (ids and vectors) to appropriate
	size and format so they can be used on model_2
	The size of the progress trace should be same as the vector_ids 
	of pt of model_2. Results are saved at model_1 path 
	'''

	model_output_1 = model_path_1 +"/outputs/"
	cost_1_file = model_path_1 + "/inputs/cost_func.txt"
	prog_1_file = model_output_1 +"progress_trace.txt"

	cost_2_file = model_path_2 + "/inputs/cost_func.txt"
	prog_2_file = model_path_2 +"/outputs/progress_trace.txt"


	cost_1 = get_formula(cost_1_file)	
	cost_1_dict = get_costDict(cost_1)
	ids_vect_1, opti_vect = get_progTrace(prog_1_file)

	cost_2 = get_formula(cost_2_file)
	cost_2_dict = get_costDict(cost_2)
	ids_vect_2, _ = get_progTrace(prog_2_file)
	
	# dict in form {'species_id' : [kCDmodel1, kCDmodel2]}
	# if not kCDmodel1 for species_id then x
	# if not kCDmodel2 for species_id then ommited, not present
	merged_dict = Odict()
	for k in cost_2_dict.keys():
		merged_dict[k] = [];
		if k in cost_1_dict:
			merged_dict[k].append(cost_1_dict[k])
		else:
			merged_dict[k].append('x')
		merged_dict[k].append(cost_2_dict[k])

	
	res_vect = []
	for vect in opti_vect:
		tail = vect[len(cost_1_dict.keys()):]
		# fill kCDs with 1.0 values with size of resulting vector
		head = [1.0 for i in range(0,len(cost_2_dict))]
		
		for k, v in merged_dict.items():
			if not v[0] == 'x':
				# pdb.set_trace()
				head[ids_vect_2.index(v[1])] = vect[ids_vect_1.index(v[0])] 
		res_vect.append(head + tail)

	res_ids = ids_vect_2[:len(cost_2_dict)] + ids_vect_1[len(cost_1_dict):]
	# pdb.set_trace()	
	mod_prog_file = model_output_1 + "swaped_progress_trace.txt"
	mod_prog = open(mod_prog_file, "w")

	mod_prog.write("VARIABLE PAR IDS: {}\n\n".format(res_ids))



	for vect in res_vect:
		mod_prog.write("VECTOR: {}\n\n".format(vect))

	mod_prog.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])