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
 - val =  which value to be assigned for add
 - ids = the pybios_ids model ids to bi handled
Output:
 - File = progress trace with modified parameters
'''
def main(model_path, funct, val, ids):

	model_output = model_path+"/outputs/"
	prog_file = model_output+"progress_trace.txt"
	opti_vect = []
	with open(prog_file) as f:	
		for line in f:
			if line.startswith("VARIABLE PAR IDS:"):
				ids_vect = eval(line[line.index("["):])
			if line.startswith("VECTOR:"):
				vect = eval(line[line.index("["):])
				assert len(ids_vect) == len(vect) 
				opti_vect.append(vect)
	
	ids_ind = []			
	for mid in ids:
		ids_ind.extend([ids_vect.index(i) for i in ids_vect if mid==i])

	if funct == 'remove':
		ids_vect = [i for j, i in enumerate(ids_vect) if j not in ids_ind]
	if funct == 'add':
		acc = 0
		for i in range(len(ids_ind)):
			name = 'kCDn'+str(ids_ind[i])
			ids_vect.insert(ids_ind[i], name)
			# acc += 1
	res = []
	for vect in opti_vect:
		if funct == 'remove':
			res_v = [i for j, i in enumerate(vect) if j not in ids_ind]
		if funct == 'add':
			res_v = vect
			acc = 0
			for i in range(len(ids_ind)):
				res_v.insert(ids_ind[i], float(val))
				# acc+=1
		else:	
			res_v = [i if j not in ids_ind else float(val) for j, i in enumerate(vect)] 
		res.append(res_v)

	# pdb.set_trace()
	mod_prog_file = model_output + str(funct) + "_progress_trace.txt"
	mod_prog = open(mod_prog_file, "w")

	mod_prog.write("VARIABLE PAR IDS:\t" + str(ids_vect) + "\n\n")



	for vect in res:
		mod_prog.write("VECTOR: {}\n\n".format(vect))

	mod_prog.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])