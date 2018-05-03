import numpy as np
import datetime
import re
import os
import sys
import math
import pdb

# model_path = "/project/V0001-1/modcell/users/a.kovachev/"
# numb_vect = 20

def main(model_path, numb_vect):

	model_output = model_path+"outputs/"
	prog_file = model_output+"simulated_vectors_file.txt"

	with open(prog_file) as f:	
		for line in f:
			if line.startswith("VARIABLE PAR IDS:"):
				ids_vector = eval(line[line.index("["):])

	n_param = len([i for i in ids_vector if "{" in i])
	n_kCD = len([i for i in ids_vector if i.startswith("k")])

	tmp_prog_file = model_output + "random_progress_trace.txt"
	tmp_prog = open(tmp_prog_file, "w")

	tmp_prog.write("VARIABLE PAR IDS:\t" + str(ids_vector) + "\n\n")
	# pdb.set_trace()
	boundsK = (-1,1)
	boundsKCD = (1.e-4, 1)
	for i in xrange(int(numb_vect)):
		vectorKCD = boundsKCD[0] + np.random.rand(int(n_kCD))	\
			*(boundsKCD[1] - boundsKCD[0]) 	
		vectorK = boundsK[0] + np.random.rand(int(n_param)) \
			*(boundsK[1] - boundsK[0])  
		tmp_prog.write("VECTOR: {}\n\n".format(np.concatenate((vectorKCD, vectorK)).tolist()))

	tmp_prog.close()
	return

if __name__ == "__main__":
	try:
		main(sys.argv[1], sys.argv[2])
	except:
		print "The function has 2 input arguments: path to a optimized model" 
		print "and number of random vectors which need to be generated as output from an optimization."
		print "Usage ex.: pyhton createRandomVectProgTrace.py /project/V0001-1/modcell/users/a.kovachev/ 20 "