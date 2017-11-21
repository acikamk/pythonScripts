import pandas
from collections import OrderedDict as Odict
import collections
import sys 
import pdb
import random
import os
'''
Shuffle optimized vectors x 5
'''
def main(input_file):


	opti_vect = []

	with open(input_file) as f:	
		for line in f:
			if line.startswith("VARIABLE PAR IDS:"):
				ids_vect = eval(line[line.index("["):])
			if line.startswith("VECTOR:"):
				vect = eval(line[line.index("["):])
				assert len(ids_vect) == len(vect) 
				opti_vect.append(vect)

	mod_prog_file = os.path.dirname(os.path.abspath(input_file)) + "/shuffled_progress_trace.txt"
	mod_prog = open(mod_prog_file, "w")
	mod_prog.write("VARIABLE PAR IDS:\t" + str(ids_vect) + "\n\n")
	for vect in opti_vect:
		for i in xrange(5):
			random.shuffle(vect)
			mod_prog.write("VECTOR: {}\n\n".format(vect))	
	mod_prog.close()

if __name__ == "__main__":
    main(sys.argv[1])