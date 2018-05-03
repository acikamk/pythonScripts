import numpy as np
import math
import pdb

''' Script that calcualtes the error and correlation between
random starting vector used to generate simulated data
and optimized vectors.
Input:
	- input_vectors_file: file where information about
	the random vector is stored - as outputed from
	runSteadyStateOpti.py, where info for the random (or optimized)
	vector that was used for generation of simulated data is given
	- progress_trace.txt file as output from simtools
	were the fitted vectors are present
	- last (int): Default(-1) parameter that tells which vectors from the 
		optimized vectors file to be taken.. possible options: -1 (last),
		0(fist) -2 (second last), else all vectors
Output:
	relative error and correlation
'''

def main(input_vectors_file, progress_trace, last):
	last = int(last)
	if not (last == -1 or last == 0 or last == -2):
		last = -3
	f = open(input_vectors_file)
	data = f.readlines()
	par_ids = eval([l for l in data if l.startswith('PAR_IDS')][0].split('=')[-1])
	var_ids = eval([l for l in data if l.startswith('VAR_PAR_IDS')][0].split("=")[-1])
	vector = eval([l for l in data if l.startswith('RND_VECTOR')][0].split('=')[-1])
	var_idx = eval([l for l in data if l.startswith('VAR_PAR_IDX')][0].split('=')[-1])
	rnd_vect = [j for i, j in enumerate(vector) if i in var_idx] 
	assert len(var_ids) == len(rnd_vect), "Something wrong with input data!"
	f.close()
	 	
	f = open(progress_trace)
	data = f.readlines()
	var_par_ids = eval([l for l in data if \
		 l.startswith('VARIABLE PAR IDS:')][0].split('\t')[-1])
	if last !=-3:
		all_vect = [eval([l for l in data if \
			l.startswith('VECTOR:')][last].split('\t')[-1])]
	else:
		all_vect = [eval(i.split('\t')[-1]) for i in \
			[l for l in data if l.startswith('VECTOR:')]]
	f.close()
	
	# calculates the relative error and prints some results
	
	for v in all_vect:
		c=0
		a = []
		b = []
		for var in var_par_ids:
			if var in var_ids:
				# print rnd_vect[var_ids.index(var)], math.pow(10, last_vector[var_par_ids.index(var)])
				# print abs(1-math.pow(10, last_vector[var_par_ids.index(var)])/rnd_vect[var_ids.index(var)])
				rel_err = abs(1-math.pow(\
					10, v[var_par_ids.index(var)])/rnd_vect[var_ids.index(var)])
				
				if rel_err < 0.03:
					print var
					# print abs(1-math.pow(10, last_vector[var_par_ids.index(var)])/rnd_vect[var_ids.index(var)])
					# print rnd_vect[var_ids.index(var)], math.pow(10, last_vector[var_par_ids.index(var)])
					c+=1
				a.append(rnd_vect[var_ids.index(var)])
				b.append(math.pow(10, v[var_par_ids.index(var)]))

		if not len(a)<4:
			print c
			print np.corrcoef(a, b)
		else:
			print [(i, j, abs(1-i/j)) for (i,j) in zip(a, b)]  


if __name__ == "__main__":
	try:
		main(sys.argv[1], sys.argv[2], sys.argv[3])
	except:
		print "Script that calcualtes the error and correlation between \
		random starting vector used to generate simulated data \
		and optimized vectors. \n\
		Input: \n\
			\t- input_vectors_file: file where information about \n\
			the random vector is stored - as outputed from \n\
			runSteadyStateOpti.py, where info for the random (or optimized) \n\
			vector that was used for generation of simulated data is given \n\
			\t- progress_trace.txt file as output from simtools \
			were the fitted vectors are present \n\
			\t- last (int): Default(-1) parameter that tells which vectors from the \n\
				optimized vectors file to be taken.. possible options: -1 (last), \n\
				0(fist) -2 (second last), else all vectors\n"