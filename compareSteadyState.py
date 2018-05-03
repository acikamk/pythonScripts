import numpy
import math
import sys
import collections 
from collections import OrderedDict as Odict
import matplotlib.pyplot as plt 

'''
Script that compares two simulated files, after optmizations
which where performed on simulated dataset.
Inputs:
	- sim files list, 
	- initialization file, where the original parameter vector 
		is present
	-x, starting offset to which values to be plotted
	-y, ending offset to which values to be plotted
	  necessary since plotting of > 30/40 parameters
	  is not comperhensive 
Usage ex.:
	python compareSteadyState.py '/project/V0001-1/modcell/users/a.kovachev/Notch_Signalling_r377624/Notch_Debug_ss_sim_test_0/outputs/simulated_vectors_file.txt /project/V0001-1/modcell/users/a.kovachev/Notch_Signalling_r377624/Notch_Debug_ss_sim_test_1/outputs/simulated_vectors_file.txt' /project/V0001-1/modcell/users/a.kovachev/data/data/init_vector_Notch_maml_abs.txt 0 5

'''

colors = ['b', 'g', 'c', 'm', 'r']
def main(sim_files, init_file, x, y):
	x=int(x)
	y=int(y)
	print "Parsing init_vector file"
	init_vector = Odict()
	with open(init_file) as file:
		for line in file:
			if line.startswith("PAR_IDS="):
				vids_vector = eval(line[line.index("["):])
			elif line.startswith("RND_VECTOR="):
				# scale >0
				rnd_vector = eval(line[line.index("["):])		
	for key, val in zip(vids_vector, rnd_vector):
		init_vector[key] = val
	sim_files =  sim_files.split()
	fig = plt.figure()
	n=0
	sim_vector = Odict()
	sim_vect_rev = Odict()
	for sim_file in sim_files:	
		print "Parsing simulated_vectors.txt file"
		sim_vector[n] = Odict()
		bounds = (-1,1)
		with open(sim_file) as file:
			for line in file:
				if line.startswith("VARIABLE PAR IDS:"):
					ids_vector = eval(line[line.index("["):])
				elif line.startswith("VECTOR:"):
					# scale log10
					opti_vector = eval(line[line.index(":")+1:])		
		opti_vector = numpy.power(10, numpy.asarray(opti_vector))
		for key,val in zip(ids_vector, opti_vector):
			sim_vector[n][key] = val
			if key not in sim_vect_rev:
				sim_vect_rev[key] = []
			sim_vect_rev[key].append(val)
		# import pdb; pdb.set_trace()
		merge_vect = {k: (init_vector[k], sim_vector[n][k], abs(init_vector[k]-sim_vector[n][k])) for k in init_vector if k in sim_vector[n]}
		
		plt.plot(range(1,y-x+1),[m[1] for m in merge_vect.values()[x:y]], 'x', color=colors[n], markeredgewidth=1.5)
		n=n+1
		plt.hold(True)
	plt.plot(range(1, y-x+1), [m[0] for m in merge_vect.values()[x:y]], 'ro', markersize='7')
	plt.xticks(range(1, y-x+1), merge_vect.keys()[x:y], rotation=90)
	ax = fig.gca()
	ax.set_xticks(numpy.arange(1, y-x+2, 1))
	plt.grid()
	
	## some additional alaysis
	## calculate  intersect between results
	# init_vector_mod = {k: init_vector[k] for k in init_vector if k in sim_vector[0]}
	# match_rand = []

	## find parameters that are far from the random one but simlar between them
	# for ids in sim_vect_rev:
	# 	if ids in init_vector_mod:
	# 		if abs(sum(numpy.diff(sim_vect_rev[ids])))<0.2 \
	# 		and abs(init_vector_mod[ids]-sim_vect_rev[ids][1])>2\
	# 		 and abs(init_vector_mod[ids]-sim_vect_rev[ids][3])>2 \
	# 		 and abs(init_vector_mod[ids]-sim_vect_rev[ids][0])>2 :
	# 			match_rand.append(ids)	
	# # import pdb; pdb.set_trace()	

	## find 4 different groups of parameters
	# match = {}
	# almost = {}
	# far = {}
	# very_far = {}

	# for n in sim_vector:
	# 	match[n]=[]
	# 	almost[n]=[]
	# 	far[n]=[]
	# 	very_far[n]=[]
	# 	for k in init_vector_mod:
	# 		if abs(init_vector_mod[k]-sim_vector[n][k])<0.05:
				
	# 			match[n].append(k)
	# 		elif abs(init_vector_mod[k]-sim_vector[n][k])<0.5 and abs(init_vector_mod[k]-sim_vector[n][k])>0.05:
				
	# 			almost[n].append(k)
	# 		elif 2<abs(init_vector_mod[k]-sim_vector[n][k])>0.5:

	# 			far[n].append(k)
	# 		else:
	# 			very_far[n].append(k)
	## import pdb; pdb.set_trace()	

	# merge_vect = {k: (init_vector[k], sim_vector[k], abs(init_vector[k]-sim_vector[k])) for k in init_vector if k in sim_vector}
	# merge_vect = Odict(sorted(merge_vect.items()))

	# fig = plt.figure()
	# plt.boxplot(merge_vect.values()[0:50])
	# plt.xticks(range(50), merge_vect.keys()[0:50], rotation=90)
	# plt.grid()
	# # fig = plt.figure()
	# plt.plot(range(1, 51),[x[0] for x in merge_vect.values()[0:50]], 'ro')
	# plt.plot(range(1, 51),[x[1] for x in merge_vect.values()[0:50]], 'bx')
	# plt.grid()
	# plt.xticks(range(1, 51), merge_vect.keys()[0:50], rotation=90)
	# # plt.boxplot(merge_vect.values()[51:100], labels=merge_vect.keys()[51:100])
	# # fig = plt.figure()
	# # plt.boxplot(merge_vect.values()[101:150], labels=merge_vect.keys()[101:150])
	# # fig = plt.figure()
	# # plt.boxplot(merge_vect.values()[151:], labels=merge_vect.keys()[151:])
	# # import pdb; pdb.set_trace()
	
	plt.show()

if __name__ == "__main__":
	try:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	except:
		print "Usage ex:"
		print "python compareSteadyState.py '/project/V0001-1/modcell/users/a.kovachev/Notch_Signalling_r377624/Notch_maml_abs_test_0/outputs/simulated_vectors_file.txt /project/V0001-1/modcell/users/a.kovachev/Notch_Signalling_r377624/Notch_maml_abs_test_1/outputs/simulated_vectors_file.txt' /project/V0001-1/modcell/users/a.kovachev/data/data/init_vector_Notch_maml_abs.txt 0 45"
		print "python compareSteadyState.py '/project/V0001-1/modcell/users/a.kovachev/Notch_Signalling_r377624/Notch_maml_abs_test_0/outputs/simulated_vectors_file.txt /project/V0001-1/modcell/users/a.kovachev/Notch_Signalling_r377624/Notch_maml_abs_test_1/outputs/simulated_vectors_file.txt' /project/V0001-1/modcell/users/a.kovachev/data/data/init_vector_Notch_maml_abs.txt 46 90"
		print "NOTE: the list of simulated should be wrap around '' so can be recognized as string an parsed as list of files."
		print "Please check the notation of the input files and the parameters size of the model for the offset."