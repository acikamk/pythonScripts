import numpy as np
import re
import os
from os import system
import sys
import time
import math
import matplotlib.pyplot as plt
import collections 
from collections import OrderedDict as Odict
import random
import cProfile
import itertools
sys.path.insert(0, "/home/H0001-1/a.kovachev/pybios/ThoughtExperiments/thoughtExperiments/scripts/testSimulations/profiler/")
from prof import profiler_start, profiler_stop

"""
Script that runs time course analysis of simple forward simulation for several models
Additionaly can run profiler on the C code simulation function

Input argument:
True - if given thatn simulations will be started with steady state starts

models = list of pybios models
ccompiled = list of ther ccompiled files
sim_time_vector  =  list of simulation times

"""
def check_steady_state(res_ccompile, cctime, not_in_steady_state):
	t_diff = cctime[-1] - cctime[-2]
	for i in range(len(res_ccompile.timeCourse[-1])):
		x_diff = abs(res_ccompile.timeCourse[-2][i]-res_ccompile.timeCourse[-1][i])                               
		diff_quot = x_diff / t_diff                               
		if diff_quot >= 0.061:                  
			not_in_steady_state = True
			# print diff_quot
	return not_in_steady_state

# model_file_3 = "/project/V0001-1/modcell/modcell_data/project_manager/32_PP20/2016.10.17_10.17.37_435_PP20/Model/"           # set path to model_file
# ccompiled_model_file_3 = "/project/V0001-1/modcell/modcell_data/models/c_pybios_models/CancerSignalling_Mutations_Drugs_v2016_Birmingham_r371850/" # set path to ccompiled model file 

# ccompiled_model_file_2 = "/home/H0001-1/a.kovachev/pybios/ThoughtExperiments/thoughtExperiments/scripts/CCompiled/WNT_Signalling_r371946/"
# model_file_2 = "/home/H0001-1/a.kovachev/pybios/ThoughtExperiments/thoughtExperiments/scripts/WNT_Signalling_r371946/"

models = ["/home/H0001-1/a.kovachev/simtools/data/models/Notch_Signalling_r377624/"]
ccompiled = ["/home/H0001-1/a.kovachev/simtools/data/models/Notch_Signalling_r377624/"]
sim_time_vector  =  [1e2, 1e3, 1e4, 1e5]


def run_simulation_test(model, sim_time_vector):

	from model import PybiosEqs
	import modeltools 

	process_times=Odict()
	real_times=Odict()	
	try:
		if (sys.argv[1] == "True" or sys.argv[1] == "true"):
			print "Running simulations with steady state as a start"
			flag_next = False
		else:	
			print "Running normal simulations"
			flag_next = True
	except IndexError:	
		print "Running normal simulations"		
		flag_next = True

	flag_first = True
	for sim_time in sim_time_vector:
		try:
		    model_raw = PybiosEqs(0.000000, sim_time, [model], modParPath = [model], filePrefix = model)
		except TypeError:# for older models                
		    model_raw = PybiosEqs(0.000000, sim_time, model)

		mp = modeltools.ModelProcessor()
		S = model_raw.S   
		F = model_raw.F   
		K = model_raw.K  

		step = sim_time/100
		cctime = np.arange(0, float(sim_time)+1 , step, dtype=np.float)     
		if flag_first: 
			start_time = time.time()
			start_ptime = time.clock()
			not_in_steady_state = False
			try:
				profiler_start("simulate_normal_"+model[model[:model.rfind("/")].rfind("/")+1:-1]+".log")
				res_ccompile = mp.simulate_time_course(S, cctime, F, K)
				S_time=res_ccompile.timeCourse[-1] 
				profiler_stop()
			except MemoryError:
				print "Memory Error"

			process_times[sim_time] = time.clock() - start_ptime
			real_times[sim_time] = time.time() - start_time
			print "T({}) Process time normal: {}".format(sim_time, process_times[sim_time])
			print "T({}) Real time normal: {}".format(sim_time, real_times[sim_time])

			flag_first = flag_next # change this if running with or without steady state (true normal, false ss)
			S_time=res_ccompile.timeCourse[-1] 
			not_in_steady_state = check_steady_state(res_ccompile, cctime, not_in_steady_state)
			if not_in_steady_state:
				print "SPECIES NOT IN STEADY STATE" 				
		else:
			start_ptime = time.clock()
			start_time = time.time()
			try:
				profiler_start("simulate_short_"+model[model[:model.rfind("/")].rfind("/")+1:-1]+".log")
				res_ccompile = mp.simulate_time_course(S_time, cctime, F, K)
				profiler_stop()
			except MemoryError:
				print "Memory Error"
			process_times[sim_time] = time.clock() - start_ptime
			real_times[sim_time] = time.time() - start_time
			print "T({}) Process time short: {}".format(sim_time, process_times[sim_time])
			print "T({}) Real time short: {}".format(sim_time, real_times[sim_time])
			
			S_time=res_ccompile.timeCourse[-1]

			not_in_steady_state = check_steady_state(res_ccompile, cctime, not_in_steady_state)
			if not_in_steady_state:
				print "SPECIES NOT IN STEADY STATE RUNNING AGAIN WITH NORMAL SIMULATIONS"
				res_ccompile = mp.simulate_time_course(S, cctime, F, K)
	
	return res_ccompile, process_times, real_times, len(S)

# time vise dict
results = {}
res_S = []
res_list = []
color_map = ['r', 'b', 'y', 'g']
width = 0.30
w = 0.0

# initialise the results dict time vise
for times in sim_time_vector:
	results[times] = []

# run simulation test for each of the models and each time
for (model, ccompiled_model, c) in zip(models, ccompiled, color_map):	

	sys.path.append(ccompiled_model)
	sys.path.append(model)

	print "Analysing model: {}".format(model)

	[res_ccompile, process_times, real_times, S_size] = run_simulation_test(model, sim_time_vector)
 
 	# remove the model from the path since all have the same name: modeltools
	sys.path.remove(ccompiled_model)		
	sys.path.remove(model)
	del sys.modules['modeltools']

	# add the model sizle to the list
	res_S = res_S + [S_size]
	for times in sim_time_vector:
		 results[times].append(real_times[times])
	
	res_list.append(real_times.values())

# normalize and standardize data	
# x_min = min(res_list)
# x_max = max(res_list)
# mean = np.mean(res_list)
# stdev = np.std(res_list)
# standardize = {}
# normalize = {}

# plot results
plt.figure(1)
process_line, = plt.semilogx(process_times.keys(), process_times.values(), 'r', label = 'Process times')
real_line, = plt.semilogx(real_times.keys(), real_times.values(), label = 'Real times')
plt.xlabel("simulation time")
plt.ylabel("time (s)")
plt.legend()
plt.draw()
figure = plt.gcf()
figure.savefig("Times_figure.pdf", dpi = figure.dpi)

plt.figure(2)
fig, ax = plt.subplots()
for (times,c) in zip(sim_time_vector, color_map):	
	# normalize = np.array([(x-x_min)/(x_max-x_min) for x in results[times]] )
	# standardize = np.array([(x-mean)/stdev for x in results[times]] )
	ax.loglog(res_S, results[times], marker='o', color=c, label=str(times))
	

lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
# ax.loglog([60,10000], [0.001,100], 'k-')
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)

plt.legend()
plt.draw()
figure = plt.gcf()
figure.savefig("Times_compare.pdf", dpi = figure.dpi)
