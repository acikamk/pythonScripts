import numpy as np
import re
import os
from os import system
import sys
import math
import random
import time
import collections 
from collections import OrderedDict as Odict
"""
Script that runs simple forward simulation with the sap and simtools c solvers. 
No experiment table or cost_function used. Model initialised with ones (random*)
data.

Input:
	ccompiled = path to the simtools c_compiled pybios model
	sap = path to the sap c_compiled pybios model
	sim_time =  time for the simulation
	x =  2 # number of species to check
	rand = true/false wether to use rand paramters or ones
	size_S = number of species in the model
	size_F = number of conditions in the model
	size_K = number of parameters in the model

"""

ccompiled = "/home/H0001-1/a.kovachev/pybios/scripts/models/CancerSignalling_Merged_2017v1_Mutations_Horst_Drugs_r388408/c_solver/"
#"/home/H0001-1/a.kovachev/pybios/scripts/models/CancerSignalling_Merged_2017v1_r378547"
#
#"/home/H0001-1/a.kovachev/simtools/data/models/Model_Optimization_Mai2017/"

sap =  "/home/H0001-1/a.kovachev/pybios/scripts/models/CancerSignalling_Merged_2017v1_Mutations_Horst_Drugs_r388408/sap_solver/"
"/home/H0001-1/a.kovachev/pybios/scripts/models/CancerSignalling_Merged_2017v1_r378547_sap"
#
sim_time = 1e9
x =  2
rand = True

# /home/H0001-1/a.kovachev/simtools/data/models/Model_Optimization_Mai2017/
# size_S = 1228
# size_F = 168
# size_K = 4686

size_S = 9217
size_F = 1235
size_K = 29797



sys.path.append(ccompiled)
sys.path.append(sap)
	
# import sap and cython libs
import modeltools 
import model_solver

mp = modeltools.ModelProcessor(abstol=1.0E-9, reltol=1.0E-5, maxStepNum=5000)
            # bdfMaxOrder=5, stabLimDet=1):
sp = model_solver.SolverWrapper()
            # 'initial_step_size' : 1.0E-4,
            # 'abs_tolerance'     : 1.0E-9,
            # 'rel_tolerance'     : 1.0E-5,
            # 'min_step_size'     : 1E-8,
            # 'max_step_size'     : 1E8, !!!!!
 # self.solver = SapOdeSolver.Solver(model_name=self.model_name)
sp.set_ode_parameters(1.0E-4, 1.0E-9, 1.49e-8, 1E-8, 10000)
# set_ode_parameters(h, abs_tol, rel_tol, min_step_size, max_step_size):
# h initial step size

# set model parameters
if rand:
	S = np.random.rand(size_S)
	F = np.random.rand(size_F)
	K = np.random.rand(size_K)
else:	
	S = np.zeros(size_S)
	F = np.ones(size_F)
	K = np.ones(size_K)

step = sim_time/100
cctime = np.arange(0, float(sim_time)+1 , step, dtype=np.float)     
#cctime = np.ndarray((100,), buffer=np.array([0,1.01010101010101,2.02020202020202,3.03030303030303,\
# 4.04040404040404,5.05050505050505,6.06060606060606,7.07070707070707,8.08080808080808,9.09090909090909,\
# 10.1010101010101,11.1111111111111,12.1212121212121,13.1313131313131,14.1414141414141,15.1515151515152,\
# 16.1616161616162,17.1717171717172,18.1818181818182,19.1919191919192,20.2020202020202,21.2121212121212,\
# 22.2222222222222,23.2323232323232,24.2424242424242,25.2525252525253,26.2626262626263,27.2727272727273,\
# 28.2828282828283,29.2929292929293,30.3030303030303,31.3131313131313,32.3232323232323,33.3333333333333,\
# 34.3434343434343,35.3535353535354,36.3636363636364,37.3737373737374,38.3838383838384,39.3939393939394,\
# 40.4040404040404,41.4141414141414,42.4242424242424,43.4343434343434,44.4444444444444,45.4545454545455,\
# 46.4646464646465,47.4747474747475,48.4848484848485,49.4949494949495,50.5050505050505,51.5151515151515,\
# 52.5252525252525,53.5353535353535,54.5454545454546,55.5555555555556,56.5656565656566,57.5757575757576,\
# 58.5858585858586,59.5959595959596,60.6060606060606,61.6161616161616,62.6262626262626,63.6363636363636,\
# 64.6464646464647,65.6565656565657,66.6666666666667,67.6767676767677,68.6868686868687,69.6969696969697,\
# 70.7070707070707,71.7171717171717,72.7272727272727,73.7373737373737,74.7474747474748,75.7575757575758,\
# 76.7676767676768,77.7777777777778,78.7878787878788,79.7979797979798,80.8080808080808,81.8181818181818,\
# 82.8282828282828,83.8383838383838,84.8484848484848,85.8585858585859,86.8686868686869,87.8787878787879,\
# 88.8888888888889,89.8989898989899,90.9090909090909,91.9191919191919,92.9292929292929,93.9393939393939,\
# 94.9494949494950,95.9595959595960,96.9696969696970,97.9797979797980,98.9898989898990,100]))

# # set and run sap simulations 
start_time = time.time()
start_ptime = time.clock()  
sp.set_model_parameters(F, K)
sp.set_ode_parameters(1.0E-4, 1.0E-9, 1.0e-5, 1E-8, 1E8)
sap_res = sp.simulate_with_time_series(S, cctime)
process_time = time.clock() - start_ptime
real_time = time.time() - start_time
print "Process time SAP: \t\t\t {}".format(process_time)
print "Real time SAP: \t\t\t {}".format(real_time)
sap_last = sap_res[-1]

start_time = time.time()
start_ptime = time.clock()  
sap_res1 = sp.simulate(S, sim_time)
process_time = time.clock() - start_ptime
real_time = time.time() - start_time
print "Process time SAP_sim: \t\t\t {}".format(process_time)
print "Real time SAP_sim: \t\t\t {}".format(real_time)
sap_last = sap_res1

# run ccompiled simulations
start_time = time.time()
start_ptime = time.clock()  
cc_results = mp.simulate_time_course(S, cctime, F, K)
process_time = time.clock() - start_ptime
real_time = time.time() - start_time
print "Process time C: \t\t\t {}".format(process_time)
print "Real time C: \t\t\t {}".format(real_time)

cc_res = cc_results.timeCourse
cc_last = cc_res[-1]

start_time = time.time()
start_ptime = time.clock()  
cc_results=mp.simulate(S, sim_time, F, K)
process_time = time.clock() - start_ptime
real_time = time.time() - start_time
print "Process time C_sim: \t\t\t {}".format(process_time)
print "Real time C_sim: \t\t\t {}".format(real_time)

#cc_last = cc_results.timeCourse[-1]
cc_last = cc_results.finalState

# calculate abs and rel tolerances
# total difference
abs_tot = abs(sum(sum(sap_res))-sum(sum(cc_res)))
rel_tot = abs(1-sum(sum(sap_res))/sum(sum(cc_res)))
print "Total concentrations error:\t\t abs_tol: {} \t rel_tol:{}".\
		format(abs_tot, rel_tot)
# last vector difference
abs_last = abs(sum(sap_last)-sum(cc_last))
rel_last = abs(1-sum(sap_last)/sum(cc_last))
print "Last t concentrations error:\t\t abs_tol: {} \t rel_tol:{}".\
		format(abs_last, rel_last)
# check x random species
ind_check = random.sample(range(1, len(S)), x)
for i in ind_check:
	sap_time= np.array(zip(*sap_res)[i])
	cc_time = np.array(zip(*cc_res)[i])
	abs_ind = sum(abs(sap_time-cc_time))
	rel_ind = sum(abs(1-sap_time/cc_time))
	print "Species id: {} concentrations error:\t\t abs_tol: {} \t rel_tol:{}".\
			format(i, abs_ind, rel_ind)
