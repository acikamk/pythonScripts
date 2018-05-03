import numpy
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
import argparse
from argparse import RawTextHelpFormatter
import itertools
import pprint
'''
Script that tests the convergence of the solver with and without 
using the steady state of the control as intiial states of the drug treatments
This script uses text file as config file with format provided and parsed using
AlacrisPython.
Example input file is tsss.cfg: 
PYBIOS_MODEL = path to the pybios files of the model
CONTROL_TREATMENT_FILE = path to a cost function file, ex. cost_func.txt
EXPERIMENT_TABLE = path to the experiment table file ex. exp_table.csv
CCOMPILED_MODEL = path to the c_compiled version of the model, 
ie.e the folder where modeltools.so is 


'''

sys.path.append("/home/H0001-1/a.kovachev/simtools/src")
from modelstuff import PyBiosModel

sys.path.append("/home/H0001-1/a.kovachev/pybios/AlacrisPython")

from AlacrisPython.Utils.AlacrisConfig import AlacrisConfig, CONFIG_GENERAL, CONFIG_DATA_INPUT, CONFIG_SIMULATION_ANALYSIS
from AlacrisPython.Utils.StringUtils import wrap_text

# initialize Alacis Python config
CONFIG = AlacrisConfig( [ CONFIG_GENERAL,
                        CONFIG_DATA_INPUT,
                        CONFIG_SIMULATION_ANALYSIS,
                        ] )

class ReadConfig(argparse.Action):
    
    def __call__(self, parser, namespace, value, option_string=None):
        # update default config by settings given in config file       
        CONFIG.load( value )

        # set attribute self.dest
        setattr(namespace, self.dest, value)


parser = argparse.ArgumentParser( description = 
            """
            The script should be run calling:
            testSSSpeed -c /path/to/config/file

            # Necessary config data:

            [DATA_INPUT]
            PYBIOS_MODEL = <PATH_TO_PYBIOS_MODEL>
            CONTROL_TREATMENT_FILE = <PATH_TO_COST_FUNC_FILE>
            EXPERIMENT_TABLE = <PATH_TO_EXP_TABLE_FILE>
            CCOMPILED_MODEL = <PATH_TO_CCOMPILED_MODEL>
            """,
                                formatter_class = RawTextHelpFormatter,
                                usage = "%(prog)s [-h --help] [OPTIONS]" )

parser.add_argument( '-c', '--config',
                    metavar="CONFIG_FILE",
                    dest = 'config_file',
                    action = ReadConfig,
                    help = wrap_text("path to a config file. \
                    You get the config file structure by calling 'testSSSpeed -h'") )

def main():

    # Take the config file from command line
    args = parser.parse_args()
      
    model_path = CONFIG.get( 'PYBIOS_MODEL', section='DATA_INPUT') 
    if model_path[-1] != "/":
        model_path += "/"
    ids_file = model_path + "identifiers.py"
    model_id = model_path[model_path[:-1].rindex("/")+1:-1] 
    # Get options     
    cost_func = CONFIG.get('CONTROL_TREATMENT_FILE', section='DATA_INPUT')
    exp_table = CONFIG.get('EXPERIMENT_TABLE', section='DATA_INPUT')
    c_model = CONFIG.get('CCOMPILED_MODEL', section='DATA_INPUT')

    # import c-compiled model
    sys.path.append(c_model)
    import modeltools

    mp = modeltools.ModelProcessor(maxStepNum=5000)
    # prepare time and range for simulation
    sim_time = 1e9
    step = sim_time/100
    cctime = numpy.arange(0, float(sim_time)+1 , step, dtype=numpy.float)  

    # initialize the model and random inital parameter vector
    model = PyBiosModel(model_id, ids_file, exp_table, cost_func)
    bounds = (-1,1)
    varKParVect =  numpy.power(10, bounds[0] + \
            numpy.random.rand(len(model.variable_indexK_arr) + \
                model.dimKCD)*(bounds[1] - bounds[0]) )
    # varKParVect = bounds[0] + \
        # numpy.random.rand(len(model.variable_indexK_arr)+ \
        # model.dimKCD)*(bounds[1] - bounds[0]) 
    process_times=Odict()
    real_times=Odict()  
    
    fixed_Ki = numpy.setdiff1d(xrange(model.dimK), model.variable_indexK_arr)
    K = numpy.zeros(model.dimK)
    if len(model.variable_indexK_arr) > 0:
        K[model.variable_indexK_arr] = varKParVect
    # sample_states = {}
    # import pdb; pdb.set_trace()
    cost_func_sample_ids = set(itertools.chain(*model.conditionDict.keys()))

    # for each sample in the cost function run formward simulation
    # print the run time
    for sampleId in cost_func_sample_ids:
        sample = model.samples[sampleId]
        if len(fixed_Ki)  > 0:
            K[fixed_Ki] = sample.K[fixed_Ki]
            start_time = time.time()
            start_ptime = time.clock()   
            res = mp.simulate(sample.S0, sim_time, sample.F, K)
            if res.success:
                process_times[sampleId] = time.clock() - start_ptime
                # real_times[sampleId] = time.time() - start_time
                print "({}) Process time normal: \t\t\t {}".\
                    format(sampleId, process_times[sampleId])
                # print "({}) Real time normal: {}".format(sampleId, real_times[sampleId])
            else:
                raise RuntimeError("Simulation time course failed with flag %s for sample %s" \
                         %( res.exitcode, sampleId))

    # import pdb; pdb.set_trace()
    print sum(process_times.values())
    # pprint.pprint(process_times)

    process_times=Odict()
    ss_dict = Odict()
    cost_func_sample_ids_sorted = sorted(list(set(itertools.chain(*model.conditionDict.keys()))))
    # import pdb; pdb.set_trace()
    # run the same simulations now with using the previous steady state as initial
    for sampleId in cost_func_sample_ids_sorted:
        sample = model.samples[sampleId]
        if len(fixed_Ki) > 0:
            K[fixed_Ki] = sample.K[fixed_Ki]
            # import pdb; pdb.set_trace() 
            if any(string in sampleId for string in ss_dict.keys()):  
                S0 = [ss_dict[i] for i in ss_dict.keys() \
                             if sampleId.startswith(i)][0]
            else:
                S0 = sample.S0
            start_time = time.time()
            start_ptime = time.clock()   
            res = mp.simulate(S0, sim_time, sample.F, K)
            # res = mp.simulate_time_course(sample.S0, cctime, sample.F, K)
            if res.success:
                process_times[sampleId] = time.clock() - start_ptime
                # real_times[sampleId] = time.time() - start_time
                print "({}) Process time ss: \t\t\t {}".\
                        format(sampleId, process_times[sampleId])
                # print "({}) Real time normal: {}".format(sampleId, real_times[sampleId])
                if "CONTROL" in sampleId or not "_" in sampleId:
                    ss_dict[sampleId] = res.finalState
                    # ss_dict[sampleId] = res.timeCourse[-1]
            else:
                raise RuntimeError("Simulation time course failed with flag %s for sample %s" \
                         %( res.exitcode, sampleId))  

    # print process_times
    print sum(process_times.values())

if __name__ == "__main__":  
    main()
