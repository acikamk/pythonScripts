import numpy as np
import os
import sys
import collections 
from collections import OrderedDict as Odict
import argparse
from argparse import RawTextHelpFormatter
import itertools
import pandas
import logging
from operator import itemgetter
import tempfile
import matplotlib.pyplot as plt

# check if the path need to be forwarded as args
sys.path.append("/home/H0001-1/a.kovachev/pybios/AlacrisPython")
sys.path.append("/home/H0001-1/a.kovachev/simtools/src")

from AlacrisPython.Utils.AlacrisConfig import AlacrisConfig 
from AlacrisPython.Utils.AlacrisConfig import CONFIG_DATA_INPUT
from AlacrisPython.Utils.AlacrisConfig import CONFIG_BIOMARKER_IDENT
from AlacrisPython.Utils.StringUtils import wrap_text
import runopt
"""
Script that runs biomarker identification using parameters swap
between responder and non-responder samples listed in the config file.
NOTE: Virtual enviroment should be activated before running the script
    - source venv/bin/activate
Input:
    -Config file with input parameters defined below.

Output:
    -Plots of the T/C values for forward simulation for all responder
    and non-responder pairs.

"""
# initialize Alacis Python config
CONFIG = AlacrisConfig([ CONFIG_DATA_INPUT,
                         CONFIG_BIOMARKER_IDENT,
                        ] )

class ReadConfig(argparse.Action):
    
    def __call__(self, parser, namespace, value, option_string=None):
        # update default config by settings given in config file       
        CONFIG.load( value )
        # set attribute self.dest
        setattr(namespace, self.dest, value)

parser = argparse.ArgumentParser(description = 
            """
            The script should be run calling:
            runBiomarkId -c /path/to/config/file

            # Necessary config data:
            [DATA_INPUT]
            PYBIOS_MODEL = <PATH_TO_PYBIOS_MODEL>
            CONTROL_TREATMENT_FILE = <PATH_TO_COST_FUNC_FILE>
            EXPERIMENT_TABLE = <PATH_TO_EXP_TABLE_FILE>
            CCOMPILED_MODEL = <PATH_TO_CCOMPILED_MODEL_FILE>
            STUDY_PATH = <PATH_TO_OPTIMIZATION_RESULTS>
            [BIOMARKER_IDENT]
            RESPONDERS = <LIST_OF_RESP_NAMES>
            NON_RESPONDERS = <LIST_OF_NON_RESP_NAMES>
            """,
                                formatter_class = RawTextHelpFormatter,
                                usage = "%(prog)s [-h --help] [OPTIONS]")

parser.add_argument('-c', '--config',
                    metavar="CONFIG_FILE",
                    dest = 'config_file',
                    action = ReadConfig,
                    help = wrap_text("path to a config file." +
                    "You get the config file structure by calling 'testSSSpeed -h'"))

def main():
    # Take the config file from command line
    args = parser.parse_args()
      
    # Get the input parameters  
    model_path = existingDir(CONFIG.get('PYBIOS_MODEL', section='DATA_INPUT'))
    cost_func_file = existingFilePath(CONFIG.get('CONTROL_TREATMENT_FILE', section='DATA_INPUT'))
    exp_table_file = existingFilePath(CONFIG.get('EXPERIMENT_TABLE', section='DATA_INPUT'))
    c_model = existingFilePath(CONFIG.get('CCOMPILED_MODEL', section='DATA_INPUT'))
    responders = CONFIG.get('RESPONDERS', section='BIOMARKER_IDENT')
    responders = CONFIG._cast_by_type(responders, list)
    non_responders = CONFIG.get('NON_RESPONDERS', section='BIOMARKER_IDENT')
    non_responders = CONFIG._cast_by_type(non_responders, list)
    simtools_results_path = existingDir(CONFIG.get('STUDY_PATH', section='DATA_INPUT'))

    progress_trace_file = existingFilePath(simtools_results_path \
                 + "/outputs/progress_trace.txt")
    tempfile.tempdir = simtools_results_path + "/outputs/"

    # parse the responders and non-responders data
    exp_table = pandas.read_table(exp_table_file, sep='\t', index_col=0)

    # create config dictionary as input for simulation with simtools
    config = initConfig(model_path, c_model, progress_trace_file)

    # turn of logging to console
    logger = logging.getLogger()
    logger.propagate = False
    fitness = Odict()

    # iterate over the responders and non-respodenrs
    for responder in responders:
        for non_responder in non_responders:
            fitness[responder] = Odict()
            fitness[non_responder] = Odict()
            tmp_exp_table = exp_table.copy()
            if len(exp_table.loc[:,  responder])!=len(exp_table.loc[:, non_responder]): 
                print "Samples with different number of parameters!"
                return
            print "Starting simulations for {} vs {}...".format(responder, non_responder)
            # iterate over each of the parameters
            for i in range(0, len(exp_table.loc[:,  responder])):
                if not "GeneSpecificScaling" in tmp_exp_table.index[i] and \
                  not (exp_table.iloc[i][responder] == exp_table.iloc[i][non_responder]):
                    tmp = exp_table.iloc[i][responder]
                    # swap parameters
                    tmp_exp_table.iloc[i][responder] = tmp_exp_table.iloc[i][non_responder]
                    tmp_exp_table.iloc[i][non_responder] = tmp
                    # run simulations
                    fitness[responder][tmp_exp_table.index[i]] = runSampleSimulation(
                                                                    config, 
                                                                    tmp_exp_table, 
                                                                    cost_func_file, 
                                                                    responder)

                    fitness[non_responder][tmp_exp_table.index[i]] = runSampleSimulation(
                                                                        config,
                                                                        tmp_exp_table,
                                                                        cost_func_file, 
                                                                        non_responder)
            # order T/C results for the run        
            if len(fitness[responder].keys())>2:                        
                sorted_responder = Odict(sorted(fitness[responder].items(), 
                                        key=itemgetter(1)))
                order = np.argsort(np.array(sorted_responder.values()))
                sorted_nonresponder = Odict(sorted(fitness[non_responder].items(), 
                                        key=itemgetter(1)))
                sorted_fitness_nrorder = [sorted_nonresponder.values()[i] for i in order]
                # import pdb; pdb.set_trace() 
                fig = plt.figure()
                ax = fig.add_subplot(111)
        
                ax.set_xlabel('Swapped parameters')
                ax.set_ylabel('Modeled T/C value')

                ax.plot(sorted_responder.values(), color = "b", 
                                label = "Responder_"+responder)
                ax.plot(sorted_fitness_nrorder, color = "r",
                                label = "Non-responder_"+non_responder)
                plt.grid()
                plt.legend()
                # plt.show()
                plt.tight_layout()
                txt = 'Responders:\n'+ '\n'.join(sorted_responder.keys()[:5])
                plt.figtext(1.00, 0.80, txt, color='black', size='x-small')
                txt = 'Non-responders:\n' + '\n'.join(sorted_nonresponder.keys()[:5])
                plt.figtext(1.00, 0.30, txt, color='black', size='x-small')
                plt.draw()
                plt.savefig(simtools_results_path + "/outputs/" + \
                            responder + "_" + non_responder + ".png",
                            dpi=200, bbox_inches='tight')
                plt.close(fig)
                
def existingFilePath(file_path):
    file_path = os.path.abspath(os.path.expanduser(file_path))
    assert os.path.isfile(file_path), "File %s was not found" % file_path
    return file_path

def existingDir(dir_path):
    dir_name = os.path.abspath(os.path.expanduser(dir_path))
    assert os.path.isdir(dir_name), "Cannot find dir: %s" % dir_path
    return dir_name

def ModifyCostFile(cost_func_file, sample):
    tmp_cost_file = tempfile.NamedTemporaryFile(delete=False)
    with open(cost_func_file) as f:
        for line in f:
            if line.startswith("ID") or sample in line:
                tmp_cost_file.write(line)
    tmp_cost_file.close()
    return tmp_cost_file

def initConfig(model_path, c_model, progress_trace_file):
    config = collections.defaultdict(dict)
    existingDir(model_path)
    ids_file = model_path + "/identifiers.py"
    model_id = model_path[model_path[:-1].rindex("/")+1:-1] 
    config["simulate"] = True
    config["optimize"] = False
    config["model_id"] = model_id
    config["auto"]["ids_file"] = ids_file
    config["vectors_file"] = progress_trace_file 
    config["log_scale_sim"] =  True
    config["get_gradient"] = False
    config["get_states"] = True
    config["take_only_last"] = True
    config["use_cluster"] = False
    config["local_cores"] = 4
    config["port"] = 47075
    config["auto"]["modeltools_so"] = c_model
    config["socket_timeout"] = 36000
    return config

def runSampleSimulation(config, tmp_exp_table, cost_func_file, sample):
    # set up the temporary config files
    config["auto"]["exp_table"] = tmp_exp_table   
    tmp_cost_file = ModifyCostFile(cost_func_file, sample)
    config["auto"]["cost_func"] = tmp_cost_file.name

    tmp_sim_file = tempfile.NamedTemporaryFile(delete=False)
    config["auto"]["simulated_vectors_file"] = tmp_sim_file.name
    tmp_sim_file.close()

    # run siumlations uisng simtools
    runopt.evaluate(config)

    # take t/c constraint from simulated file
    with open(config["auto"]["simulated_vectors_file"]) as f:
        for line in f:
            if line.startswith("FITNESS:"):
                fitness = eval(line.split('FITNESS:')[-1])
            if line.startswith("CONSTRAINTS:"):  
                constr_dict_raw = eval(line.split('CONSTRAINTS:')[-1])
                cell_division = float([dd.values()[0] for k, dd in \
                         constr_dict_raw.items() if len(k) == 2][0])
        # import pdb; pdb.set_trace()        
        os.remove(config["auto"]["cost_func"])
    os.remove(config["auto"]["simulated_vectors_file"])
    return cell_division



if __name__ == "__main__":
   main()


# generate the combinations of possible swaps 
# regarding the number of swaps 
#for i in range(0, comb_length):
#    index_combi = list(itertools.combinations([range(0,len(responder))],comb_length))
            