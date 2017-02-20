'''General job parameters.
'''
# A model folder with this name should present in the "data/models" folder
# of the simtools package
model_id = "HORST_July2016_Drugs_r377479"
# A folder with this name will be created in the "results/<model_id>"
# folder. All results produced by the job will be stored insite this folder.
job_id = "OncoTrack_sim_23_12_16_test"
# Log level for the console logger. File logger will be initialized with
# "DEBUG" level anyway, so all log messages will be availabel in the
# "results/<model_id>/<job_id>/joblog.txt" file.
loglevel = "INFO"

results_dir = "/project/V0001-1/modcell/users/a.kovachev/"
''' Tasks.
Boolean flags indicating which tasks should be performed.
'''
# If True, an experimental table is assumed to be in the format produced by
# the project manager and will be converted to the format acceptable by the
# simtools package. Else the table is supposed to be already in the simtools
# compatible format and will be taken directly without any processing.
process_exp_table = False   
# If False, a cost function will be taken directly as it is. Else it will be
# processed taking into account shuffle_conditions and n_training_conditions
# parameters.
process_cost_function = False
# If True a c code for the integration will be generated and compiled de
# novo. Else a precompiled modeltools.so file will be taken from the
# location specified by the modeltools_so parameter (specified below).
generate_code = False
# If True, optimization will be started with parameters from the
# "Optimization" section below.
optimize = False
# If True, simulation will be started with parameters from the "Simulation"
# section below.
simulate = True


''' Experimental table.
'''
# A path to the file with experimental data for model initialization (RPKM
# values, mutations etc.). This file should be either in the format produced
# by the project manager or in the format required by the simtools package.
# Both formats are described below.
exp_table = ("/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/OncoTrack_sim_23_12_16/inputs/exp_table.csv")
if process_exp_table:
    # If True, gene scaling parameters for all samples will be fixed to the
    # value specified in the identifiers.py file of the model.
    fixed_gene_scaling = True
    # If True, drug target interaction parameters will be fixed to the value
    # specified in the identifiers.py file of the model.
    fixed_drug_target_interaction = True
    # If True, drug translocation parameters will be fixed to the value
    # specified in the identifiers.py file of the model.
    fixed_drug_translocation = True


''' Cost function.
'''
# A path to the cost function file. Required format for the cost function
# is described below.
cost_func = ("/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/OncoTrack_sim_23_12_16/inputs/tmp_cost_func.txt")
if process_cost_function:
    # If True, conditions in the cost function will be shuffled. Else the
    # order of conditions will remain as in the original file.
    shuffle_conditions = True
    # first n_training_conditions conditions will be taken for training
    # if None all conditions will be used for training
    n_training_conditions = None


''' Compiled model.
'''
if not generate_code:
    # A path to a dynamically loaded library with precompiled code for the
    # integration of the model.
    modeltools_so = ("/project/V0001-1/modcell/users/a.kovachev/HORST_July2016_Drugs_r377479/OncoTrack_sim_23_12_16/inputs/modeltools.so")

init_with_prev_final =  True

''' Computation.
These parameters are use only if "optimize" or "simulate" task flag is True.
'''
# number of cores to use on the local computer
local_cores = 4
# socket timeout in seconds
socket_timeout = 60
# if True cluster nodes will be used else only local cores will be active
use_cluster = False
# if use_cluster:
#     port = 47051
#     cluster_nodes = { "compute01" : 12 }

''' Simulation.
These parameters are use only if "simulate" task flag is True.
'''
# A file with vectors to simulate in an appropriate format
vectors_file = ("VectorsFile")
# take only the last vector in the vectors_file
take_only_last = True
# evaluate gradient
get_gradient = False
# write 'steady' states
get_states = True
# use all except kCD parameters on the log10 scale
log_scale_sim = True
