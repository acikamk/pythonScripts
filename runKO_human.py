import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse, ast
from argparse import RawTextHelpFormatter
from os.path import isfile, join
import re
import pandas as pd
import numpy as np
import itertools
import pdb
sys.path.append('/home/H0001-1/a.kovachev/simtools/src')
from modelstuff import PyBiosModel

global exp_table
global cost_file 
global identifiers_file 
global cmodel_file
global all_ids 

def existingFilePath(file_path):
	file_path = os.path.abspath(os.path.expanduser(file_path))
	assert os.path.isfile(file_path), \
		'File %s was not found' % file_path
	return file_path
def existingDir(dir_path):
    dir_name = os.path.abspath(os.path.expanduser(dir_path))
    assert os.path.isdir(dir_name), \
    	'Cannot find dir: %s' % dir_path
    return dir_name

def newFileInExistingDir(file_path):
   dir_name = existingDir(os.path.dirname(file_path))
   file_name = os.path.basename(file_path)
   assert file_name, \
   		'Empty file name for dir %s' % dir_name
   return os.path.join(dir_name, file_name)

def read_data(model):
	print 'INFO Reading input data.'
	global exp_table
	global cost_file 
	global identifiers_file 
	global cmodel_file
	global all_ids

	exp_table = model + '/exp_table_mini.csv'
	cost_file = model + '/cost_func_empty.txt'
	identifiers_file = model + '/identifiers.py'
	outname = '/'+ exp_table[exp_table.rindex('_')+1:exp_table.rindex('.')] + '_'

	exp_df = pd.read_table(exp_table, sep = '\t', 
		low_memory = False, index_col = 0)
	exec('global diff_var_ids, fixed_var_ids, par_ids\n'
		 + open(identifiers_file).read())

	all_ids = diff_var_ids.copy()
	all_ids.update(fixed_var_ids)

	with open(cost_file) as f:
		lines = f.readlines()
		cost_data = [(l.split('\t')[0].split(':')[1].strip(),
			l.split('\t')[1].strip()) \
		 	for l in lines if l.startswith('ID')] 

	return exp_df, cost_data, outname

def parseOptiFile(opti_file):
	kcds_dict = Odict()
	opti_vector = Odict()
	kcds_dict = Odict()
	with open(opti_file) as file:
		for l in file:
			if 'VARIABLE PAR IDS:' in l:
				ids_vector = eval(l[l.index('['):])
			elif 'FITNESS:' in l:
				fitness = eval(l[l.index(':')+1:])
			elif 'VECTOR:' in l:
				opti_vector = eval(l[l.index(':')+1:])		
	for var in ids_vector:
		if var.startswith('k'):
			kcds_dict[var] = opti_vector[ids_vector.index(var)]			
	return ids_vector, opti_vector, kcds_dict

def modifyExpTable(exp_df,  args, outname):

	print 'INFO Generating experiment table.'
	exp_file_name = args.out_dir + '/' + \
		 outname +'_exp_table.csv'
	if not args.sp:
		not_in_exp = [i for i in args.sp 
		 if i not in exp_df.index.tolist()]

		assert not not_in_exp, \
			('Error in input data. Species: ' + not_in_exp + 
			'are missing in the experiment table definition.')
	splist = args.sp
	plist = []
	prtvals = args.vals
	if args.genes:
		splist += [i for (i,j) in fixed_var_ids.iteritems() 
			if fixed_var_ids[i][5] == 'gene' ]
	if args.rpkms:
		plist += [i for (i,j) in par_ids.iteritems() 
			if 'RPKM' in i]

	exp_df.columns = [i.replace('_','-') for i in exp_df.columns.tolist()]
	exp_df_copy = exp_df.copy()
	
	for s in splist:
		for column in exp_df:
			ofs = 1
			for p in prtvals:
				col_name = column.replace('_','-') + '_' + \
				fixed_var_ids[s][1].replace(' ','-')
				col_copy = exp_df[column].copy()
				if args.perc:
					col_copy[s] += col_copy[s]*p/100
					col_name += '_percKO'+ str(p)
				else:
					if col_copy[s] == 1.0 and args.genes:
						p = 0.0
					col_copy[s] = p
					col_name += '_KO'+ str(col_copy[s])
				if exp_df[column][s] == col_copy[s]:
					print 'Continuing for ' + all_ids[s][1] + ' : ' + str(col_copy[s])
					continue
				
				col_copy.name = col_name
				col_loc = exp_df_copy.columns.get_loc(column)
				exp_df_copy.insert(col_loc+ofs, col_name, col_copy)
				ofs+=1
	
	for par in plist:
		for column in exp_df:
			ofs = 1
			for p in prtvals:
				col_name = column.replace('_','-') + '_' + par.replace('_','-')
				col_copy = exp_df[column].copy()
				if args.perc:
					col_copy[par]+=col_copy[par]*p/100
					col_name += '_percKO'+ str(p)
				else:
					col_copy[par] += p
					col_name += '_KO'+ str(p)
				if exp_df[column][par] == col_copy[par]:
					print 'Continuing for ' + par + ':' + str(p)
					continue
				
				col_copy.name = col_name
				col_loc = exp_df_copy.columns.get_loc(column)
				exp_df_copy.insert(col_loc+ofs, col_name, col_copy)
				ofs+=1

	exp_df_copy.to_csv(exp_file_name, sep = '\t')

	return exp_df_copy, exp_file_name

def modifyCostFunct(samples, cost_data, args, outname):	

	cost_file_name = args.out_dir + '/'+ \
			outname + '_cost_func.txt'
	f  = open(cost_file_name, 'w')
	for cost_pair in cost_data:
		f.write('ID:{}\t{}\n'.format(cost_pair[0], cost_pair[1]))
		for sample in samples:
			f.write('{}\t1.0\n'.format(sample))
	f.close()

	return cost_file_name

def setUpK(args, model):
	bounds = (-1,1)
	if not args.opti_vect:
		varKParVect =  np.power(10, bounds[0] + \
			np.random.rand(len(model.variable_indexK_arr) + \
			model.dimKCD)*(bounds[1] - bounds[0]) ) 
	else:
		ids_vector, varKParVect, kcds_dict = \
							parseOptiFile(args.opti_vect)
		varKParVect = np.power(10,varKParVect)
		varIds = sorted([par_ids[i][0] for i in ids_vector 
				if not i.startswith('k')])
		if varIds  != model.variable_indexK_arr.tolist():
			print 'Problem with optmized vector!'
			pdb.set_trace()


	K = np.zeros(model.dimK)
	if len(model.variable_indexK_arr) > 0:
		K[model.variable_indexK_arr] = varKParVect

	return K

def runKOexperiments(exp_df_KO_f, cost_file_KO_f, args):

    print 'INFO Running knockout simulations.'
    control_conc = Odict()
    ko_conc = Odict()
    all_conc = Odict()
    sys.path.append(args.dir)
    import modeltools
    model_id = args.dir[args.dir[:-1].rindex('/')+1:-1] 
    mp = modeltools.ModelProcessor(maxStepNum=5000)

    model = PyBiosModel(model_id, 
    		identifiers_file, 
    		exp_df_KO_f, 
    		cost_file_KO_f)

    K = setUpK(args, model)
    fixed_Ki = np.setdiff1d(xrange(model.dimK), model.variable_indexK_arr)
    cost_sample_ids = set(itertools.chain(*model.conditionDict.keys()))

    for sampleId in sorted(cost_sample_ids):
        sample = model.samples[sampleId]
        if len(fixed_Ki)  > 0:
            K[fixed_Ki] = sample.K[fixed_Ki]   
            res = mp.simulate(sample.S0, 1e8, sample.F, K)
            if not res.success:
               print ('Simulation failed: flag= %s, sample= %s'
                  %( res.exitcode, sampleId))
               continue	
            if not 'KO' in sampleId:
                control_conc[sampleId] = res.finalState    
            else:
                ko_conc[sampleId] = res.finalState
            all_conc[sampleId] = res.finalState

    return control_conc, ko_conc, all_conc

def analyseDataReadOut(read_out, ko_conc, th):
	significantRO = []
	print 'INFO Analyse read out data.'
	for k in ko_conc.keys():
			ko_n = k
			ko_v = read_out[k]
			control_n = k.split('_')[0]
			control_v = read_out[control_n]
			if float(control_v)==0.0:
				fold = np.inf
			else:
				fold = float(ko_v)/control_v
			if fold > th or fold < 1/float(th):
				significantRO.append(
					{'Sample': control_n,
					 'Knockout of': k.split('_')[1], 
					 'ko_value': k.split('_')[2].replace('KO',''),
					 'read_out_contrl_val': control_v, 
					 'read_out_exper_val':ko_v, 
					 'fold_change': fold})
				print control_n, k, control_v,  ko_v, fold
	significantRO_df = pd.DataFrame(significantRO, 
				columns = ['Sample', 'Knockout of', 'ko_value',
		 					'read_out_contrl_val', 'read_out_exper_val','fold_change'])
	return significantRO_df

def analyseDataSA(all_conc, ko_conc, th):
	significantSA = []
	print 'INFO Analyse sensitivity data.'
	diff_var_rev = {v[0]: v[1] for k, v in diff_var_ids.iteritems()}
	for k,v in ko_conc.iteritems():
			control_n = k.split('_')[0]
			control_v = all_conc[control_n]
			fold = np.zeros(len(v))
			for i,j in enumerate(v):
				if control_v[i] == 0.0: 
					fold[i] = np.inf if v[i]!= 0.0 else 1
				else:
					fold[i] = v[i]/float(control_v[i])
			for i,j in enumerate(fold):
			 	if j > th or j < 1/float(th):
					significantSA.append(
						{'Sample': control_n,
						'Knockout of': k.split('_')[1], 
						'ko_value': k.split('_')[2].replace('KO',''),
          				'species-affected': diff_var_rev[i], 
          				'control_val': control_v[i], 
          				'experiment_val': v[i], 
          				'fold_change': j})
	significantSA_df = pd.DataFrame(significantSA, 
			columns = ['Sample', 'Knockout of', 'ko_value',
					   'species-affected', 'control_val', 
					   'experiment_val', 'fold_change'])

	return significantSA_df

def getReadOutsVals(all_conc, cost_data):
	print 'INFO Get read outs data.'
	read_out = Odict()
	for cost_pair in cost_data:
		pybios_ids = re.findall('(\[\d+\])', cost_pair[1])
		list_ids = [all_ids[i][0] for i in pybios_ids]
		for c, v in all_conc.iteritems():
			if not c in read_out:
				read_out[c] = 0.0
			formula =  re.sub('(\[\d+\])', lambda x: 'v['+str(all_ids[x.group()][0])+']', cost_pair[1])
			read_out[c] += eval(formula)

	return read_out

def runWorkflow(exp_df_KO_f, cost_func_KO_f, args):
	control_conc, ko_conc, all_conc = \
			runKOexperiments(exp_df_KO_f, cost_func_KO_f, args)
	if args.repeat > 1:
		control_conc_df = pd.DataFrame(control_in)
		ko_conc_df = pd.DataFrame(ko_in)
		all_conc_df = pd.DataFrame(all_in)
		for i in xrange(1, args.repeat):
			print "Starting "+ i + "-th repeath..."
			control_conc, ko_conc, all_conc = \
				runKOexperiments(exp_df_KO_f, cost_func_KO_f, args)
			if (len(control_conc[control_conc.keys()[0]])>1):
				control_conc = {(key, index): v 
								for key, val in control_conc.iteritems()
								for index, v in enumerate(val) }
				ko_conc = {(key, index): v 
							for key, val in ko_conc.iteritems()
							for index, v in enumerate(val) }
				all_conc = {(key, index): v 
							for key, val in all_conc.iteritems() 
							for index, v in enumerate(val) }
			control_conc_df.add(pd.DataFrame(control_conc))
			ko_conc_df.add(pd.DataFrame(ko_conc_df))
			all_conc_df.add(pd.DataFrame(all_conc_df))

		control_conc_df.apply(lambda x: x/args.repeat)
		ko_conc_df.apply(lambda x: x/args.repeat)
		all_conc_df.apply(lambda x: x/args.repeat)
		all_conc = Odict(sorted(all_conc_df.to_dict(orient='list').iteritems()))
		control_conc = Odict(sorted(control_conc_df.to_dict(orient='list').iteritems()))
		ko_conc = Odict(sorted(ko_conc_df.to_dict(orient='list').iteritems()))

	return all_conc, control_conc, ko_conc

arg_parser = argparse.ArgumentParser(description='Run knockout experiments.') 

arg_parser.add_argument('--dir',
            type = existingDir,
            default = './',
            help = 'Path to the model directory were'\
            'experiment table and cost file are present.')
arg_parser.add_argument('--sp',
            nargs = '+',
            type = str,
            default = [],
            help = 'List of species to be pertrubed.')
arg_parser.add_argument('--genes',
			action = 'store_true',
			default = False,
            help = 'Whether to pertrube all genes.'\
           	'Default: False.')
arg_parser.add_argument('--rpkms',
			action = 'store_true',
			default = False,
            help = 'Whether to pertrube all RPKM values."\
            " Default: False.')
arg_parser.add_argument('--vals',
            nargs = '+',
            type = float,
            default = [],
            help = 'List of perturbation values.')	
arg_parser.add_argument('--perc',
			action = 'store_true',
			default = False,
	        help = 'Whether to calculate percentage of the values.'\
	        'Default: False.')
arg_parser.add_argument('--opti-vect',
			type = existingFilePath,
	        default = None,
	        help = 'Path to optimized vector file. If none '\
	        'random values for parameters will be used.' )
arg_parser.add_argument('--threshold',
			nargs='?',
	        type = float,
	        default = 1,
	        help = 'Set the threshold'\
	        'for output results. Default 1 - print all.' )
arg_parser.add_argument('--repeat',
			nargs='?',
	        type = int,
	        default = 1,
	        help = 'How many times to repeat the study '\
	        'with new random initial parameters.' )
arg_parser.add_argument('--only-sim',
			action = 'store_true',
	        default = False,
	        help = 'Whether to run forward simulations'\
	        ' on the given exp_table.' )
arg_parser.add_argument('--out-dir',
			nargs='?',
	        type = existingDir,
	        help = 'Directory to store the new experiment'
	        'table and cost function.' )


args = arg_parser.parse_args()

if not args.out_dir:
 	args.out_dir = args.dir

exp_df, cost_data, outname = read_data(args.dir)

if args.genes:
        outname+='genes_'
if args.rpkms:
        outname+='rpkms_'
if args.perc:
        outname+='perc_'
if args.opti_vect:
        outname+='opti_vect_'
if al.only_sim:
        outname+='only_sim_'

exp_df_KO, exp_df_KO_f = modifyExpTable(exp_df, args, outname)
cost_func_KO_f = modifyCostFunct(exp_df_KO.columns.tolist(), cost_data, args, outname)

all_conc, control_conc, ko_conc = runWorkflow(exp_df_KO_f, cost_func_KO_f, args)
read_out = getReadOutsVals(all_conc, cost_data)
significantRO = analyseDataReadOut(read_out, ko_conc, args.threshold)
significantSA = analyseDataSA(all_conc, ko_conc, args.threshold)

print 'INFO Saving data to output files with affix:' + outname
for c in set(significantSA['ko_value'].tolist()):
	tmp = significantSA[significantSA['ko_value']==c]
	tmp.to_csv(args.dir + outname + c + '_SA.csv', sep = '\t')
outname+=('_').join([str(v) for v in args.vals])+'_'

# significantSA.to_csv(args.dir + outname + 'SA.csv', sep = '\t')
significantRO.to_csv(args.dir + outname + 'RO.csv', sep = '\t')


