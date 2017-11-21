import collections 
from collections import OrderedDict as Odict
from collections import deque 
from collections import Counter
import os
import sys
import argparse, ast
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
import re
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import pdb
sys.path.append('/home/H0001-1/a.kovachev/simtools/src')
from modelstuff import PyBiosModel
global exp_table
global cost_file 
global identifiers_file 
global cmodel_file

def existingFilePath(file_path):
	file_path = os.path.abspath(os.path.expanduser(file_path))
	assert os.path.isfile(file_path), 'File %s was not found' % file_path
	return file_path
def existingDir(dir_path):
    dir_name = os.path.abspath(os.path.expanduser(dir_path))
    assert os.path.isdir(dir_name), 'Cannot find dir: %s' % dir_path
    return dir_name

def newFileInExistingDir(file_path):
   dir_name = existingDir(os.path.dirname(file_path))
   file_name = os.path.basename(file_path)
   assert file_name, 'Empty file name for dir %s' % dir_name
   return os.path.join(dir_name, file_name)

def read_data(dir_model):
	print 'INFO Reading input data.'
	global exp_table
	global cost_file 
	global identifiers_file 
	global cmodel_file
	exp_table = dir_model + '/exp_table_mean.csv'
	cost_file = dir_model + '/cost_func_empty.txt'
	identifiers_file = dir_model + '/identifiers.py'

	exp_df = pd.read_table(exp_table, sep = '\t', low_memory = False, index_col = 0)
	exec('global diff_var_ids, fixed_var_ids, par_ids\n' + open(identifiers_file).read())
	
	all_ids = diff_var_ids.copy()
	all_ids.update(fixed_var_ids)

	with open(cost_file) as f:
		lines = f.readlines()
		cost_data = [(l.split('\t')[0].split(':')[1].strip(), l.split('\t')[1].strip()) \
			for l in lines if l.startswith('ID')] 

	return exp_df, cost_data, all_ids 

def parseOptiFile(opti_file):
	kcds_dict = Odict()
	opti_vector = Odict()
	kcds_dict = Odict()
	with open(opti_file) as file:
		for line in file:
			if 'VARIABLE PAR IDS:' in line:
				ids_vector = eval(line[line.index('['):])
			elif 'FITNESS:' in line:
				fitness = eval(line[line.index(':')+1:])
			elif 'VECTOR:' in line:
				opti_vector = eval(line[line.index(':')+1:])		
	for var in ids_vector:
		if var.startswith('k'):
			kcds_dict[var] = opti_vector[ids_vector.index(var)]			
	return ids_vector, opti_vector, kcds_dict

def modifyExpTable(exp_df, all_ids, args, outname):

	print 'INFO Generating experiment table.'
	exp_file_name = args.out_dir + '/' + outname +'_exp_table.csv'
	if not args.sp:
		not_in_exp = [i for i in args.sp if i not in exp_df.index.tolist()]

		assert not not_in_exp, ('Error in input data. Species: ' + not_in_exp + 
			'are missing in the experiment table definition.')
	splist = args.sp
	plist = []
	prtvals = args.vals
	if args.genes:
		splist += [i for (i,j) in all_ids.iteritems() if all_ids[i][5] == 'gene' ]
	if args.proteins:
		splist += [i for (i,j) in all_ids.iteritems() if all_ids[i][5] == 'protein' ]
	if args.rpkms:
		plist += [i for (i,j) in par_ids.iteritems() if 'RPKM' in i]
	# pdb.set_trace()
	exp_df_copy = exp_df.copy()
	for s in splist:
		for column in exp_df:
			offset = 1
			for p in prtvals:
				col_copy = exp_df[column].copy()
				if args.perc:
					col_copy[s] = col_copy[s]+col_copy[s]*p/100
					col_name = column + '_' + all_ids[s][1].replace(' ','-')+ '_percKO'+ str(p)
				else:
					if col_copy[s] == 0.0:
						col_copy[s] = 1.0
					elif col_copy[s] == 1.0:
						col_copy[s]=0.0
					else:	
						col_copy[s] = p
					col_name = column + '_' + all_ids[s][1].replace(' ','-')+ '_KO'+ str(col_copy[s])
				if exp_df[column][s] == col_copy[s]:
					print 'Continue for ' + s + ':' + str(p)
					continue
				
				col_copy.name = col_name
				exp_df_copy.insert(exp_df_copy.columns.get_loc(column)+offset, col_name, col_copy)
				offset+=1
	
	for par in plist:
		for column in exp_df:
			offset = 1
			for p in prtvals:
				col_copy = exp_df[column].copy()
				if args.perc:
					col_copy[par] =col_copy[par]+col_copy[par]*p/100
					col_name = column + '_' + par + '_percKO_'+ str(p)
				else:
					col_copy[par] = p
					col_name = column + '_' + par+ '_KO_'+ str(p)
				if exp_df[column][par] == col_copy[par]:
					print 'Continue for ' + par + ':' + str(p)
					continue
				
				col_copy.name = col_name
				exp_df_copy.insert(exp_df_copy.columns.get_loc(column)+offset, col_name, col_copy)
				offset+=1

	exp_df_copy.to_csv(exp_file_name, sep = '\t')
	#pdb.set_trace()
	return exp_df_copy, exp_file_name

def modifyCostFunct(samples, cost_data, args, outname):	

	cost_file_name = args.out_dir + '/'+ outname+ '_cost_func.txt'
	f  = open(cost_file_name, 'w')
	for cost_pair in cost_data:
		f.write('ID:{}\t{}\n'.format(cost_pair[0], cost_pair[1]))
		for sample in samples:
			f.write('{}\t1.0\n'.format(sample))
	f.close()

	return cost_file_name

def runKOexperiments(exp_df_KO_f, cost_file_KO_f, args):

    print 'INFO Running knockout simulations.'
    control_conc = Odict()
    ko_conc = Odict()
    all_conc = Odict()
    sys.path.append(args.dir)
    import modeltools
    model_id = args.dir[args.dir[:-1].rindex('/')+1:-1] 
    mp = modeltools.ModelProcessor(maxStepNum=5000)
    time_vect = np.arange(0, float(1e9)+1 , 1e9/100.0, dtype=np.float)  

    cmodel = PyBiosModel(model_id, identifiers_file, exp_df_KO_f, cost_file_KO_f)
    bounds = (-1,1)
    if not al.opti_vect:
        varKParVect =  np.power(10, bounds[0] + \
        np.random.rand(len(cmodel.variable_indexK_arr) + \
            cmodel.dimKCD)*(bounds[1] - bounds[0]) ) 
    else:
    	ids_vector, varKParVect, kcds_dict = parseOptiFile(al.opti_vect)
    	varKParVect = np.power(10,varKParVect)
    	varIds = sorted([par_ids[i][0] for i in ids_vector if not i.startswith('k')])
    	if varIds  != cmodel.variable_indexK_arr.tolist():
    		print 'Problem with optmized vector!'
    		pdb.set_trace()

    fixed_Ki = np.setdiff1d(xrange(cmodel.dimK), cmodel.variable_indexK_arr)
    K = np.zeros(cmodel.dimK)
    if len(cmodel.variable_indexK_arr) > 0:
        K[cmodel.variable_indexK_arr] = varKParVect
    # print K[0:10]    
    cost_func_sample_ids = set(itertools.chain(*cmodel.conditionDict.keys()))
    # pdb.set_trace()
    for sampleId in sorted(cost_func_sample_ids):
        sample = cmodel.samples[sampleId]
        if len(fixed_Ki)  > 0:
            K[fixed_Ki] = sample.K[fixed_Ki] 
            if not args.only_sim:  
            	# pdb.set_trace()
                res = mp.simulate(sample.S0, 1e8, sample.F, K)
                if not res.success:
                   print ('Simulation time course failed with flag %s for sample %s'\
                      %( res.exitcode, sampleId))
                   continue	
		if not 'KO' in sampleId:
                    control_conc[sampleId] = res.finalState    
                else:
                    ko_conc[sampleId] = res.finalState
                all_conc[sampleId] = res.finalState
            else:
                res = mp.simulate_time_course(sample.S0, time_vect, sample.F, K)
                if not res.success:
                    print ('Simulation time course failed with flag %s for sample %s'\
                      %( res.exitcode, sampleId))
 		    continue
                if not 'KO' in sampleId:
                    control_conc[sampleId] = res.timeCourse    
                else:
                    ko_conc[sampleId] = res.timeCourse
                all_conc[sampleId] = res.timeCourse
    # pdb.set_trace()
    return control_conc, ko_conc, all_conc

def analyseDataReadOut(read_out):
	significantRO = []
	print 'INFO Analyse read out data.'
	for k,v in read_out.iteritems():
		if not 'KO' in k:
			control_n = k
			control_v = v
			# significantRO[control_n] = []
			continue
		elif control_n in k:
			fold = control_v/float(v)
			if fold > 2 or fold < 0.5:
				significantRO.append({'Sample': control_n, 'Knockout': k.split('_')[1], 'ko_value': k.split('_')[2],
          			'org_val': control_v, 'exp_val':v, 'fold_change': fold})
				# significantRO[control_n].append((k, control_v, v, fold))
				print control_n, k, control_v,  v, fold
		else:
			print 'Problem'
			pdb.set_trace()
	significantRO_df = pd.DataFrame(significantRO, columns = ['Sample',  'Knockout', 'ko_value',\
			 'org_val', 'exp_val', 'fold_change'])
	return significantRO_df

def analyseDataSA(all_conc, all_ids):
	significantSA = []
	print 'INFO Analyse sensitivity data.'
	all_ids_rev = {v[0]: v[1] for k, v in all_ids.iteritems()}
	for k,v in all_conc.iteritems():
		if not 'KO' in k:
			control_n = k
			control_v = v
			# significantSA[control_n] = Odict()
			continue
		elif control_n in k:
			fold = [(i+1e-8)/(j+1e-8) for (i,j) in zip(control_v, v)]
			# s_fold = [(all_ids_rev[i], control_v[i], v[i], j) for i,j in enumerate(fold) if j>2 or j<0.5]
			# significantSA[control_n][k] = s_fold
			for i,j in enumerate(fold):
			 	if j>2 or j<0.5:
					significantSA.append({'Sample': control_n, 'Knockout': k.split('_')[1], 'ko_value': k.split('_')[2],
          			'species-affected': all_ids_rev[i], 'org_val': control_v[i], 'exp_val': v[i], 'fold_change': j})

	significantSA_df = pd.DataFrame(significantSA, columns = ['Sample',  'Knockout', 'ko_value','species-affected',\
			 'org_val', 'exp_val', 'fold_change'])
	return significantSA_df

def getReadOutsVals(all_conc, cost_data, all_ids):
	print 'INFO Get read outs data.'
	read_out = Odict()
	for cost_pair in cost_data:
		pybios_ids = re.findall('(\[\d+\])', cost_pair[1])
		list_ids = [all_ids[i][0] for i in pybios_ids ]
		for c, v in all_conc.iteritems():
			read_out[c] = sum([v[i] for i,j in enumerate(v) if i in list_ids]) 
	return read_out

def runRepeat(control_in, ko_in, all_in, exp_df_KO_f, cost_func_KO_f, al):
	control_conc_df = pd.DataFrame(control_in)
	ko_conc_df = pd.DataFrame(ko_in)
	all_conc_df = pd.DataFrame(all_in)
	for i in xrange(1, al.repeat):
		print i
		control_conc, ko_conc, all_conc = runKOexperiments(exp_df_KO_f, cost_func_KO_f, al)
		# pdb.set_trace()	
		if (len(control_conc[control_conc.keys()[0]])>1):
			control_conc = {(key, index): v for key, val in control_conc.iteritems() for index, v in enumerate(val) }
			ko_conc = {(key, index): v for key, val in ko_conc.iteritems() for index, v in enumerate(val) }
			all_conc = {(key, index): v for key, val in all_conc.iteritems() for index, v in enumerate(val) }
		control_conc_df.add(pd.DataFrame(control_conc))
		ko_conc_df.add(pd.DataFrame(ko_conc_df))
		all_conc_df.add(pd.DataFrame(all_conc_df))

	control_conc_df.apply(lambda x: x/al.repeat)
	ko_conc_df.apply(lambda x: x/al.repeat)
	all_conc_df.apply(lambda x: x/al.repeat)
	all_conc = Odict(sorted(all_conc_df.to_dict(orient='list').iteritems()))
	control_conc = Odict(sorted(control_conc_df.to_dict(orient='list').iteritems()))
	ko_conc = Odict(sorted(ko_conc_df.to_dict(orient='list').iteritems()))

	return all_conc, control_conc, ko_conc

arg_parser = argparse.ArgumentParser(description='Run knockout experiments.') 

arg_parser.add_argument('--dir',
            type = existingDir,
            default = './',
            help = 'Path to the model directory were experiment table and cost file are present.')
arg_parser.add_argument('--sp',
            nargs = '+',
            type = str,
            default = [],
            help = 'List of species to be pertrubed.')
arg_parser.add_argument('--genes',
			action = 'store_true',
			default = False,
            help = 'Whether to pertrube all genes. Default: False.' )
arg_parser.add_argument('--proteins',
			action = 'store_true',
			default = False,
            help = 'Whether to pertrube all proteins. Default: False.')
arg_parser.add_argument('--rpkms',
			action = 'store_true',
			default = False,
            help = 'Whether to pertrube all RPKM values. Default: False.')
arg_parser.add_argument('--vals',
            nargs = '+',
            type = float,
            default = [],
            help = 'List of pertrubation values.')	
arg_parser.add_argument('--perc',
			action = 'store_true',
			default = False,
	        help = 'Whether to callulate percentage or absoulte values. Default: False.')
arg_parser.add_argument('--opti-vect',
			type = existingFilePath,
	        default = None,
	        help = 'Whether to use optimized vector as input or random values.' )
arg_parser.add_argument('--repeat',
			nargs='?',
	        type = int,
	        default = 1,
	        help = 'How many times to repeat the study with new random initial parameters.' )
arg_parser.add_argument('--out-dir',
			nargs='?',
	        type = existingDir,
	        help = 'Directory to store the new experiment table and cost function.' )
arg_parser.add_argument('--only-sim',
			action = 'store_true',
	        default = False,
	        help = 'Whether to run forward simulations on the given exp_table.' )

al = arg_parser.parse_args()

if not al.out_dir:
 	al.out_dir = al.dir

exp_df, cost_data, all_ids = read_data(al.dir)
# pdb.set_trace()

outname = '/KO_mean_'
if al.genes:
        outname+='genes_'
if al.proteins:
        outname+='proteins_'
outname+=('_').join([str(v) for v in al.vals])+'_'
if al.rpkms:
        outname+='rpkms_'
if al.perc:
        outname+='perc_'
if al.opti_vect:
        outname+='opti_vect_'
if al.only_sim:
        outname+='only_sim_'

if al.only_sim:
	exp_df_f = al.dir + '/exp_table.csv'
	cost_func_f = modifyCostFunct(exp_df.columns.tolist(), cost_data, al)
	control_conc, ko_conc, all_conc = runKOexperiments(exp_df_f, cost_func_f, al)
	control_conc_r = {(key, index): v for key, val in control_conc.iteritems() for index, v in enumerate(val) }
	ko_conc_r = {(key, index): v for key, val in ko_conc.iteritems() for index, v in enumerate(val) }
	all_conc_r = {(key, index): v for key, val in all_conc.iteritems() for index, v in enumerate(val) }
	# pdb.set_trace()
	if al.repeat > 1:
		all_conc, control_conc, ko_conc = runRepeat(control_conc_r, ko_conc_r, all_conc_r, exp_df_f, cost_func_f, al)
	all_conc_b = Odict()
	pdb.set_trace()
	for k in list(set(zip(*all_conc.keys())[0])):
		all_conc_b[k] = np.zeros((max(zip(*all_conc.keys())[1]), len(all_conc[(k, 0)])))
		for i in range(0, max(zip(*all_conc.keys())[1])):
			all_conc_b[k][i] = all_conc[(k,i)]
	#pdb.set_trace()
	set([all_conc_b[j][0:100][99].tolist().index(i) for j in all_conc_b.keys() for i in all_conc_b[j][0:100][99] if i>2000])
else:	
	# exp_df_KO, exp_df_KO_f = modifyExpTable(exp_df, all_ids, al, outname)
	# cost_func_KO_f = modifyCostFunct(exp_df_KO.columns.tolist(), cost_data, al, outname)
	exp_df_KO_f = al.dir +'/KO_mean_genes_1.0_opti_vect__exp_table.csv'
	cost_func_KO_f = al.dir + '/KO_mean_genes_1.0_opti_vect__cost_func.txt'
	control_conc, ko_conc, all_conc = runKOexperiments(exp_df_KO_f, cost_func_KO_f, al)
	#pdb.set_trace()
	if al.repeat > 1:
		all_conc, control_conc, ko_conc=runRepeat(control_conc, ko_conc, all_conc, exp_df_KO_f, cost_func_KO_f, al)

	read_out = getReadOutsVals(all_conc, cost_data, all_ids)
	significantRO = analyseDataReadOut(read_out)
	significantSA = analyseDataSA(all_conc, all_ids)

	significantSA.to_csv(al.dir + outname + 'SA.csv', sep = '\t')
	significantRO.to_csv(al.dir + outname + 'RO.csv', sep = '\t')
#pdb.set_trace()
# diff_var_ids['[374]'] = (11, '''Pik3ca-001:Pik3r1-201 [complex in cytoplasm]''', 0.0, 53332, {}, 'complex', []) 72
# sum_ko = [sum(v) for k,v in ko_conc.iteritems()]
# z=np.histogram(sum_ko)
# # plt.hist(sum_ko)
# # plt.show()
# sum_mean_ko = np.mean(sum_ko)
# sum_std_ko = np.std(sum_ko)


# sum_cont = [sum(v) for k,v in control_conc.iteritems()]
# t=np.histogram(sum_cont)
# # plt.hist(sum_cont)
# # plt.show()
# sum_mean = np.mean(sum_cont)
# sum_std = np.std(sum_cont)
