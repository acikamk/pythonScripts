import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import expit
import pdb
import time
import os
from os import listdir
from os.path import isfile, join
import sys

'''
Script that fits generalized logistic function to real drug response curve data.
Normalization of the activation across each cellline is performed before the fitting.
The first low drug concentrations are removed and standard deviations values
are not used in the predictions since these settings gave the best results.

Developed for CCLE dataset

Input: 
	- path to the folder that has the pharmacology data ex. for CCLE data 
		/project/V0001-1/modcell/modcell_data/modcellexpdb/CCLE/pharmacology/
	- cost formula (now is fixed value, can be set as input!!!!) 
		now defined for the Hasenauer, ERBB_RAS_AKT_DRUGS model
Example: python fitSigmoid.py /project/V0001-1/modcell/modcell_data/modcellexpdb/CCLE/pharmacology/
Output:
	- Separate folder for each drug is created. In each folder plots for the real and 
	fitted drug response curves are drawn. Failed optimizations (fittings) are marked with 
	'Fail-' in the name. Results .csv file is given for each drug and cost function file is generated
	only including the drug cellline combinations which have rmse smaller than the median value.

NOTE: I try to authomatically transform the sample names as outputed from the Project Manager,
however this needs to be checked and maybe addapted.
	Drugs are not mapped correclty here!!!! probably the drug name should come form the folder name
	or some mapping dictionary!

'''

# create class that will holds all the data:
# from file: sample, name, drug name, list of dosages, list of activation values, stdeviation
# generated: predicted-fitted activiity, rmse error and mae error of fitting
# simulated: drug doses and simulated activations

punctuations = '''!()-[]{};:'"\,<>./?@#$%^&*_~ '''

class CelllineData:
	def __init__(self, sample, drug, dose, actv, stdev, \
			pred_actv=np.inf, rmse=np.inf, mae=np.inf, \
			sim_dose = [], sim_actv=[], nactv=np.inf):

		self.sample = ''.join(ch.upper() for ch in sample if ch not in punctuations)
		self.drug = drug
		self.dose = dose
		self.actv = actv
		self.stdev = stdev

	def setPredActv(self, pred_actv):
		self.pred_actv = pred_actv
		self.rmse = np.sqrt(np.square(np.subtract(self.nactv, \
			self.pred_actv)).mean())
		self.mae = np.abs(np.subtract(self.nactv, self.pred_actv)).mean()

	def setSimData(self, sim_dose, sim_actv):
		self.sim_actv = sim_actv
		self.sim_dose = sim_dose

# the used sigmoid like function for fitting
def sigmoid(x, k, q, b):
	return 1-((k-1)/(1+q*np.exp(-b*x)))

def existingFilePath(file_path):
	file_path = os.path.abspath(os.path.expanduser(file_path))
	assert os.path.isfile(file_path), "File %s was not found" % file_path
	return file_path

def existingDir(dir_path):
	dir_name = os.path.abspath(os.path.expanduser(dir_path))
	assert os.path.isdir(dir_name), "Cannot find dir: %s" % dir_path
	return dir_name

def import_data(filename):
	data = []
	df=pd.read_table(filename, sep = '\t', low_memory=False)
	# df=df_org[['DrugName', 'BioSource', 'Doses','Activities','STDEVs']].copy()
	for index, row in df.iterrows(): 
		try:
			# remove the fist 2 values (low concentrations from the data)
			dose_list = np.array(row['Doses'].split(';')[2:], dtype=np.float_)
			actv_list = np.array(row['Activities'].split(';')[2:], dtype=np.float_)
			stdev_list = np.array(row['STDEVs'].split(';')[2:], dtype=np.float_)
		except:
			print('Error in parsing data!')
			pdb.set_trace()
		try:	
			data.append(CelllineData(row['BioSource'], row['DrugName'],\
			 dose_list, actv_list, stdev_list))
		except:
			print('Error in object generation!')
			pdb.set_trace()
	return data

def runAnalysis(data, path, drug):
	
	# data = import_data(file)
	try:
		os.mkdir(path)
	except OSError:
		print('ERROR: Directory with same name exists!')
		pdb.set_trace()

	results = []
	failed = []
	sim_dose = np.array([0.0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5,\
		 0.75, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])

	for el in data:
		try:
			# fit the function to the data and predict the parameters
			# predict activiity for the original doses
			# generate simulated data for the defined doses
			popt, pcov = curve_fit(sigmoid, el.dose, el.nactv)
			pred_actv = sigmoid(el.dose, *popt)
			sim_actv = sigmoid(sim_dose, *popt)
			# set the results
			el.setPredActv(pred_actv)
			el.setSimData(sim_dose, sim_actv)
			results.append(el)
			
			# check for failed fitting, and plot the real data	
		except RuntimeError:
			print("Optimization failed for sample {} and drug {}.".\
				format(el.sample, el.drug))
			plt.errorbar(el.dose*1000, el.nactv, yerr=el.stdev, label='data')
			ax = plt.gca()
			ax.set_xscale('log', nonposy="clip")
			plt.ylim(0,100)
			plt.xlabel('Drug concentration')
			plt.ylabel('Normalized activity')
			title = '{}-{}'.format(el.sample, el.drug).replace('/', '_')
			plt.title(title)
			plt.legend()
			plt.grid(True)
			plt.savefig(path +'/Fail-'+title+'.png', bbox_inches='tight')
			plt.close()
			failed.append(el)
			continue

		# plot results
		plt.errorbar(el.dose*1000, el.nactv, yerr=el.stdev, label='data')
		ax = plt.gca()
		ax.set_xscale('log', nonposy="clip")
		plt.semilogx(el.dose*1000, pred_actv, 'r-', \
			label='fit: k=%5.3f, q=%5.3f, b=%5.3f' % tuple(popt))
		plt.semilogx(el.sim_dose[2:]*1000, el.sim_actv[2:], 'r+')
		plt.ylim(0,100)
		title = '{}-{}'.format(el.sample, el.drug).replace('/', '_')
		plt.xlabel('Drug concentration')
		plt.ylabel('Normalized activity')
		plt.title(title)
		plt.grid(True)
		plt.legend()
		# pdb.set_trace()
		try:
			plt.savefig(path+'/'+title+'.png', bbox_inches='tight')
		except:
			pdb.set_trace()
		plt.close()
	
	# generate DataFrame to save results
	df_res=pd.DataFrame([vars(f) for f in results])
	df = df_res[['sample','drug', 'dose', 'actv', 'stdev', 'nactv', \
	'pred_actv','mae','rmse', 'sim_dose', 'sim_actv']]
	df.sort_values(by='rmse', inplace=True)
	df_cost = df[df['rmse']<df['rmse'].median()]
	df_cost.sort_values(by='sample', inplace=True)

	# create cost function file
	cost_file = open(path + '/cost_func_' + drug + '.txt', 'w')
	cost_formula = '(kCD0*[131]+kCD1*[1293]+kCD2*[130]+kCD3*[326]+kCD4*[902]'\
		'+kCD5*[115]+kCD6*[325]+kCD7*[323]+kCD8*[829]+kCD9*[127]+kCD10*[919]'\
		'+kCD11*[129])/(1+kCD12*[18]+kCD13*[19]+kCD14*[20]+kCD15*[21])\n'
	cost_file.write('ID:CellDivisionFunction\t' +  cost_formula)
	for index, row in df_cost.iterrows():
		# TUMOR-SW1417-cellline-01-01, example of the sample format for control
		control_name = 'TUMOR-' + row['sample'] + '-cellline-01-01'
		# control_value = round(row['sim_actv'][0], 5)
		# cost_file.write('{}\t{}\t{}\n'.format(control_name, control_value,'1.0'))
		# for i, dose in enumerate(row['sim_dose'][1:],1):
		for i, dose in enumerate(row['sim_dose']):
			# TUMOR-SW1417-cellline-01-01_PLX-4720_Conc2.5nM, example of the sample format for drug
			tumor_name = control_name + '_' + row['drug'] + '_Conc' + str(dose*1000) + 'nM'
			tumor_value = round(row['sim_actv'][i], 5)
			cost_file.write('{}\t{}\t{}\n'.format(tumor_name, tumor_value, '1.0'))
			# pdb.set_trace()
	cost_file.close()

	# prepare the DataFrame for nice output format
	df['dose'] = df['dose'].apply(\
		np.array_str).apply(lambda x: x.replace('\n', ''))
	df['actv'] = df['actv'].apply(\
		np.array_str).apply(lambda x: x.replace('\n', ''))
	df['stdev'] = df['stdev'].apply(\
		np.array_str).apply(lambda x: x.replace('\n', ''))
	df['nactv'] = df['nactv'].apply(\
		np.array_str).apply(lambda x: x.replace('\n', ''))
	df['pred_actv'] = df['pred_actv'].apply(\
		np.array_str).apply(lambda x: x.replace('\n', ''))
	df['sim_dose'] = df['sim_dose'].apply(\
		np.array_str).apply(lambda x: x.replace('\n', ''))
	df['sim_actv'] = df['sim_actv'].apply(\
		np.array_str).apply(lambda x: x.replace('\n', ''))

	df.to_csv(path+'/'+drug+'.csv', sep  = '\t', index=False)
	# df_cost.to_csv(path+'/median_rmse_'+drug+'.csv', sep  = '\t', index=False)

	return

if __name__ == '__main__':
	
	dir_path = existingDir(sys.argv[1])
	files = [f for f in listdir(dir_path) if isfile(join(dir_path, f)) and '.csv' in f]

	sample_dict = {}
	data_dict = {}
	for file in files:
		# for each file, read data and generate the class list
		sample_df= pd.read_table(join(dir_path, file), sep = '\t', \
		 usecols=['BioSource', 'Activities']) 
		file_name = join(dir_path, file)
		data = import_data(file_name)
		data_dict[file] = data

		# retrieve the data for normalization 
		# for each cellline (biosource) append all possible activity values
		for index, row in sample_df.iterrows():
			if not row['BioSource'] in sample_dict:
				sample_dict[row['BioSource']] = []
			sample_dict[row['BioSource']].extend(map(\
				float,  row['Activities'].split(';')[2:]))

	# normalize the activity values from 0-100 (this can be addapted ex:0-1)
	# set the normalized activity feature (nactv) which will be used for fitting
	# it is done separately for each activity value, might not be very efficient
	a = 0
	b = 100
	for file, data in data_dict.iteritems():
		drug = file.split('/')[-1].split('_')[1]
		print('Analysing drug:\t' + drug)
		path = '/project/V0001-1/modcell/users/a.kovachev/ThomasProject/'+ drug
		for el in data:
			# find the min and max for this sample and calc the normalized value
			x_min = min(sample_dict[el.sample])
			x_max = max(sample_dict[el.sample])
			el.nactv = np.array([round(a + (y - x_min)*(b-a)/(x_max-x_min),5) \
				for y in el.actv])

		runAnalysis(data, path, drug)
