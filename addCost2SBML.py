import collections 
from collections import OrderedDict as Odict
from collections import deque 
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
from libsbml import *
import re
import pdb

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def recParam(math):
    leafs = []
    def _recParam(m):
        if m is not None:
            if m.getNumChildren() == 0 and m.isName(): 
                leafs.append(m.getName())
            for i in range(0, m.getNumChildren()):
                _recParam(m.getChild(i))
    _recParam(math)
    return leafs

		

def main(identifiers_file, sbml_file, cost_function_file):
	compartments = {'cellular component]' : '1',
				'extracellular]': '2',
				'golgi apparatus]': '11',
				'endoplasmic reticulum]':'12',
				'peroxisome]': '13',
				'plasma membrane]': '3',
				'intracellular]':'4',
				'cytoplasm]':'5',
				'nucleus]':'6',
				'nucleolus]':'7',
				'ribosome]':'8',
				'mitochondrion]':'9',
				'lysosome]':'10'}
	# identifiers file anaysis
	try:
		f = open(identifiers_file)
		identifiers_code = compile(f.read(), identifiers_file, 'exec')
		f.close()
		eval(identifiers_code)

		diff_sp_dict = locals()['diff_var_ids']
		fixed_sp_dict = locals()['fixed_var_ids']

		sbml_dict_diff =  {k:"SP_"+ str(j[3]) + "_" + \
		compartments[j[1].split(" in ")[1]] for k,j in diff_sp_dict.iteritems()}      
		
		sbml_dict_fixed =  {k:"SP_"+ str(j[3]) + "_" + \
		compartments[j[1].split(" in ")[1]] for k,j in fixed_sp_dict.iteritems()}  

		intersect_keys = [i for i in sbml_dict_fixed.keys() if i in sbml_dict_diff.keys()]
		if intersect_keys:
			print intersect_keys

		# sbml_dict_all = sbml_dict_fixed.copy()
		# sbml_dict_all.update(sbml_dict_diff) 
	except:
		print  bcolors.FAIL \
		 + "The identifiers.py file has incompatile structure!" \
		 + bcolors.ENDC
		return

	# sbml file analysis
	try:
		document = SBMLReader().readSBML(sbml_file)
		if document.getNumErrors():
			print  "{}Warning: The SBML model contains {} errors!{}"\
			.format(bcolors.WARNING, document.getNumErrors(), bcolors.ENDC)

		model = document.getModel()

		diff_sb_dict = {i.getId(): i.getName() \
				 for i in model.getListOfSpecies() if not i.getConstant()}
		fixed_sb_dict = {i.getId(): i.getName() \
				for i in model.getListOfSpecies() if i.getConstant()}
		pars_sb_dict = {i.getId(): i.getName() \
				for i in model.getListOfParameters()}


		# if not pars_sb_dict:
		# 	pars_lists = [ i.getKineticLaw().getListOfParameters()\
		# 			for i in model.getListOfReactions()]
		# 	list_pars = [j.getListOfAllElements()\
		# 			for j in pars_lists]
		# 	pars_sb_dict = {i.getId(): i.getName()\
		# 			for l in list_pars for i in l}
	except:
		print bcolors.FAIL \
		+ "The SBML file was not processed correctly!"\
		+ bcolors.ENDC
		return

	# parse cost function file and formula
	try:
		f = open(cost_function_file)
		lines = f.readlines()
		cost_lines = [l for l in lines if l.startswith("ID:")]
		for c in cost_lines:
			c_id = c.split()[0].split(':')[1].replace('|', '_')
			c_id_sbml = 'observable_' + c_id
			c_id_sigma = c_id_sbml + '_sigma'
			c_formula = c.split()[1].strip().replace('kCD','k')
			c_all_params = [c_id_sbml, c_id_sigma]
			try:
				if 'k' in c_formula and '[' in c_formula:
					c_sbml = re.sub('(k\d+)\*(\[\d+\])', lambda x: 'weight_' + c_id + '_' + sbml_dict_diff[x.group(2)] + '*' + sbml_dict_diff[x.group(2)] , c_formula)
				elif '[' in c_formula:
					c_sbml = re.sub('(\[\d+\])', lambda x: sbml_dict_diff[x.group()] , c_formula)
			except:
				print "Error in sub"
				pybios_ids = re.findall('(\[\d+\])', c_formula)
				print pybios_ids
			if not ("scaling_" in c_sbml or '/' in c_sbml):
				c_sbml = "scaling_"+ c_id + "_genotypespecific*(" + c_sbml + ")"
			elif '/' in c_sbml:
				c_sbml = "scaling_"+ c_id + "_genotypespecific*(" + c_sbml.split('/')[0] + ")/" +  c_sbml.split('/')[1]
			# pdb.set_trace()
			try: 
				sbml_math = parseL3FormulaWithModel(c_sbml, model) 
				params = []
				params = recParam(sbml_math)
				params = [e for e in params if e not in diff_sb_dict.keys() + fixed_sb_dict.keys() + pars_sb_dict.keys()]	
				c_all_params+=params
			except:
				print "Error in recParam or Model"
				pdb.set_trace()

			try:
				for par in c_all_params: 
					observable = model.createParameter()
					observable.setId(par)
					observable.setConstant(False)
					model.addParameter(observable)
				
				rule = model.createAssignmentRule()
				rule.setMath(sbml_math)
				rule.setVariable(c_id_sbml)
				model.addRule(rule)
			except:
				print "Error in rules"
				pdb.set_trace()
		try:
			document.checkL2v4Compatibility()
		except:	
			pdb.set_trace()
		# document.checkConsistency()

		writeSBMLToFile(document, sbml_file[:sbml_file.rindex('.')] + '_wCostwSP.sbml') 
	except:
		print bcolors.FAIL \
		+ "The cost function file was not processed correctly!"\
		+ bcolors.ENDC
		return

	return

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])