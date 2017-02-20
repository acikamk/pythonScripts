import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
from libsbml import *
"""
Script that checks the consistency between PyBios identifiers.py file and 
sbml model file. This version only checks the matching number of species,
constants and parameters and reports which values are missmatched. 
#### TO DO: Report which values are missing/missmatched. #####
"""
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def main():

	try:
		identifiers_file = sys.argv[1]
		sbml_file = sys.argv[2]
	except IndexError:	
		print bcolors.WARNING\
		 + "Plese give identifiers.py ans sbml file as input parameters!"\
		 + bcolors.ENDC
		return
	# identifiers file anaysis
	f = open(identifiers_file)
	try:
		identifiers_code = compile(f.read(), identifiers_file, 'exec')
		f.close()
		eval(identifiers_code)

		diff_species = locals()['diff_var_ids']
		fixed_species = locals()['fixed_var_ids']

		fixed_sp_dict = {j[1]: k for k,j in fixed_species.iteritems()}
		diff_sp_dict =  {j[1]: k for k,j in diff_species.iteritems()}        
		pars_sp_dict = locals()['par_ids']
	except:
		print  bcolors.FAIL \
		 + "The identifiers.py file has incompatile structure!" \
		 + bcolors.ENDC
		return

	# sbml file analysis
	try:
		document = SBMLReader().readSBML(sbml_file)
		if document.getNumErrors():
			print  "{}Warning: The SBML model contains {} errors! {}"\
			.format(bcolors.WARNING, document.getNumErrors(), bcolors.ENDC)

		model = document.getModel()


		diff_sb_dict = {i.getId(): i.getName() \
				 for i in model.getListOfSpecies() if not i.getConstant()}
		fixed_sb_dict = {i.getId(): i.getName() \
				for i in model.getListOfSpecies() if i.getConstant()}
		pars_sb_dict = {i.getId(): i.getName() \
				for i in model.getListOfParameters()}

		if not pars_sb_dict:
			pars_lists = [ i.getKineticLaw().getListOfParameters()\
					for i in model.getListOfReactions()]
			list_pars = [j.getListOfAllElements()\
					for j in pars_lists]
			pars_sb_dict = {i.getId(): i.getName()\
					for l in list_pars for i in l}
	except:
		print bcolors.FAIL \
		+ "The SBML file was not processed correctly!"\
		+ bcolors.ENDC
		return

	# preform check
	if len(diff_sb_dict.keys())!=len(diff_sp_dict.keys()):
		if len(fixed_sb_dict.keys())!=len(fixed_sp_dict.keys()):
			if len(pars_sb_dict.keys())!=len(pars_sp_dict.keys()):
				print bcolors.WARNING \
				+ "PROBLEM! The identifiers.py and model.sbml file have different number of species, constants and parameters!" \
				+ bcolors.ENDC
			else:
				print bcolors.WARNING\
				 + "PROBLEM! The identifiers.py and model.sbml file have different number of species and constants!"\
				 + bcolors.ENDC
		else:
			print bcolors.WARNING \
			 + "PROBLEM! The identifiers.py and model.sbml file have different number of species!"\
			 + bcolors.ENDC
	elif len(fixed_sb_dict.keys())!=len(fixed_sp_dict.keys()):
		if len(pars_sb_dict.keys())!=len(pars_sp_dict.keys()):
			print bcolors.WARNING \
			+ "PROBLEM! The identifiers.py and model.sbml file have different number of constants and parameters!"\
			+ bcolors.ENDC
		else:
			print bcolors.WARNING \
			+ "PROBLEM! The identifiers.py and model.sbml file have different number of constants!"\
			+ bcolors.ENDC
	elif len(pars_sb_dict.keys())!=len(pars_sp_dict.keys()):
		print  bcolors.WARNING \
		+ "PROBLEM! The identifiers.py and model.sbml file have different number of parameters!"\
		+ bcolors.ENDC
	else:
		print  bcolors.OKGREEN \
		+ "OK! The identifiers.py and model.sbml file have the SAME number of species, constants and parameters!"\
		+ bcolors.ENDC
	return		

if __name__ == "__main__":
    main()