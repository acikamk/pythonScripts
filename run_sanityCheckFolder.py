import sanityCheck
import sys
import os
"""
Runs sanity check on the sbml and pybios models folders
"""
file_list=[]
sbml_path = "/project/V0001-1/modcell/modcell_data/models/sbml_pybios_models/"
pybios_path = "/project/V0001-1/modcell/modcell_data/models/python_pybios_models/"
for file in os.listdir(sbml_path):
    if file.endswith(".sbml"):
        sbml_file = sbml_path + file
        id_file = pybios_path + file[:file.index(".")] + "/identifiers.py"
        try:
        	os.path.isfile(sbml_file)
        	os.path.isfile(id_file)
        	print file
        	sanityCheck.main(id_file, sbml_file)
        except:
        	print "ERROR:" + file
        	continue