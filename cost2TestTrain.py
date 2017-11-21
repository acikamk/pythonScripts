import collections 
from collections import OrderedDict as Odict
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from os import listdir
from os.path import isfile, join
import compiler
import re
import random
import pdb

# cost_function_file = "/home/H0001-1/a.kovachev/simtools/data/models/Human_Model_0417/cost_function.txt"
# output_name = "cost_func"
# n_tests = 5
# CellDivisionFunction	(k0*[131]+k1*[1293]+k2*[130]+k3*[326]+k4*[902]+k5*[115]+k6*[325]+k7*[323]+k8*[829]+k9*[127]+k10*[919]+k11*[129])/(1+k12*[18]+k13*[19]+k14*[20]+k15*[21])
def getSampleName(sample):
	# start = sample.split("/")[0].index('-')
	# end = sample.split("/")[0].rindex('-01')
	# control = sample[start+1:end+2]
	control = sample[0:sample.split("/")[0].index('_')]
	# import pdb; pdb.set_trace()
	return control

def saveData(output_name, cost_func_dict, i, training_list, testing_list):

	f_train = open("./data/" + output_name + "_train_"+ str(i) + ".txt", 'w') 
	f_test = open("./data/" + output_name + "_test_"+ str(i) + ".txt", 'w') 
	for k in cost_func_dict.keys():

		f_train.write("ID:" + k.split(':')[0] + '\t' + k.split(':')[1] + "\n")
		f_test.write("ID:" + k.split(':')[0] + '\t' + k.split(':')[1] + "\n")
		for key, value in cost_func_dict[k].items():
			if getSampleName(key) in training_list:
				f_train.write(key + "\t" + value )
			elif getSampleName(key) in testing_list:
				f_test.write(key + "\t" + value )
			else:
				import pdb; pdb.set_trace()
				print "Key Error :)! " +key +"\n"
	f_train.close() 
	f_test.close()


def main(cost_function_file, output_name, n_tests, cv=True, shuffle_one=False):
	n_tests = int(n_tests)
	controls = set()
	# conc_dict_samples = Odict()
	cost_func_dict = Odict()
	f = open(cost_function_file)
	for line in f:
		if line.startswith("ID"):
			cost_id = line.split("\t")[0].replace("ID:","")
			cost_val = line.split("\t")[1].strip()
			cost_func_dict[cost_id+":"+cost_val] = Odict()
		else:
			sample = line.split("\t")[0]
			value =  "\t".join(line.split("\t")[1:])
			if "/" in sample:
				drug = sample.split("/")[0]
				# pdb.set_trace()
				# control =  getSampleName(sample.split("/")[1])
				control = sample.split("/")[1]
			elif "CONTROL" in sample:
				control = getSampleName(sample)
			elif not '_' in sample:
				control = getSampleName(sample)
			controls.add(control)
			cost_func_dict[cost_id+":"+cost_val][sample] = value
	f.close()
	print len(controls)
	controls_list = list(controls)
	bucket_size = len(controls)

	if cv:
		testing_size = int(round(bucket_size/n_tests)) if not shuffle_one else 1
		ids = range(0, n_tests*testing_size, testing_size)
		for i, start in enumerate(ids):
			end = len(controls) if i+1>=len(ids) else ids[i+1]
			testing_list = controls_list[start:end]
			training_list = list(controls-set(testing_list))
			# print results to file
			saveData(output_name, cost_func_dict, i, training_list, testing_list)

	else:
		training_size = int(round(bucket_size*80/100)) if not shuffle_one else bucket_size-1
		for i in xrange(n_tests):
			training_list = random.sample(controls, training_size)
			testing_list = list(controls-set(training_list))
			# print results to file
			saveData(output_name, cost_func_dict, i, training_list, testing_list)

	return

if __name__ == "__main__":
	if len(sys.argv) == 5:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	if len(sys.argv) == 6:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
	else:
		main(sys.argv[1], sys.argv[2], sys.argv[3])