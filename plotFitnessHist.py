import pandas as pd
import numpy as np
from collections import OrderedDict as Odict
import matplotlib.pyplot as plt
import collections
import sys 
import pdb
import random
import os
import seaborn as sns
# plt.style.use('seaborn-deep')
'''
Plots fitness distribution of given random simulation and
optimized vectors file

'''

def mad_based_outlier(points, thresh=3.5):
	if len(points.shape) == 1:
		points = points[:,None]
	median = np.median(points, axis=0)
	diff = np.sum((points - median)**2, axis=-1)
	diff = np.sqrt(diff)
	med_abs_deviation = np.median(diff)

	modified_z_score = 0.6745 * diff / med_abs_deviation

	return modified_z_score > thresh

def parse_file(input_file, only_last=False):

	
	fitness_list = []
	
	with open(input_file) as file:
		for line in file:
			if "FITNESS:" in line:
				fitness_list.append(eval(line[line.index(":")+1:]))
	if only_last:
		return np.array(fitness_list[-1])			
	else:
		return np.array(fitness_list)


def main(argv):

	random_data_file = argv[0]
	optimized_data_file = argv[1]

	if argv[2] and argv[2].lower() in ('yes', 'true', 't', 'y', '1'):
		only_last = True
	else:
		only_last = False

	random_fitness = parse_file(random_data_file, False)
	optimized_fitness = parse_file(optimized_data_file, only_last)

	outliers = random_fitness[mad_based_outlier(random_fitness)]

	modified_rand = random_fitness


	for i in outliers:
		if i > 1000:
			modified_rand = np.delete(modified_rand, np.where(modified_rand==i))

	fig, ax = plt.subplots(figsize=(1800/100, 1200/100), dpi=100)
	fig.suptitle("Distribution of the optimized and random fitness values.")

	# pdb.set_trace()
	if only_last:

		sns.distplot(modified_rand, bins=range(1, 1000, 10), \
			ax=ax, kde=False, \
			label = 'random_fitness', \
			hist_kws = {'edgecolor':'black', 'alpha': 0.7})
		plt.plot(optimized_fitness, 0.1, color='r', \
				marker = 'o', label='optimized fitness final')
	else:
		for (x, l) in [(modified_rand, 'random_fitness' ), (optimized_fitness, 'optimized_fitness')]:
			sns.distplot(x, bins=range(1, 1000, 10), \
				ax=ax, kde=False, \
			 	label = l, hist_kws = {'edgecolor':'black', 'alpha': 0.7})
	ax.set_xlim([0, 1000])
	plt.legend(loc="upper right")
	plt.draw()
	plt.show()
	
	return


if __name__ == '__main__':
	try:
		main(sys.argv[1:])
	except:
		print "Prints plots of fitness distribution from random simulated file and simtools optmized file"
		print "Input arguments: random simulated file (progress_trace like struct), progres_trace.txt and flag whether to analyse only the last vector (yes, y, t, true or 1)"
		print "Usage ex.: python plotFitnessHist.py simulated_results.txt progress_trace.txt 1"