import pdb
import os
import sys

'''
Script that modifes cost function data based on a 
given list of celllines which can be used
All the celllines not in the list will be excluded 
from the data
Input:
	- cost function file
	- celllines which should be used
Example:
	python modifyCost_Cellines.py cost_func.txt 669V HC1996 BT20 C32
IMPROVEMENT: celllines can be input file, smth like
	f = open('celllines.txt',  'r')
	celllines = [el.strip() for el in f.readlines()]
'''


def main(cost_file, celllines):

	with open(cost_file, 'r') as fileIn:
		with open(cost_file[:-4]+'_modi_cellFab.txt', 'w') as fileOut:
			for line in fileIn:
				# modify if input are to ommit or to include
				if any(s for s in celllines if s in line):
					fileOut.write(line)
	return

if __name__ == "__main__":
	try:
		main(sys.argv[1], sys.argv[2:])
	except:
		print "Input:\n\
			\t- cost function file\n\
			\t - celllines which should be used\n\
		Example:\n\
			\t python modifyCost_Cellines.py cost_func.txt 669V HC1996 BT20 C32.