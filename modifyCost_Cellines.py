import pdb
import os
import sys


def main(cost_file, celllines):

	with open(cost_file, 'r') as fileIn:
		with open(cost_file[:-4]+'_modi_cellFab.txt', 'w') as fileOut:
			for line in fileIn:
				# modify if input are to ommit or to include
				if any(s for s in celllines if s in line):
					fileOut.write(line)
	return

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2:])